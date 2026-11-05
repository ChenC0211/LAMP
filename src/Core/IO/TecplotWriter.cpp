#include "Core/IO/TecplotWriter.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cmath>
#include <iterator>
#include <system_error>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>

namespace TecplotWriter
{
namespace
{
void require(bool condition, const char* reason)
{
    if (!condition) throw std::invalid_argument(std::string("TecplotWriter: ") + reason);
}

std::string quoted(const std::string& input)
{
    std::string result = "\"";
    for (unsigned char c : input)
    {
        require(c >= 32 && c != 127, "strings cannot contain control characters");
        if (c == '"' || c == '\\') result += '\\';
        result += static_cast<char>(c);
    }
    result += '"';
    // Each string is emitted on a single line; quoted strings cannot continue.
    require(result.size() <= 31900, "quoted string exceeds the ASCII line limit");
    return result;
}

bool identifier(const std::string& name)
{
    const auto alpha = [](unsigned char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
    };
    if (name.empty() || !alpha(static_cast<unsigned char>(name.front()))) return false;
    for (unsigned char c : name)
        if (!alpha(c) && !(c >= '0' && c <= '9')) return false;
    return true;
}

struct Header
{
    std::string title;
    std::vector<std::string> variables;
    std::vector<std::size_t> cellVariables;
    std::vector<std::string> auxiliary;
    const char* zoneType = nullptr;
};

Header validate(const Zone& zone)
{
    Header header;
    header.title = quoted(zone.title);
    require(zone.title.size() <= 128, "zone title exceeds 128 bytes");
    require(std::isfinite(zone.solutionTime), "solution time must be finite");
    require(zone.strandId > 0, "strand ID must be positive");
    const auto maxInt = static_cast<std::size_t>(std::numeric_limits<int>::max());
    require(!zone.points.empty() && zone.points.size() <= maxInt, "invalid node count");
    require(zone.connectivity.size() <= maxInt, "invalid element count");
    require(zone.fields.size() <= maxInt - 2, "too many variables");
    std::size_t nodesPerElement = 0;
    switch (zone.element)
    {
    case Element::Line: nodesPerElement = 2; header.zoneType = "FELINESEG"; break;
    case Element::Triangle: nodesPerElement = 3; header.zoneType = "FETRIANGLE"; break;
    case Element::Quadrilateral: nodesPerElement = 4; header.zoneType = "FEQUADRILATERAL"; break;
    case Element::OrderedPoints: header.zoneType = "ORDERED"; break;
    default: require(false, "unsupported element type");
    }
    if (zone.element == Element::OrderedPoints)
        require(zone.connectivity.empty(), "ordered points cannot carry FE connectivity");
    else
        require(!zone.connectivity.empty(), "invalid element count");
    for (const auto& point : zone.points)
        require(std::isfinite(point.x) && std::isfinite(point.y), "nonfinite coordinate");
    for (const auto& element : zone.connectivity)
    {
        require(element.size() == nodesPerElement, "incorrect connectivity size");
        for (std::size_t k = 0; k < element.size(); ++k)
        {
            require(element[k] >= 0 && static_cast<std::size_t>(element[k]) < zone.points.size(),
                "connectivity node is out of range");
            for (std::size_t j = 0; j < k; ++j)
                require(element[j] != element[k], "element repeats a node");
        }
    }
    std::set<std::string> names;
    const auto addName = [&](const std::string& name) {
        require(!name.empty(), "variable name cannot be empty");
        require(names.insert(name).second, "variable names must be unique, including coordinates");
        header.variables.push_back(quoted(name));
    };
    for (const auto& name : zone.coordinateNames) addName(name);
    for (const auto& field : zone.fields)
    {
        addName(field.name);
        std::size_t count = 0;
        switch (field.location)
        {
        case Location::Nodal: count = zone.points.size(); break;
        case Location::CellCentered:
            require(zone.element != Element::OrderedPoints, "ordered points require nodal fields");
            count = zone.connectivity.size();
            header.cellVariables.push_back(header.variables.size());
            break;
        default: require(false, "unsupported field location");
        }
        require(field.values.size() == count, "field value count does not match its location");
        for (double value : field.values) require(std::isfinite(value), "nonfinite field value");
    }
    for (const auto& entry : zone.auxiliary)
    {
        require(identifier(entry.first), "invalid AUXDATA identifier");
        const std::string line = "AUXDATA " + entry.first + "=" + quoted(entry.second);
        require(line.size() <= 32000, "AUXDATA exceeds the ASCII line limit");
        header.auxiliary.push_back(line);
    }
    return header;
}

// Bulk numeric serialization avoids one locale/ostream sentry per scalar.
// general/17 preserves the public DOUBLE precision and existing text layout.
// Every Zone is fully validated before this path can open/truncate its file.
class NumericBuffer {
public:
    explicit NumericBuffer(std::ostream& out):out_(out) {buffer_.reserve(65536);}
    void number(double value) {
        char text[64];
        const auto result=std::to_chars(text,text+sizeof(text),value,std::chars_format::general,17);
        require(result.ec==std::errc{},"DOUBLE text conversion failed");
        buffer_.append(text,static_cast<std::size_t>(result.ptr-text));
    }
    void number(int value) {
        char text[32];
        const auto result=std::to_chars(text,text+sizeof(text),value);
        require(result.ec==std::errc{},"connectivity text conversion failed");
        buffer_.append(text,static_cast<std::size_t>(result.ptr-text));
    }
    void separator(char value) {
        buffer_+=value;
        if(buffer_.size()>=65000)flush();
    }
    void flush() {
        out_.write(buffer_.data(),static_cast<std::streamsize>(buffer_.size()));buffer_.clear();
    }
private:
    std::ostream& out_;
    std::string buffer_;
};

template<class Value>
void block(std::ostream& out, std::size_t count, Value value)
{
    NumericBuffer buffer(out);
    for (std::size_t k = 0; k < count; ++k) {
        buffer.number(value(k));
        buffer.separator((k + 1) % 6 == 0 || k + 1 == count ? '\n' : ' ');
    }
    buffer.flush();
}
}

History::History(const std::filesystem::path& path, const std::string& title,
    const std::vector<std::string>& columns,
    const std::map<std::string, std::string>& auxiliary)
    : path_(path), columns_(columns.size())
{
    require(!path.empty() && !std::filesystem::exists(path), "history requires a new output path");
    require(!columns.empty() && title.size() <= 128, "invalid history columns or title");
    std::ostringstream header;
    header.imbue(std::locale::classic());
    header << "VARIABLES=\n";
    std::set<std::string> unique;
    for (const auto& column : columns)
    {
        require(!column.empty() && unique.insert(column).second, "empty or duplicate history column");
        header << quoted(column) << '\n';
    }
    header << "ZONE T=" << quoted(title) << ", DATAPACKING=POINT\nDT=(\n";
    for (std::size_t k = 0; k < columns_; ++k) header << "DOUBLE\n";
    header << ")\n";
    for (const auto& entry : auxiliary)
    {
        require(identifier(entry.first), "invalid history AUXDATA identifier");
        const std::string line = "AUXDATA " + entry.first + "=" + quoted(entry.second);
        require(line.size() <= 32000, "history AUXDATA exceeds the ASCII line limit");
        header << line << '\n';
    }
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(path_, std::ios::binary);
    out << header.str();
    out.close();
    bytes_ = std::filesystem::file_size(path_);
}

void History::append(const std::vector<long double>& values)
{
    require(!failed_, "cannot append after history I/O failure");
    require(values.size() == columns_, "history row has incorrect column count");
    std::ostringstream record;
    record.imbue(std::locale::classic());
    record << std::setprecision(17) << std::defaultfloat;
    for (std::size_t k = 0; k < values.size(); ++k)
    {
        require(std::isfinite(values[k]), "history values must be finite; use an explicit validity flag");
        record << (k ? " " : "") << values[k];
    }
    record << '\n';
    try
    {
        require(std::filesystem::file_size(path_) == bytes_, "history file was modified outside its owner");
        std::ofstream out;
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out.open(path_, std::ios::binary | std::ios::app);
        out << record.str();
        out.close();
        bytes_ = std::filesystem::file_size(path_);
    }
    catch (const std::exception& error)
    {
        failed_ = true;
        throw std::runtime_error("TecplotWriter: history append failed for " + path_.string() + ": " + error.what());
    }
}

namespace
{
void writeZones(const std::filesystem::path& path, const std::vector<const Zone*>& zones)
{
    require(!zones.empty(), "output requires at least one zone");
    std::vector<Header> headers;
    headers.reserve(zones.size());
    for (const Zone* zone : zones)
    {
        headers.push_back(validate(*zone));
        require(headers.back().variables == headers.front().variables,
            "zone coordinate names and field names/order must match");
    }
    require(!path.empty(), "output path cannot be empty");
    std::ofstream out;
    out.imbue(std::locale::classic());
    out.exceptions(std::ios::failbit | std::ios::badbit);
    try
    {
        out.open(path, std::ios::out | std::ios::trunc);
        out << std::setprecision(17) << std::defaultfloat;
        out << "VARIABLES=\n";
        for (const auto& name : headers.front().variables) out << name << '\n';
        for (std::size_t index = 0; index < zones.size(); ++index)
        {
            const Zone& zone = *zones[index];
            const Header& header = headers[index];
            out << "ZONE T=" << header.title << '\n';
            if (zone.element == Element::OrderedPoints)
                out << "I=" << zone.points.size() << ", J=1, K=1";
            else
                out << "NODES=" << zone.points.size() << ", ELEMENTS=" << zone.connectivity.size();
            out << ", ZONETYPE=" << header.zoneType << ", DATAPACKING=BLOCK\n"
                << "SOLUTIONTIME=" << zone.solutionTime << ", STRANDID=" << zone.strandId << '\n';
            out << "DT=(\n";
            for (std::size_t k = 0; k < header.variables.size(); ++k) out << "DOUBLE\n";
            out << ")\n";
            if (!header.cellVariables.empty())
            {
                out << "VARLOCATION=([\n";
                for (std::size_t k = 0; k < header.cellVariables.size(); ++k)
                    out << (k ? "," : "") << header.cellVariables[k] << '\n';
                out << "]=CELLCENTERED)\n";
            }
            for (const auto& entry : header.auxiliary) out << entry << '\n';
            block(out, zone.points.size(), [&](std::size_t k) { return zone.points[k].x; });
            block(out, zone.points.size(), [&](std::size_t k) { return zone.points[k].y; });
            for (const auto& field : zone.fields)
                block(out, field.values.size(), [&](std::size_t k) { return field.values[k]; });
            NumericBuffer connectivity(out);
            for (const auto& element : zone.connectivity)
            {
                for (std::size_t k = 0; k < element.size(); ++k) {
                    connectivity.number(element[k]+1);
                    connectivity.separator(k + 1 == element.size() ? '\n' : ' ');
                }
            }
            connectivity.flush();
        }
        out.close();
    }
    catch (const std::ios_base::failure& error)
    {
        throw std::runtime_error("TecplotWriter: I/O failure for " + path.string() + ": " + error.what());
    }
}
}

void write(const std::filesystem::path& path, const Zone& zone)
{
    writeZones(path, {&zone});
}

void write(const std::filesystem::path& path, const std::vector<Zone>& zones)
{
    std::vector<const Zone*> views;
    views.reserve(zones.size());
    for (const auto& zone : zones) views.push_back(&zone);
    writeZones(path, views);
}

namespace
{
using CsvRows = std::vector<std::vector<std::string>>;
std::string trimmed(const std::string& value)
{
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return {};
    return value.substr(first,value.find_last_not_of(" \t\r\n")-first+1);
}
std::string lower(std::string value)
{
    for(auto& c:value)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return value;
}
CsvRows parseCsvRecords(std::string text)
{
    if(text.compare(0,3,"\xef\xbb\xbf")==0)text.erase(0,3);
    CsvRows rows;std::vector<std::string> row;std::string field;
    enum State { Start, Plain, Quoted, Closed } state=Start;
    const auto finishField=[&]{row.push_back(std::move(field));field.clear();state=Start;};
    const auto finishRow=[&]{finishField();rows.push_back(std::move(row));row.clear();};
    for(std::size_t i=0;i<text.size();++i){const char c=text[i];
        require(c!='\0',"CSV contains a NUL byte");
        if(state==Quoted){
            if(c=='"'){
                if(i+1<text.size()&&text[i+1]=='"'){field+='"';++i;}else state=Closed;
            }else field+=c;
        }else if(c==',')finishField();
        else if(c=='\r'||c=='\n'){
            finishRow();if(c=='\r'&&i+1<text.size()&&text[i+1]=='\n')++i;
        }else if(c=='"'){
            require(state==Start,"CSV quote inside an unquoted field or after closing quote");state=Quoted;
        }else{
            require(state!=Closed,"CSV characters after closing quote");field+=c;state=Plain;
        }
    }
    require(state!=Quoted,"CSV has an unterminated quoted field");
    if(state!=Start||!field.empty()||!row.empty())finishRow();
    return rows;
}
CsvRows readCsv(const std::filesystem::path& path)
{
    std::ifstream input(path,std::ios::binary);
    if(!input)throw std::runtime_error("TecplotWriter: cannot read CSV "+path.string());
    std::string text((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    require(!input.bad(),"CSV read failed");
    auto rows=parseCsvRecords(std::move(text));
    if(rows.empty())return rows;
    std::set<std::string> names;
    for(auto& name:rows.front()){
        name=trimmed(name);require(!name.empty()&&names.insert(name).second,"CSV has empty or duplicate column names");
    }
    for(std::size_t i=1;i<rows.size();++i)
        require(rows[i].size()==rows.front().size(),"CSV row has incorrect column count");
    return rows;
}
bool decimal(const std::string& value)
{
    std::size_t i=0,digits=0;
    if(i<value.size()&&(value[i]=='+'||value[i]=='-'))++i;
    while(i<value.size()&&value[i]>='0'&&value[i]<='9'){++i;++digits;}
    if(i<value.size()&&value[i]=='.'){
        ++i;while(i<value.size()&&value[i]>='0'&&value[i]<='9'){++i;++digits;}
    }
    if(!digits)return false;
    if(i<value.size()&&(value[i]=='e'||value[i]=='E')){
        ++i;if(i<value.size()&&(value[i]=='+'||value[i]=='-'))++i;
        const auto first=i;while(i<value.size()&&value[i]>='0'&&value[i]<='9')++i;
        if(first==i)return false;
    }
    return i==value.size();
}
struct CsvNumber {
    enum Kind { Missing, Numeric, Text } kind=Missing;
    std::string token;long double value=0;
};
CsvNumber csvNumber(const std::string& raw)
{
    const auto token=trimmed(raw);auto key=lower(token);
    if(key.empty()||key=="na")return {};
    if(key=="true")return {CsvNumber::Numeric,"1",1};
    if(key=="false")return {CsvNumber::Numeric,"0",0};
    if(key.front()=='+'||key.front()=='-')key.erase(0,1);
    if(key=="inf"||key=="infinity"||key=="nan"||(key.size()>4&&key.compare(0,4,"nan(")==0&&key.back()==')'))return {};
    if(!decimal(token))return {CsvNumber::Text,token,0};
    std::istringstream parser(token);parser.imbue(std::locale::classic());long double value=0;
    parser>>value;
    // Finite validation only: emit the original token, never the parsed value.
    if(!parser||!std::isfinite(value))return {};
    return {CsvNumber::Numeric,token,value};
}
std::string visible(const std::string& value)
{
    std::string result;const char* hex="0123456789abcdef";
    for(unsigned char c:value){
        if(c=='\n')result+="\\n";else if(c=='\r')result+="\\r";else if(c=='\t')result+="\\t";
        else if(c<32||c==127){result+="\\x";result+=hex[c>>4];result+=hex[c&15];}
        else result+=static_cast<char>(c);
    }
    return result;
}
bool categoryName(const std::string& name)
{
    static const std::set<std::string> keys{"stage","status","module","event","reason","rejected_reason","kind",
        "phase","species","side","boundary","boundary_name","boundary_id","species_id","phase_id","target","type",
        "pol_solver_mode","top1_name","top2_name","top3_name","dominant_source","dominant_boundary",
        "dominant_reaction_channel","curvature_measure","force_discretization","area_update"};
    return keys.count(lower(name))!=0;
}
std::size_t clockColumn(const std::vector<std::string>& names,std::initializer_list<const char*> choices)
{
    for(const auto* choice:choices)for(std::size_t i=0;i<names.size();++i)if(lower(names[i])==choice)return i;
    return names.size();
}
struct CsvCurve {
    std::string column;std::vector<std::string> category;
    std::vector<std::array<std::string,3>> samples;
};
}

bool exportCsvHistory(const std::filesystem::path& csvPath,const std::filesystem::path& datPath)
{
    require(!datPath.empty(),"CSV history output path cannot be empty");
    std::error_code error;
    require(!std::filesystem::equivalent(csvPath,datPath,error),"CSV history cannot overwrite its source");
    const auto table=readCsv(csvPath);if(table.size()<2)return false;
    const auto& names=table.front();const auto none=names.size();
    const auto time=clockColumn(names,{"time","accepted_time","time_s","time_code","t_si_s","time_si_s","physical_time","t"});
    const auto step=clockColumn(names,{"step","accepted_step","timestep","attempted_step"});
    const auto axis=time!=none?time:step;
    std::vector<std::size_t> categories,metrics;
    for(std::size_t column=0;column<names.size();++column){
        bool numeric=false,text=false;
        for(std::size_t row=1;row<table.size();++row){const auto value=csvNumber(table[row][column]);
            numeric|=value.kind==CsvNumber::Numeric;text|=value.kind==CsvNumber::Text;}
        if(column==axis||column==step){require(!text,"CSV clock column contains text");continue;}
        if(categoryName(names[column])||text){
            require(categoryName(names[column])||!numeric,"CSV column mixes numeric and categorical values");categories.push_back(column);
        }else if(numeric)metrics.push_back(column);
    }
    std::vector<CsvCurve> curves;
    for(const auto metric:metrics){
        struct State { std::size_t curve=0;long double x=0;bool active=false; };
        std::map<std::vector<std::string>,State> groups;
        for(std::size_t row=1;row<table.size();++row){
            std::vector<std::string> group;for(const auto column:categories)group.push_back(table[row][column]);
            auto& state=groups[group];const auto value=csvNumber(table[row][metric]);
            const auto record=std::to_string(row);
            const auto x=axis==none?CsvNumber{CsvNumber::Numeric,record,static_cast<long double>(row)}:csvNumber(table[row][axis]);
            const auto s=step==none?CsvNumber{CsvNumber::Numeric,record,static_cast<long double>(row)}:csvNumber(table[row][step]);
            if(value.kind!=CsvNumber::Numeric||x.kind!=CsvNumber::Numeric||s.kind!=CsvNumber::Numeric){state.active=false;continue;}
            if(!state.active||x.value<=state.x){curves.push_back({names[metric],group,{}});state.curve=curves.size()-1;}
            curves[state.curve].samples.push_back({x.token,value.token,s.token});state.x=x.value;state.active=true;
        }
    }
    if(curves.empty())return false;
    const std::string axisName=axis==none?"Record":names[axis];
    std::string stepName=step==none?"Record":"Step";
    if(lower(stepName)==lower(axisName))stepName="Source"+stepName;
    const std::string meaning=axis==none?"source logical data record (1-based); no physical clock column":
        (time==none?"actual source step column; no physical time column":"actual source time column; no resampling");
    const auto source=std::filesystem::weakly_canonical(csvPath).generic_string();
    const std::string generator="LAMP_CsvHistory_v1";
    const auto generatorLine="DATASETAUXDATA Generator="+quoted(generator);
    const auto sourceLine="DATASETAUXDATA SourceCsv="+quoted(visible(source));
    if(std::filesystem::exists(datPath)){
        require(!std::filesystem::is_symlink(datPath),"CSV history cannot replace a symlink");
        std::ifstream old(datPath,std::ios::binary);require(static_cast<bool>(old),"cannot inspect existing CSV companion ownership");
        bool owned=false,sameSource=false;std::string line;
        while(std::getline(old,line)){
            if(line.compare(0,4,"ZONE")==0)break;
            line=trimmed(line);owned|=line==generatorLine;sameSource|=line==sourceLine;
        }
        require(!old.bad()&&owned&&sameSource,"existing DAT is not a matching generated CSV companion");
    }
    std::ostringstream output;output.imbue(std::locale::classic());
    output<<"VARIABLES="<<quoted(axisName)<<' '<<quoted("Value")<<' '<<quoted(stepName)<<'\n'
        <<generatorLine<<'\n'<<sourceLine<<'\n';
    for(const auto& curve:curves){
        auto title=visible(curve.column);for(std::size_t i=0;i<categories.size();++i)title+=" | "+visible(names[categories[i]])+'='+visible(curve.category[i]);
        if(title.size()>128){std::size_t end=125;while(end>0&&(static_cast<unsigned char>(title[end])&0xc0)==0x80)--end;title.resize(end);title+="...";}
        output<<"ZONE T="<<quoted(title)<<", I="<<curve.samples.size()<<", ZONETYPE=ORDERED, DATAPACKING=POINT\nDT=(DOUBLE DOUBLE DOUBLE)\n";
        const auto auxiliary=[&](const std::string& key,const std::string& value){output<<"AUXDATA "<<key<<'='<<quoted(visible(value))<<'\n';};
        auxiliary("SourceCsv",source);auxiliary("Column",curve.column);
        auxiliary("AxisColumn",axis==none?"Record":names[axis]);auxiliary("XMeaning",meaning);
        auxiliary("StepColumn",step==none?"Record":names[step]);
        if(step==none)auxiliary("RecordMeaning","source logical data record (1-based); not a physical step");
        auxiliary("MissingPolicy","NA, empty and nonfinite omitted; gaps and nonincreasing X start new zones");
        for(std::size_t i=0;i<categories.size();++i){auxiliary("Category"+std::to_string(i)+"Name",names[categories[i]]);auxiliary("Category"+std::to_string(i)+"Value",curve.category[i]);}
        for(const auto& sample:curve.samples)output<<sample[0]<<' '<<sample[1]<<' '<<sample[2]<<'\n';
    }
    // Finish every parse/validation/serialization before touching the target.
    static std::atomic<unsigned long long> sequence{0};
    const auto suffix=std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count())+'.'+std::to_string(sequence++);
    const std::filesystem::path temporary=datPath.string()+".tmp."+suffix;
    try{
        std::ofstream file;file.exceptions(std::ios::failbit|std::ios::badbit);file.open(temporary,std::ios::binary|std::ios::trunc);
        file<<output.str();file.close();std::filesystem::rename(temporary,datPath);
    }catch(const std::exception& failure){
        std::error_code ignored;std::filesystem::remove(temporary,ignored);
        throw std::runtime_error("TecplotWriter: CSV history publication "+datPath.string()+": "+failure.what());
    }
    return true;
}
std::size_t exportCsvHistories(const std::filesystem::path& directory)
{
    require(std::filesystem::is_directory(directory),"CSV histories require an existing series directory");
    std::vector<std::filesystem::path> paths;
    for(const auto& entry:std::filesystem::recursive_directory_iterator(directory))
        if(!entry.is_symlink()&&entry.is_regular_file()&&lower(entry.path().extension().string())==".csv")paths.push_back(entry.path());
    std::sort(paths.begin(),paths.end());std::size_t exported=0;
    for(const auto& csv:paths){auto dat=csv;dat.replace_extension(".dat");exported+=exportCsvHistory(csv,dat)?1:0;}
    return exported;
}

namespace
{
void requireNewTablePath(const std::filesystem::path& path)
{
    require(!path.empty(),"table history output path cannot be empty");
    std::error_code error;
    const auto status=std::filesystem::symlink_status(path,error);
    if(error&&error!=std::errc::no_such_file_or_directory)
        throw std::runtime_error("TecplotWriter: cannot inspect table history path "+path.string()+": "+error.message());
    require(!std::filesystem::exists(status)&&!std::filesystem::is_symlink(status),
        "table history requires a new output path");
}
std::string byteHex(const std::string& value)
{
    static const char digits[]="0123456789abcdef";
    std::string result;result.reserve(value.size()*2);
    for(unsigned char c:value){result+=digits[c>>4];result+=digits[c&15];}
    return result;
}
std::string tableAuxiliary(const std::string& scope,const std::string& name,const std::string& value)
{
    require(identifier(name),"invalid table history AUXDATA identifier");
    const auto line=scope+' '+name+'='+quoted(value);
    require(line.size()<=32000,"table history AUXDATA exceeds the ASCII line limit");
    return line+'\n';
}
// Long diagnostic categories retain all bytes without violating Tecplot's
// single-line string bound. Ordinary short category keys stay byte-identical.
std::string tableCategoryAuxiliary(const std::string& name,const std::string& value)
{
    constexpr std::size_t chunkBytes=12000;
    if(value.size()<=chunkBytes)return tableAuxiliary("AUXDATA",name,value);
    std::vector<std::string> chunks;
    for(std::size_t begin=0;begin<value.size();) {
        auto end=std::min(value.size(),begin+chunkBytes);
        // Keep a UTF-8 continuation sequence in one readable display chunk.
        if(end<value.size())while(end>begin&&(static_cast<unsigned char>(value[end])&0xc0)==0x80)--end;
        require(end>begin,"invalid long category chunk boundary");
        chunks.push_back(value.substr(begin,end-begin));begin=end;
    }
    auto result=tableAuxiliary("AUXDATA",name+"ChunkCount",std::to_string(chunks.size()));
    result+=tableAuxiliary("AUXDATA",name+"ByteCount",std::to_string(value.size()));
    for(std::size_t i=0;i<chunks.size();++i)
        result+=tableAuxiliary("AUXDATA",name+"Chunk"+std::to_string(i),chunks[i]);
    return result;
}
}

struct TableHistory::Data
{
    struct Position
    {
        bool active=false,hasClock=false;
        long double clock=0;
        std::vector<std::string> category;
        std::size_t records=0,zones=0;
    } position;
    std::filesystem::path path;
    std::string title,header;
    std::vector<TableColumn> columns;
    std::vector<std::size_t> numeric,categories;
    std::size_t clock=0;
    std::ofstream stream;
    std::uintmax_t expectedBytes=0,flushedBytes=0;
    bool opened=false,closed=false,failed=false;

    void inspectSize() const
    {
        require(std::filesystem::is_regular_file(path)&&!std::filesystem::is_symlink(path),
            "table history file was removed or replaced outside its owner");
        const auto size=std::filesystem::file_size(path);
        // The persistent stream may already have emitted some buffered bytes.
        require(size>=flushedBytes&&size<=expectedBytes,
            "table history file size changed outside its owner");
    }
};

TableHistory::TableHistory(const std::filesystem::path& path,const std::string& title,
    const std::vector<TableColumn>& columns,const std::map<std::string,std::string>& auxiliary)
    :data_(std::make_unique<Data>())
{
    requireNewTablePath(path);
    require(title.size()<=128&&!columns.empty(),"invalid table history title or columns");
    quoted(title);
    auto& data=*data_;data.path=path;data.title=title;data.columns=columns;
    std::vector<std::string> names;
    std::set<std::string> unique;
    for(std::size_t i=0;i<columns.size();++i){const auto& column=columns[i];
        require(!trimmed(column.name).empty()&&unique.insert(lower(column.name)).second,
            "empty or duplicate table history column");
        quoted(column.name);
        switch(column.kind){
        case TableColumnKind::Numeric:data.numeric.push_back(i);names.push_back(column.name);break;
        case TableColumnKind::Category:data.categories.push_back(i);break;
        default:require(false,"unsupported table history column kind");
        }
    }
    require(!data.numeric.empty(),"table history requires at least one numeric column");
    for(const auto column:data.numeric)
        require(unique.insert(lower(columns[column].name+"Available")).second,
            "table history Available column collides with a source column");
    const auto time=clockColumn(names,{"time","accepted_time","time_s","time_code","t_si_s","time_si_s","physical_time","t"});
    const auto step=clockColumn(names,{"step","accepted_step","timestep","attempted_step"});
    data.clock=time!=names.size()?time:step;
    std::ostringstream header;header.imbue(std::locale::classic());header<<"VARIABLES=\n";
    for(const auto& name:names)header<<quoted(name)<<'\n';
    for(const auto& name:names)header<<quoted(name+"Available")<<'\n';
    std::set<std::string> reserved;
    const auto dataset=[&](const std::string& name,const std::string& value){
        reserved.insert(lower(name));header<<tableAuxiliary("DATASETAUXDATA",name,value);
    };
    dataset("Generator","LAMP_TableHistory_v1");
    dataset("MissingPolicy","NA, empty and nonfinite numeric values use 0 placeholders with <column>Available=0; use matching Available variables for XY value blanking; 0 with Available=1 is real data");
    dataset("CategoryEncoding","CategoryNValue is readable escaped display; CategoryNValueHex is the exact original CSV field byte sequence in lowercase hexadecimal");
    dataset("LongCategoryEncoding","Long CategoryNValue/ValueHex use <key>ChunkCount, <key>ByteCount and consecutive <key>Chunk0... keys; concatenate unescaped chunks in order before hex decoding");
    dataset("ClockColumn",data.clock==names.size()?"none":names[data.clock]);
    dataset("ZonePolicy","category byte changes or nonincreasing finite clock start a zone; missing values remain flagged rows; repeated records are retained");
    for(std::size_t i=0;i<columns.size();++i){
        dataset("Column"+std::to_string(i)+"Name",columns[i].name);
        dataset("Column"+std::to_string(i)+"Kind",columns[i].kind==TableColumnKind::Numeric?"Numeric":"Category");
    }
    for(const auto& entry:auxiliary){
        require(reserved.insert(lower(entry.first)).second,"table history AUXDATA collides with generated metadata");
        header<<tableAuxiliary("DATASETAUXDATA",entry.first,entry.second);
    }
    data.header=header.str();
}

TableHistory::~TableHistory()=default;

void TableHistory::append(const std::string& completeCsvRecords)
{
    auto& data=*data_;
    require(!data.closed,"cannot append after table history close");
    require(!data.failed,"cannot append after table history I/O failure");
    const auto rows=parseCsvRecords(completeCsvRecords);
    if(rows.empty())return;
    auto next=data.position;
    std::ostringstream output;output.imbue(std::locale::classic());
    for(const auto& row:rows){
        require(row.size()==data.columns.size(),"table history row has incorrect column count");
        std::vector<CsvNumber> numeric;numeric.reserve(data.numeric.size());
        for(const auto column:data.numeric){
            const auto value=csvNumber(row[column]);
            require(value.kind!=CsvNumber::Text,"table history numeric column contains text");
            numeric.push_back(value);
        }
        std::vector<std::string> category;category.reserve(data.categories.size());
        for(const auto column:data.categories)category.push_back(row[column]);
        const bool finiteClock=data.clock<numeric.size()&&numeric[data.clock].kind==CsvNumber::Numeric;
        const bool restart=finiteClock&&next.hasClock&&numeric[data.clock].value<=next.clock;
        if(!next.active||category!=next.category||restart){
            ++next.zones;
            auto title=data.title;
            for(std::size_t i=0;i<data.categories.size();++i)
                title+=" | "+visible(data.columns[data.categories[i]].name)+'='+visible(category[i]);
            if(title.size()>128){std::size_t end=125;while(end>0&&(static_cast<unsigned char>(title[end])&0xc0)==0x80)--end;title.resize(end);title+="...";}
            output<<"ZONE T="<<quoted(title)<<", ZONETYPE=ORDERED, DATAPACKING=POINT\nDT=(\n";
            for(std::size_t i=0;i<data.numeric.size()*2;++i)output<<"DOUBLE\n";
            output<<")\n";
            output<<tableAuxiliary("AUXDATA","FirstRecord",std::to_string(next.records+1));
            for(std::size_t i=0;i<data.categories.size();++i){
                const auto prefix="Category"+std::to_string(i);
                output<<tableAuxiliary("AUXDATA",prefix+"Name",data.columns[data.categories[i]].name);
                output<<tableCategoryAuxiliary(prefix+"Value",visible(category[i]));
                output<<tableCategoryAuxiliary(prefix+"ValueHex",byteHex(category[i]));
            }
            next.active=true;next.category=std::move(category);
        }
        for(std::size_t i=0;i<numeric.size();++i)
            output<<(i?" ":"")<<(numeric[i].kind==CsvNumber::Numeric?numeric[i].token:"0");
        for(const auto& value:numeric)output<<' '<<(value.kind==CsvNumber::Numeric?'1':'0');
        output<<'\n';++next.records;
        if(finiteClock){next.hasClock=true;next.clock=numeric[data.clock].value;}
    }
    // No stream, file, or accepted-position mutation until the entire batch has
    // passed parsing, schema checks, finite classification and AUX serialization.
    auto text=output.str();if(!data.opened)text=data.header+text;
    require(text.size()<=static_cast<std::uintmax_t>(std::numeric_limits<std::streamsize>::max()),
        "table history batch exceeds stream size");
    try{
        if(!data.opened){
            requireNewTablePath(data.path);
            data.stream.imbue(std::locale::classic());
            data.stream.exceptions(std::ios::failbit|std::ios::badbit);
            data.stream.open(data.path,std::ios::binary|std::ios::app);
            require(std::filesystem::file_size(data.path)==0,"table history target appeared before first write");
            data.opened=true;
        }
        data.stream.write(text.data(),static_cast<std::streamsize>(text.size()));
        data.expectedBytes+=static_cast<std::uintmax_t>(text.size());
        data.position=std::move(next);
    }catch(const std::exception& failure){
        data.failed=true;
        throw std::runtime_error("TecplotWriter: table history append "+data.path.string()+": "+failure.what());
    }
}

void TableHistory::flush()
{
    auto& data=*data_;
    require(!data.failed,"cannot flush after table history I/O failure");
    if(data.closed||!data.opened)return;
    try{
        data.inspectSize();data.stream.flush();
        require(std::filesystem::file_size(data.path)==data.expectedBytes,
            "table history flushed file size differs from accepted bytes");
        data.flushedBytes=data.expectedBytes;
    }catch(const std::exception& failure){
        data.failed=true;
        throw std::runtime_error("TecplotWriter: table history flush "+data.path.string()+": "+failure.what());
    }
}

void TableHistory::close()
{
    auto& data=*data_;if(data.closed)return;
    flush();
    try{if(data.opened)data.stream.close();data.closed=true;}
    catch(const std::exception& failure){
        data.failed=true;
        throw std::runtime_error("TecplotWriter: table history close "+data.path.string()+": "+failure.what());
    }
}
}
