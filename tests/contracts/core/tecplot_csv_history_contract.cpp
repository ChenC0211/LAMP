#include "Core/IO/TecplotWriter.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace fs=std::filesystem;
static void check(bool value,const char* reason){if(!value)throw std::runtime_error(reason);}
static void put(const fs::path& path,const std::string& text){std::ofstream out(path,std::ios::binary);out<<text;check(static_cast<bool>(out),"fixture write failed");}
static std::string get(const fs::path& path){std::ifstream in(path,std::ios::binary);return std::string(std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>());}
static std::size_t count(const std::string& text,const std::string& token){std::size_t result=0,pos=0;while((pos=text.find(token,pos))!=std::string::npos){++result;pos+=token.size();}return result;}
template<class Action> static void rejects(Action action){bool caught=false;try{action();}catch(const std::exception&){caught=true;}check(caught,"malformed CSV or foreign target accepted");}

static std::string csvField(const std::string& value)
{
    std::string result="\"";
    for(const char c:value){if(c=='"')result+='"';result+=c;}
    return result+'"';
}
static std::string hexCategory(const std::string& text,const std::string& key,std::size_t from=0)
{
    const auto marker="AUXDATA "+key+"=\"";const auto found=text.find(marker,from);
    check(found!=std::string::npos,"category hex metadata absent");
    const auto begin=found+marker.size(),end=text.find('"',begin);
    check(end!=std::string::npos&&(end-begin)%2==0,"category hex metadata malformed");
    const auto nibble=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::runtime_error("noncanonical category hex");};
    std::string result;
    for(auto i=begin;i<end;i+=2)result+=static_cast<char>((nibble(text[i])<<4)|nibble(text[i+1]));
    return result;
}
static void tableHistoryContract(const fs::path& root)
{
    using Kind=TecplotWriter::TableColumnKind;
    const std::vector<TecplotWriter::TableColumn> schema{
        {"step",Kind::Numeric},{"time",Kind::Numeric},{"stage",Kind::Category},
        {"q",Kind::Numeric},{"flag",Kind::Numeric},{"missing",Kind::Numeric},{"reason",Kind::Category}};
    const auto lazy=root/"table_lazy.dat";
    {TecplotWriter::TableHistory history(lazy,"lazy",schema);history.append("");history.flush();history.close();history.close();}
    check(!fs::exists(lazy),"empty table created an output file");
    const auto foreign=root/"table_foreign.dat";put(foreign,"foreign table\n");
    rejects([&]{TecplotWriter::TableHistory history(foreign,"foreign",schema);});
    check(get(foreign)=="foreign table\n","table constructor changed a foreign file");
    const auto collision=root/"table_collision.dat";
    rejects([&]{TecplotWriter::TableHistory history(collision,"collision",{{"q",Kind::Numeric},{"qAvailable",Kind::Category}});});
    rejects([&]{TecplotWriter::TableHistory history(collision,"collision",{{"q",Kind::Numeric},{"Q",Kind::Numeric}});});
    rejects([&]{TecplotWriter::TableHistory history(collision,"collision",schema,{{"Generator","foreign"}});});
    check(!fs::exists(collision),"invalid schema or metadata created output");
    const auto invalid=root/"table_invalid.dat";
    {TecplotWriter::TableHistory history(invalid,"invalid",schema);
        rejects([&]{history.append("0,0,A,1,true,NA,ok\n1,.1,A,text,false,NA,ok\n");});
        rejects([&]{history.append("0,0,A,1,true,NA,ok\n1,.1,A,2,false,NA,\"unclosed\n");});
        rejects([&]{history.append("0,0,A,1,true,NA,ok\n1,.1,A,2,false\n");});
        check(!fs::exists(invalid),"invalid batch published its valid prefix");
        history.flush();history.close();}
    const auto dat=root/"table.dat";
    const std::string precise="1.2345678901234567890123456789012345e-18";
    const std::string reason="line one\nliteral \\n, \"quoted\"\tend\r"+std::string("\xe4\xb8\xad");
    const auto record=[&](const std::string& step,const std::string& time,const std::string& stage,
        const std::string& q,const std::string& flag,const std::string& missing,const std::string& why){
        return step+','+time+','+csvField(stage)+','+q+','+flag+','+missing+','+csvField(why)+'\n';
    };
    TecplotWriter::TableHistory history(dat,"wide table",schema,{{"Units","native diagnostics"}});
    history.append(record("0","0","predictor, A",precise,"true","NA",reason)
        +record("1",".1","predictor, A","0","false","inf",reason));
    history.append(record("2",".2","predictor, A","-2","true","",reason));
    history.flush();const auto first=get(dat);
    const std::string variables="VARIABLES=\n\"step\"\n\"time\"\n\"q\"\n\"flag\"\n\"missing\"\n"
        "\"stepAvailable\"\n\"timeAvailable\"\n\"qAvailable\"\n\"flagAvailable\"\n\"missingAvailable\"\n";
    check(first.compare(0,variables.size(),variables)==0,"table numeric names/order or availability schema changed");
    check(count(first,"ZONE T=")==1&&first.find(", I=")==std::string::npos,"three increasing table points did not share one live zone");
    check(first.find("wide table | stage=predictor, A")!=std::string::npos,"table zone list omits readable category labels");
    check(first.find("DATAPACKING=POINT")!=std::string::npos&&count(first,"DOUBLE\n")==10,"wide table is not all-DOUBLE POINT data");
    check(first.find("0 0 "+precise+" 1 0 1 1 1 1 0\n")!=std::string::npos
        &&first.find("1 .1 0 0 0 1 1 1 1 0\n")!=std::string::npos
        &&first.find("2 .2 -2 1 0 1 1 1 1 0\n")!=std::string::npos,
        "table narrowed tokens or confused real zero/false with missing values");
    check(first.find("DATASETAUXDATA Generator=\"LAMP_TableHistory_v1\"")!=std::string::npos
        &&first.find("Available=0")!=std::string::npos&&first.find("XY value blanking")!=std::string::npos,
        "table placeholder/blanking metadata absent");
    check(first.find("Category0Name=\"stage\"")!=std::string::npos
        &&first.find("Category1Name=\"reason\"")!=std::string::npos
        &&hexCategory(first,"Category0ValueHex")=="predictor, A"
        &&hexCategory(first,"Category1ValueHex")==reason,"table category bytes were not reversible");
    check(first.find("Category1Value=\"line one\\\\nliteral \\\\n, \\\"quoted\\\"\\\\tend\\\\r")!=std::string::npos,
        "table category display is not readable escaped text");
    rejects([&]{history.append(record("3",".3","predictor, A","3","true","0",reason)
        +record("4",".4","predictor, A","bad","true","0",reason));});
    check(get(dat)==first,"invalid table batch changed accepted prefix");
    history.append(record("3",".3","predictor, A","3","true","0",reason));
    history.flush();const auto afterRecovery=get(dat);
    check(count(afterRecovery,"ZONE T=")==1&&afterRecovery.find("3 .3 3 1 0 1 1 1 1 1\n")!=std::string::npos,
        "table validation failure advanced or poisoned accepted state");
    auto otherReason=reason;otherReason.replace(otherReason.find('\n'),1,"\\n");
    history.append(record("4",".4","corrector","4","false","nan",otherReason));
    history.append(record("5",".4","corrector","5","false","-inf",otherReason));
    history.append(record("6",".2","corrector","6","false","+inf",otherReason));
    history.append(record("7",".5","corrector","7","false","1",otherReason));
    history.append(record("8",".6","corrector","8","false","1",reason));
    history.close();history.close();const auto complete=get(dat);
    check(count(complete,"ZONE T=")==5,"category/reason/repeated/backward clock segmentation changed");
    const auto secondZone=complete.find("ZONE T=",complete.find("ZONE T=")+1);
    check(hexCategory(complete,"Category1ValueHex",secondZone)==otherReason&&otherReason!=reason,
        "actual newline and literal backslash-n category bytes collapsed");
    check(complete.find("AUXDATA FirstRecord=\"5\"")!=std::string::npos
        &&complete.find("AUXDATA FirstRecord=\"6\"")!=std::string::npos
        &&complete.find("AUXDATA FirstRecord=\"7\"")!=std::string::npos
        &&complete.find("AUXDATA FirstRecord=\"9\"")!=std::string::npos,
        "segmentation dropped repeated records or changed logical record numbers");
    rejects([&]{history.append(record("9",".7","corrector","9","true","0",otherReason));});
    check(get(dat)==complete&&!fs::exists(root/"table.csv"),"closed table changed or created a CSV");
    // A real failure diagnostic can exceed 32 KiB, and hexadecimal bytes
    // require twice that space. Reconstruct every original byte independently.
    const auto largePath=root/"table_long_failure.dat";
    std::string largeReason;
    for(int i=0;i<4000;++i)largeReason+="detail \"quoted\" \\ path\n"+std::string("\xe4\xb8\xad");
    {TecplotWriter::TableHistory large(largePath,"physical failure",schema);
        large.append(record("1202",".17471699090331874","FAILED",precise,"false","NA",largeReason));large.close();}
    const auto largeText=get(largePath);std::istringstream largeLines(largeText);
    for(std::string line;std::getline(largeLines,line);)check(line.size()<=32000,"long category exceeded Tecplot ASCII line limit");
    const auto attribute=[&](const std::string& key){
        const std::string marker="AUXDATA "+key+"=\"";const auto pos=largeText.find(marker);
        check(pos!=std::string::npos,"long category metadata absent");const auto first=pos+marker.size(),last=largeText.find('"',first);
        check(last!=std::string::npos,"long category metadata unterminated");return largeText.substr(first,last-first);
    };
    const auto chunks=std::stoull(attribute("Category1ValueHexChunkCount"));std::string hex;
    for(std::size_t i=0;i<chunks;++i)hex+=attribute("Category1ValueHexChunk"+std::to_string(i));
    check(hex.size()==std::stoull(attribute("Category1ValueHexByteCount"))&&hex.size()==2*largeReason.size(),"long category chunk bytes were dropped");
    const auto nibble=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::runtime_error("bad long category hex");};
    std::string recovered;for(std::size_t i=0;i<hex.size();i+=2)recovered+=static_cast<char>(16*nibble(hex[i])+nibble(hex[i+1]));
    check(recovered==largeReason&&largeText.find("1202 .17471699090331874 "+precise+" 0 0 1 1 1 1 0\n")!=std::string::npos,
        "long failure category lost exact bytes or numeric record");
    check(largeText.find("LongCategoryEncoding")!=std::string::npos&&count(largeText,"ZONE T=")==1,"long failure metadata changed actual zone count");
    const auto raced=root/"table_raced.dat";
    {TecplotWriter::TableHistory delayed(raced,"late foreign",schema);put(raced,"late foreign owner\n");
        rejects([&]{delayed.append(record("0","0","A","1","true","NA","reason"));});}
    check(get(raced)=="late foreign owner\n","lazy table overwrote an existing first-write target");
    const auto changed=root/"table_changed.dat";
    {TecplotWriter::TableHistory owned(changed,"owned",schema);owned.append(record("0","0","A","1","true","NA","reason"));owned.flush();
        {std::ofstream external(changed,std::ios::binary|std::ios::app);external<<"foreign append\n";}
        const auto outside=get(changed);rejects([&]{owned.flush();});rejects([&]{owned.append(record("1",".1","A","2","true","NA","reason"));});
        check(get(changed)==outside,"table I/O failure changed a foreign-modified prefix");}
}
int main()
{
    try{
        const auto root=fs::current_path()/("tecplot_csv_history_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        fs::create_directories(root);const auto csv=root/"mixed.csv",dat=root/"mixed.dat";
        const std::string precise="1.2345678901234567890123456789012345e-18";
        put(csv,"\xef\xbb\xbfstep,time,stage,q,flag,unavailable\r\n"
            "0,0,\"predictor, \"\"A\"\"\","+precise+",true,NA\r\n"
            "1,.1,\"predictor, \"\"A\"\"\",2.25,false,\r\n"
            "2,.2,\"predictor, \"\"A\"\"\",NA,true,inf\r\n"
            "3,.3,\"predictor, \"\"A\"\"\",4.5,true,nan\r\n"
            "4,.4,corrector,6,false,\r\n"
            "5,.5,\"predictor, \"\"A\"\"\",7.5,false,\r\n");
        const auto original=get(csv);check(TecplotWriter::exportCsvHistory(csv,dat),"mixed history skipped");const auto text=get(dat);
        check(text.find("VARIABLES=\"time\" \"Value\" \"Step\"")!=std::string::npos,"default XY variables differ");
        check(text.find(precise)!=std::string::npos,"numeric token precision narrowed");
        check(count(text,"AUXDATA Column=\"q\"")==3,"category grouping or NA segment split differs");
        check(text.find("Category0Value=\"predictor, \\\"A\\\"\"")!=std::string::npos,"quoted CSV category lost comma/quote");
        check(text.find("0 "+precise+" 0\n.1 2.25 1\n")!=std::string::npos,"first finite run changed");
        check(text.find(".3 4.5 3\n.5 7.5 5\n")!=std::string::npos,"interleaved category history not retained");
        check(text.find(".1 0 1\n")!=std::string::npos,"boolean false was not explicit numeric flag");
        check(text.find("AUXDATA Column=\"unavailable\"")==std::string::npos,"all-missing field forged a curve");
        check(get(csv)==original,"CSV source changed");
        check(TecplotWriter::exportCsvHistory(csv,dat)&&get(dat)==text,"owned companion cannot be replaced deterministically");
        put(root/"empty.csv","step,time,q\n");check(!TecplotWriter::exportCsvHistory(root/"empty.csv",root/"empty.dat")&&!fs::exists(root/"empty.dat"),"header-only DAT generated");
        put(root/"missing.csv","step,time,q\n0,0,NA\n1,.1,\n2,.2,-inf\n3,.3,NaN\n");
        check(!TecplotWriter::exportCsvHistory(root/"missing.csv",root/"missing.dat")&&!fs::exists(root/"missing.dat"),"all-missing DAT generated");
        put(root/"gaps.csv","step,time,q\n0,0,1\n1,.1,nan\n2,.2,2\n3,.3,+inf\n4,.4,3\n4,.4,4\n");
        TecplotWriter::exportCsvHistory(root/"gaps.csv",root/"gaps.dat");check(count(get(root/"gaps.dat"),"ZONE T=")==4,"nonfinite/repeated-X segments connected");
        put(root/"step.csv","timestep,q\n7,3\n8,4\n");TecplotWriter::exportCsvHistory(root/"step.csv",root/"step.dat");
        const auto stepped=get(root/"step.dat");
        check(stepped.find("VARIABLES=\"timestep\" \"Value\" \"Step\"")!=std::string::npos
            &&stepped.find("AUXDATA AxisColumn=\"timestep\"")!=std::string::npos
            &&stepped.find("actual source step column; no physical time column")!=std::string::npos
            &&stepped.find("7 3 7\n8 4 8\n")!=std::string::npos,"step-only abscissa mislabeled or source values changed");
        put(root/"only_step.csv","step,q\n7,3\n8,4\n");TecplotWriter::exportCsvHistory(root/"only_step.csv",root/"only_step.dat");
        check(get(root/"only_step.dat").find("VARIABLES=\"step\" \"Value\" \"SourceStep\"")!=std::string::npos,
            "step-only variables collide when source abscissa is named step");
        put(root/"time.csv","time_s,q\n.1,3\n.2,4\n");TecplotWriter::exportCsvHistory(root/"time.csv",root/"time.dat");
        const auto time=get(root/"time.dat");check(time.find("\"Record\"")!=std::string::npos&&time.find("not a physical step")!=std::string::npos&&time.find(".1 3 1\n.2 4 2\n")!=std::string::npos,"record index forged a step");
        put(root/"multiline.csv","step,time,reason,q\n0,0,\"line one\r\nline two\",4\n");
        check(TecplotWriter::exportCsvHistory(root/"multiline.csv",root/"multiline.dat"),"quoted multiline CSV rejected");
        const auto bad=root/"bad.csv",guard=root/"guard.dat";put(guard,"foreign data\n");
        for(const auto& content:{"step,time,q\n0,0\n","step,time,q\n0,0,\"open\n","step,time,q\n0,0,un\"quoted\n","step,time,q\n0,0,\"1\"x\n","step,time,q\n0,0,1\n1,.1,bad\n"}){
            put(bad,content);rejects([&]{TecplotWriter::exportCsvHistory(bad,guard);});check(get(guard)=="foreign data\n","invalid input touched target");
        }
        rejects([&]{TecplotWriter::exportCsvHistory(csv,guard);});check(get(guard)=="foreign data\n","foreign target truncated");
        put(root/"other.csv",original);rejects([&]{TecplotWriter::exportCsvHistory(root/"other.csv",dat);});check(get(dat)==text,"different SourceCsv target overwritten");
        const auto nested=root/"nested";fs::create_directories(nested/"diagnostics");put(nested/"one.csv","step,q\n1,2\n");put(nested/"diagnostics/two.csv","step,q\n1,3\n");put(nested/"skip.csv","step,q\n1,NA\n");
        check(TecplotWriter::exportCsvHistories(nested)==2&&fs::exists(nested/"diagnostics/two.dat"),"recursive series export differs");
        tableHistoryContract(root);
        for(const auto& entry:fs::recursive_directory_iterator(root))check(entry.path().filename().string().find(".tmp.")==std::string::npos,"temporary publication leaked");
        std::cout<<"Tecplot CSV finite segmented history contract PASS\n";return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
