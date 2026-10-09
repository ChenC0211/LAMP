#include "Core/IO/Dictionary.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/DictionaryBlockText.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

struct Reader : Dictionary {
    using Dictionary::parseIndex;
    using Dictionary::parseIndex_brace;
};
template<class Error, class F>
void readerRejects(F action, const std::filesystem::path& source, const char* reason)
{
    bool caught=false;
    try { action(); }
    catch (const Error& error) {
        const std::string message=error.what();
        caught=message.find(source.string())!=std::string::npos && message.find(reason)!=std::string::npos;
    }
    if (!caught) throw std::runtime_error("dictionary failure lost its exception type, source or reason");
}
void readerFailures()
{
    namespace fs=std::filesystem;
    const auto cwd=fs::current_path();
    const auto folder=cwd/("dictionary_failures_"+
        std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    fs::create_directories(folder);
    const auto geometry=folder/"shape.dict", mesh=folder/"grid.dict";
    const auto absent=folder/"absent.dict";
    const auto write=[](const fs::path& path,const std::string& text) {
        std::ofstream stream(path); stream<<text;
        if (!stream) throw std::runtime_error("cannot write dictionary fixture");
    };
    const std::string validGeometry="Domain\n{\n min (0.2, -0.5)\n max (1.2, 0.5)\n}\n";
    const std::string validMesh="adaptiveMesh false\nBackgroundMesh\n{\n divide (8, 10)\n}\n";
    const auto readMesh=[&] { DictionaryMesh value; value.readDict(mesh.string(),geometry.string()); };
    readerRejects<std::runtime_error>([&] { DictionaryGeometry value; value.readDict(absent.string()); },
        absent,"cannot open geometry dictionary");
    readerRejects<std::runtime_error>([&] { DictionaryMesh value; value.readDict(absent.string(),geometry.string()); },
        absent,"cannot open mesh dictionary");
    write(mesh,validMesh);
    readerRejects<std::runtime_error>(readMesh,geometry,"cannot open geometry dictionary");
    write(geometry,validGeometry);
    write(mesh,"BackgroundMesh\n{\n divide (0, 10)\n}\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"invalid grid division");
    write(mesh,"adaptiveMesh false\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"missing divide");
    write(mesh,validMesh);
    write(geometry,"Domain\n{\n min (0.2, -0.5)\n}\n");
    readerRejects<std::invalid_argument>(readMesh,geometry,"missing Domain/min/max");
    write(geometry,"Domain\n{\n min (0.2, -0.5)\n max (0.1, 0.5)\n}\n");
    readerRejects<std::invalid_argument>(readMesh,geometry,"invalid domain bounds");
    write(geometry,validGeometry);
    const std::string region="GeometryRefinement\n{\n refineRegion\n{\n entry\n{\n";
    write(mesh,validMesh+region+" name focus\n shape circle\n}\n}\n}\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"entry 'focus' uses unsupported shape 'circle'");
    write(mesh,validMesh+region+" shape box\n}\n}\n}\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"<unnamed> is missing required field 'name'");
    write(mesh,validMesh+region+" name focus\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"unterminated GeometryRefinement.refineRegion entry");
    // A new record remains usable in the same process after every failure.
    write(mesh,validMesh);
    DictionaryGeometry parsedGeometry; parsedGeometry.readDict(geometry.string());
    DictionaryMesh parsedMesh; parsedMesh.readDict(mesh.string(),geometry.string());
    if (!parsedGeometry.hasDomain || parsedGeometry.domainMin.x!=.2 || parsedMesh.Nx!=8 || parsedMesh.Ny!=10
        || parsedMesh.minCoord.y!=-.5 || parsedMesh.maxCoord.x!=1.2 || fs::current_path()!=cwd
        || fs::exists(folder/"output") || fs::exists(folder/"cache"))
        throw std::runtime_error("dictionary failure changed process state or later valid input");
}
void extensions()
{
    struct Restore {
        std::filesystem::path original=std::filesystem::current_path();
        ~Restore() { std::filesystem::current_path(original); }
    } restore;
    const auto folder=restore.original/("dictionary_extensions_"+
        std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(folder/"system"); std::filesystem::current_path(folder);
    {
        std::ofstream out("system/geometry.dict");
        out<<"ExtraPhysics\n{\n enabled off\n name \"literal {\"\n Domain\n{\n min (99, 99)\n max (100, 100)\n enabled true\n}\n}\n"
              "Domain\n{\n min (0, 0)\n max (1, 1)\n}\n";
    }
    DictionaryGeometry geometry; geometry.readDict();
    if (!geometry.hasDomain || geometry.domainMin.x!=0. || geometry.domainMax.x!=1.
        || !geometry.extensionBlock("ExtraPhysics") || geometry.extensionEnabled("ExtraPhysics")
        || geometry.extensionEnabled("Absent"))
        throw std::runtime_error("opaque extension nesting/capability/source isolation");
}
void blockScope()
{
    DictionaryBlockText reader;
    using Scope=DictionaryBlockText::Scope;
    const auto check=[](bool condition,const char* message) {
        if(!condition)throw std::runtime_error(message);
    };
    std::string body="stale";
    check(reader.extractBlockText("\xEF\xBB\xBFReactions { bom 1 }","Reactions",body,Scope::TopLevel)
        && body==" bom 1 ","UTF-8 BOM must not hide the first block");
    const std::string text=R"(Output {
 Diagnostics { Reactions { enabled true } }
 label "literal } Reactions {"
}
reactionNames (Reactions { ignored })
name Reactions
Other { value 9 }
// Reactions { ignored }
Reactions # ignored { }
{
 physical true
 nested { value "escaped \" } {" }
}
)";
    check(reader.extractBlockText(text,"Reactions",body,Scope::TopLevel) &&
        body.find("physical true")!=std::string::npos && body.find("escaped")!=std::string::npos,
        "root block must skip diagnostic, list, scalar, quoted and commented names/braces");
    check(reader.extractBlockText(text,"Reactions",body) && body==" enabled true ",
        "nested solver block extraction remains available explicitly or by default");
    for(const auto& input:{std::string("Reactions false\nOther { x 1 }"),
        std::string("Output { Reactions { enabled true } }"),std::string("ReactionsExtra { x 1 }"),
        std::string("Reactions { missing"),std::string("\"Reactions\" { x 1 }")})
        check(!reader.extractBlockText(input,"Reactions",body,Scope::TopLevel) && body.empty(),
              "false block match or stale body on failed extraction");
    check(reader.extractFirstBlockText("gas { second 2 }\nGas { first 1 }\nGas { later 3 }",
        {"Gas","gas"},body,Scope::TopLevel) && body==" first 1 ","caller alias priority or first occurrence changed");
    check(!reader.extractFirstBlockText("Gas {}",{},body) && body.empty(),"empty aliases retain stale body");
    check(reader.extractBlockText("A {}; Reactions { x 1 }","Reactions",body,Scope::TopLevel)
        && body==" x 1 ","adjacent top-level blocks");
}
int main() try {
    readerFailures();
    extensions();
    blockScope();
    Reader reader;
    for (bool brace : {false,true}) {
        const auto parse=[&](const std::string& text, std::vector<SignInt>& result) {
            return brace ? reader.parseIndex_brace("{"+text+"}",result)
                         : reader.parseIndex("("+text+")",result);
        };
        for (const auto text : {"bad","1; 2","1, -","(1,2)","1.5","1,2 suffix"}) {
            std::vector<SignInt> values;
            if (parse(text,values)) throw std::runtime_error("malformed indices accepted");
        }
        std::vector<SignInt> values;
        if (!parse("1, -12  3",values) || values.size()!=3 || values[0].value!=1
            || values[0].isNegative || values[1].value!=12 || !values[1].isNegative
            || values[2].value!=3 || values[2].isNegative)
            throw std::runtime_error("valid signed index list changed");
    }
    std::cout<<"Malformed index lists return instead of hanging; signed lists preserved\n";
    return 0;
} catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
