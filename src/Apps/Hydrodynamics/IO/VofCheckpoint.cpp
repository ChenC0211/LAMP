#include "Hydrodynamics/IO/VofCheckpoint.h"
#include "Core/IO/BinaryArchive.h"
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#endif

namespace {
namespace fs=std::filesystem;
namespace C=LAMP::IO::Checkpoint;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error("[VOF restart] "+why);}
std::string readFile(const fs::path& path){
    std::ifstream in(path,std::ios::binary);check(bool(in),"cannot read "+path.u8string());
    const auto size=fs::file_size(path);check(size<=std::numeric_limits<std::uint32_t>::max(),"state/manifest file exceeds supported byte range");
    std::string result(static_cast<std::size_t>(size),'\0');
    if(size)in.read(result.data(),static_cast<std::streamsize>(size));
    check(bool(in)&&in.peek()==std::char_traits<char>::eof(),"truncated or changing file "+path.u8string());return result;
}
void writeFile(const fs::path& path,const std::string& bytes){
    std::ofstream out(path,std::ios::binary|std::ios::trunc);check(bool(out),"cannot create "+path.u8string());
    out.write(bytes.data(),static_cast<std::streamsize>(bytes.size()));out.flush();check(bool(out),"failed writing "+path.u8string());
    out.close();check(!out.fail(),"failed closing "+path.u8string());
}
// Streaming checksum avoids retaining a capacity-sized native mesh file in
// memory. The sealed state and its exact mesh-identity words use SHA256 too.
std::string meshChecksum(const fs::path& path){
    std::ifstream in(path,std::ios::binary);check(bool(in),"cannot read mesh cache");
    std::array<char,65536> buffer{};std::uint64_t hash=14695981039346656037ull;
    while(in){in.read(buffer.data(),buffer.size());const auto count=in.gcount();
        for(std::streamsize i=0;i<count;++i){hash^=static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]);hash*=1099511628211ull;}}
    check(in.eof()&&!in.bad(),"failed hashing mesh cache");std::ostringstream out;out<<std::hex<<std::setfill('0')<<std::setw(16)<<hash;return out.str();
}
void replaceFile(const fs::path& source,const fs::path& target){
#ifdef _WIN32
    check(MoveFileExW(source.c_str(),target.c_str(),MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)!=0,
        "cannot atomically publish "+target.u8string());
#else
    fs::rename(source,target);
#endif
}
std::string compatibility(const HydrodynamicsInput& input,const std::string& app){
    C::Writer out("LAMP/VOF/configuration/1");out.text(app);
    // Effective non-run dictionaries contain every actual geometry, material,
    // boundary, discretization and AMR choice. Paths/comments are irrelevant.
    for(const auto& doc:input.documents)if(fs::u8path(doc.file).filename()!="run.dict"){
        out.text(fs::u8path(doc.file).filename().u8string());out.text(doc.effectiveText);}
    out.scalar(input.run.dt);out.boolean(input.run.adaptive.has_value());
    if(input.run.adaptive){const auto& a=*input.run.adaptive;out.scalar(a.minDt);out.scalar(a.maxDt);out.scalar(a.targetCourant);out.scalar(a.safety);out.scalar(a.maxGrowth);}
    out.boolean(input.run.stepRetry.has_value());if(input.run.stepRetry){const auto& r=*input.run.stepRetry;
        out.scalar(r.minDt);out.scalar(r.maxRetries);out.scalar(r.recoverAfter);}
    // endTime/output/log/thread/checkpoint cadence may change in a new segment.
    return LAMP::IO::sha256(out.finish());
}
std::map<std::string,std::string> manifest(const fs::path& directory){
    const auto text=readFile(directory/"manifest.txt");std::istringstream in(text);std::string line;
    check(static_cast<bool>(std::getline(in,line))&&line=="LAMP_VOF_CHECKPOINT 1","unsupported or incomplete checkpoint manifest");
    std::map<std::string,std::string> result;
    while(std::getline(in,line)){const auto p=line.find('=');check(p!=std::string::npos&&p>0,"invalid checkpoint manifest record");
        check(result.emplace(line.substr(0,p),line.substr(p+1)).second,"duplicate checkpoint manifest record");}
    for(const char* key:{"app","compatibility","state_bytes","state_sha256","mesh_bytes","mesh_fnv1a64","step","time"})
        check(result.count(key)==1,std::string("missing manifest field ")+key);
    check(result.size()==8,"unknown checkpoint manifest fields");return result;
}
void verifyFiles(const fs::path& directory,const std::map<std::string,std::string>& m,const std::string& state){
    check(std::to_string(state.size())==m.at("state_bytes")&&LAMP::IO::sha256(state)==m.at("state_sha256"),"state size/SHA256 mismatch");
    const auto mesh=directory/"mesh_cache.bin";
    check(fs::is_regular_file(mesh)&&std::to_string(fs::file_size(mesh))==m.at("mesh_bytes")&&meshChecksum(mesh)==m.at("mesh_fnv1a64"),"mesh size/checksum mismatch");
}
}
VofCheckpointStore::VofCheckpointStore(const HydrodynamicsInput& input,std::string app)
    :input_(input),app_(std::move(app)),
    compatibility_(input.run.checkpointEverySteps>0||!input.run.restartDirectory.empty()?compatibility(input_,app_):std::string{}),
    root_(fs::absolute("output/checkpoints")){}
std::optional<VofCheckpointState> VofCheckpointStore::read() const {
    if(input_.run.restartDirectory.empty())return std::nullopt;
    VofCheckpointState result;result.directory=fs::absolute(fs::u8path(input_.run.restartDirectory));
    check(fs::is_directory(result.directory)&&!fs::is_symlink(result.directory),"checkpoint directory missing or linked");
    const auto m=manifest(result.directory);check(m.at("app")==app_,"checkpoint belongs to a different App");
    check(m.at("compatibility")==compatibility_,"physical/discretization/time-policy inputs differ from checkpoint");
    const auto state=readFile(result.directory/"state.bin");verifyFiles(result.directory,m,state);
    C::Reader in(state,"LAMP/VOF/accepted-checkpoint/1");check(in.text()==app_&&in.text()==compatibility_,"state configuration differs from manifest");
    result.clock.step=in.scalar<std::int64_t>();result.clock.begin=in.scalar<double>();result.clock.end=in.scalar<double>();
    const auto position=in.scalar<std::uint8_t>();check(position==0||position==2,"checkpoint is not accepted Initial/End");
    result.clock.position=static_cast<LAMP::TimeLayer::Position>(position);(void)result.clock.time();
    check(result.clock.step>=0&&((position==0&&result.clock.step==0&&result.clock.begin==0&&result.clock.end==0)
        ||(position==2&&result.clock.step>0&&result.clock.end>result.clock.begin)),"invalid accepted clock");
    result.nextDt=in.scalar<double>();result.withoutRetry=in.scalar<int>();result.transaction=in.text();in.finish();
    check(result.clock.step<=1000000&&std::isfinite(result.nextDt)&&result.nextDt>0
        &&result.withoutRetry>=0&&result.withoutRetry<=result.clock.step,"invalid adaptive/retry continuation state");
    check(result.clock.time()<input_.run.endTime,"checkpoint already reached endTime; choose a later endTime for another segment");
    check(m.at("step")==std::to_string(result.clock.step),"manifest/source step mismatch");
    std::ostringstream acceptedTime;acceptedTime<<std::setprecision(17)<<result.clock.time();
    check(m.at("time")==acceptedTime.str(),"manifest/source time mismatch");
    return result;
}
void VofCheckpointStore::restoreMesh(QuadtreeMesh& mesh,const VofCheckpointState& state) const {
    mesh.readCache(state.directory.u8string()+"/");
}
void VofCheckpointStore::write(const QuadtreeMesh& mesh,const VofCheckpointState& state){
    const auto start=std::chrono::steady_clock::now();
    check(state.clock.position!=LAMP::TimeLayer::Position::Midpoint&&state.nextDt>0&&std::isfinite(state.nextDt)
        &&state.clock.step>=0&&state.clock.step<=1000000&&state.withoutRetry>=0
        &&state.withoutRetry<=state.clock.step&&!state.transaction.empty(),"cannot save incomplete accepted checkpoint");
    fs::create_directories(root_);std::ostringstream name;name<<"checkpoint_"<<std::setfill('0')<<std::setw(12)<<state.clock.step;
    const auto destination=root_/name.str(),pending=root_/(".pending_"+name.str());
    check(!fs::exists(destination)&&!fs::exists(pending),"checkpoint generation already exists; use a fresh run directory");
    fs::create_directory(pending);
    C::Writer out("LAMP/VOF/accepted-checkpoint/1");out.text(app_);out.text(compatibility_);
    out.scalar<std::int64_t>(state.clock.step);out.scalar(state.clock.begin);out.scalar(state.clock.end);
    out.scalar<std::uint8_t>(static_cast<std::uint8_t>(state.clock.position));out.scalar(state.nextDt);
    out.scalar(state.withoutRetry);out.text(state.transaction);const auto bytes=out.finish();
    check(bytes.size()<=std::numeric_limits<std::uint32_t>::max(),"checkpoint state exceeds supported byte range");
    mesh.writeCache(pending.u8string()+"/");writeFile(pending/"state.bin",bytes);
    const auto cache=pending/"mesh_cache.bin";std::ostringstream meta;meta<<std::setprecision(17)
        <<"LAMP_VOF_CHECKPOINT 1\napp="<<app_<<"\ncompatibility="<<compatibility_
        <<"\nstate_bytes="<<bytes.size()<<"\nstate_sha256="<<LAMP::IO::sha256(bytes)
        <<"\nmesh_bytes="<<fs::file_size(cache)<<"\nmesh_fnv1a64="<<meshChecksum(cache)
        <<"\nstep="<<state.clock.step<<"\ntime="<<state.clock.time()<<'\n';
    writeFile(pending/"manifest.txt",meta.str());fs::rename(pending,destination);
    writeFile(root_/"latest.tmp",name.str()+"\n");replaceFile(root_/"latest.tmp",root_/"latest.txt");
    generations_.push_back(destination);
    // Only remove paths created by THIS store in this new run. The parent run
    // and incomplete generations are never enumerated or silently repaired.
    if(generations_.size()>2){const auto obsolete=generations_.front();
        check(obsolete.parent_path()==root_&&obsolete.filename().u8string().find("checkpoint_")==0,"checkpoint cleanup escaped owned root");
        fs::remove_all(obsolete);generations_.erase(generations_.begin());}
    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"[CHECKPOINT] step="<<state.clock.step<<" time="<<state.clock.time()<<" wall="<<seconds<<"s "<<destination.u8string()<<std::endl;
}
