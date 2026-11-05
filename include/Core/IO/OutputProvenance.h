#pragma once
#include "Core/IO/ConfigurationTrace.h"
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// Passive serialization and exact-byte SHA-256 for App output provenance.
// No configuration reading, numerical state or filesystem policy.
namespace LAMP::IO {
inline std::ostringstream stream(){std::ostringstream s;s.imbue(std::locale::classic());s<<std::scientific<<std::setprecision(16);return s;}
inline std::string quoted(const std::string& text)
{
    std::ostringstream s;s<<'"';for(unsigned char c:text){if(c=='"'||c=='\\')s<<'\\'<<c;
        else if(c<32)s<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<static_cast<unsigned>(c);else s<<c;}s<<'"';return s.str();
}
inline std::string csv(const std::string& text)
{std::string s="\"";for(char c:text){if(c=='"')s+='"';s+=c;}return s+'"';}
inline std::string value(const DictionaryInput::ConfigurationValue& v)
{
    auto s=stream();std::visit([&](const auto& x){using T=std::decay_t<decltype(x)>;
        if constexpr(std::is_same_v<T,bool>)s<<(x?"true":"false");
        else if constexpr(std::is_same_v<T,std::vector<int>>){s<<'(';for(std::size_t i=0;i<x.size();++i)s<<(i?" ":"")<<x[i];s<<')';}
        else s<<x;},v);return s.str();
}
inline std::uint32_t rotate(std::uint32_t v,unsigned n){return (v>>n)|(v<<(32-n));}
// Standard SHA-256 over the exact file bytes. Independently checked against
// the platform SHA-256 in the App output contract, not used numerically.
inline std::string sha256(const std::string& text)
{
    constexpr std::uint32_t k[64]={0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
        0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
        0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::array<std::uint32_t,8> h{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::vector<unsigned char> bytes(text.begin(),text.end());const auto bits=static_cast<std::uint64_t>(bytes.size())*8;
    bytes.push_back(0x80);while(bytes.size()%64!=56)bytes.push_back(0);
    for(int i=7;i>=0;--i)bytes.push_back(static_cast<unsigned char>(bits>>(8*i)));
    for(std::size_t p=0;p<bytes.size();p+=64){std::array<std::uint32_t,64> w{};
        for(int i=0;i<16;++i)for(int j=0;j<4;++j)w[i]=(w[i]<<8)|bytes[p+4*i+j];
        for(int i=16;i<64;++i){const auto a=w[i-15],b=w[i-2];w[i]=w[i-16]+(rotate(a,7)^rotate(a,18)^(a>>3))+w[i-7]+(rotate(b,17)^rotate(b,19)^(b>>10));}
        auto a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],v=h[7];
        for(int i=0;i<64;++i){const auto t1=v+(rotate(e,6)^rotate(e,11)^rotate(e,25))+((e&f)^(~e&g))+k[i]+w[i];
            const auto t2=(rotate(a,2)^rotate(a,13)^rotate(a,22))+((a&b)^(a&c)^(b&c));
            v=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;}
        const std::array<std::uint32_t,8> out{a,b,c,d,e,f,g,v};for(int i=0;i<8;++i)h[i]+=out[i];
    }
    std::ostringstream s;s<<std::hex<<std::setfill('0');for(auto v:h)s<<std::setw(8)<<v;return s.str();
}
inline std::string utc()
{const auto t=std::time(nullptr);std::tm tm{};
#ifdef _WIN32
 gmtime_s(&tm,&t);
#else
 gmtime_r(&t,&tm);
#endif
 std::ostringstream s;s<<std::put_time(&tm,"%Y-%m-%dT%H:%M:%SZ");return s.str();}
inline std::string compiler()
{
#ifdef _MSC_VER
 return "MSVC "+std::to_string(_MSC_VER);
#else
 return __VERSION__;
#endif
}
inline const char* buildType()
{
#ifdef NDEBUG
 return "Release";
#else
 return "Debug";
#endif
}
}
