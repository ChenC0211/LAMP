#pragma once
#include <array>
#include "Core/IO/OutputProvenance.h"
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

// Native scalar persistence only. Physics owners define their own versioned
// records and validate restored invariants; Apps own files and atomic publish.
namespace LAMP::IO::Checkpoint {
inline void require(bool ok,const char* why){if(!ok)throw std::runtime_error(std::string("[Checkpoint] ")+why);}
class Writer {
    std::string bytes_;
public:
    explicit Writer(const std::string& tag){
        require(sizeof(long double)==sizeof(double)
            &&std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits,
            "this checkpoint version requires an unpadded binary64 long double ABI");
        text(tag);scalar<std::uint32_t>(0x01020304u);
        scalar<std::uint32_t>(sizeof(int));scalar<std::uint32_t>(sizeof(unsigned));scalar<std::uint32_t>(sizeof(std::size_t));
        scalar<std::uint32_t>(sizeof(double));scalar<std::uint32_t>(sizeof(long double));
        scalar<std::uint32_t>(std::numeric_limits<long double>::digits);
        scalar<std::int32_t>(std::numeric_limits<long double>::min_exponent);
        scalar<std::int32_t>(std::numeric_limits<long double>::max_exponent);
    }
    template<class T> void scalar(T value){
        static_assert((std::is_arithmetic_v<T>||std::is_enum_v<T>)&&!std::is_same_v<T,bool>,"scalar fields only; encode booleans explicitly");
        bytes_.append(reinterpret_cast<const char*>(&value),sizeof value);
    }
    void boolean(bool value){scalar<std::uint8_t>(value?1:0);}
    void count(std::size_t size){scalar<std::uint64_t>(static_cast<std::uint64_t>(size));}
    void text(const std::string& value){count(value.size());bytes_.append(value);}
    template<class T> void values(const std::vector<T>& value){
        static_assert(std::is_arithmetic_v<T>&&!std::is_same_v<T,bool>,"native scalar vectors only");
        count(value.size());require(value.size()<=std::numeric_limits<std::size_t>::max()/sizeof(T),"vector byte overflow");
        if(!value.empty())bytes_.append(reinterpret_cast<const char*>(value.data()),value.size()*sizeof(T));
    }
    std::string finish() const {return bytes_+LAMP::IO::sha256(bytes_);}
};
class Reader {
    const std::string& bytes_;std::size_t at_=0,end_=0;
public:
    Reader(const std::string& bytes,const std::string& tag):bytes_(bytes){
        require(sizeof(long double)==sizeof(double)
            &&std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits,
            "this checkpoint version requires an unpadded binary64 long double ABI");
        require(bytes.size()>=64,"truncated checksum");end_=bytes.size()-64;
        require(bytes.substr(end_)==LAMP::IO::sha256(bytes.substr(0,end_)),"payload SHA256 mismatch");
        require(text()==tag,"unsupported payload kind/version");
        require(scalar<std::uint32_t>()==0x01020304u,"native byte order differs");
        require(scalar<std::uint32_t>()==sizeof(int)&&scalar<std::uint32_t>()==sizeof(unsigned)
            &&scalar<std::uint32_t>()==sizeof(std::size_t),"native integer ABI differs");
        require(scalar<std::uint32_t>()==sizeof(double)&&scalar<std::uint32_t>()==sizeof(long double)
            &&scalar<std::uint32_t>()==std::numeric_limits<long double>::digits
            &&scalar<std::int32_t>()==std::numeric_limits<long double>::min_exponent
            &&scalar<std::int32_t>()==std::numeric_limits<long double>::max_exponent,"native floating format differs");
    }
    template<class T> T scalar(){
        static_assert((std::is_arithmetic_v<T>||std::is_enum_v<T>)&&!std::is_same_v<T,bool>,"scalar fields only");
        require(at_<=end_&&sizeof(T)<=end_-at_,"truncated scalar");T value;
        std::memcpy(&value,bytes_.data()+at_,sizeof value);at_+=sizeof value;return value;
    }
    bool boolean(){const auto value=scalar<std::uint8_t>();require(value<=1,"invalid boolean");return value!=0;}
    std::size_t count(std::size_t minimumBytes=1){
        const auto value=scalar<std::uint64_t>();
        require(minimumBytes>0&&at_<=end_&&value<=static_cast<std::uint64_t>((end_-at_)/minimumBytes),"invalid or truncated record count");
        return static_cast<std::size_t>(value);
    }
    std::string text(){const auto size=count();std::string value=bytes_.substr(at_,size);at_+=size;return value;}
    template<class T> std::vector<T> values(){
        static_assert(std::is_arithmetic_v<T>&&!std::is_same_v<T,bool>,"native scalar vectors only");
        const auto size=count(sizeof(T));std::vector<T> value(size);
        if(size)std::memcpy(value.data(),bytes_.data()+at_,size*sizeof(T));at_+=size*sizeof(T);return value;
    }
    void finish() const {require(at_==end_,"trailing payload bytes");}
};
}
