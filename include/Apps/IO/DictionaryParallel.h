#pragma once

#include "Core/IO/ConfigurationTrace.h"
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

// Process resources are owned by AppRuntime, never by a numerical module.
namespace AppRuntime
{
inline int readParallelThreads(const DictionaryInput::Node& root,
    DictionaryInput::ConfigurationTrace* trace=nullptr,std::optional<int> override={})
{
    const auto* parallel=root.find("Parallel");
    const DictionaryInput::Node* field=nullptr;
    int requested=2;
    if(parallel){
        parallel->allow({"threads"});
        field=&parallel->get("threads");
        requested=field->integer();
        if(requested<1 || requested>256)field->fail("an integer thread count in [1,256] for AppRuntime.OpenMP");
    }
    if(override && (*override<1 || *override>256))
        throw std::invalid_argument("--threads requires an integer in [1,256]");
    const int effective=override.value_or(requested);
    if(trace){
        const std::string consumer=override
            ?"AppRuntime.OpenMP; explicit --threads overrides the dictionary/default"
            :"AppRuntime.OpenMP";
        if(field)trace->record(*field,effective,"threads",consumer);
        else{
            auto owner=root;owner.path="Parallel";
            trace->defaultValue(owner,"threads",effective,"threads",consumer);
        }
    }
    return effective;
}

inline void configureParallelThreads(int threads,const char* source)
{
    const auto set=[](const char* name,const std::string& value){
#ifdef _WIN32
        if(_putenv_s(name,value.c_str())!=0)
#else
        if(setenv(name,value.c_str(),1)!=0)
#endif
            throw std::runtime_error(std::string("Cannot set process resource: ")+name);
    };
    set("OMP_NUM_THREADS",std::to_string(threads));
    set("OMP_THREAD_LIMIT",std::to_string(threads));
    set("OMP_DYNAMIC","false");
    set("OPENBLAS_NUM_THREADS","1");
    set("MKL_NUM_THREADS","1");
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#if _OPENMP >= 200805
    omp_set_max_active_levels(1);
#endif
    const int actual=omp_get_max_threads();
#else
    const int actual=1;
#endif
    std::cout<<"[Runtime] threads="<<actual<<" requested_threads="<<threads
        <<" source="<<source<<std::endl;
}
}
