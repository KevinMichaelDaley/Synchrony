#pragma once
#include <cstdint>
#include "Config.h"
#include "Allocate.h"
static_assert(Config::ThreadLocalTemporaryHeapSizeGPU * Config::MaxGPUThreads <= CUDA_MAX_STATIC_ARRAY_SIZE);
static DEVICE char d_linear_mem_pool[Config::MaxGPUThreads*Config::ThreadLocalTemporaryHeapSizeGPU]={0};
static DEVICE size_t d_alloc_head[Config::MaxGPUThreads]={0};

template<typename T> inline DEVICE T* AllocateTemporaryBuffer(size_t N, int tidx){
    size_t sz=N*sizeof(T);
    d_alloc_head[tidx]+=sz;
    if(__builtin_expect(d_alloc_head[tidx]>=Config::ThreadLocalTemporaryHeapSize,false)){
        //we've run out of memory in this heap (on the gpu) so fail.
        return nullptr;
    }
    size_t head0=d_alloc_head[tidx]-sz;
    char* res=&(d_linear_mem_pool[ tidx*Config::ThreadLocalTemporaryHeapSizeGPU + head0 ]);
    return reinterpret_cast<T*>(res);
    
}

template<typename T> inline DEVICE  int FreeTemporaryBuffer(const T* ptr, int tidx){
    size_t head0=d_alloc_head[tidx];
    std::ptrdiff_t head=reinterpret_cast<std::ptrdiff_t>(ptr) - reinterpret_cast<std::ptrdiff_t>(&d_linear_mem_pool[tidx* Config::ThreadLocalTemporaryHeapSizeGPU]);
    if(__builtin_expect(head>=0 && head0<Config::ThreadLocalTemporaryHeapSize && static_cast<size_t>(head)<head0,true)){
            d_alloc_head[tidx]=head;
            return 0;
    }
    return -1;
}
