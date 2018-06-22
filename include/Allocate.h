#pragma once
#include <cassert>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include <memory>
#include "Config.h"

static char* linear_mem_pool[Config::MaxThreads];
static size_t alloc_head[Config::MaxThreads]={0};

template<typename T> inline T* AllocateTemporaryBuffer(size_t N, int tidx){
    if(__builtin_expect(linear_mem_pool[tidx]==nullptr, false)){
        linear_mem_pool[tidx]=new char[Config::ThreadLocalTemporaryHeapSize];
        std::memset(linear_mem_pool[tidx], 0xba, Config::ThreadLocalTemporaryHeapSize);
    }
    size_t sz=N*sizeof(T);
    alloc_head[tidx]+=sz;
    if(__builtin_expect(alloc_head[tidx]>=Config::ThreadLocalTemporaryHeapSize,false)){
        //we've run out of memory in this heap so just use malloc.
        return static_cast<T*>(std::malloc(sz));
    }
    size_t head0=alloc_head[tidx]-sz;
    char* res=&(linear_mem_pool[tidx][head0]);
    return reinterpret_cast<T*>(res);
    
}

template<typename T> inline void FreeTemporaryBuffer(const T* ptr, int tidx){    //temporary buffers must be freed in reverse order of allocation,
                                                                    //so that we always know where the head is.  these are generally used for local per-function storage that "goes out of scope" when we are finished with it,
                                                                    //so they can be treated mathematically like stack allocations, even though they are too big to live on the actual stack.
                                                                    //alternatively, the oldest one can be freed, which automatically frees all the newer ones.
                                                                    //given the above definition of allocatetemporaryBuffer(), this assurance is guaranteed by requiring that the new allocation head be before the old one.
                                                                    //we can also free e.g. the second half of a temporary buffer; it is perfectly valid.
    size_t head0=alloc_head[tidx];
    std::ptrdiff_t head=reinterpret_cast<std::ptrdiff_t>(ptr) - reinterpret_cast<std::ptrdiff_t>(linear_mem_pool[tidx]);
    if(__builtin_expect(head>=0 && head0<Config::ThreadLocalTemporaryHeapSize && static_cast<size_t>(head)<head0,true)){
            alloc_head[tidx]=head;
            return;
    }
    free(const_cast<T*>(ptr));//we are outside the range so we must have allocated this using malloc.
}
