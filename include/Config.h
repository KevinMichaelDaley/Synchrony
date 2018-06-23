#pragma once
#include <cstdio>
#include <cstdint>
namespace Config{
#if defined(__CUDA__)
    #define HOSTDEVICE __host__ __device__
    #define DEVICE __device__
    #define HOST __host__
    #define SHARED __shared__
#else
    #define DEVICE static_assert(false, "using device code outside of cuda is not allowed!");
    #define HOSTDEVICE 
    #define HOST
    #define SHARED 
#endif
#define CUDA_MAX_STATIC_ARRAY_SIZE 2048*1024*1024uL
    static constexpr size_t MaxThreads=256;
    static constexpr size_t MaxGPUThreads=512;
    static constexpr size_t MaxStackBufferSize=4096;
    static constexpr size_t ThreadLocalTemporaryHeapSize=8*1024*1024;
    static constexpr size_t ThreadLocalTemporaryHeapSizeGPU=1024*1024;

}
