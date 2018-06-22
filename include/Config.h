#include <cstdio>
#include <cstdint>
namespace Config{
    static constexpr size_t MaxThreads=256;
    static constexpr size_t MaxStackBufferSize=4096;
    static constexpr size_t ThreadLocalTemporaryHeapSize=8*1024*1024;
}
