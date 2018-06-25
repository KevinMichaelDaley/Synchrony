#pragma once
#include <cstdint>
#include <cstdio>
#include "Config.h"
class FixedSizeRealBuffer{
private:
    real_t* _buf;
    size_t _len;
public:
    HOSTDEVICE FixedSizeRealBuffer(): _buf(nullptr), _len(0uL)
                                                    {}//otherwise we cant declare an array of these the easy way.
    HOSTDEVICE FixedSizeRealBuffer(real_t* data, size_t N): _buf(data), 
                                                   _len(N)     
                                                    {}
    //note: deleting this container doesn't delete the underlying storage,
    //since it might be allocated elsewhere (or static).
    HOSTDEVICE ~FixedSizeRealBuffer(){}
    HOSTDEVICE void BoundsAssert(size_t index) const;
    HOSTDEVICE size_t GetLength() const;
    HOSTDEVICE real_t Get(size_t index) const;
    HOSTDEVICE void Set(size_t index, real_t val);
    HOSTDEVICE void AccumulateWeighted(real_t weight, const FixedSizeRealBuffer& term, long start=0, long end=-1, size_t o_start=0);
    HOSTDEVICE void FillWith(real_t val);
    HOSTDEVICE void CopyFrom(const FixedSizeRealBuffer& other, long start=0, long end=-1, size_t o_start=0);
};
    
