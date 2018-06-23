#pragma once
#include <cstdint>
#include <cstdio>
#include "Config.h"
class FixedSizeDoubleBuffer{
private:
    double* _buf;
    size_t _len;
public:
    HOSTDEVICE FixedSizeDoubleBuffer(): _buf(nullptr), _len(0uL)
                                                    {}//otherwise we cant declare an array of these the easy way.
    HOSTDEVICE FixedSizeDoubleBuffer(double* data, size_t N): _buf(data), 
                                                   _len(N)     
                                                    {}
    //note: deleting this container doesn't delete the underlying storage,
    //since it might be allocated elsewhere (or static).
    HOSTDEVICE ~FixedSizeDoubleBuffer(){}
    HOSTDEVICE void BoundsAssert(size_t index) const;
    HOSTDEVICE size_t GetLength() const;
    HOSTDEVICE double Get(size_t index) const;
    HOSTDEVICE void Set(size_t index, double val);
    HOSTDEVICE void AccumulateWeighted(double weight, const FixedSizeDoubleBuffer& term, long start=0, long end=-1, size_t o_start=0);
    HOSTDEVICE void FillWith(double val);
    HOSTDEVICE void CopyFrom(const FixedSizeDoubleBuffer& other, long start=0, long end=-1, size_t o_start=0);
};
    
