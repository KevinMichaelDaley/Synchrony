#pragma once
#include <cstdint>
#include <cstdio>
class FixedSizeDoubleBuffer{
private:
    double* _buf;
    size_t _len;
public:
    FixedSizeDoubleBuffer(): _buf(nullptr), _len(0uL)
                                                    {}//otherwise we cant declare an array of these the easy way.
    FixedSizeDoubleBuffer(double* data, size_t N): _buf(data), 
                                                   _len(N)     
                                                    {}
    //note: deleting this container doesn't delete the underlying storage,
    //since it might be allocated elsewhere (or static).
    ~FixedSizeDoubleBuffer(){}
    void BoundsAssert(size_t index) const;
    size_t GetLength() const;
    double Get(size_t index) const;
    void Set(size_t index, double val);
    void AccumulateWeighted(double weight, const FixedSizeDoubleBuffer& term);
    void FillWith(double val);
    void CopyFrom(const FixedSizeDoubleBuffer& other);
};
    
