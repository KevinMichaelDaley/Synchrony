#include "FixedSizeBuffer.h"
#include <cassert>
#include <cmath>
inline void FixedSizeDoubleBuffer::BoundsAssert(size_t index) const{
    assert(index<_len);
}

size_t FixedSizeDoubleBuffer::GetLength() const{
    return _len;
}

double FixedSizeDoubleBuffer::Get(size_t index) const{
    BoundsAssert(index);
    return _buf[index];
}

void FixedSizeDoubleBuffer::Set(size_t index, double val){
    BoundsAssert(index);
    assert(!std::isnan(val));
    _buf[index]=val;
}

void FixedSizeDoubleBuffer::AccumulateWeighted(double weight, const FixedSizeDoubleBuffer& term){
    //*this += weight .* term
    BoundsAssert(term.GetLength()-1);
    for(size_t entry_i=0; entry_i<_len; ++entry_i){
        _buf[entry_i]+=weight*term.Get(entry_i);
    }
}
void FixedSizeDoubleBuffer::FillWith(double val){
    for(size_t entry_i=0; entry_i<_len; ++entry_i){
        _buf[entry_i]=val;
    }
}
void FixedSizeDoubleBuffer::CopyFrom(const FixedSizeDoubleBuffer& other){
    BoundsAssert(other.GetLength()-1);
    for(size_t entry_i=0; entry_i<_len; ++entry_i){
        _buf[entry_i]=other.Get(entry_i);
    }
}
     
