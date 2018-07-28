#include "FixedSizeBuffer.h"
#if defined(__CUDA__)
#else
    #include <cassert>
#endif
#include <cmath>
inline HOSTDEVICE void FixedSizeRealBuffer::BoundsAssert(size_t index) const{
    assert(index<_len);
}

HOSTDEVICE size_t FixedSizeRealBuffer::GetLength() const{
    return _len;
}

HOSTDEVICE real_t FixedSizeRealBuffer::Get(size_t index) const{
    BoundsAssert(index);
    assert(!std::isnan(_buf[index]));
    return _buf[index];
}

HOSTDEVICE void FixedSizeRealBuffer::Set(size_t index, real_t val){
    BoundsAssert(index);
    assert(!std::isnan(val));
    _buf[index]=val;
}

HOSTDEVICE void FixedSizeRealBuffer::AccumulateWeighted(real_t weight, const FixedSizeRealBuffer& term, long start, long end, size_t o){
    if(end<0){
        end=term.GetLength();
    }
    //*this += weight .* term
    
    size_t s=static_cast<size_t>(start);
    size_t e=static_cast<size_t>(end);
    BoundsAssert(e-1-s+o);
    for(size_t entry_i=s; entry_i<e; ++entry_i){
        _buf[entry_i-s+o]+=weight*term.Get(entry_i);
    }
}
HOSTDEVICE void FixedSizeRealBuffer::FillWith(real_t val){
    for(size_t entry_i=0; entry_i<_len; ++entry_i){
        _buf[entry_i]=val;
    }
}
HOSTDEVICE void FixedSizeRealBuffer::CopyFrom(const FixedSizeRealBuffer& other, long start, long end, size_t o){
    if(end<0){
        end=other.GetLength();
    }
    
    size_t s=static_cast<size_t>(start);
    size_t e=static_cast<size_t>(end);
    BoundsAssert(e-1-s+o);
    for(size_t entry_i=s; entry_i<e; ++entry_i){
        _buf[entry_i-s+o]=other.Get(entry_i);
    }
}
     
