#pragma once
class ThreadInfo{
private:
    const size_t _local_tidx;
    ThreadInfo()=delete;
public:
    static HOSTDEVICE const ThreadInfo& Serial(){
        static const ThreadInfo s={0};
        return s;
    }
    HOSTDEVICE ThreadInfo(size_t tidx): _local_tidx(tidx)  {}
    HOSTDEVICE size_t GetIndex() const{
        return _local_tidx;
    }
};
