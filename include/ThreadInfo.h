#pragma once
class ThreadInfo{
private:
    const size_t _local_tidx;
    ThreadInfo()=delete;
public:
    static const ThreadInfo& Serial(){
        static const ThreadInfo s={0};
        return s;
    }
    ThreadInfo(size_t tidx): _local_tidx(tidx)  {}
    size_t GetIndex() const{
        return _local_tidx;
    }
};
