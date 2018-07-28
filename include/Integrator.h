#pragma once

#if defined(__CUDA__)
#include "AllocateGPU.h"
#else
#include "Allocate.h"
#endif

#include "FixedSizeBuffer.h"
#include "ThreadInfo.h"
#include "RK.h"
#include <type_traits>
#include <cassert>

template<typename SystemT> class System {
    //static polymorphism so all types that are used as arguments to template need to derive from it. 
    //that way we can use the below virtual function call without actually dereferencing anything.
protected:
    
#if defined(__CUDA__)
    DEVICE void Fn( FixedSizeRealBuffer& out, const FixedSizeRealBuffer& in, real_t t, const ThreadInfo& thread_info) const{
        static_cast<const SystemT*>(this)->FnGPU( out, in, t, thread_info );
    }
#endif
    
    HOST void Fn( FixedSizeRealBuffer& out, const FixedSizeRealBuffer& in, real_t t, const ThreadInfo& thread_info) const{
        static_cast<const SystemT*>(this)->FnHost( out, in, t, thread_info );
    }
    HOSTDEVICE size_t GetDimension() const{
        return static_cast<const SystemT*>(this)->GetDimension();
    }
    HOSTDEVICE size_t GetIndexOffset(const ThreadInfo& thread_info) const{
           return static_cast<const SystemT*>(this)->GetIndexOffset(thread_info);
    }
    HOSTDEVICE real_t Integrate(FixedSizeRealBuffer& in_out_state, 
                           real_t t,
                           real_t step_size, 
                           const ThreadInfo& thread_info, 
                           bool adaptive){//override this virtual method with the default integration routine you want to use (if not rk4/ode45).
                                                          //for discrete maps this can just be Fn(in_out_state, in_out_state, k)
        if(!adaptive){
            IntegrateRKExplicitFixed(in_out_state, t, step_size, thread_info);
            return 0;
        } else{
            return IntegrateRKExplicitAdaptive(in_out_state, t, step_size, thread_info);
        }
    }
    //add new integration routines to this template so they can be reused.
    //parameters fixed over the whole trajectory are set in the constructor/whatever get-setter methods of the derived class
    //state should always be modified in-place; the caller is responsible for storing previous values and the associated allocations,
    //since they may not be necessary and to store two state arrays here when we might not need to is terribly inefficient for large dimensions.
    
public:
    HOST real_t Solve(FixedSizeRealBuffer& in_out_state,
               real_t t,
               real_t step_size, 
               const ThreadInfo& thread_info        =   ThreadInfo::Serial(),
                bool adaptive=true
               
         )
    {
        return static_cast<SystemT*>(this)->Integrate(in_out_state, t, step_size, thread_info, adaptive);
    }
#if defined(__CUDA__)
    DEVICE real_t Solve(FixedSizeRealBuffer& in_out_state,
                      real_t t,
                      real_t step_size,
                      bool adaptive=true//if we are not using an adaptive method, then we always return 0
                       ){
        size_t global_id = blockIdx.x *blockDim.x + threadIdx.x;
        return static_cast<SystemT*>(this)->Integrate(in_out_state, t, step_size, ThreadInfo(global_id), adaptive);
    }
#endif
                      
private:
    void HOSTDEVICE IntegrateRKExplicitFixed ( FixedSizeRealBuffer& in_out_state, 
                                    real_t t,
                                    real_t step_size, 
                                    const ThreadInfo& thread_info,
                                    const RKButcherTableau& table = RKButcherTableau::RK4Classic()
                                  )
    {//explicit runge-kutta method of arbitrary order and coefficients.
        
    static_assert(std::is_base_of_v<System<SystemT>, SystemT>);
        long offset = GetIndexOffset(thread_info);
        size_t ord=table.GetNumStages();
        size_t dimension=GetDimension();
        
        FixedSizeRealBuffer intermediate_step[10];
        bool needs_free=false;
        real_t lmem[32]={0};
        real_t* tmp_buffer=&lmem[0];
        if((ord+1)*dimension>=32){
            needs_free=true;
            tmp_buffer = AllocateTemporaryBuffer<real_t>( (ord+1)*dimension , thread_info.GetIndex());
        }
            
        
        FixedSizeRealBuffer wsum( &(tmp_buffer[ord*dimension]) , dimension);//this stores the intermediate offset
    
        for(size_t step_i=0; step_i<ord; ++step_i){
            
            real_t* tmp_buffer_i = &tmp_buffer[step_i*dimension];
            FixedSizeRealBuffer& k_i = intermediate_step[step_i];
            k_i = FixedSizeRealBuffer( tmp_buffer_i, dimension);
            wsum.CopyFrom(in_out_state, offset, offset+static_cast<long>(dimension), 0);
            for(size_t prev_j=0; prev_j<step_i; ++prev_j){                    
                real_t aij = table.GetXCoeff(step_i, prev_j);
                const FixedSizeRealBuffer& k_j = intermediate_step[prev_j];
                wsum.AccumulateWeighted( (aij * step_size) , k_j );
            }//for prev_j
            static_cast<SystemT*>(this)->Fn(k_i, wsum, t + step_i * step_size * table.GetTCoeff(step_i), thread_info);
        }//for step_i
        
        
        for(size_t step_i=0; step_i<ord; ++step_i){
            const FixedSizeRealBuffer& k_i = intermediate_step[step_i];
            double bij = table.GetYCoeff(step_i);
            in_out_state.AccumulateWeighted( (bij * step_size) , k_i, 0, dimension, offset);
        }
        if(needs_free){
            FreeTemporaryBuffer(tmp_buffer, thread_info.GetIndex());
        }
        
        
    }//method
    real_t HOSTDEVICE IntegrateRKExplicitAdaptive ( FixedSizeRealBuffer& in_out_state, 
                                    real_t t,
                                    real_t step_size, 
                                    const ThreadInfo& thread_info,
                                    const RKButcherTableau& table_lo = RKButcherTableau::RKDormandPrince<false>(),
                                    const RKButcherTableau& table_hi = RKButcherTableau::RKDormandPrince<true>(),
                                    bool continue_from_hi=true
                                  )
    {//explicit runge-kutta method of arbitrary order and coefficients.
        const size_t dimension=GetDimension();
        const long offset = GetIndexOffset(thread_info);
        bool needs_free=false;
        real_t lmem[16]={0};
        real_t* tmp_buffer=&lmem[0];
        if(2*dimension>=16){
            needs_free=true;
            tmp_buffer = AllocateTemporaryBuffer<real_t>( 2*dimension , thread_info.GetIndex());
        }
        FixedSizeRealBuffer s1(&tmp_buffer[0], dimension);
        FixedSizeRealBuffer s2(&tmp_buffer[dimension], dimension);
        s1.CopyFrom(in_out_state, offset, offset+static_cast<long>(dimension), 0);
        s2.CopyFrom(in_out_state, offset, offset+static_cast<long>(dimension), 0);
        IntegrateRKExplicitFixed ( s1, t, step_size, thread_info, table_lo );
        IntegrateRKExplicitFixed ( s2, t, step_size, thread_info, table_hi );
        real_t normsq=0.0;
        for(size_t i=0; i<dimension; ++i){
            real_t d=s1.Get(i)-s2.Get(i);
            normsq+=d*d;
        }
        if(continue_from_hi){
            in_out_state.CopyFrom(s2, 0, static_cast<long>(dimension), offset);
        }
        else{
            in_out_state.CopyFrom(s1, 0, static_cast<long>(dimension), offset);
        }
        if(needs_free){
            FreeTemporaryBuffer(tmp_buffer, thread_info.GetIndex());
        }
        return normsq;
    }//method
};//class
        
