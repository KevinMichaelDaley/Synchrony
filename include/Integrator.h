#pragma once
#include "AllocateGPU.h"
#include "FixedSizeBuffer.h"
#include "ThreadInfo.h"
#include <type_traits>
#include <cassert>
template<size_t N> class DoubleArray{
private:
    double arr[N];
public:
    HOSTDEVICE DoubleArray(double a[]){
        for(size_t i=0; i<N; ++i){
            arr[i]=a[i];
        }
    }
    HOSTDEVICE DoubleArray(double a[], size_t subN){
        for(size_t i=0; i<subN; ++i){
            arr[i]=a[i];
        }
    }
    HOSTDEVICE double operator[](int i) const{
        return arr[i];
    }
};

class RKButcherTableau {
private:
    const size_t s;
    const DoubleArray<81> A;
    const DoubleArray<9> B;
    const DoubleArray<9> C;
    HOSTDEVICE RKButcherTableau()=delete;
    HOSTDEVICE RKButcherTableau(size_t p_s, 
                     DoubleArray<81> p_A,
                     DoubleArray<9> p_B,
                     DoubleArray<9> p_C):
                    s(p_s),
                    A(p_A), B(p_B), C(p_C)
                     
    {
        
    }
public:
    HOSTDEVICE static RKButcherTableau RK4Classic() {//TODO: add more rk methods including adaptive ones (ode45 and etc.)
        double pA[16]={0,   0, 0,   0,
                    0.5, 0, 0,   0,
                    0,   0.5, 0, 0,
                    0,   0,   1.0, 0 };
        double pB[9]={1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0};
        double pC[9]={0,0.5,0.5,0,0};
        const RKButcherTableau rk4(4,
                DoubleArray<81>(pA,16),
                DoubleArray<9>(pB),
                DoubleArray<9>(pC)
                );
            return rk4;
    }   
    inline HOSTDEVICE size_t GetNumStages() const{
        assert(s<=9);
        return s;
    }
    inline HOSTDEVICE double GetXCoeff(size_t step_i, size_t prev_j) const{ 
        assert(step_i<GetNumStages() && prev_j<GetNumStages());
        return A[step_i * s + prev_j]; 
    }
    inline HOSTDEVICE double GetYCoeff(size_t step) const{ 
        assert(step<GetNumStages());
        return B[step];
    }
    inline HOSTDEVICE double GetTCoeff(size_t step) const {
        assert(step<GetNumStages());
        return C[step];
    }
};

template<typename SystemT> class System {
    //static polymorphism so all types that are used as arguments to template need to derive from it. 
    //that way we can use the below virtual function call without actually dereferencing anything.
protected:
    HOSTDEVICE void Fn( FixedSizeDoubleBuffer& out, const FixedSizeDoubleBuffer& in, double t, const ThreadInfo& thread_info) const{
        static_cast<const SystemT*>(this)->Fn( out, in, t, thread_info );
    }
    HOSTDEVICE size_t GetDimension() const{
        return static_cast<const SystemT*>(this)->GetDimension();
    }
    HOSTDEVICE size_t GetIndexOffset(const ThreadInfo& thread_info) const{
           return static_cast<const SystemT*>(this)->GetIndexOffset(thread_info);
    }
    HOSTDEVICE void Integrate(FixedSizeDoubleBuffer& in_out_state, 
                           double t,
                           double step_size, 
                           const ThreadInfo& thread_info){//override this virtual method with the default integration routine you want to use (if not rk4).
                                                          //for discrete maps this can just be Fn(in_out_state, in_out_state, k)
                                                          
        IntegrateRKExplicitFixed(in_out_state, t, step_size, thread_info);
    }
    //add new integration routines to this template so they can be reused.
    //parameters fixed over the whole trajectory are set in the constructor/whatever get-setter methods of the derived class
    //state should always be modified in-place; the caller is responsible for storing previous values and the associated allocations,
    //since they may not be necessary and to store two state arrays here when we might not need to is terribly inefficient for large dimensions.
    
public:
    HOST void Solve(FixedSizeDoubleBuffer& in_out_state,
               double t,
               double step_size = 0.001, 
               const ThreadInfo& thread_info        =   ThreadInfo::Serial()
         )
    {
        static_cast<SystemT*>(this)->Integrate(in_out_state, t, step_size, thread_info);
    }
#if defined(__CUDA__)
    DEVICE void Solve(FixedSizeDoubleBuffer& in_out_state,
                      double t,
                      double step_size = 0.001){
        size_t global_id = blockIdx.x *blockDim.x + threadIdx.x;
        static_cast<SystemT*>(this)->Integrate(in_out_state, t, step_size, ThreadInfo(global_id));
    }
#endif
                      
private:
    void HOSTDEVICE IntegrateRKExplicitFixed ( FixedSizeDoubleBuffer& in_out_state, 
                                    double t,
                                    double step_size, 
                                    const ThreadInfo& thread_info,
                                    const RKButcherTableau& table = RKButcherTableau::RK4Classic()
                                  )
    {//explicit runge-kutta method of arbitrary order and coefficients.
        
    static_assert(std::is_base_of_v<System<SystemT>, SystemT>);
        long offset = GetIndexOffset(thread_info);
        size_t ord=table.GetNumStages();
        size_t dimension=GetDimension();
        
        FixedSizeDoubleBuffer intermediate_step[256];
        
        double* const tmp_buffer = AllocateTemporaryBuffer<double>( (ord+1)*dimension , thread_info.GetIndex());
                
        
        FixedSizeDoubleBuffer wsum( &(tmp_buffer[ord*dimension]) , dimension);//this stores the intermediate offset
    
        for(size_t step_i=0; step_i<ord; ++step_i){
            
            double* tmp_buffer_i = &tmp_buffer[step_i*dimension];
            FixedSizeDoubleBuffer& k_i = intermediate_step[step_i];
            k_i = FixedSizeDoubleBuffer( tmp_buffer_i, dimension);
            wsum.CopyFrom(in_out_state, offset, offset+static_cast<long>(dimension), 0);
            for(size_t prev_j=0; prev_j<step_i; ++prev_j){                    
                double aij = table.GetXCoeff(step_i, prev_j);
                const FixedSizeDoubleBuffer& k_j = intermediate_step[prev_j];
                wsum.AccumulateWeighted( (aij * step_size) , k_j );
            }//for prev_j
            static_cast<SystemT*>(this)->Fn(k_i, wsum, t + step_i * step_size * table.GetTCoeff(step_i), thread_info);
        }//for step_i
        
        
        for(size_t step_i=0; step_i<ord; ++step_i){
            const FixedSizeDoubleBuffer& k_i = intermediate_step[step_i];
            double bij = table.GetYCoeff(step_i);
            in_out_state.AccumulateWeighted( (bij * step_size) , k_i, 0, dimension, offset);
        }
        
        FreeTemporaryBuffer(tmp_buffer, thread_info.GetIndex());
        
        
    }//method
};//class
        
