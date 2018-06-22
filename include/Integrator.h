#pragma once
#include "Allocate.h"
#include "FixedSizeBuffer.h"
#include "ThreadInfo.h"
#include <type_traits>
#include <cassert>

class RKButcherTableau {
private:
    const size_t s;
    const std::array<double,256> A;
    const std::array<double,16> B;
    const std::array<double,16> C;
    RKButcherTableau()=delete;
    RKButcherTableau(size_t p_s, 
                     std::array<double, 256> p_A,
                     std::array<double, 16> p_B,
                     std::array<double, 16> p_C):
                    s(p_s),
                    A(p_A), B(p_B), C(p_C)
                     
    {
        
    }
public:
    static const RKButcherTableau& RK4Classic() {//TODO: add more rk methods including adaptive ones (ode45 and etc.)
        static const RKButcherTableau rk4(4,
                   {0,   0, 0,   0,
                    0.5, 0, 0,   0,
                    0,   0.5, 0, 0,
                    0,   0,   1.0, 0 } ,
                    
                    {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0} ,
                   
                    {0, 0.5, 0.5, 0}
                );
            return rk4;
    }   
    inline size_t GetNumStages() const{
        assert(s<16);
        return s;
    }
    inline double GetXCoeff(size_t step_i, size_t prev_j) const{ 
        assert(step_i<GetNumStages() && prev_j<GetNumStages());
        return A[step_i * s + prev_j]; 
    }
    inline double GetYCoeff(size_t step) const{ 
        assert(step<GetNumStages());
        return B[step];
    }
    inline double GetTCoeff(size_t step) const {
        assert(step<GetNumStages());
        return C[step];
    }
};

template<typename SystemT> class System {
    //static polymorphism so all types that are used as arguments to template need to derive from it. 
    //that way we can use the below virtual function call without actually dereferencing anything.
protected:
    virtual void Fn( FixedSizeDoubleBuffer& out, const FixedSizeDoubleBuffer& in, double t) const=0;
    virtual void Integrate(FixedSizeDoubleBuffer& in_out_state, 
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
    void Solve(FixedSizeDoubleBuffer& in_out_state,
               double t,
               double step_size = 0.001, 
               const ThreadInfo& thread_info        =   ThreadInfo::Serial()
         ){
        static_cast<SystemT*>(this)->Integrate(in_out_state, t, step_size, thread_info);
    }
private:
    void IntegrateRKExplicitFixed ( FixedSizeDoubleBuffer& in_out_state, 
                                    double t,
                                    double step_size, 
                                    const ThreadInfo& thread_info,
                                    const RKButcherTableau& table = RKButcherTableau::RK4Classic()
                                  )
    {//explicit runge-kutta method of arbitrary order and coefficients.
        
    static_assert(std::is_base_of_v<System<SystemT>, SystemT>);
        size_t ord=table.GetNumStages();
        size_t dimension=in_out_state.GetLength();
        
        FixedSizeDoubleBuffer intermediate_step[ord];
        bool needs_alloc=true;
        size_t sz=(Config::MaxStackBufferSize/sizeof(double))-1;
        int loop_c=0;
        do{//this is kind of a tricky construction but here is what it does:
                    //it tries to allocate a vla on the stack with the maximum allowed size
                    //if that size isn't large enough, it frees (almost all of) that array by going out of scope,
                    //setting the allocation variable to size 1, and reallocating that array with the new size.
                    //then it overwrites the resulting pointer with the pointer to a heap-allocated variable.
            double stack_pool[sz+1];
            double* tmp_buffer=&stack_pool[0];
            if(__builtin_expect (sz < (ord+1)*dimension, false) )
            {
                if(sz){
                    sz=0;
                    continue;
                } else {
                    tmp_buffer = AllocateTemporaryBuffer<double>( (ord+1)*dimension , thread_info.GetIndex());
                    needs_alloc=false;
                }
            }
            FixedSizeDoubleBuffer wsum( &(tmp_buffer[ord*dimension]) , dimension);//this stores the intermediate offset
            
            for(size_t step_i=0; step_i<ord; ++step_i){
                
                double* tmp_buffer_i = &tmp_buffer[step_i*dimension];
                FixedSizeDoubleBuffer& k_i = intermediate_step[step_i];
                k_i = FixedSizeDoubleBuffer( tmp_buffer_i, dimension );
                wsum.CopyFrom(in_out_state);
                for(size_t prev_j=0; prev_j<step_i; ++prev_j){                    
                    double aij = table.GetXCoeff(step_i, prev_j);
                    const FixedSizeDoubleBuffer& k_j = intermediate_step[prev_j];
                    wsum.AccumulateWeighted( (aij * step_size) , k_j );
                }//for prev_j
                static_cast<SystemT*>(this)->Fn(k_i, wsum, t + step_i * step_size * table.GetTCoeff(step_i));
            }//for step_i
            
            
            for(size_t step_i=0; step_i<ord; ++step_i){
                const FixedSizeDoubleBuffer& k_i = intermediate_step[step_i];
                double bij = table.GetYCoeff(step_i);
                in_out_state.AccumulateWeighted( (bij * step_size) , k_i);
            }
            if(__builtin_expect(tmp_buffer!=&stack_pool[0],false)){
                FreeTemporaryBuffer(tmp_buffer, thread_info.GetIndex());
                break;
            }
            loop_c++;
        }  while(__builtin_expect(needs_alloc,false) && loop_c<2);  //we optimize the branch prediction for the case where only stack allocation must occur,
                                                        //because heap allocation is way slower anyway (especially in the worst case, since it might call malloc() from a thread).  
                                                        //also, this branch prediction works because the compiler
                                                        //knows the body of a do/while statement must be entered at least once.
        
    }//method
};//class
        
