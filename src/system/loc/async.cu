#include <cmath>
#include <array>
#include <iostream>
#include <random>

#include <cuda.h>
#include <curand_kernel.h>
#include <sm_60_atomic_functions.h>

#include "Integrator.h"
#include "Config.h"


#include "../../../src/FixedSizeBuffer.cpp"//FIXME: 


    
#if defined(__CUDA__)
	DEVICE unsigned Id(){
	    return blockIdx.x*blockDim.x+threadIdx.x;
	}
#endif


DEVICE real_t random_normal(curandState_t* seed, real_t mu=0.0, real_t sigma=1.0){
    return curand_normal (seed)*sigma+mu;
}


static constexpr int Nb=1;

	
    constexpr float grav=9.81;
 
	class Async : public System<Async> {
	friend class System<Async>;
	private:
	int Np;
	real_t m,L,fp,bmin, u, step_next;
	
	public:
	DEVICE Async (unsigned _Np, curandState_t* seed){
            
                Np      =  _Np;
                
                
                m       =  random_normal(seed, 76.976, 16.1631);   // Mass
                L       =  random_normal(seed, 1.167, 0.092);    // Inverted pendulum length
                fp      =  random_normal(seed, 0.6,   0.1);        // Walking frequency
                bmin    =  random_normal(seed, 0.0157, 0.007);   // Margin of stability
                
                
                
                real_t Omegap = std::sqrt(grav/L);
            
                u=bmin/(1.0-std::tanh(Omegap*0.25/fp));
                step_next=0.5/fp;
                bmin=-bmin;
            
                
	}
        
    DEVICE void Init(FixedSizeRealBuffer S){
            int id=Id();
            
            S.Set(0,0);
	    real_t b=static_cast<real_t>(id%(Nb+Np)>=Nb);
            real_t Omegap = std::sqrt(grav/L);
            S.Set(0,u*Omegap*std::tanh(Omegap*0.25/fp)*b);
            S.Set(1,0);
    }
	private:
	    
	    DEVICE size_t GetIndexOffset(const ThreadInfo& thread_info) const{
            return 0;
	    }
        DEVICE size_t GetDimension() const{
            return 2;
	    }
	    
	    HOST void FnHost(FixedSizeRealBuffer& out, const FixedSizeRealBuffer& in, real_t t, const ThreadInfo& thread_info) const{
            throw;
        }
	    DEVICE void FnGPU(FixedSizeRealBuffer& out, const FixedSizeRealBuffer& in, real_t t, const ThreadInfo& thread_info) const{
            int i=thread_info.GetIndex();
            int id=i;
            SHARED real_t bridge_accel[1];
            SHARED real_t H[1];
            H[0]=0.0;
            
            bridge_accel[0]=0.0;
            __syncthreads();
            
            real_t Omegap = std::sqrt(grav/L);
            real_t F = Omegap*Omegap*(u-in.Get(1));
            
            
            
            
            
                
            
                
            
            
            
            __syncthreads();    
            atomicAdd(&H[0], F*m*(id%(Nb+Np)>=Nb));
            
            __syncthreads();
            
            
            constexpr real_t M[2]       =  {1150.0*1000.0, 1150.0*1000.0};   
            constexpr real_t fb[2]      =  {0.524, 0.746};   
            constexpr real_t zeta[2] = {0.58/100.0, 0.68/100.0};
    
            int j=id%(Nb);
            
            real_t omegab = 2*M_PI*fb[j];
            real_t a=(id<Nb)*(H[0]/M[j] - 2*omegab*zeta[j]*in.Get(0) - omegab*omegab*in.Get(1));
            __syncthreads();
            atomicAdd(&bridge_accel[0],a);
            
            
            
            __syncthreads();
            
            
            a+=((id>=Nb))*(-F-bridge_accel[0]);
            
            out.Set(0, a);
            
            out.Set(1, in.Get(0));
            
	    }
    public:
	    DEVICE real_t GetTimeUntilNextStep(){
            int i=Id()%(Nb+Np)-Nb;
            real_t tau=step_next;
            return tau+100000000*(i<0);
        }
        DEVICE void OnStep( curandState_t* seed, FixedSizeRealBuffer state, real_t T){
            
            real_t Omegap = std::sqrt(grav/L);
            if(std::abs(GetTimeUntilNextStep()-T)<1e-10){
	            step_next=T+0.5/fp*random_normal(seed,1.0,0.05);
        	    u=state.Get(1)+state.Get(0)/Omegap+bmin*(random_normal(seed,1.0,0.2));
	            bmin=-bmin;
	    }
        }
	};

    HOSTDEVICE static inline uint32_t FloatFlip(uint32_t f)
    {
        uint32_t mask = -int32_t(f >> 31) | 0x80000000;
        return f ^ mask;
    }

    HOSTDEVICE static inline uint32_t IFloatFlip(uint32_t f)
    {
        uint32_t mask = ((f >> 31) - 1) | 0x80000000;
        return f ^ mask;
    }
    
	__global__ void AsyncSweep(int* random_seed, int Np,  real_t* traces, real_t dt){
		
	    
	    int o=Id();
	    if(o>=(Nb+Np)){
	            return;
        }
        	
        __shared__ uint32_t min_next_step[1];
#ifdef ADAPTIVE
	    real_t threshold=1e-7;
        __shared__ uint32_t max_err[1];
#endif
	    real_t state[2];
	    FixedSizeRealBuffer S(state,2);
	    curandState_t seed_state;
	    curand_init (random_seed[o], Id(), 0, &seed_state);
	    Async tn=Async(Np, &seed_state);
	    tn.Init(S);
	    real_t T=0.0;
	    real_t h=0.001, hmax=std::min(dt/2.0,0.01);
	    real_t hmin=0.0001;
	    unsigned t=0;
        while(T<180.0){
                    
                    
                    min_next_step[0]=std::numeric_limits<uint32_t>::max();
                    __syncthreads();
                    real_t tnext=tn.GetTimeUntilNextStep();
                    union{
                            uint32_t u;
                            float f;
                        } cvt;
                    cvt.f=static_cast<float>(tnext);
                        
                    __syncthreads();
                    atomicMin(&min_next_step[0], FloatFlip(cvt.u));
                    __syncthreads();
                    cvt.u=IFloatFlip(min_next_step[0]);
                    float min_next=cvt.f;
                    while(T<min_next-h){
                    	if(o<Nb+5){
                        	if(T>=t*dt){
                            
                           	 traces[o+t*(5+Nb)]=S.Get(1);
			   	 ++t;
                        	}

                    	}
#ifdef ADAPTIVE
                                while(true){
                                    max_err[0]=0;
                                    __syncthreads();
                                    real_t S0[2]={S.Get(0), S.Get(1)};
                                    cvt.f =tn.Solve(S, T, h);
                                    __syncthreads();
                                    atomicMax(&max_err[0], FloatFlip(cvt.u));
                                    __syncthreads();
                                    
                                    cvt.u=IFloatFlip(max_err[0]);
                                    real_t err=cvt.f;
                                    if(err>0.0){
                                        h*=pow(threshold/(err),0.2);
                                        
                                    	h=std::max(std::min(h,hmax),hmin);
                                    }

                                    if(err<threshold || h==hmin){
                                        break;
                                    }
                                    state[0]=S0[0];
                                    state[1]=S0[1];
                                }
#else
                                tn.Solve(S,T,h,false);
     
#endif
                           T+=h;
                           
                    }
                    tn.Solve(S,T,min_next-T);
                    T=min_next;
                    tn.OnStep(&seed_state,S,T);
                    tn.Solve(S,T,h);
                    T=T+h;
                    
                        
        }
                
}
        




void Check(){

    cudaError_t error = cudaGetLastError();
    if(error!=cudaSuccess)
    {
        fprintf(stderr,"ERROR: %s\n", cudaGetErrorString(error) );
        exit(-1);
    }
}

void Finish(){
    Check();
    cudaDeviceSynchronize();
    Check();
}

HOST int main(int argc, char** argv){

    real_t dt=0.01;
    int Np=atoi(argv[1]);
    int N=Np+Nb;
    int bsize=int(N/32)+((N%32)!=0);
    bsize*=32;
    
    printf("%i \n", bsize);
    int bcount=std::max(1,(Np+Nb)/bsize);
    int* random_seed; 
    real_t* traces;
    cudaMallocManaged(&traces, (Nb+5)*int(180/dt)*sizeof(real_t));
    Check();
    assert(traces!=nullptr);
    cudaMallocManaged(&random_seed, (Nb+Np)*sizeof(int));
    Check();
    assert(random_seed!=nullptr);
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_int_distribution<int> uniform_dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
    for(int i=0; i<Np+Nb; ++i){
        random_seed[i]=uniform_dist(e1);
    }
    AsyncSweep<<<bcount,bsize>>>(random_seed, Np, traces, dt); Finish();
    for(int t=0; t<int(180.0/dt); ++t){
        std::cout << t*dt<< " ";
        for(int n=0; n<Nb+5; ++n){
           std::cout<<  traces[t*(Nb+5)+n] << " ";
        }
       	std::cout << "\n";
    }
    cudaFree(traces);
    cudaFree(random_seed);
        
    return 0;
}
