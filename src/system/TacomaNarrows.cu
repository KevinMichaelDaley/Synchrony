#include <cmath>
#include <array>

#include "Integrator.h"
#include "Config.h"
#include <cuda.h>
#include "../../src/FixedSizeBuffer.cpp"//FIXME: 
inline HOSTDEVICE real_t sgn(real_t n){
	return (n > 0.0) - (n < 0.0);
}
struct TacomaParameters{

    real_t bridge_damping;
    real_t support_damping;
    real_t support_intrinsic_frequency;
    real_t bridge_intrinsic_frequency;
    real_t bridge_forcing_frequency;
    real_t bridge_forcing_amplitude;
    real_t vortex_shedding_strength;
    real_t support_to_bridge_coupling;
};
HOSTDEVICE real_t characteristic(TacomaParameters params,real_t x){//characteristic equation to find the normal modes of the linear coupled system
    real_t lambda=params.support_damping;
    real_t w=params.support_intrinsic_frequency;
    real_t mu=params.support_to_bridge_coupling;
    real_t A=params.bridge_forcing_amplitude;
    real_t h=params.bridge_damping;
    real_t Omega=params.bridge_intrinsic_frequency;
    
    return (x*x+lambda*x+w*w)*(x*x+h*x+Omega*Omega)-2*mu*mu*x*x*x*x;
}
HOSTDEVICE real_t d_characteristic(TacomaParameters params,real_t x){//analytic derivative of the above for use with newton's method
    real_t lambda=params.support_damping;
    real_t w=params.support_intrinsic_frequency;
    real_t mu=params.support_to_bridge_coupling;
    real_t A=params.bridge_forcing_amplitude;
    real_t h=params.bridge_damping;
    real_t Omega=params.bridge_intrinsic_frequency;
    
    real_t derivative_product=   (2*x+lambda) * (x*x+h*x+Omega*Omega)//x'y
                                +(x*x+lambda*x+w*w) * (2*x+h);//xy'
    return derivative_product-8*mu*mu   *x*x*x;//power rule
}
HOSTDEVICE real_t find_zero(TacomaParameters params){
    real_t b=params.support_intrinsic_frequency;
    real_t a=params.bridge_intrinsic_frequency;
    if(b>a){
	    real_t tmp=b;
	    b=a;
	    a=tmp;
    }
    a+=0.25;
    b-=0.25;
    if(b<0){
	    b=0;
    }
    //bisection method as a first pass
    for(int n=0; n<20; ++n){
        if(std::abs(b-a)<0.01){
            break;
        }
        real_t fa=characteristic(params,a);
        real_t fb=characteristic(params,b);
        real_t fc=characteristic(params, (b+a)/2.0);
        if(std::abs(fc)<0.0001){
            break;
        }
        if(sgn(fa)==sgn(fc)){
            a=(a+b)/2.0;
        } else{
            b=(b+a)/2.0;
        }
        
    }
    real_t x=(b+a)/2.0;
    real_t threshold=0.01;
    //newton's method with early exit on success
    for(int i=0; i<100; ++i){
        real_t f2=characteristic(params,x);
        if(std::abs(f2)<threshold){
            return x;
        }
        real_t fp2=d_characteristic(params,x);
        x-=f2/fp2;
    }
    return (b+a)/2.0;
}
        
        
HOSTDEVICE TacomaParameters make_default_params(real_t bridge_intrinsic_frequency, 
                                                real_t forcing_frequency,
                                                real_t bridge_forcing_amplitude){
                    TacomaParameters p= {0.05,
                            0.1,
                            0.97,
                            bridge_intrinsic_frequency,
                            forcing_frequency,
                            bridge_forcing_amplitude,
                            1,
                            0.01
                    };
                    return p;
}

HOSTDEVICE void support_system_derivatives(TacomaParameters params, 
                                RealArray<2> support_vel,
                                RealArray<2> support_pos,
                                real_t bridge_accel,
                                RealArray<2>* support_out_accel,
                                bool linear,
                                real_t mismatch,
                                real_t eig
                               )
{
    real_t lambda=params.support_damping;
    real_t w=(!linear)? params.support_intrinsic_frequency : sqrt(eig);
    real_t mu=params.support_to_bridge_coupling;
    real_t gamma=(!linear) ? params.vortex_shedding_strength : 0.0;
    for(int i=0; i<2; ++i){
        real_t dx=support_vel[i];
        real_t x=support_pos[i];
        real_t ddv=bridge_accel;
        w+=mismatch*i;
        ((*(support_out_accel)))[i]=-lambda*dx-w*w*x-mu*ddv+gamma*sgn(dx);
    }
}

HOSTDEVICE real_t bridge_accel(TacomaParameters params,
                    RealArray<2> support_vel,
                    RealArray<2> support_pos,
                    real_t bridge_vel,
                    real_t bridge_pos,
                    real_t t,
                    bool linear,
                    real_t mismatch,
                    real_t eig
                   )
{
    real_t v=bridge_pos;
    real_t dv=bridge_vel;
    RealArray<2> ddx({0,0});
    support_system_derivatives(params,support_vel,support_pos,0,&ddx, linear,mismatch,eig);
    real_t sum=ddx[0]+ddx[1];
    real_t A=params.bridge_forcing_amplitude;
    real_t h=params.bridge_damping;
    real_t Omega=!linear ? params.bridge_intrinsic_frequency : sqrt(eig);
    real_t gamma=(!linear) ? params.vortex_shedding_strength : 0;
    real_t mu=params.support_to_bridge_coupling;
    real_t r=1.0;
    return ((-2*h*dv-Omega*Omega*v-sum)+gamma*sgn(dv)+A*sin(t*params.bridge_forcing_frequency))/(1-mu*2*r);
}

class ParamRange{
private:
    const real_t s,e;
public:
    HOSTDEVICE ParamRange(real_t p_s, real_t p_e): s(p_s), e(p_e){}
    HOSTDEVICE ParamRange(real_t p_v): s(p_v), e(p_v){}
    HOSTDEVICE real_t Get(int i, int N) const{
        return s+(e-s)*i/real_t(N);
    }
};

struct AmplitudeState{
    real_t max, min;
    HOSTDEVICE AmplitudeState(): max(-99999999), min(999999999){};
    HOSTDEVICE void Accum(real_t pos){
        max=std::max(pos,max);
        min=std::min(pos,max);
    }
    HOSTDEVICE real_t Get() const{
        return max-min;
    }
};

struct PhaseState{
    int zero_crossing0;
    int delta_avg;
	    HOSTDEVICE PhaseState(): zero_crossing0(-1), delta_avg(0)
	    {
	    }
	    HOSTDEVICE void Accum(int t, real_t x, real_t x0, real_t y, real_t y0){
		if(x>=0 && x0<0){
		    zero_crossing0=t;
		}
		if(y>=0 && y0<0 && zero_crossing0>0){
		    delta_avg+=t-zero_crossing0;
		}
	    }
	    HOSTDEVICE real_t Get() const{
		return delta_avg;
	    }
	};
#if defined(__CUDA__)
	DEVICE unsigned Id(){
	    return blockIdx.x*blockDim.x+threadIdx.x;
	}
#endif

	HOST unsigned Id(){
	    return 0;
	}


	class TacomaNarrows : public System<TacomaNarrows> {
	friend class System<TacomaNarrows>;
	public:
	HOSTDEVICE TacomaNarrows(){
	}
	HOSTDEVICE TacomaNarrows (unsigned N, const ParamRange& amp_forcing, const ParamRange& freq_forcing, const ParamRange& freq_intrin, bool is_linear, const ParamRange& mismatch):
		_params( make_default_params(freq_intrin.Get(Id()%N,N), freq_forcing.Get(Id()%N,N), amp_forcing.Get(Id()/N,N)) ),
		_is_linear(is_linear),
		_eigenvalue(is_linear?find_zero(_params):0),
		_mismatch(mismatch.Get(Id()/N, N))
	    {}
	public:
	    HOSTDEVICE TacomaParameters GetParams() const{
		return _params;
	    }
	    HOSTDEVICE real_t GetMismatch() const{
		return _mismatch;
	    }
	    HOSTDEVICE real_t GetEigenvalue() const{
		return _eigenvalue;
	    }
	private:
	    TacomaParameters _params;
	    bool _is_linear;
	    real_t _eigenvalue;
	    real_t _mismatch;
	    
	    HOSTDEVICE size_t GetIndexOffset(const ThreadInfo& thread_info) const{
		return thread_info.GetIndex()*6;
	    }
	    HOSTDEVICE size_t GetDimension() const{
		return 6;
	    }
	    HOSTDEVICE void Fn(FixedSizeRealBuffer& out, FixedSizeRealBuffer& in, real_t t, const ThreadInfo& thread_info) const{
		int param_index=thread_info.GetIndex();
		
		real_t eig=_eigenvalue;
		real_t mismatch=_mismatch;
		bool linear=_is_linear;
		TacomaParameters params=_params;
		real_t svel[2]={in.Get(0),in.Get(2)};
		real_t spos[2]={in.Get(1),in.Get(3)};
		real_t zero[2]={0,0};
		RealArray<2> support_vel(svel);
		RealArray<2> support_pos(spos);
		real_t bridge_vel=in.Get(4);
		real_t bridge_pos=in.Get(5);
		RealArray<2> support_new_accel(zero);
		real_t bridge_new_accel=bridge_accel(params,support_vel,support_pos,bridge_vel,bridge_pos,t,linear, mismatch,eig);
		support_system_derivatives(params, support_vel, support_pos, bridge_new_accel, &support_new_accel,linear, mismatch,eig);
		out.Set(0,support_new_accel[0]);
		out.Set(1,support_vel[0]);
		out.Set(2,support_new_accel[1]);
		out.Set(3,support_vel[1]);
		out.Set(4,bridge_new_accel);
		out.Set(5,bridge_vel);
	    }
	};


	constexpr real_t h=0.001;
	constexpr int tmax=int(1000.0/h);
	__global__ void TacomaNarrowsSweep(real_t* p_Ay, real_t* p_phi, real_t* state, bool sweep_forcing_freq, bool sweep_intrin_freq, bool sweep_amplitude=true, bool is_linear=false, bool mismatch=false){
		
	    FixedSizeRealBuffer S;
	    TacomaNarrows tn;
	    {
		    int N=gridDim.x*blockDim.x;
		    S=FixedSizeRealBuffer(state,N*6);
		    tn=TacomaNarrows(std::sqrt(N), sweep_amplitude? ParamRange(0,100) : ParamRange(10),
				     sweep_forcing_freq?  ParamRange(0.4,1.4) : ParamRange(0.65),
				     sweep_intrin_freq?  ParamRange(0.4, 1.4) : ParamRange(0.97),
				     is_linear, mismatch? ParamRange(0,1):ParamRange(0));
	    }
	    int o=Id()*6;
	    PhaseState phi;
	    AmplitudeState Ay;

	    
	    for(int t=0; t<tmax; ++t){
		real_t T=t*h;
		real_t s0[2];
		for(int i=0; i<2; ++i){
		    s0[i]=S.Get(i*2+1+o);
		}
		tn.Solve(S, T, h);
		if(t>15/(2.0*h)){
		    phi.Accum(t, S.Get(3+o), s0[1], S.Get(1+o), s0[0]);
        }
        Ay.Accum(S.Get(5+o));
    }
    p_Ay[o/6]=Ay.Get();
    p_phi[o/6]=phi.Get();
    {
    state[o]=tn.GetParams().bridge_intrinsic_frequency/tn.GetParams().bridge_forcing_frequency;
    state[o+1]=tn.GetParams().bridge_forcing_amplitude;
    state[o+3]=tn.GetMismatch();
    }
}

#include <fstream>
void init_default(real_t* state, int N){
    
    for(int i=0; i<N*N; ++i){
        state[i*6]=0.03;
        state[i*6+1]=0.03;
        state[i*6+2]=0.04;
        state[i*6+3]=0.04;
        state[i*6+4]=0.0;
        state[i*6+5]=-0.0;
    }
}
void init_basin(real_t* state, int N){
    
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            int n=i*N+j;
            state[n*6]=i/real_t(N)*0.1-0.05;
            state[n*6+1]=i/real_t(N)*0.1-0.05;
            state[n*6+2]=j/real_t(N)*0.1-0.05;
            state[n*6+3]=j/real_t(N)*0.1-0.05;
            state[n*6+4]=0.0;
            state[n*6+5]=0.0;
        }
    }
}

void write_sweep(std::ofstream& ofs, real_t* Ay, real_t* phi, real_t* state, int N=32){
    ofs << "# i j x1 x0  Omega As Omega0 mismatch eigv dx       Ay  dphi\n";
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            real_t n=real_t(N);
            ofs << (j/n) << " " << (i/n)<< "    0 0     ";
            for(int k=0; k<6; ++k){
                ofs << state[(i*N+j)*6+k]<<" ";
            }
            ofs<< "         " << Ay[i*N+j] << "    " << phi[i*N+j] << "\n";
        }
    }
}

void write_basin(std::ofstream& ofs, real_t* Ay, real_t* phi, real_t* state, int N=32){
    ofs << "# i j x1 x0   forcing_freq forcing_ampl intrinsic_freq mismatch eigv dx       Ay  dphi\n";
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            real_t n=real_t(N);
            ofs << (j/n) << " " << (i/n)<< "    "<<i/n-0.5<< " " << j/n-0.5 << "    ";
            for(int k=0; k<6; ++k){
                ofs << state[(i*N+j)*6+k]<<" ";
            }
            ofs<< "         " << Ay[i*N+j] << "    " << phi[i*N+j] << "\n";
        }
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
/*
HOST void Hysteresis(init_x, init_y){
	real_t s=0.0;//forcing frequency sweep offset
	//always continue from the last state even when changing parameters
	double state[6]={init_x, init_x, init_y, init_y, 0, 0};
	double state0[6];
	FixedSizeDoubleBuffer S(&state[0], 6);
	while(true){//sweep over the parameter
	     TacomaNarrows tn(1, ParamRange(20), ParamRange(s+0.4), ParamRange(0.97), false, ParamRange(0));
	     PhaseState phi;
	     for(int t=0; t<25/0.01; ++t){
		     for(int i=0; i<6; ++i){
			     state0[i]=state[i];
		     }
		     tn.Solve(S, t*0.01, 0.01);
		     if(t>=1250){
			     phi.Accum(t, S[3], state0[3], S[1], state0[1]);
			}
		     
		}
			    
		     
	}
	
}
*/
HOST int main(int argc, char** argv){
    real_t* Ay;
    real_t* phi;
    real_t* state;
    int N=96;
    const int bsize=768;
    int bcount=N*N/bsize;
    cudaMallocManaged(&Ay, N*N*sizeof(real_t));
    Check();
    cudaMallocManaged(&phi, N*N*sizeof(real_t));
    Check();
    cudaMallocManaged(&state, 6*N*N*sizeof(real_t));
    Check();
    {
    init_default(state,N);
    TacomaNarrowsSweep<<<bcount,bsize>>>(Ay, phi, state, true, false, true, false); Finish();
    std::ofstream ofs("sweep_as_forcing.txt");
    write_sweep(ofs, Ay, phi,state, N);
    }
    {
    init_default(state,N);
    TacomaNarrowsSweep<<<bcount,bsize>>>(Ay, phi, state, true, false, true, true); Finish();
    std::ofstream ofs("eig_as_forcing.txt");
    write_sweep(ofs, Ay, phi,state, N);
    }
    {
    init_default(state,N);
    TacomaNarrowsSweep<<<bcount,bsize>>>(Ay, phi, state, true, false, false, false, true); Finish();
    std::ofstream ofs("mismatch_as_forcing.txt");
    write_sweep(ofs, Ay, phi, state, N);
    }
    {
    std::ofstream ofs("basin.txt");
    init_basin(state,N);
    TacomaNarrowsSweep<<<bcount,bsize>>>(Ay, phi, state, false, false, false, false, false); Finish();
    write_basin(ofs, Ay, phi, state, N);
    }
//    Hysteresis(0.01,0.01);
    cudaFree(Ay);
    cudaFree(phi);
    cudaFree(state);
    return 0;
}
