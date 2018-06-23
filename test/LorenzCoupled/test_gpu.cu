#include "Integrator.h"
#include <cuda.h>
#include "../../src/FixedSizeBuffer.cpp"
class Lorenz : public System<Lorenz> {
friend class System<Lorenz>;
protected:
    HOSTDEVICE size_t GetIndexOffset(const ThreadInfo& thread_info) const{
        return thread_info.GetIndex()*6;
    }
    HOSTDEVICE size_t GetDimension() const{
        return 6;
    }
    HOSTDEVICE void Fn(FixedSizeDoubleBuffer& out, FixedSizeDoubleBuffer& in, double t, const ThreadInfo& thread_info) const{
        int global_id=thread_info.GetIndex();
        double coupling_strength = 8/256.0 * global_id;
        for(int i=0; i<2; ++i){
            
                
            double x=in.Get(i*3);
            double y=in.Get(i*3+1);
            double z=in.Get(i*3+2);
            size_t i2=1-i;
            double diffus=coupling_strength*(in.Get((i2)*3)-x);
            double dx=10*(y-x);
            double dy=x*(28-z)-y;
            double dz=x*y-(8/3.0)*z;
            out.Set(i*3, dx+diffus);
            out.Set(i*3+1, dy);
            out.Set(i*3+2, dz);
        }
    }
};

__global__ void integrate(double* s, size_t N){
    FixedSizeDoubleBuffer inout_state(s,N);
    Lorenz lz;
    for(int t=0; t<10; ++t){
        lz.Solve(inout_state, t*0.001, 0.001);
    }
}

__host__ int main(int argc, char** argv) {
    double* s;
    cudaMallocManaged(&s, 6*256*sizeof(double));
    for(int i=0; i<6*256; ++i){
        s[i]=(rand()%1000)/10000.0;
    }

    
    for(int t=0; t<1000000/200; ++t){
            integrate<<<1,256>>>(s, 6*256);
            cudaDeviceSynchronize();
            cudaError_t error = cudaGetLastError();
            if(error!=cudaSuccess)
            {
                fprintf(stderr,"ERROR: %s\n", cudaGetErrorString(error) );
                exit(-1);
            }
            printf(" %lf ", t*0.0001);
            
            for(int j=0; j<6*256; ++j){
                printf(" %lf ", s[j]);
            }
            
            printf("\n");
    }//for t
    cudaFree(s);
    return 0;
}//main
