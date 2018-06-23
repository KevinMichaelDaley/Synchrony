#include "Integrator.h"


class Lorenz : public System<Lorenz> {
friend class System<Lorenz>;
protected:
    virtual void Fn(FixedSizeDoubleBuffer& out, const FixedSizeDoubleBuffer& in, double t) const override{
        for(int i=0; i<100; ++i){
            
                
            double x=in.Get(i*3);
            double y=in.Get(i*3+1);
            double z=in.Get(i*3+2);
            double i2=(i+1)%100;
            double i1=i-1+(i==0)*99;//nearest neighbor cycle graph
            double diffus=800.0*(in.Get((i2)*3)-x);
            diffus+=800.0*(in.Get((i1)*3)-x);    
            double dx=10*(y-x);
            double dy=x*(28-z)-y;
            double dz=x*y-(8/3.0)*z;
            out.Set(i*3, dx+diffus);
            out.Set(i*3+1, dy);
            out.Set(i*3+2, dz);
        }
    }
};
int main(int argc, char** argv) {
    double s[300];
    for(int i=0; i<300; ++i){
        s[i]=(rand()%1000)/10000.0;
    }
    FixedSizeDoubleBuffer inout_state(s,300);
    Lorenz lz;
    for(int t=0; t<1000000; ++t){
        lz.Solve(inout_state, t*0.0001, 0.0001);
        if(t%200==0){
            printf(" %lf ", t*0.001);
        
            for(int j=0; j<300; ++j){
                printf(" %lf ", s[j]);
            }
            
            printf("\n");
        }//if sample
    }//for t
    return 0;
}//main
