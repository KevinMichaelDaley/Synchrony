#include "Integrator.h"


class Fitzhugh : public System<Fitzhugh> {
friend class System<Fitzhugh>;
protected:
    virtual void Fn(FixedSizeDoubleBuffer& out, const FixedSizeDoubleBuffer& in, double t) const override{
        double x=in.Get(0);
        double y=in.Get(1);
        out.Set(0, x-(x*x*x/3.0)-y+0.5);
        out.Set(1, (x+0.7-0.8*y)/12.5);
    }
};
int main(int argc, char** argv) {
    double s[2]={-0.1,0.1};
    FixedSizeDoubleBuffer inout_state(s,2);
    Fitzhugh fh;
    for(int t=0; t<10000; ++t){
        fh.Solve(inout_state, t*0.001);
        printf("%lf %lf %lf\n", t*0.001, s[0], s[1]);
    }
}
