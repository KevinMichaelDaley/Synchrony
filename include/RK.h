
#include <cassert>
template<size_t N> class RealArray{
private:
    real_t arr[N];
public:
    HOSTDEVICE RealArray(real_t a[]){
        for(size_t i=0; i<N; ++i){
            arr[i]=a[i];
        }
    }
    HOSTDEVICE RealArray(real_t a[], size_t subN){
        for(size_t i=0; i<subN; ++i){
            arr[i]=a[i];
        }
    }
    HOSTDEVICE real_t operator[](int i) const{
        return arr[i];
    }
    HOSTDEVICE real_t& operator[](int i){
        return arr[i];
    }
};

class RKButcherTableau {
private:
    const size_t s;
    const RealArray<81> A;
    const RealArray<9> B;
    const RealArray<9> C;
    HOSTDEVICE RKButcherTableau()=delete;
    HOSTDEVICE RKButcherTableau(size_t p_s, 
                     RealArray<81> p_A,
                     RealArray<9> p_B,
                     RealArray<9> p_C):
                    s(p_s),
                    A(p_A), B(p_B), C(p_C)
                     
    {
        
    }
public:
    HOSTDEVICE static RKButcherTableau RK4Classic() {//TODO: add more rk methods including adaptive ones (ode45 and etc.)
        real_t pA[16]={0,   0, 0,   0,
                    0.5, 0, 0,   0,
                    0,   0.5, 0, 0,
                    0,   0,   1.0, 0 };
        real_t pB[9]={1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0};
        real_t pC[9]={0,0.5,0.5,0,0};
        const RKButcherTableau rk4(4,
                RealArray<81>(pA,16),
                RealArray<9>(pB),
                RealArray<9>(pC)
                );
            return rk4;
    }   
    template<bool high_ord> HOSTDEVICE static RKButcherTableau RKDormandPrince() {//ode45
        real_t pA[81]={0,   0, 0,   0, 0, 0, 0,
                    0.2, 0, 0,   0, 0, 0, 0,
                    3/40.0, 9/40.0, 0,0,0,0,0,
                    44/45.0, -56/15.0, 32/9.0, 0, 0, 0, 0,
                    19372.0/6561.0, -25360/2187.0, 64448.0/6561.0, -212.0/729.0, 0,0,0,
                    35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11/84.0, 0
        };
        real_t pB0[9]={35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11/84.0, 0};
        real_t pB1[9]={5179/57600.0, 0.0, 7571/16695.0, 393/640.0,-92097.0/339200.0, 187/2100.0, 1/40.0, 0};
        real_t pC[9]={0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0, 0};
        const RKButcherTableau rkhi(7,
                RealArray<81>(pA,16),
                RealArray<9>(pB0),
                RealArray<9>(pC)
                );
        const RKButcherTableau rklo(7,
                RealArray<81>(pA,16),
                RealArray<9>(pB1),
                RealArray<9>(pC)
                );
            return high_ord? rkhi : rklo;
    }  
    inline HOSTDEVICE size_t GetNumStages() const{
        assert(s<=9);
        return s;
    }
    inline HOSTDEVICE real_t GetXCoeff(size_t step_i, size_t prev_j) const{ 
        assert(step_i<GetNumStages() && prev_j<GetNumStages());
        return A[step_i * s + prev_j]; 
    }
    inline HOSTDEVICE real_t GetYCoeff(size_t step) const{ 
        assert(step<GetNumStages());
        return B[step];
    }
    inline HOSTDEVICE real_t GetTCoeff(size_t step) const {
        assert(step<GetNumStages());
        return C[step];
    }
};
