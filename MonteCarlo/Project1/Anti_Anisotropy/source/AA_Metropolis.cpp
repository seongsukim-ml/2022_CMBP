// #include "AA_CONf.hpp"
#include "AA_Metropolis.hpp"
// #include ""

namespace model::AA{
    void Metropolis::zerostep(){
        this->total_step  = 0;
        this->fliped_step = 0;
    }

    Metropolis::Metropolis(model::AA::MODEL_CONF &PROFILE, ewald_ND &ewd, myrnd &rand)
    : AA_MODEL(PROFILE,ewd,rand){
        zerostep();
    }

    Metropolis::Metropolis(model::AA::MODEL_CONF &PROFILE)
    : AA_MODEL(PROFILE){
        zerostep();
    }

    void Metropolis::SetTemp(FLOAT1 cur_beta){

    }
    void Metropolis::do_step(bool wipeRandom = false){
        INT1 k, n;
        INT1 XNN = 1;
        for(int i = 0; i < N; i++){
            if(wipeRandom){
                k = (this->N)*(this->rand.unif())
            } else if(N % 2 == 0){
                k = 2*i;
                if(k < N) k = (int(k/Lx))%2 == 0 ? k+1 : k;
                else k = (int(k/Lx))%2 == 0 ? k-N : k-N+1;
            } else {
                k = 2*i >= N ? 2*i-N: 2*i;
            }

            FLOAT2 delta = 0;
            for(int jj = k%Lx; jj < N; jj += Lx)
                delta += Jy*e2d.pi_ij_1D(jj,k)*ss(sc[jj]);

            INT1 nn;

            if(((nn = k - XNN)+1)%Lx == 0) nn += this->Lx;
            delta += Jx*ss(sc[nn]);
            if(((nn = k + XNN)-1)%Lx == Lx-1) nn -= this->Lx;
            delta += Jx*ss(sc[nn]);

            delta *= 2*ss(sc[k]);
            this->total_step++;

            if((delta <= 0) || (dis(gen) < Prob(delta))){
                this->fliped_step++;
                this->ss(sc[k]) *= -1;
                this->total_spin += 2*ss(sc[k]);
                this->staggered_spin +=2*ss(!(sc[k]^sign[k]));
                this->HH += delta;
            }
        }
    }
    void Metropolis::Measure_fast(){return;}
    void Metropolis::IterateUntilEquilibrium(INT1 equil_time, bool wipeRandom){
        for(int j = 0; j < equil_time; j++)
            Calculate(0,random);
    }
}