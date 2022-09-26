#include "AA_MODEL.hpp"

namespace model::AA{
    void BIT_MODEL::Initialize(){
        this->Initialize(ewald_ND(1,vector<int>({Lx,Ly}),alpha));
    }

    void BIT_MODEL::Initialize(ewald_ND e1d){
        this->sc = boost::dynamic_bitset<>(N);
        this->sb = boost::dynamic_bitset<>(N);
        // Random Initialize;
        if(isTinf){
            for(int i = 0; i < N; i++){
                sc[i] = rand.bern();
            }
        }

        for(int i = 0; i < N; i++)
            this->sb[i] = (i/Lx)%2;

        this->e1d = e1d;
        this->Measure();
    }

    void BIT_MODEL::Measure(){
        this->total_spin = 2*sc.count()-N;
        this->staggered_spin = 2*N-(sc^sb).count()-N;

        double result = 0;

        for(int i = 0; i < N; i++){
            double sum = 0;
            for(int j = i + Lx; j < N; j+= Lx)
                sum += e1d.pi_ij_1D(i,j)*ss(sc[j]);
            cout << sum << '\n';
            sum *= Jy;
            if(i%Lx != Lx-1)
                sum += Jx*ss(sc[i+1]);
            else
                sum += Jx*ss(sc[i+1-Lx]);
            result += sum*ss(sc[i]);
        }

        this->HH = -result-B*total_spin;
    }
    void BIT_MODEL::Measure_fast(){
        Measure();
        // Need to be implemented
        return;
    }
}
