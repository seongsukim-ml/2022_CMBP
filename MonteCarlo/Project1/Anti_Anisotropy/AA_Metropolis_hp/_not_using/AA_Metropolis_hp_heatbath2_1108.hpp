#ifndef ____AA_Metropolis____
#define ____AA_Metropolis____

// File & IO System
#include <iostream>

// Data Structure
#include <vector>
#include <tuple>
#include <string>
#include <tuple>
// Mathmatics
#include <math.h>
#include <random>
// Etc.
#include <stdlib.h>
#include <ctime>
// Ewald sum
#include "../../Ewald_sum.cpp"
#include <boost/dynamic_bitset.hpp>
#include <cblas.h>

#ifndef __BIT_CONVERSION__
#define __BIT_CONVERSION__
#define bs(bit) (2*bit-1)
#define sb(spin) ((spin+1)/2)
#endif

#ifndef __TYPE_DEFINE__
#define __TYPE_DEFINE__
    typedef int32_t INT1;
    typedef int64_t INT2;
    typedef int16_t INT8;
    typedef double FLOAT1;
    // typedef long double FLOAT1;
    typedef long double FLOAT2;
#endif
// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

using namespace std;

static random_device rd;  // Will be used to obtain a seed for the random number engine
static long unsigned int seed = rd();
// vector<int> seeds(16);
// static long unsigned int seed = static_cast<long unsigned int>(time(0));

static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<long double> dis(0.0, 1.0);
static bernoulli_distribution bern(0.5);

//In 2D Ising model => 2/log(1+sqrt(2))
// const double T_CRIT = 2.269185;

class AA_Metropolis{

    public:
        const INT1 Lx; // Short range (Time domain)
        const INT1 Ly; // Long range  (Spatial domain)
        const INT1 N;
        const INT1 Bin;
        const FLOAT1 B;
        const FLOAT1 Jx;
        const FLOAT1 Jy;
        const FLOAT1 alpha; // dimension = 1
        bool isTinf;

        FLOAT1 cur_beta;

        const INT1 XNN = 1;
        const INT1 YNN;

        FLOAT1 HH;
        INT1 sigma;
        INT1 staggered;
        clock_t __start__, __finish__;

        ewald_ND e2d;

        vector<FLOAT1> MV;
        vector<FLOAT1> CV;
        vector<FLOAT1> BV;
        vector<FLOAT1> TV;
        vector<FLOAT1> BetaV;
        vector<FLOAT1> res;

        // short* sc; // Square lattice configuration of 2D Ising model
        vector<INT8> sc;
        vector<FLOAT1> heatbath;
        // N size of correation result vector: cor[n] is correlation btw 0th site spin and nth spin;
        float* cor_long;
        float* cor_short;
        float *cor_short_res, *cor_long_res;


        vector<INT8> sign;
        vector<INT8> linked_cb; // linked checkerboard;
        vector<INT8> linked_rn; // linked checkerboard;
        vector<INT8> linked_ln; // linked checkerboard;
        // FLOAT1 prob[5];
        // bool* sign;

        INT2 Fliped_Step = 0;
        INT2 Total_Step  = 0;
        INT2 Calc_call = 0;

        static string Name(){return "AA_Metropolis";}
        AA_Metropolis(INT1 Lx, INT1 Ly, INT1 bin, FLOAT1 B, FLOAT1 Jx, FLOAT1 Jy, FLOAT1 alpha, FLOAT1 Tsrt, FLOAT1 Tfin, bool isTinf);
        AA_Metropolis(vector<FLOAT1> &args);
        ~AA_Metropolis(){
            // delete e2d;
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (FLOAT1)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }
        FLOAT2 Prob(FLOAT1 delta);
        void Initialize(FLOAT1 beta);
        // void Initialzie(int idx);
        void Measure();
        void Measure_fast();
        // void Calculate(INT1 _n = 0,bool Random = true);
        void Calculate(INT1 _n = 0,bool Random = false);
        void IterateUntilEquilibrium(INT1 equil_time,bool random = true);
        void CalculateCorrelation_reset();
        void CalculateCorrelation_update();
        void CalculateCorrelation();
        void Heatbath_Measure();
        // void CalculateCorrelationStaggered();
};

AA_Metropolis::AA_Metropolis(INT1 Lx, INT1 Ly, INT1 bin, FLOAT1 B, FLOAT1 Jx, FLOAT1 Jy, FLOAT1 alpha, FLOAT1 Tsrt, FLOAT1 Tfin, bool isTinf)
:Lx(Lx), Ly(Ly), N(Lx*Ly), Bin(bin), B(B), Jx(Jx), Jy(Jy), YNN(Lx), alpha(alpha) {
    this-> isTinf = isTinf;
    // this-> sc = new short[N];
    this-> sc        = vector<INT8>(N,1);
    this-> heatbath  = vector<FLOAT1>(N);

    this-> cor_long  = new float[N];
    this-> cor_short = new float[N];
    this-> sign      = vector<INT8>(N);
    this-> linked_cb = vector<INT8>(N);
    this-> linked_ln = vector<INT8>(N);
    this-> linked_rn = vector<INT8>(N);

    cor_short_res = new float[Lx];
    cor_long_res = new float[Ly];

    // Checkboard style
    // if(Lx%2 == 1)
    //     for(INT1 i = 0; i < N; i++)
    //         sign[i] = i%2;
    // else
    //     for(INT1 i = 0; i < N; i++)
    //         sign[i] = (i+(i/Lx))%2;
    // for(INT1 i = 0; i < N; i++){
    //     cout << -sign[i]*2+1;
    // }

    for(INT1 i = 0; i < N; i++)
        sign[i] = -2*((i/Lx)%2)+1;

    for(int i = 0; i < N; i++){
        int k;
        if(N % 2 == 0){ // N % 2 == 0
            k = 2*i;
            if(k < N) k = (INT1(k/Lx))%2 == 0 ? k+1 : k;
            else k = (INT1(k/Lx))%2 == 0 ? k-N : k-N+1;
        } else{
            k = 2*i >= N ? 2*i-N: 2*i;
        }
        // cout << i << ' ' << k << '\n';
        linked_cb[i] = k;
        // cout << i << linked_cb[i] << '\n';
    }

    for(int i = 0; i < N; i++){
        int nn;
        if(((nn = i - XNN)+1)%Lx == 0) nn += this->Lx;
        linked_ln[i] = nn;
        if(((nn = i + XNN)-1)%Lx == Lx-1) nn -= this->Lx;
        linked_rn[i] = nn;
    }

    this-> MV = vector<FLOAT1>(Bin);
    this-> CV = vector<FLOAT1>(Bin);
    this-> BV = vector<FLOAT1>(Bin);
    this-> TV = vector<FLOAT1>(Bin);
    this-> BetaV = vector<FLOAT1>(Bin);

    __start__ = clock();

    for(INT1 i = 0; i < bin; i++){
        if(!Tsrt){
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(FLOAT1)(bin))*(i+1);
        } else if(bin == 1){
            this->TV[i] = Tsrt;
        } else {
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(FLOAT1)(bin-1))*(i);
        }
        this->BetaV[i] = 1/TV[i];
    }
    // Ewald sum
    e2d = ewald_ND(2,vector<INT1>({Lx,Ly}),alpha);
}

AA_Metropolis::AA_Metropolis(vector<FLOAT1> &args):
AA_Metropolis(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9]){}

FLOAT2 AA_Metropolis::Prob(FLOAT1 delta){
    // Look Up table
    return expl(-this->cur_beta*delta);
}

void AA_Metropolis::Initialize(FLOAT1 beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    // Initialize the spin configuration
    sc = vector<INT8>(N,1);
    // if(Random)
    for(INT1 i = 0; i < N; i++){
        // T = \inf start
        if(this->isTinf) this->sc[i] = bern(gen);
    }
    // this->ProbCalc(beta);
    cur_beta = beta;

    this->Measure();
    Heatbath_Measure();
}

void AA_Metropolis::Measure(){ //O(N^2)
    INT1 i, j;
    INT1 sigma = 0;
    INT1 staggered = 0;
    FLOAT1 HH, result = 0;
    // cout << e2d.pi_ij(0,56) << '\n';
    for(i = 0; i < N; i++){
        FLOAT1 partial_sum  = 0;
        for(j = i+Lx; j < N; j+=Lx) // y 방향 long range 계산
            partial_sum += e2d.pi_ij_1D(i,j)*sc[j];
        partial_sum *= Jy;
        if(i%Lx != Lx-1) // x 방향 short ragne 계산
            partial_sum += Jx*sc[i+1];
        else
            partial_sum += Jx*sc[i+1-Lx];
        sigma += sc[i];
        staggered += sc[i]*sign[i];
        result += partial_sum*sc[i];
    }
    // result*= 0.5;   // 전체의 절반만 계산했음. 해밀토니안에 2를 안곱해도 계수 조정해서 찾을 수 있지 않을까?
                    // i=0, j=i
    HH = -result-B*sigma;
    this->HH = HH;
    this->sigma = sigma;
    this->staggered = staggered;
}

void AA_Metropolis::Heatbath_Measure(){ //O(N^2)
    INT1 i, j;
    // cout << e2d.pi_ij(0,56) << '\n';
    for(i = 0; i < N; i++){
        FLOAT1 partial_sum  = 0;

        for(int jj = i%Lx; jj < N; jj+=Lx)  // Long range diff
            partial_sum += e2d.pi_ij_1D(jj,i)*sc[jj];
        partial_sum *= Jy;
        // cout << partial_sum << " ";
        partial_sum += Jx*(sc[linked_ln[i]] + sc[linked_rn[i]]);
        // partial_sum *= 2*sc[i];
        heatbath[i] = partial_sum;
        // cout << heatbath[i] << '\n';
    }
}

void AA_Metropolis::Measure_fast(){
    // FLOAT1 res[] = {HH,sigma,staggered};
    // return vector<FLOAT1>(res,res +sizeof(res)/sizeof(res[0]));
    return; // do nothing
}

// void AA_Metropolis::Calculate(INT1 _n, bool Random){ //O(N^2)
//     INT1 i, k, n;
//     // n = !_n ? (this->N) : _n;
//     n = N;
//     for(i = 0; i < n; i++){
//         // // Sweep Randomly
//         // if(Random){ // 골고루 분포되어있는 랜덤을 찾아보면 쓸 수 있음
//         //     k = (this->N)*dis(gen);
//         // // Sweep Sequential
//         // } else {
//         //     k = linked_cb[i];
//         // }
//         k = linked_cb[i];


//         FLOAT1 delta = 0;
//         // FIXME: Linked list
//         for(INT1 jj = k%Lx; jj < N; jj+=Lx)  // Long range diff
//             delta += e2d.pi_ij_1D(jj,k)*sc[jj]; // 여기에 오류가 있었음 pi_ij(jj,k) 함수를 호출하고 있었음
//         delta *= Jy;
//                                              // Short range diff
//         delta += Jx*(sc[linked_ln[k]]+sc[linked_rn[k]]);

//         delta *= 2*sc[k];
//         this->Total_Step++;

//         // if((delta <= 0) || (dis(gen) < Prob(delta))){
//         if((delta <= 0) || (logl(dis(gen)) < -cur_beta*delta)){
//             this->Fliped_Step++;
//             sc[k] = -sc[k];
//             this->sigma     += 2*(sc[k]);
//             this->staggered += 2*(sc[k]*sign[k]);
//             this->HH += delta;
//         }
//     }
// }

void AA_Metropolis::Calculate(INT1 _n, bool Random){ //O(N^2)
    INT1 i, k, n;
    // n = !_n ? (this->N) : _n;
    n = N;
    for(i = 0; i < n; i++){
        if(Random){ // 골고루 분포되어있는 랜덤을 찾아보면 쓸 수 있음
            k = (this->N)*dis(gen);
        // Sweep Sequential
        } else {
            k = linked_cb[i];
        }
        FLOAT1 delta = 2 * heatbath[k] * sc[k];
        Total_Step++;
        // cout << delta << '\n';

        // if((delta <= 0) || (dis(gen) < Prob(delta))){
        if((delta <= 0) || (logl(dis(gen)) < -cur_beta*delta)){
            this->Fliped_Step++;
            sc[k] = -sc[k];
            this->sigma     += 2*(sc[k]);
            this->staggered += 2*(sc[k]*sign[k]);
            this->HH += delta;
            // O(sqrt(N))
            for(INT1 jj = k%Lx; jj < N; jj+=Lx){
                heatbath[jj] += Jy*e2d.pi_ij_1D(jj,k)*2*sc[k];
            }  // Long range diff
            heatbath[linked_ln[k]] += Jx*2*sc[k];
            heatbath[linked_rn[k]] += Jx*2*sc[k];
        }
    }
}

void AA_Metropolis::IterateUntilEquilibrium(INT1 equil_time,bool random){
    for(INT1 j = 0; j < equil_time; j++)
        Calculate(0,random);
}

void AA_Metropolis::CalculateCorrelation_reset(){
    for(int i = 0; i < N; i++)
        cor_short[i] = 0;
    for(int i = 0; i < N; i++)
        cor_long[i] = 0;
}

void AA_Metropolis::CalculateCorrelation_update(){
    for(int i = 0; i < 1; i++){
        for(int j = 0; j < Lx; j++){
            cor_short_res[j] = sc[i*Lx + j];
        }
        cblas_sspr(CblasRowMajor,CblasUpper,Lx,1,cor_short_res,1,cor_short);
    }
    for(int i = 0; i < 1; i++){
        for(int j = 0; j < Ly; j++){
            cor_long_res[j] = sc[j*Lx + i];
        }
        cblas_sspr(CblasRowMajor,CblasUpper,Ly,1,cor_long_res,1,cor_long);
    }
}

void AA_Metropolis::CalculateCorrelation(){
    CalculateCorrelation_reset();
    CalculateCorrelation_update();
    // memset(cor_short,0,sizeof(cor_short));
    // memset(cor_long,0,sizeof(cor_long));

    // for(int i = 0; i < N; i++)
    //     cor_short[i] = 0;
    // for(int i = 0; i < N; i++)
    //     cor_long[i] = 0;

    // for(int i = 0; i < 1; i++){
    //     for(int j = 0; j < Lx; j++){
    //         cor_short_res[j] = sc[i*Lx + j];
    //     }
    //     cblas_sspr(CblasRowMajor,CblasUpper,Lx,1,cor_short_res,1,cor_short);
    // }
    // for(int i = 0; i < 1; i++){
    //     for(int j = 0; j < Ly; j++){
    //         cor_long_res[j] = sc[j*Lx + i];
    //     }
    //     cblas_sspr(CblasRowMajor,CblasUpper,Ly,1,cor_long_res,1,cor_long);
    // }
}
#endif // ____AA_Metropolis____