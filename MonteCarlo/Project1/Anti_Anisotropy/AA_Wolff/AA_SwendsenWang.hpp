#ifndef ____AA_SWENDSENWANG____
#define ____AA_SWENDSENWANG____ 

// File & IO System
#include <iostream>

// Data Structure
#include <vector>
#include <tuple>
#include <string>
#include <tuple>
#include <set>
#include <algorithm>
#include <queue>

// Mathmatics
#include <math.h>
#include <random>
// Etc.
#include <stdlib.h>
#include <ctime>
// Ewald sum
#include "../../Ewald_sum.cpp"

using namespace std;

// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

//In 2D Ising model => 2/log(1+sqrt(2))
// const double T_CRIT = 2.269185;

typedef tuple<double,int> duo;

class AA_SwendsenWang{

    public:
        const int Lx; // Short range (Time domain)
        const int Ly; // Long range  (Spatial domain)
        const int N;
        const int Bin;
        const double B;
        const double Jx;
        const double Jy;
        double Jx_sign = 1;
        double Jy_sign = 1;
        const double alpha; // dimension = 1
        bool isTinf;

        double cur_beta;

        const int XNN = 1;
        const int YNN;

        double HH;
        int sigma;
        int staggered;
        clock_t __start__, __finish__;

        ewald_ND e2d;

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        short* sc; // Square lattice configuration of 2D Ising model
        // double Prob[5];
        bool* sign;

        unsigned long long Fliped_Step = 0;
        unsigned long long Total_Step  = 0;
        long Calc_call = 0;

        static string Name(){return "AA_SwendsenWang";}
        AA_SwendsenWang(int Lx, int Ly, int bin, double B, double Jx, double Jy, double alpha, double Tsrt, double Tfin, bool isTinf);
        AA_SwendsenWang(vector<double> args);
        ~AA_SwendsenWang(){
            // delete e2d;
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }
        long double Prob(double delta);

        void Initialize(double beta);
        // void Initialzie(int idx);
        void Measure();
        void Measure_fast();
        void Calculate(int _n = 0,bool Random = true);
        void IterateUntilEquilibrium(int equil_time,bool random = true);
};

AA_SwendsenWang::AA_SwendsenWang(int Lx, int Ly, int bin, double B, double Jx, double Jy, double alpha, double Tsrt, double Tfin, bool isTinf)
:Lx(Lx), Ly(Ly), N(Lx*Ly), Bin(bin), B(B), Jx(Jx), Jy(Jy), YNN(Lx), alpha(alpha) {
    this-> isTinf = isTinf;
    this-> sc = new short[N];
    this-> sign = new bool[N];
    // Checkboard style
    // if(Lx%2 == 1)
    //     for(int i = 0; i < N; i++)
    //         sign[i] = i%2;
    // else
    //     for(int i = 0; i < N; i++)
    //         sign[i] = (i+(i/Lx))%2;
    // for(int i = 0; i < N; i++){
    //     cout << -sign[i]*2+1;
    // }
    for(int i = 0; i < N; i++)
        sign[i] = (i/Lx)%2;

    this-> MV = vector<double>(Bin);
    this-> CV = vector<double>(Bin);
    this-> TV = vector<double>(Bin);
    this-> BetaV = vector<double>(Bin);
    if(Jx < 0) Jx_sign = -1;
    if(Jy < 0) Jy_sign = -1;


    __start__ = clock();

    for(int i = 0; i < bin; i++){ 
        if(!Tsrt){
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin))*(i+1);
        } else if(bin == 1){
            this->TV[i] = Tsrt;
        } else {
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin-1))*(i);
        }
        this->BetaV[i] = 1/TV[i];
    }
    // Ewald sum
    e2d = ewald_ND(2,vector<int>({Lx,Ly}),alpha);
}

AA_SwendsenWang::AA_SwendsenWang(vector<double> args):
AA_SwendsenWang(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9]){}

long double AA_SwendsenWang::Prob(double delta){
    return expl(-this->cur_beta*delta);
}

void AA_SwendsenWang::Initialize(double beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    for(int i = 0; i < N; i++){
        // T = 0 start
        sc[i] = 1;

        // T = \inf start
        if(this->isTinf) this->sc[i] -= int(dis(gen)*2)*2;
    }
    // this->ProbCalc(beta);
    cur_beta = beta;

    this->Measure();
}

void AA_SwendsenWang::Measure(){ //O(N^2)
    int i, j;
    int sigma = 0;
    int staggered = 0;
    double HH, result = 0;
    // cout << e2d.pi_ij(0,56) << '\n';
    for(i = 0; i < N; i++){
        double temp = 0;
        for(j = i+Lx; j < N; j+=Lx) // y 방향 long range 계산
            temp += e2d.pi_ij_1D(i,j)*sc[j];
        temp *= Jy;
        if(i%Lx != Lx-1) // x 방향 short ragne 계산
            temp += Jx*sc[i+1];
        else
            temp += Jx*sc[i+1-Lx];
        sigma += sc[i];
        staggered += sc[i]*(-sign[i]*2+1);
        result += temp*sc[i];
    }
    HH = -result-B*sigma;
    this->HH = HH;
    this->sigma = sigma;
    this->staggered = staggered;
}

void AA_SwendsenWang::Measure_fast(){
    // double res[] = {HH,sigma,staggered};
    // return vector<double>(res,res +sizeof(res)/sizeof(res[0]));
    return; // do nothing
}

void AA_SwendsenWang::Calculate(int _n, bool Random){ //O(N^2)
    // 1) 모든 i,j site에 대해서 bond를 계산함
    // 1-1) 만약에 s(i) != s(j) 이면 skip
    // 1-2) s(i) == s(j)이면 P = 1-exp(-2*beta*Jij)의 확률로 bond 생성
    // 2) 모든 site에 대해서 bond 생성 후, BFS로 cluster 탐색. 1/2의 확률로 flip

    double delta;
    // bool bond[N][N];
    bool visited[N];
    memset(visited,0,sizeof(visited));
    vector<set<int>> adj = vector<set<int>>(N);
    // memset(bond,0,sizeof(bond));

    // step 1.
    for(int i = 0; i < N; i++){
        // a) long range (y dir 계산)
        int j = i + Lx;
        while(j < N){
            if(sc[i]*(-sign[i]*2+1) == sc[j]*(-sign[j]*2+1)){
                delta = Jy*e2d.pi_ij_1D(i,j)*sc[i]*sc[j];
                if(dis(gen) < (1-Prob(2*delta))){
                    adj[i].insert(j);
                    adj[j].insert(i);
                }
                // cout << i << " " << j << " " << delta << '\n';
            }
            j += Lx;
        }
        // b) short range (x dir 계산)
        if(((j = i + XNN)-XNN)%Lx == Lx-1) j -= this->Lx;
        if(sc[i] == sc[j]*Jx_sign){
            delta = Jx;
            if(dis(gen) < (1-Prob(2*delta))){
                adj[i].insert(j);
                adj[j].insert(i);
            }
            // cout << i << " " << j << " " << delta << '\n';
        }
    }

    // step 2.
    queue<int> que;
    for(int i = 0; i < N; i++){ //O(N+lambda) = O(N+E)
        // cout << "i is " << i << "\n";
        if(visited[i] == true)
            continue;
        if(adj[i].empty()){
            // cout << "i is empty" << '\n';
            visited[i] = true;
            continue;
        }

        int mul = 2*(int)(dis(gen)*2) -1;

        que.push(i);
        visited[i] = true;
        // BFS
        while(!que.empty()){
            int s = que.front();
            que.pop();
            
            // cout << s << endl;

            sc[s] *= mul; // Flip
            this-> Total_Step++;
            if(mul == -1)
                this->Fliped_Step++;

            auto itr = adj[s].begin();
            while(itr != adj[s].end()){
                int n = *itr++;
                if(!visited[n]){
                    visited[n] = true;
                    que.push(n);
                }
            }
        }
    }
}

void AA_SwendsenWang::IterateUntilEquilibrium(int equil_time,bool random){
    for(int j = 0; j < equil_time; j++)
        Calculate(0,random);
}

#endif // ____AA_SWENDSENWANG____ 