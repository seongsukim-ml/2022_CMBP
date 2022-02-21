#ifndef ____SwendsenWang_HK____
#define ____SwendsenWang_HK____ 

#define cor(x,y) ((y)*L+x+N)%N

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


using namespace std;

// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

//In 2D Ising model => T_CRIT = 2/log(1+sqrt(2))
const double T_CRIT = 2.269185; 

typedef tuple<int,int> duo;

class SwendsenWang_2D{ 
    public:
        const int L;
        const int N;
        const int Bin;
        const int B;
        const int J;
        bool isTinf;

        const static int init = 1.01;
        const int XNN = 1;
        const int YNN; 

        clock_t __start__, __finish__;

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        // short* sc;                          // Square lattice configuration of 2D Ising model
        vector<short> sc;
        double prob;

        vector<short> iBondLeft;
        vector<short> jBondAbove;

        vector<int> label;
        vector<int> labels;


        unsigned long Fliped_Step = 0;
        unsigned long Total_Step  = 0;
        unsigned long Calc_call = 0;

        // Return name of Model
        static string Name(){return "SwendsenWang_HK";}

        // Constructor of Model
        SwendsenWang_2D(int L, int bin, double B, double J, double Tsrt, double Tfin, bool isTinf);
        SwendsenWang_2D(vector<double> args);
        
        // Destructor of Model
        ~SwendsenWang_2D(){
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }

        // Initialize the Ising model according to beta = 1/kT
        void Initialize(double beta);
        // Initialize the Ising model according to beta = 1/kT
        void ProbCalc(double beta);

        // Boundary condition used for Measurement (calculate right and below interaction)
        int BoundaryHelical(int i);
        // Boundary condition on Sweep (calculate)
        int SweepHelical(int i);

        // Measure sum of magnetizaiton and Hamiltonian of Spin configuration on current situation
        duo Measure();

        // Sub step of Monte Carlo step.
        int Calculate();
        int Calculate(int site);
        void IterateUntilEquilibrium(int equil_time);

        // HK algorithms parts for SwendsenWang
        void MakeBonds();
        void ClusterAndLabel();
        void Union(int x, int y);
        int Find(int x);
        void FlipCluster();
};

SwendsenWang_2D::SwendsenWang_2D(int L, int bin, double B, double J, double Tsrt, double Tfin, bool isTinf) :L(L), N(L*L), Bin(bin), B(B), J(J), YNN(L){
    this-> isTinf = isTinf;
    // this-> sc = new short[N];
    this-> sc = vector<short>(N);

    this-> MV = vector<double>(Bin);
    this-> CV = vector<double>(Bin);
    this-> TV = vector<double>(Bin);
    this-> BetaV = vector<double>(Bin);

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
}

SwendsenWang_2D::SwendsenWang_2D(vector<double> args): SwendsenWang_2D(args[0],args[1],args[2],args[3],args[4],args[5],args[6]){}

void SwendsenWang_2D::ProbCalc(double beta){
    this->prob = 1-exp(-2*beta*J);
}

void SwendsenWang_2D::Initialize(double beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    for(int i = 0; i < N; i++){
        // T = 0 start
        sc[i] = 1;

        // T = \inf start
        if(this->isTinf) this->sc[i] -= int(dis(gen)*2)*2;
    }
    this->ProbCalc(beta);

}

int SwendsenWang_2D::SweepHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;

    if((nn = i - XNN) < 0) nn += this->N;
    sum += this->sc[nn];
    if((nn = i + XNN) >= this->N) nn -= this->N;
    sum += this->sc[nn];
    if((nn = i - YNN) < 0) nn += this->N;
    sum += this->sc[nn];
    if((nn = i + YNN) >= this->N) nn -= this->N;
    sum += this->sc[nn];

    return sum;
}

int SwendsenWang_2D::BoundaryHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;
    
    if((nn = i + XNN) == N) nn = 0;
    sum += this->sc[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += this->sc[nn];

    return sum;
}

duo SwendsenWang_2D::Measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = this->BoundaryHelical(i);
        res += J*sum*sc[i];
        sigma += sc[i];
    }
    
    HH = -res-B*sigma;

    return make_tuple(HH,sigma);
}

// void SwendsenWang_2D::Calculate(int _n){
//     int i, k, delta, n;
//     n = !_n ? (this->N) : _n;
//     double a;
//     for(i = 0; i < n; i++){
//         // Sweep Randomly
//         k = (this->N)*dis(gen);
//         // Sweep Sequential
//         // k = n

//         delta = (this->SweepHelical(k))*(this->sc[k]);
        
//         a = dis(gen);
//         this->Total_Step++;

//         if(delta <= 0){
//             this->Fliped_Step++;
//             this->sc[k] = -(this->sc[k]);
//         } else if(a < prob[delta]){
//             this->Fliped_Step++;
//             this->sc[k] = -(this->sc[k]);
//         }
//     }
// }

int SwendsenWang_2D::Calculate(){
    /*In SwendsenWang algorithm, each clustering and flipping process is one mcs*/
    int size = Fliped_Step;
    
    MakeBonds();
    ClusterAndLabel();
    FlipCluster();
    
    Total_Step++;
    return Fliped_Step-size;
}

int SwendsenWang_2D::Calculate(int site){
    /*In SwendsenWang algorithm, each clustering and flipping process is one mcs*/
    int i, k, l, sp, cur, size = 0;
    
    MakeBonds();
    ClusterAndLabel();
    FlipCluster();

    Total_Step++;
    return size;
}


void SwendsenWang_2D::IterateUntilEquilibrium(int equil_time){
    double k = 0;
    while(++k < equil_time){
        Calculate();
    }
}

void SwendsenWang_2D::MakeBonds(){
    this->iBondLeft  = vector<short>(N);
    this->jBondAbove = vector<short>(N);
    
    int nn;
    for(int i = 0; i < N; i++){
        nn = (i-XNN+N)%N;
        if(sc[nn] == sc[i] && dis(gen) < prob){
            iBondLeft[nn] = true;
        }
        nn = (i-YNN+N)%N;
        if(sc[nn] == sc[i] && dis(gen) < prob){
            jBondAbove[nn] = true;
        }
    }
}

// void SwendsenWang_2D::ClusterAndLabel(){
//     int largest_label = 0;
//     vector<int> label = vector<int>(N);
    
//     delete labels;
//     labels = new int[N];
//     for(int cnt = 0; cnt < N; cnt++) labels[cnt] = cnt;

//     for(int i = 0; i < L; i++){
//         for(int j = 0; j < L; j++){
//             int left  = iBondLeft[cor(i,j)]  ? label[cor(i-1,j)] : 0; // label = 0 means unclusterd coordinate
//             int above = jBondAbove[cor(i,j)] ? label[cor(i,j-1)] : 0;
//             if(left == 0 && above == 0){
//                 largest_label += 1; // Unused new cluster label
//                 label[cor(i,j)] = largest_label;
//             } else if(left != 0 && above == 0){
//                 label[cor(i,j)] = Find(left);
//             } else if(left == 0 && above != 0){
//                 label[cor(i,j)] = Find(above);
//             } else {
//                 Union(left,above);
//                 label[cor(i,j)] = Find(left);
//             }
//         }
//     }
// }

string A = "·╴╶─╷┐┌┬╵┘└┴│┤├┼";

void SwendsenWang_2D::ClusterAndLabel(){
    int largest_label = 0;
    label = vector<int>(N);    
    labels = vector<int>(N+1);
    for(int cnt = 0; cnt < N+1; cnt++) labels[cnt] = cnt;

    int bond[4];
    int bonds_num;

    string res = "";

    for(int j = 0; j < L; j++){
        for(int i = 0; i < L; i++){
            res += to_string(iBondLeft[cor(i,j)]);
        }
        res += "\n";
    }
    cout << "iBondLeft\n" << res;

    res = "";

    for(int j = 0; j < L; j++){
        for(int i = 0; i < L; i++){
            res += to_string(jBondAbove[cor(i,j)]);
        }
        res += "\n";
    }
    cout << "iBondAbove\n" << res;
    
    res = "";
    for(int j = 0; j < L; j++){
        for(int i = 0; i < L; i++){
            res += A.at(iBondLeft[cor(i,j)] + 2*iBondLeft[cor(i+1,j)] + \
                        4*jBondAbove[cor(i,j-1)]+ 8*jBondAbove[cor(i,j)]);
        }
        res += "\n";
    }
    cout << "Box Draw\n" << res;

    for(int j = 0; j < L; j++){
        for(int i = 0; i < L; i++){
            bonds_num = 0;
            
            // PBC(periodic boundary condition)
            // 0,0 always skip
            // i,0 skip above
            // 0,j pbc
            // i,Lx-1 below
            // L-1,L-1 4 direction

            if(!(i == 0 && j == 0) && iBondLeft[cor(i,j)])      bond[bonds_num++] = cor(i-1,j);
            if(j > 0 && jBondAbove[cor(i,j)])                   bond[bonds_num++] = cor(i,j-1);
            if(j == L-1 && jBondAbove[cor(i,j+1)])              bond[bonds_num++] = cor(i,j+1);
            if((i == L-1 && j == L-1) && iBondLeft[0])          bond[bonds_num++] = cor(i+1,j);

            if(bonds_num == 0){
                label[cor(i,j)] = ++largest_label;
            } else {
                // if only one bond
                int min_label = label[bond[0]];
                if (bonds_num != 1){
                    // Find minimum label among bonds
                    min_label = largest_label;
                    for(int ll = 0; ll < bonds_num; ll++){
                        min_label = min_label > label[bond[ll]] ? label[bond[ll]] : min_label;
                    }
                    
                    // Union near label to minimum label
                    for(int ll = 0; ll < bonds_num; ll++){
                        Union(label[bond[ll]],min_label);
                    }
                }
                    label[cor(i,j)] = Find(min_label);
            }
        }
    }
}

// void SwendsenWang_2D::ClusterAndLabel(){
//     int largest_label = 0;
//     label = vector<int>(N);    
//     labels = vector<int>(N);
//     for(int cnt = 0; cnt < N; cnt++) labels[cnt] = cnt;

//     int bond[2];
//     int bonds_num = 0;

//     for(int j = 0; j < L; j++){
//         for(int i = 0; i < L; i++){
//             bonds_num = 0;

//             if(iBondLeft[cor(i,j)])  bond[bonds_num++] = cor(i-1,j);
//             if(jBondAbove[cor(i,j)]) bond[bonds_num++] = cor(i,j-1);

//             if(bonds_num == 0){
//                 label[cor(i,j)] = ++largest_label;
//             } else if (bonds_num == 1) {
//                 if(label[cor(i,j)] == 0){
//                     if(label[bond[0]] == 0)
//                         label[bond[0]] = ++largest_label;
//                     label[cor(i,j)] = Find(label[bond[0]]);
//                 } else {
//                     if(label[bond[0]] == 0)
//                         label[bond[0]] = label[cor(i,j)];
//                     else
//                         Union(label[bond[0]],label[cor(i,j)]);
//                 }
//             } else {
//                 Union(label[bond[0]],label[bond[1]]);
//                 if(label[cor(i,j)] == 0){
//                     label[cor(i,j)] = Find(label[bond[0]]);           
//                 }
//                 else{
//                     Union(label[bond[0]],label[cor(i,j)]);
//                 }
//             }
//         }
//     }
// }

void SwendsenWang_2D::FlipCluster(){
    bool sNewChosen[N] = {};
    short sNew[N];
    for(int i = 0; i < N; i++){
        label[i] = Find(label[i]);
    }

    for(int i = 0; i < N; i++){
        if (!sNewChosen[label[i]]) {    
            sNew[label[i]] = dis(gen) < 0.5 ? +1 : -1;
            sNewChosen[label[i]] = true;
        }
        if(this->sc[i] != sNew[label[i]]){
            this->sc[i] = sNew[label[i]];
            Fliped_Step++;
        }
    }
}

void SwendsenWang_2D::Union(int x, int y){
    // int xx = Find(x);
    // int yy = Find(y);
    // if (xx > yy)
    //   labels[xx] = yy;
    // else
    //   labels[yy] = xx;

    // maybe y should be smaller I guess
    labels[Find(x)] = Find(y);
}

int SwendsenWang_2D::Find(int x){
    int y = x;
    while(labels[y] != y)
        y = labels[y];
    while(labels[x] != x){
        int z = labels[x];
        labels[x] = y;
        x = z;
    }
    return y;
}

#endif // ____SwendsenWang_HK____ 