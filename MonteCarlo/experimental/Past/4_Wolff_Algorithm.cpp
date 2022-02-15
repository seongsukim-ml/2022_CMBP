// * beta   = inverse temperature
// * prob[] = array of acceptance probability
// * s[]    = lattice of spins with helical boundary conditions
// * L      = const ant edge length of lattice

// File & IO System
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef linux
#include <sys/stat.h>
#endif
// Data Structure
#include <vector>
#include <tuple>
#include <string>
// Mathmatics
#include <math.h>
#include <random>
// Etc.
#include <stdlib.h>
#include <ctime>
// #include <format>

using namespace std;

#define L 100
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 25

#ifdef _WIN32
static string Filename = ".\\Result\\Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
#endif
#ifdef linux
static string Filename = "./Result/Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
#endif


#define _init 1.01
#define Tstart 2.2
#define Tfin 2.4


static short s[N]; // Square lattice configuration of 2D Ising model

static string c; // Cluster cache of Wolff algorithm step
static double prob; // Probability to add cluster

static bool isTinf; // Indicator of whether T starts at inf or not (if T is inf, then intial s[N] start with random distribution)
static int init;

/* Random number generator */
static random_device rd;  // Will be used to obtain a seed for the random number engine
static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
static uniform_real_distribution<> dis(0.0, 1.0);

static int Fliped_Step;
static int Total_Step;

void prob_calc(double beta){
    prob = 1-exp(-2*beta*J);
}

void initialize(double beta){
    int i;
    // T = 0 start
    for(i = 0; i < N; i++){
        // T = 0 start
        isTinf = false; 
        s[i] = 1; 

        // T = \inf start
        isTinf = true;
        if(isTinf) s[i] -= int(dis(gen)*2)*2;
    }
    prob_calc(beta);
    
}

int sweep_helical(int i){
    int nn, sum = 0;

    if((nn = i - XNN) < 0) nn += N;
    sum += s[nn];
    if((nn = i + XNN) >= N) nn -= N;
    sum += s[nn];
    if((nn = i - YNN) < 0) nn += N;
    sum += s[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += s[nn];

    return sum;
}

int helical(int i){
    int nn, sum;
    
    if((nn = i + XNN) == N) nn = 0;
    sum = s[nn];

    if((nn = i + YNN) >= N) nn -= N;
    sum += s[nn];

    return sum;
}

vector<int> measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = helical(i);
        res += J*sum*s[i];
        sigma += s[i];
    }
    
    HH = -res-B*sigma;

    return vector<int>({HH,sigma});
}

void flush_clsuster(){
    c = string(N,'0');
}

int flip_cluster(){
    int i, size = 0;
    for(i = 0; i < N; i++){
        if(c.at(i) == '1'){
            Fliped_Step++;
            size++;
            s[i] = -s[i];
        }
    }
    return size;
}

void clustering(int k, double padd = prob){
    if(c.at(k) == '1') return;
    else{
        if((s[k] == init) && (dis(gen) < padd)) c.at(k) = '1';
        else return;
    }
    int nn;
    if((nn = k - XNN) < 0) nn += N;
    clustering(nn);
    if((nn = k + XNN) >= N) nn -= N;
    clustering(nn);
    if((nn = k - YNN) < 0) nn += N;
    clustering(nn);
    if((nn = k + YNN) >= N) nn -= N;
    clustering(nn);
}

double calculate(){
    int i, k, size = 0;
    /*In Wolff algorithm, each clustering and flipping process is one mcs*/
    
    /*k: Choosing random initial site of Wolff step */
    k = N*dis(gen);
    
    flush_clsuster();
    init = s[k];        
    clustering(k,_init);
    size += flip_cluster();
    
    Total_Step++;

    return size;
}

int main(){
    vector<int> value;
    vector<double> m(bin);
    vector<double> c(bin);
    vector<double> size(bin);
    vector<vector<double>> res(bin,vector<double>(4,0));
    string Tat = isTinf ? "inf" : "0";

    cout << "Wolff Algorithm\n";
    cout << "L = " << L << ", " << "bin = " << bin << ", Start T at " << Tat << "\n";
    cout << "--------------------------------------------------------------------" << "\n";
    cout << "magnetization-----specific heat-----Fliped Step-----Total Step------" << "\n";
    cout << "--------------------------------------------------------------------" << "\n";

    clock_t start, finish;
    double duration;
    start = clock();

    int sigma, HH, equil_time, epoch;
    double T, beta;

    for(int i = 0; i < bin; i++){ // i is bin
        Fliped_Step = 0;
        Total_Step = 0; 
        T = Tstart + ((Tfin-Tstart)/bin)*(i+1);
        beta = 1/T;

        initialize(beta);

        // Basic equlibrium time is almost 1000 so that twice of that is enough
        equil_time = 2000;

        // Equilibrium time is larger when T is near Tc
        if(T<2.4 || T>2.0) equil_time =10000;

        for(int j = 0; j < equil_time; j++){ // j is epoch
            calculate();
        }
        value = measure();
        HH = value[0];
        sigma = value[1];
        cout <<"idx: " << i << " || T/J=" << T << " |equi time passed| sig=" << sigma << " || H=" << HH << "\n";

        // Calculating result
        epoch = 18000;
        for(int j = 0; j < epoch; j++){ // j is epoch
            size[i] += calculate()/epoch;
            
            value = measure();
            HH = value[0];
            sigma = value[1];
           
            res[i][0] += abs(sigma)/double(epoch);
            res[i][1] += (sigma*sigma)/double(epoch);
            res[i][2] += HH/double(epoch);
            res[i][3] += (HH*HH)/double(epoch);
        }
        m[i] = res[i][0]/N;
        c[i] = (beta*beta)/N*(res[i][3]-(res[i][2]*res[i][2]));
        cout << m[i] << " " << c[i] << " " << Fliped_Step << " " << Total_Step << '\n';
    }
    

    // Save the data (It's NOT real time saving process)
    bool save = true;
    if(save){
        #ifdef _WIN32
        CreateDirectory("Result", NULL);
        #endif

        #ifdef linux
        mkdir("Result", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        #endif

        ofstream myfile;
        string NewFilename = Filename + ".csv";
        ifstream f(NewFilename);
        int filenum = 1;

        if(f.good()){
            NewFilename = Filename + "_" + to_string(filenum++) + ".csv";
            ifstream f(NewFilename);
        }
        while(true){
            NewFilename = Filename + "_" + to_string(filenum++) + ".csv";
            ifstream f(NewFilename);
            if(!f.good()) break;
        }

        myfile.open(NewFilename);

        myfile << "idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,HH,HH**2,size\n";
        for(int i = 0; i <bin; i++){
            string temp = to_string(i) + "," + to_string(Tstart + ((Tfin-Tstart)/double(bin))*(i+1)) + "," + to_string(m[i]) + "," + to_string(c[i]) + ",";
            temp = temp + to_string(res[i][0]) + "," + to_string(res[i][1]) + "," + to_string(res[i][2]) + "," + to_string(res[i][3]) + ",";
            temp = temp + to_string(res[i][3]) + "\n";
            myfile << temp;
        }
        myfile.close();
        cout << "Save Completed: " << NewFilename << "\n";
    }

    finish = clock();
    cout << "Program Exit. Spent time: " << (double)(finish-start)/CLOCKS_PER_SEC << "\n";
}