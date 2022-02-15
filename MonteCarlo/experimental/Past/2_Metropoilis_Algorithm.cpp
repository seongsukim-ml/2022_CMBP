// * beta   = inverse temperature
// * prob[] = array of acceptance probability
// * s[]    = lattice of spins with helical boundary conditions
// * L      = const ant edge length of lattice

// File & IO System
#include <iostream>
#include <fstream>

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

using namespace std;

#define L 5 /*Parameter: lattice size*/
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 40 /*Parametr: Change binning of temperature*/

#define Tstart 0
#define Tfin 5

#ifdef _WIN32
#endif
#ifdef linux
static string Filename = "./Result/Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
#endif


static short s[N]; // Square lattice configuration of 2D Ising model
static double prob[5]; // 1 1 exp 1 exp

static random_device rd;  // Will be used to obtain a seed for the random number engine
static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
static uniform_real_distribution<> dis(0.0, 1.0);

static long Fliped_Step;
static long Total_Step;
static bool isTinf;

void prob_calc(double beta){
    int i;
    for(i = 2; i < 5; i += 2){
        prob[i] = exp(-2*beta*i);
    }
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

void calculate(int n = N){
    int i, k, delta;
    double a;
    for(i = 0; i < n; i++){
        // Sweep Randomly
        k = N*dis(gen);
        // Sweep Sequential
        // k = n

        delta = sweep_helical(k)*s[k];
        
        a = dis(gen);
        Total_Step++;

        if(delta <= 0){
            Fliped_Step++;
            s[k] = -s[k];
        } else if(a < prob[delta]){
            Fliped_Step++;
            s[k] = -s[k];
        }
    }
}

int main(){
    vector<int> value;
    vector<double> m(bin);
    vector<double> c(bin);
    vector<double> T(bin);
    vector<vector<double>> res(bin,vector<double>(5,0));
    string Tat = isTinf ? "inf" : "0";

    cout << "Metropolis Algorithm\n";
    cout << "L = " << L << ", " << "bin = " << bin << ", Start T at " << Tat << "\n";
    cout << "-------------------------------------------------------------------------------------------" << "\n";
    cout << "magnetization-----specific heat-----Fliped Step-----Total Step-----------------------------" << "\n";
    cout << "-------------------------------------------------------------------------------------------" << "\n";

    clock_t start, finish;
    double duration;
    start = clock();

    int sigma, HH, equil_time, epoch;
    double beta;
    for(int i = 0; i < bin; i++){
        if(!Tstart){
            T[i] = Tstart + ((Tfin-Tstart)/(double)(bin))*(i+1);
        } else{
            T[i] = Tstart + ((Tfin-Tstart)/(double)(bin-1))*(i);
        }
    }

    epoch = 18000;
    for(int i = 0; i < bin; i++){ // i is bin
        Fliped_Step = 0;
        Total_Step = 0;

        beta = 1/T[i];

        initialize(beta);

        // Basic equlibrium time is almost 1000 so that twice of that is enough
        equil_time = 2000;

        // Equilibrium time is larger when T is near Tc
        if(T[i]<2.4 || T[i]>2.0) equil_time =10000;

        for(int j = 0; j < equil_time; j++){ // j is epoch
            calculate();
        }
        value = measure();
        HH = value[0];
        sigma = value[1];

        cout <<"idx: " << i << "||" << sigma << " " << HH << "\n";

        // Calculating result
        epoch = 20000;
        for(int j = 0; j < epoch; j++){ // j is epoch
            calculate();
            
            value = measure();
            HH = value[0];
            sigma = value[1];
    
            res[i][0] += abs(sigma)/double(epoch);
            res[i][1] += (sigma/double(epoch)*sigma);
            res[i][2] += (sigma/double(epoch)*sigma)*(sigma*sigma);
            res[i][3] += HH/double(epoch);
            res[i][4] += (HH/double(epoch)*HH);
        }
        m[i] = res[i][0]/N;
        c[i] = (beta*beta)/N*(res[i][4]-(res[i][3]*res[i][3]));

        // cout << m[i] << " " << c[i] << " " << " "<< Fliped_Step << " " << Total_Step << '\n';

    }
    

    // Save the data (Not real time process)
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

        myfile << "idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2,c_dev(bootstrap)\n";
        for(int i = 0; i <bin; i++){
            string temp = to_string(i) + "," + to_string(T[i]) + "," + to_string(m[i]) + "," + to_string(c[i]) + ",";
            temp = temp + to_string(res[i][0]) + "," + to_string(res[i][1]) + "," + to_string(res[i][2]) + ",";
            temp = temp + to_string(res[i][3]) + "," + to_string(res[i][4]) + "\n";

            myfile << temp;
        }
        myfile.close();
        cout << "Save Completed: " << NewFilename << "\n";
    }

    finish = clock();
    cout << "Program Exit. Spent time: " << (double)(finish-start)/CLOCKS_PER_SEC << "\n";
}