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
#define bin 40

#ifdef _WIN32
static string Filename = ".\\Result\\Wolff_c_"+to_string(L)+"_int"+to_string(bin);
#endif
#ifdef linux
static string Filename = "./Result/Wolff_c_"+to_string(L)+"_int"+to_string(bin);
#endif

#define _init 1.01
#define Tstart 0
#define Tfin 5

static short s[N]; // Square lattice configuration of 2D Ising model

static string c; // Cluster cache of Wolff algorithm step
static double prob; // Probability to add cluster

static bool isTinf; // Indicator of whether T starts at inf or not (if T is inf, then intial s[N] start with random distribution)

/* Random number generator */
// 이게 진짜 random이 아니었어..? seed가 항상 같았던거? -> 같았음! rd()는 system의 hardware에 따라서 제대로 작동 안할 수 있음..
// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

static unsigned long long int Fliped_Step;
static unsigned long long int Total_Step;

// 최적화 옵션

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

double calculate(){
    /*In Wolff algorithm, each clustering and flipping process is one mcs*/
    int i, k, l, sp, cur, size = 0;
    short snew, sold;
    vector<int> stack(N);
    vector<int> nn(4);

    /*k: Choosing random initial site of Wolff step */
    k = N*dis(gen);
    
    stack[0] = k;
    sp = 1; // stack pointer
    
    sold = s[k];
    snew = -s[k];
    
    s[k] = snew;

    while(sp){
        cur = stack[--sp];
        if((nn[0] = cur - XNN) <  0) nn[0] += N;
        if((nn[1] = cur + XNN) >= N) nn[1] -= N;
        if((nn[2] = cur - YNN) <  0) nn[2] += N;
        if((nn[3] = cur + YNN) >= N) nn[3] -= N;

        for(l = 0; l < 4; l++){
            if(s[nn[l]] == sold){
                if(dis(gen) < prob){
                    stack[sp++] = nn[l];
                    s[nn[l]] = snew;
                    Fliped_Step++;
                    size++;
                }
            }
        }
    }
    
    Total_Step++;
    return size;
}

int main(){
    vector<int> value;
    vector<double> m(bin);
    vector<double> c(bin);
    vector<double> size(bin);
    vector<double> T(bin);
    vector<vector<double>> res(bin,vector<double>(5,0));
    string Tat = isTinf ? "inf" : "0";

    cout << "Wolff Algorithm\n";
    cout << "Radnomness test(seed): " << seed << endl;
    cout << "L = " << L << ", " << "bin = " << bin << ", Start T at " << Tat << "\n";
    cout << "--------------------------------------------------------------------" << "\n";
    cout << "magnetization-----specific heat-----Fliped Step-----Total Step------" << "\n";
    cout << "--------------------------------------------------------------------" << "\n";

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
    
    for(int i = 0; i < bin; i++){ // i is bin
        Fliped_Step = 0;
        Total_Step = 0; 
        beta = 1/T[i];

        initialize(beta);

        // Basic equlibrium time is almost 1000 so that twice of that is enough
        equil_time = 1000;

        // // Equilibrium time is larger when T is near Tc
        // if(T<2.4 || T>2.0) equil_time =10000;

        for(int j = 0; j < equil_time; j++){ // j is epoch
            calculate();
        }

        value = measure();
        HH = value[0];
        sigma = value[1]; 
        cout <<"idx: " << i << " || T/J=" << T[i] << " |equi time passed| sig=";
        cout << sigma << " (" <<sigma/double(N)<< ")" << " || H=" << HH << " (" <<HH/double(N)<< ")" "\n";

        // Calculating result
        epoch = 14000;
        for(int j = 0; j < epoch; j++){ // j is epoch
            size[i] += calculate()/double(epoch);
            
            value = measure();
            HH = value[0];
            sigma = value[1];
           
            res[i][0] += abs(sigma)/double(epoch);
            res[i][1] += (sigma/double(epoch)*sigma);
            res[i][2] += (sigma/double(epoch)*sigma)*(sigma*sigma);
            res[i][3] += HH/double(epoch)/L;
            res[i][4] += (HH/double(epoch)*HH)/N;
        }
        // 1step = 1mcs * <n>/L^d << 이걸 고려해야하나?
        m[i] = res[i][0]/N;
        c[i] = (beta*beta)*(res[i][4]-(res[i][3]*res[i][3])); // over N을 위에 코드에 넣음
        cout << m[i] << " " << c[i] << " " << Fliped_Step << " " << Total_Step << '\n';
    }
    
    cout.copyfmt(ios(NULL));

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

        myfile << "idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2\n";
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