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

#define L 5
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 40

#ifdef _WIN32
static string Filename = ".\\Result\\Exact_cpp_"+to_string(L)+"_int"+to_string(bin);
#endif
#ifdef linux
static string Filename = "./Result/Exact_cpp_"+to_string(L)+"_int"+to_string(bin);
#endif

#define Tstart 0
#define Tfin 5

static short s[N];
static long count;

static vector<double> ZHH2(bin,0);
static vector<double> ZHH(bin,0);
static vector<double> ZM4(bin,0);
static vector<double> ZM2(bin,0);
static vector<double> ZM(bin,0);
static vector<double> Z(bin,0);
static vector<double> T(bin);

// Function for prints spin configuration
void print_s(){
    string res = "";
    for(int i = 0; i < N; i++){
        res += char(s[i]+'0');
        res += " ";
        if(i%L == L-1){
            res.pop_back();
            res += "\n";
        }
    }
    cout << res << '\n';
}

// Fucntion that makes all spins up configuration
void initialize(){
    int i, k;

    for(i = 0; i < N; i++)
        s[i] = 1;

    for(int i = 0; i < bin; i++){
        if(!Tstart){
            T[i] = Tstart + ((Tfin-Tstart)/(double)(bin))*(i+1);
        } else{
            T[i] = Tstart + ((Tfin-Tstart)/(double)(bin-1))*(i);
        }
    }
}

// Function for calculate Hamiltonian
void sweep(){
    int i, k;
    int nn, sum, res = 0, sigma = 0;
    int HH;
    double Temp,exp_HH;
    for(i = 0; i < N; i++){

        if((nn = i + XNN) == N) nn = 0;
        sum = s[nn];
        if((nn = i + YNN) >= N) nn -= N;
        sum += s[nn];

        sigma += s[i];
        res += sum*s[i];
    }
    HH = -J*res + B*sigma;
    
    for(k = 0; k < bin; k++){
        exp_HH = exp(-HH/T[k]);
        Z[k] += exp_HH;
        ZM[k] += abs(sigma)*exp_HH;
        ZM2[k] += sigma*sigma*exp_HH;
        ZM4[k] += sigma*sigma*sigma*sigma*exp_HH;
        ZHH[k] += HH*exp_HH;
        ZHH2[k] += HH*exp_HH*HH;
    }
}

void generator_sym(int i){
    if(i == N) return;
    if(i == 1) sweep();
    s[i] = -1;
    sweep();
    count++;
    generator_sym(i+1);
    s[i] = 1;
    generator_sym(i+1);
}

void print_all(double* a, int len){
    string res = "";
    for(int i = 0; i < len; i++){
        res += to_string(a[i]) + " ";
    }
    res += "\n";
    cout << res << endl;
}

void print_all(int* a, int len){
    string res = "";
    for(int i = 0; i < len; i++){
        res += to_string(a[i]) + " ";
    }
    res += "\n";
    cout << res << endl;
}

int main(){
    vector<double> m(bin);
    vector<double> c(bin);
    vector<double> size(bin);
    vector<vector<double>> res(bin,vector<double>(5,0));

    cout << "Exact calculation\n";
    cout << "L = " << L << ", " << "bin = " << bin << "\n";

    clock_t start, finish;
    double duration;
    start = clock();

    // Actual calculation

    initialize();
    generator_sym(1);

    for(int i = 0; i < bin; i++){
        m[i] = ZM[i]/Z[i]/N;
        c[i] = (1/T[i]/T[i]/N)*(ZHH2[i]/Z[i]-(ZHH[i]/Z[i])*(ZHH[i]/Z[i]));
        cout << m[i] << " " << c[i] << endl;
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

        myfile << "idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2\n";
        for(int i = 0; i <bin; i++){
            string temp = to_string(i) + "," + to_string(T[i]) + "," + to_string(m[i]) + "," + to_string(c[i]) + ",";
            temp = temp + to_string(ZM[i]/Z[i]) + "," + to_string(ZM2[i]/Z[i]) + "," + to_string(ZM4[i]/Z[i]) + ",";
            temp = temp + to_string(ZHH[i]/Z[i]) + "," + to_string(ZHH2[i]/Z[i]) + "\n";
            // temp = temp + to_string(ZM[i]) + "," + to_string(Z[i]) + "\n";

            myfile << temp;
        }
        myfile.close();
        cout << "Save Completed: " << NewFilename << "\n";
    }

    finish = clock();
    cout << "Total step of calculation: " << count << "\n";
    cout << "Program Exit. Spent time: " << (double)(finish-start)/CLOCKS_PER_SEC << "\n";
}