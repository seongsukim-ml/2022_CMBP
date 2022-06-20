#include "Cluster_Long_range.hpp"
#include "../../headers/Writer.hpp"
#include <iostream>
#include <iomanip>

/***************** Parameters 1 *****************/
const int kL          = 8;          /*Parameter: lattice size*/
const int kN          = kL*kL;
const int kBin        = 20;          /*Parameter: Change binning of temperature*/
const int kB          = 0;
const int kJ          = 1;
const double alpha    = 2+1;

// const double Tsrt = T_CRIT*(1-0.08);
// const double Tfin = T_CRIT*(1+0.08);

const double Tsrt = 5.0;
const double Tfin = 5.5;

const double isTinf = false;
const bool Random = true;

int equil_time = 1000;
const int equil_time_base = 200;
int mcs = 1e4;
/***************** Parameters 1 *****************/

typedef Cluster_LR_2D Model;

// Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
#ifdef _WIN32
static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs);
#endif
#ifdef linux
static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs);
#endif

// clock used to measure time
clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << Model::Name() + "Algorithm\n";
    cout << "Radnomness test(seed): " << seed << '\n';
    cout << "L = " << kL << ", " << "bin = " << kBin << ", Start T at " << Tat << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << "\n";
    cout << "--index--||---Temp----||EQ:sig------HH----------||magnetization---specific heat||Fliped Step------Total Step------" << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << fixed <<setprecision(6);
    
    // Show +/- sign
    cout << showpos;


    __start__ = clock();
}

void Farewell(int N = 0){
    __finish__ = clock();
    if(N)
        cout << "\nProgram Abonormally Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    else
        cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    cout << "-------------------------------------------------------------------------------------------\n";
}

// Event handler that notify the spent time when program abonormally stop
void handler(int A){
    cout << endl;
    Farewell(1);
    exit(A);
}

// arguments list that helps to pass the args to model
vector<double> args = {kL,kBin,kB,kJ,alpha,Tsrt,Tfin,isTinf};

int main(){
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    Greetings();

    for(int gg = 0; gg < 1; gg++){
        Model model = Model(args);
        cout << model.J_tot << '\n';

        double MM, HH;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        cout << "blah" << '\n';
        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            double t = 0;

            for(int j = 0; j < kL; j++){
            for(int k = 0; k < kL; k++){
                cout << model.sc[j*kL+k];
            }
            cout << '\n';
            }

            for(int mm = 0; mm < kN*(kN-1)/2; mm++){
                if(model.Walker_Table_A[mm] != -1)
                    t += 1;
                else
                    t += model.Walker_Table_P[mm];
            }
            cout << "total Prob is " << t << endl;
            model.Calculate();

            for(int j = 0; j < kL; j++){
            for(int k = 0; k < kL; k++){
                cout << model.sc[j*kL+k];
            }
            cout << '\n';
            }
        }
    }
    Farewell();
}