#include "AA_Metropolis_hp_heatbath2.hpp"
#include "../../../headers/Writer.hpp"
#include <iostream>
#include <iomanip>
#define base_L24

#ifdef base_L128
    /***************** (Test) Parameters 1 *****************/
    int kLx         = 128;          /*Parameter: lattice size*/
    int kLy         = 128;

    int kN          = kLx*kLy;
    int kBin        = 11;          /*Parameter: Change binning of temperature*/

    double kB       = 0;
    double kJx      = 1;
    double kJy      = -1;

    double alpha    = 1.5;

    double Tsrt = 1.78;
    double Tfin = 1.80;

    double isTinf = false;
    bool Random = false;

    int equil_time_base = 1e3;
    int equil_time = equil_time_base;
    int mcs = 1e4;
    /***************** (Test) Parameters 1 *****************/
#endif

#ifdef base_L24
    /***************** (Test) Parameters 1 *****************/
    int kLx         = 2;          /*Parameter: lattice size*/
    int kLy         = 4;

    int kN          = kLx*kLy;
    int kBin        = 2;          /*Parameter: Change binning of temperature*/

    double kB       = 0;
    double kJx      = 1;
    double kJy      = -1;

    double alpha    = 2.0;
    // double alpha    = 100;

    double Tsrt = 1.91;
    double Tfin = 1.92;
    // double Tsrt = 2.26918;
    // double Tfin = 2.26918;

    double isTinf = false;
    bool Random = false;

    int equil_time_base = 1e3;
    int equil_time = equil_time_base;
    int mcs = 1e5;
    /***************** (Test) Parameters 1 *****************/
#endif

typedef AA_Metropolis Model;

// clock used to measure time
clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << Model::Name() + "Algorithm\n";
    cout << "Radnomness test(seed): " << seed << '\n';
    cout << "L = " << kLx << "," << kLy << ", " << "bin = " << kBin << ", Start T at " << Tat <<  "\n";
    cout << "N = " << kN << ", " << "alpha = " << alpha << ", mcs " << mcs <<  "\n";
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