#include "../AA_CONF.hpp"
#include "../AA_MODEL.hpp"
#include <iostream>

using namespace std;

namespace PARAM {
    static int Lx          = 24;          /*Parameter: lattice size*/
    static int Ly          = 24;

    static int N           = Lx*Ly;       /*Automatically setting*/
    static int Bin         = 40;          /*Parameter: Change binning of temperature*/

    static double B        = 0;
    static double Jx       = 1;
    static double Jy       = -1;
    static double alpha    = 3;

    static double Tsrt = 1.5;
    static double Tfin = 4;

    static double isTinf = false;
    static bool Random = false;

    static int equil_time_base = 1e4;
    static int equil_time = equil_time_base;
    static int mcs = 1e4;
}

void handler(){};

int main(){
    // signal(SIGSEGV, &greet::handler); // SIGSEGV : Segmentation Fault
    // signal(SIGINT , &greet::handler); // SIGINT  : Signal interupt

    model::AA::MODEL_CONF PROFILE = {
        PARAM::Lx,PARAM::Ly,PARAM::Bin,
        PARAM::B,PARAM::Jx,PARAM::Jy,PARAM::alpha,
        PARAM::Tsrt,PARAM::Tfin,PARAM::isTinf
        };
    model::AA::BIT_MODEL model(PROFILE);
    cout << model.total_spin << '\n';
    cout << model.HH << '\n';
}