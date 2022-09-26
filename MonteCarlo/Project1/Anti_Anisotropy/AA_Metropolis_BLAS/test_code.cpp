#include "AA_Metropolis_bit.hpp"
#include <iostream>

using namespace std;

namespace PARAM {
    static int Lx          = 24;          /*Parameter: lattice size*/
    static int Ly          = 24;

    static int N           = Lx*Ly;
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

extern std::clock_t greet::__start__,greet::__finish__;

int main(){
    signal(SIGSEGV, &greet::handler); // SIGSEGV : Segmentation Fault
    signal(SIGINT , &greet::handler); // SIGINT  : Signal interupt

}