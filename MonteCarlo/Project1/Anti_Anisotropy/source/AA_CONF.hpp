#pragma once
// File & IO System
#include <iostream>
#include <bitset>
#include <vector>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <string>

#include <signal.h>
#include <cstdlib>

#include "../../../headers/my_random.hpp"
// Random number generator from my_random

// #define ss(spin) (2*spin-1) // convert boolean spin data to -1, 1 data
// #define bb(bit)  (bit+1)/2

/* Note */
// In 2D Ising model, the critical temperture is 2/log(1+sqrt(2)), that is about 2.269.
// const double T_CRIT = 2.269185;

namespace model::AA{
    struct MODEL_CONF{
        // The configuration of Anti-Anistorpy model
        public:
            const int Lx; // Short range interaction (Time domain)
            const int Ly; // Long range interaction (Spatioal domain)
            const int N;  // Total number of spin, that is (N = Lx * Ly).

            const int Bin;

            const double B;  // the z directional external magnetic field. Usually 0
            const double Jx; // the x directional interaction coeffcient. Usually 1
            const double Jy; // the y directional interaction coeffcient. Usually -1i

            const double alpha; // the exponential coefficient of long range interaction. (z = 1 + alpha)

            double Tsrt; // Set the starting temperature of the model calculation. It can be 0, but in that case it skip the 0 temp.
            double Tfin; // Set the finishing temperature of the model calculation. It cannot be 0.

            bool isTinf; // Determin the strating spin configuration as random(Infinite T) or not(0 T).
        MODEL_CONF(int _Lx, int _Ly, int _Bin, double _B, double _Jx, double _Jy, double _alpha, double _Tsrt, double _Tfin, bool _isTinf)
        : Lx(_Lx), Ly(_Ly), N(_Lx*_Ly), Bin(_Bin), B(_B), Jx(_Jx), Jy(_Jy), alpha(_alpha), isTinf(_isTinf), Tsrt(_Tsrt), Tfin(_Tfin)
        {};

        MODEL_CONF(const MODEL_CONF &PROFILE)
        : Lx(PROFILE.Lx), Ly(PROFILE.Ly), N(PROFILE.N), Bin(PROFILE.Bin), B(PROFILE.B), Jx(PROFILE.Jx), Jy(PROFILE.Jy), alpha(PROFILE.alpha),
        Tsrt(PROFILE.Tsrt), Tfin(PROFILE.Tfin), isTinf(PROFILE.isTinf)
        {}

        void print_config(){
            std::string Tat = isTinf ? "inf" : "0";

            std::cout << "L = " << Lx << "," << Ly << ", " << "bin = " << Bin << ", Start T at " << Tat <<  "\n";
            std::cout << "N = " << N << ", " << "alpha = " << alpha << "\n";
        }
    };
}

