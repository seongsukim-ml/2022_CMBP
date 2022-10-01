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

#ifndef TYPE_DEFINE
#define TYPE_DEFINE
typedef INT1 int32_t
typedef INT2 int64_t
typedef INT8 int8_t
typedef FLOAT1 long double
typedef FLOAT2 long double
#endif

/* Note */
// In 2D Ising model, the critical temperture is 2/log(1+sqrt(2)), that is about 2.269.
// const double T_CRIT = 2.269185;

namespace model::AA{
    struct MODEL_CONF{
        // The configuration of Anti-Anistorpy model
        public:
            const INT1 Lx; // Short range lattice length (Time domain)
            const INT1 Ly; // Long range lattice length (Spatioal domain)
            const INT1 N;  // Total number of spin, that is (N = Lx * Ly).

            const INT1 Bin;

            const FLOAT1 B;  // the z directional external magnetic field. Usually 0
            const FLOAT1 Jx; // the x directional interaction coeffcient. Usually 1
            const FLOAT1 Jy; // the y directional interaction coeffcient. Usually -1i

            const FLOAT1 alpha; // the exponential coefficient of long range interaction. (z = 1 + alpha)

            FLOAT1 Tsrt; // Set the starting temperature of the model calculation. It can be 0, but in that case it skip the 0 temp.
            FLOAT1 Tfin; // Set the finishing temperature of the model calculation. It cannot be 0.

            bool isTinf; // Determin the strating spin configuration as random(Infinite T) or not(0 T).
        MODEL_CONF(INT1 _Lx, INT1 _Ly, INT1 _Bin, FLOAT1 _B, FLOAT1 _Jx, FLOAT1 _Jy, FLOAT1 _alpha, FLOAT1 _Tsrt, FLOAT1 _Tfin, bool _isTinf)
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

