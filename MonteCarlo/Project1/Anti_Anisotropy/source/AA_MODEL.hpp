#pragma once

#include "AA_CONF.hpp"
#include <boost/dynamic_bitset.hpp>
#include <string>
// #include "../../../headers/my_random.hpp"
// #include "../../../headers/Ewald_sum.cpp"
#include <my_random.hpp>
#include <My_Ewald_sum.hpp>

#define ss(bit) (2*bit-1) // convert boolean spin data to -1, 1 data
#define bb(spin)  (spin+1)/2

namespace model::AA{
    class BIT_MODEL : public MODEL_CONF{
        public:
            boost::dynamic_bitset<> sc;  // spin configuration
            boost::dynamic_bitset<> sb;  // background configuration of staggered board
            boost::dynamic_bitset<> bit0;
            boost::dynamic_bitset<> bit1;

            boost::dynamic_bitset<> correation;

            int total_spin;
            int staggered_spin;
            double HH; // Total Energy
            vector<double> correation_sum;

            ewald_ND e1d;
            myrnd rand;

            BIT_MODEL(model::AA::MODEL_CONF &PROFILE, ewald_ND &ewd, myrnd &rand): MODEL_CONF(PROFILE){
                // this->e1d = ewd;
                this->rand = rand;
                Initialize(ewd);
            };

            BIT_MODEL(model::AA::MODEL_CONF &PROFILE): MODEL_CONF(PROFILE){
                this->rand = myrnd();
                Initialize();
            }

            ~BIT_MODEL(){
            }

            void Initialize(ewald_ND e1d);
            void Initialize();
            void Measure();
            void Measure_fast();
            // void Correation_Measure();
    };
};