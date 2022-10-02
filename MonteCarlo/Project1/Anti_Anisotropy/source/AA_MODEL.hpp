#pragma once

#include "AA_CONF.hpp"
#include <boost/dynamic_bitset.hpp>
#include <string>
// #include "../../../headers/my_random.hpp"
// #include "../../../headers/Ewald_sum.cpp"
#include <my_random.hpp>
#include <my_ewald_sum.hpp>

#define ss(bit) (2*bit-1) // convert boolean spin data to -1, 1 data
#define bb(spin)  (spin+1)/2

namespace model::AA{
    class BIT_MODEL : public MODEL_CONF{
        public:
            boost::dynamic_bitset<> sc;  // spin configuration
            boost::dynamic_bitset<> sb;  // background configuration of staggered board
            FLOAT1 cur_beta;

            INT1 total_spin;
            INT1 staggered_spin;
            FLOAT2 HH; // Total Energy

            ewald_ND e1d;
            myrnd rand;

            BIT_MODEL(const model::AA::MODEL_CONF &PROFILE, const ewald_ND &ewd, const myrnd &rand);
            BIT_MODEL(const model::AA::MODEL_CONF &PROFILE);
            // ~BIT_MODEL();

            void Initialize();
            void Measure();
            void SetTemp(FLOAT2 beta);
            virtual void Measure_fast();
            FLOAT2 HH1();
            INT1 MM1();
    };
};