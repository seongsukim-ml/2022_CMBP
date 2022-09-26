#include <iostream>
#include "../Model_Driver.hpp"
#include "../my_random.hpp"
#include "../../Project1/Anti_Anisotropy/AA_Metropolis_BLAS/AA_CONF.hpp"

int main(){
    model::AA::MODEL_CONF PROFILE = {1,1,1,1,1,1,1,1,1,1};
    myrnd rnd;
    int a = 1;

    model::Model_Driver<model::AA::MODEL_CONF,int,myrnd> MD(PROFILE,a,rnd);
}