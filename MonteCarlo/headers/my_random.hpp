#pragma once
#include <random>

class myrnd{
    public:
        // If you cannot use random device, use below code.
        long unsigned int seed;
        std::mt19937 gen; // Standard mersenne_twister_engine seeded with time()

        // If you can use random device, use below code.
        // std::random_device rd;  // Will be used to obtain a seed for the random number engine
        // std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

        // Random number generator
        std::uniform_real_distribution<> udis; // uniform random number generator, ex) udis(gen) -> [0,1]
        std::bernoulli_distribution bdis;           // bernoilli random number generator of probability 0.5, ex) bdis(gen)
    public:
        myrnd(){
            // If you cannot use random device, use below code
            this-> seed = static_cast<long unsigned int>(time(0));
            this-> gen  = std::mt19937(seed); // Standard mersenne_twister_engine seeded with time()

            this-> udis = std::uniform_real_distribution<> (0.0, 1.0); // uniform random number generator, ex) udis(gen) -> [0,1]
            this-> bdis = std::bernoulli_distribution (0.5);           // bernoilli random number generator of probability 0.5, ex) bdis(gen)
        };
        double bern(){
            return this->bdis(this->gen);
        };
        double unif(){
            return this->udis(this->gen);
        };
};