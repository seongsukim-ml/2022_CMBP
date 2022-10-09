#include <random>
#include <iostream>
#include <math.h>
#include <time.h>

long double ls[100000000];

int main(){
    std::ios::sync_with_stdio(false); std::cin.tie(NULL); std::cout.tie(NULL);
    clock_t start = clock();
    long double k = 0, kk;
    for(int i = 0; i < 10e7; i++){
        k -= 0.000001;
        ls[i] = expl(k);
        // std::cout << k << ' ' << expl(k) << '\n';
    }
    clock_t finish = clock();
    std::cout << (double)(finish - start)/CLOCKS_PER_SEC << '\n';
    start = clock();
    k = 0;
    for(int i = 0; i < 10e7; i++){
        k -= 0.000001;
        kk = ls[i];
        // std::cout << k << ' ' << expl(k) << '\n';
    }
    finish = clock();
    std::cout << (double)(finish - start)/CLOCKS_PER_SEC << '\n';

}