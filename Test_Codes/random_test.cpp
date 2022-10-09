#include <random>
#include <iostream>

int main(){
    std::random_device rd;
    std::cout << int(rd()) << '\n';
}