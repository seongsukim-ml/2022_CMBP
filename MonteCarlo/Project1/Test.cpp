#include <iostream>
#include <boost/random.hpp>
#include <boost/math/special_functions/gamma.hpp>
// #include "C:/Users/Seongsu/header/boost_1_78_0/boost/random.hpp"
 
int main()
{
     boost::random::mt19937 dice;
     boost::random::uniform_int_distribution<> sides(1, 6);
     for (int i = 0; i < 10; i++)
     {
         std::cout << sides(dice) << " ";
     } 
     std::cout << "\n";\
     std::cout << floor(-3) << '\n';
     return 0;
}