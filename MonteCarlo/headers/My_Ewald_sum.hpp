/**
 * @file Ewald_sum.cpp
 * @author Seongsu, Kim (you@domain.com)
 * @brief Ewald method header file that is made for caculation of Ising model
 * @version 0.1
 * @date 2022-03-22
 * 
 * @copyright Copyright (c) 2022
 * @see
 * Gamma : https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/sf_gamma/tgamma.html
 * Incomplete Gamma : https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
 * Exponent integral : https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/expint/expint_i.html
 */

#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <iostream>

using namespace std;

// Periodic boundary condition
class ewald_ND{
    
    public:
        int d; // Dimension of the Ewald model.. but it looks useless
        double alpha; // alpha is "d + sigma" (overall inverse power)
        double kappa;
        
        vector<int> L; // Lattice size of each edge L[0] = Lx, L[1] = Ly

        ewald_ND(){};
        ewald_ND(int kd, vector<int> kL, double kalpha);

        int xof(int i){return i%L[0];}
        int yof(int i){return i/L[0];}

        // For general dimension Ndof = N dimensional of, Ndofh = N dimensional of helper
        int Ndof(int i, int d){return Ndofh(i,d)%L[d];}
        int Ndofh(int i,int d){
            if(d == 0) return i;
            else return Ndofh(i,d-1)/L[d-1];
        }

        // 2D function
        double dist_2pt(double x1, double y1, double x2, double y2){return pow((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2),0.5);};
        double dist_ij(int i, int j){return dist_2pt(xof(i),yof(i),xof(j),yof(j));};
        double dist_ij_min(int i, int j){
            double diffx = min(abs(xof(i)-xof(j)),L[0]-abs(xof(i)-xof(j)));
            // cout << diffx;
            double diffy = min(abs(yof(i)-yof(j)),L[1]-abs(yof(i)-yof(j)));
            // cout << diffy;
            return pow(diffx*diffx + diffy*diffy,0.5);
        };
        double vec_mul(double x1, double y1, double x2, double y2){return x1*x2 + y1* y2;};
        // int cache_point_2D(int i, int j){return (N)*i-i*(i+1)/2+(j-i-1);};
        int cache_point_2D(int i, int j){
            int nx = abs(xof(i)-xof(j)); nx = 2*nx > L[0] ? L[0] - nx : nx;
            int ny = abs(yof(i)-yof(j)); ny = 2*ny > L[0] ? L[0] - ny : ny;
            return nx+ny*((L[0])/2+1);
            // return (abs(abs(xof(i)-xof(j))-(L[0])/2))+(abs(abs(yof(i)-yof(j))-(L[1])/2))*(L[0]/2);
        };

        double pi_1(int i, int j, double kappa);
        double pi_2(int i, int j, double kappa);
        double pi_ij(int i, int j);

        double pi_1_1D(int i, int j, double kappa);
        double pi_2_1D(int i, int j, double kappa);
        double pi_ij_1D(int i, int j, bool consider_parallel = true);
        double periodic_sum_comparison_1D(int i, int j, int idx);

        double gam(double z);
        double gaminc_up(double z, double alpha);

        vector<double> dist_cache_2D; // O(N^4) 메모리 크기는 좀 힘들지 않을까-> v.max_size = 1e19 나오는디..?
                                   // N= 1000까지도 가능할듯? -> 실제로 해보려니까 뻗어버리네..
                                   // dist_cache_2D 크기는 n/4만 해도 가능 한 듯
        vector<double> dist_cache_1D;
};