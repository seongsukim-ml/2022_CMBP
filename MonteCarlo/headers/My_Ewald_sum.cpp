#include "My_Ewald_sum.hpp"

double pi = 3.1415926535;

ewald_ND::ewald_ND(int kd, vector<int> kL, double kalpha): d(kd){ //kL: x,y,z
    alpha = kalpha;
    L = kL;
    kappa = 2/(double)L[0]; // Multi dimension 에서는 다르게 해야할 듯..?

    // Memoization
    unsigned long long size = ((L[0])/2+1)*((L[1])/2+1);
    // cout << dist_cache_2D.max_size() << ' ' << N << ' ' << size << '\n';
    dist_cache_2D = vector<double>(size,-1);
    dist_cache_1D = vector<double>(L[1],-1);
}

// Gamma function
double ewald_ND::gam(double z){
    return boost::math::tgamma(z);
}

double ewald_ND::gaminc_up(double alpha, double z){ // not normalized
    if(alpha > 0) return boost::math::tgamma(alpha,z);
    else if(alpha == 0) return -(boost::math::expint(-z)); // if(alpha == 0) return -Ei(-z)
    else{ // negative alpha is not defined in cpp boost
        double sum_res = 0;
        double mul_res = 1;

        double s = alpha;
        int n = -floor(s);
        // cout << s << ' ' << n << ' ' << z << '\n';

        for(int i = 1; i <= n; i++){
            double temp = pow(z,s+i-1)*exp(-z);
            for(int j = 1; j <= i; j++) temp *= 1/(double)(s+j-1);
            sum_res += temp;

            mul_res *= 1/(s+i-1);
        }
        mul_res *= gaminc_up(s+n,z);
        // cout << -sum_res << ' ' << mul_res << '\n';
        return (-sum_res) + mul_res;
    }
}

double ewald_ND::pi_1(int i,int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    for(int idx = 0; idx <= 5; idx++){
        for(int ii = 0; ii <= idx; ii++){ // nx, ny are nu_x and nu_y
        int nx = ii, ny = idx - ii; // if idx = 2, (nx, ny) = (0, 2) (1, 2) (2, 2)
            for(int kk = 0; kk < 2-(nx==0); kk++)
            for(int ll = 0; ll < 2-(ny==0); ll++){
                // dis = |r_i - r_j + L|
                double dis = dist_2pt(xof(i)-xof(j),yof(i)-yof(j),nx*L[0]*(1-2*kk),ny*L[1]*(1-2*ll));
                double inv_dist = pow(dis,-alpha);
                double temp = inv_dist*gaminc_up(alpha/2,(kappa*dis)*(kappa*dis));
                res += temp;
            }
        }
    }
    return 1/gam(alpha/2)*res;
}

double ewald_ND::pi_2(int i,int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    for(int idx = 1; idx <= 5; idx++){ // idx = 0 should be skipped? Maybe.. in thesis, h=0 should be vanished
        for(int ii = 0; ii <= idx; ii++){ // nx, ny are nu_x and nu_y
        double nx = (double)ii/L[0], ny = (double)(idx - ii)/L[1]; // if idx = 2, (nx, ny) = (0, 2) (1, 2) (2, 2)
            for(int kk = 0; kk < 2-(nx==0); kk++)
            for(int ll = 0; ll < 2-(ny==0); ll++){
                double hx = nx*(1-2*kk), hy = ny*(1-2*ll);
                double k_res = cos(2*pi*vec_mul(hx,hy,xof(i)-xof(j),yof(i)-yof(j)));
                double pi_k = pi*dist_2pt(hx,hy,0,0);
                // k_res *= pow(0.5*pi_k,alpha-d);
                k_res *= 0.5*pow(pi_k,alpha-d);

                k_res *= gaminc_up(-0.5*(alpha-d),(pi_k/kappa)*(pi_k/kappa));
                // cout << idx << ' ' << hx << ' ' << hy << ' ' << k_res <<'\n';
                res += k_res;
            }
        }
    }
    return 2*pow(pi,d/2)/gam(alpha/2)/(L[0]*L[1])*res;
}

double ewald_ND::pi_ij(int i,int j){
    if(i == j) return 0;
    // else if(i > j) swap(i,j); // i < j always true
    double &res = dist_cache_2D[cache_point_2D(i,j)];
    // cout << cache_point_2D(i,j) << '\n';
    if(res == -1){
        res = pi_1(i,j,kappa) + pi_2(i,j,kappa);
    }
    return res;
    // return pi_1(i,j,kappa) + pi_2(i,j,kappa);
}

// 만약에 parallel한 성분도 고려해야한다면 p1_1_1D, pi_2_2D 모두 다 pararell 방향의 mirror 이미지도 고려할 수 있게 수정이 필요함
double ewald_ND::pi_1_1D(int i, int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    for(int idx = 0; idx <= 5; idx++){
        int nx = idx;
        for(int kk = 0; kk < 2-(nx==0); kk++){ // This loop is used for convert sign of nx
            nx *= -1;
            double dis = dist_2pt(yof(i)-yof(j),0,-nx*L[1],0);
            // double dis = pow(abs(yof(j)-yof(i)+nx*L[1]),0.5);
            // cout << dis1 << " " << dis << '\n';
            double inv_dist = pow(dis,-alpha);
            double temp = inv_dist*gaminc_up(alpha/2,(kappa*dis)*(kappa*dis));
            res += temp;
            // cout << idx <<" " << yof(i)-yof(j) <<  " " << nx << " " << temp << '\n';
        }
    }
    return 1/gam(alpha/2)*res;
}

double ewald_ND::pi_2_1D(int i,int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    int d = 1;
    for(int idx = 1; idx <= 5; idx++){ // idx = 0 should be skipped? Maybe.. in thesis, h=0 should be vanished
        double nx = (double)idx/L[1];
        for(int kk = 0; kk < 2; kk++){ // This loop is used for convert sign of nx
            nx *= -1;
            double k_res = cos(2*pi*vec_mul(nx,0,yof(i)-yof(j),0));
            // double pi_k = pi*dist_2pt(nx,0,0,0);
            double pi_k = pi*abs(nx);

            k_res *= 0.5*pow(pi_k,alpha-d);
            k_res *= gaminc_up(-0.5*(alpha-d),(pi_k/kappa)*(pi_k/kappa));
            // cout << idx << ' ' << hx << ' ' << hy << ' ' << k_res <<'\n';
            res += k_res;
            // cout << idx << " " << nx << " " << k_res << '\n';
        }
    }
    return 2*pow(pi,d/2)/gam(alpha/2)/L[1]*res;
}

double ewald_ND::pi_ij_1D(int i, int j, bool consider_parallel){
    if(i == j) return 0;
    if(consider_parallel && (xof(i) != xof(j))) return 0;
    double &res = dist_cache_1D[abs(yof(i)-yof(j))];
    if(res == -1){
        res = pi_1_1D(i,j,kappa) + pi_2_1D(i,j,kappa);

    }
    return res;
}

double ewald_ND::periodic_sum_comparison_1D(int i, int j, int idx){
    double res = 0;
    int dis = yof(j)-yof(i); // just difference
    res += pow(abs(dis),-alpha);
    for(int k = 1; k < idx; k++){
        res += pow(abs(dis+L[0]*k*1),-alpha);
        res += pow(abs(dis+L[0]*k*-1),-alpha);
    }
    return res;
}

// 1) print_pi_ij
double print_pi_ij(int lx, int ly, double alpha, int i, int j){
    ewald_ND e2d(2,vector<int>({lx,ly}),alpha);
    double res = e2d.pi_ij(i,j);
    cout << res << '\n';
    return res;
    // cout << e2d.pi_1(i,j,e2d.kappa)+e2d.pi_2(i,j,e2d.kappa) << '\n';
}
