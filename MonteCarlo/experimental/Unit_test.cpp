#include<iostream>
#include<vector>
#include<bitset>
#include<string>
#include<math.h>

using namespace std;

#define L 3 // Number of Lattice
#define XLL L
#define XLL L*L

int partition_function(int H,int T,int J=1){
    // k = 1.38E-23 # boltzman constant
    int k = 1;
    return exp(-1*J*H/(k*T));
}

int main(){
    unsigned int s = 1 << (L*L-1); // symmetry 때문에 마지막 칸 고정
    double p =0, a = 1;
    int i, sum;
    for(i = 0; i < ((1 << (L*L-1))); i++){
        sum = 0;
        for(int j = 0; j < L*L; j++){
            sum += (s >> j & 1);
        }
        s += 1;
        cout << s << " " << sum << endl;
    }
}