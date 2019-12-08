// 균등 분포 및 가우시안 분포의 난수 생성기를 구현해 본다.

#include <iostream>
#include <algorithm>
#include <time.h>
#include <vector>
#include "../recipes_cpp/ran1.cpp"
#include "../recipes_cpp/gasdev.cpp"

using namespace std;

vector<DP> x_unif1, x_unif2, x_unif3, x_unif4;
vector<DP> x_gauss1, x_gauss2, x_gauss3, x_gauss4;


int main() {
    DP a = -3.0, b = 2.0, m = -0.5, s = 1.5;
    int x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=1000; i++) {
        x_unif1.push_back(NR::ran1(x)*(b-a)+a); // 1000 Samples, Uniform distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=100; i++) {
        x_unif2.push_back(NR::ran1(x)*(b-a)+a); // 100 Samples, Uniform distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=10000; i++) {
        x_unif3.push_back(NR::ran1(x)*(b-a)+a); // 10000 Samples, Uniform distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=30000; i++) {
        x_unif4.push_back(NR::ran1(x)*(b-a)+a); // 100000 Samples, Uniform distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=1000; i++) {
        DP temp = NR::gasdev(x)*s+m;
        if(temp <= -3.0*s || temp >= 3.0*s) {
            i--;
            continue;
        }
        x_gauss1.push_back(temp); // 1000 Samples, Gaussian distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=100; i++) {
        DP temp = NR::gasdev(x)*s+m;
        if(temp <= -3.0*s || temp >= 3.0*s) {
            i--;
            continue;
        }
        x_gauss2.push_back(temp); // 100 Samples, Gaussian distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=10000; i++) {
        DP temp = NR::gasdev(x)*s+m;
        if(temp <= -3.0*s || temp >= 3.0*s) {
            i--;
            continue;
        }
        x_gauss3.push_back(temp); // 10000 Samples, Gaussian distribution
    }

    x=clock()*clock()%10000007%100000; // 시간을 이용한 난수 초기값 생성
    for(int i=1; i<=30000; i++) {
        DP temp = NR::gasdev(x)*s+m;
        if(temp <= -3.0*s || temp >= 3.0*s) {
            i--;
            continue;
        }
        x_gauss4.push_back(temp); // 100000 Samples, Gaussian distribution
    }

    for(int i=0; i<x_unif1.size(); i++) {
        cout << x_unif1[i] << endl;
    }
    return 0;
}