// Gaussian Probability Density Fuction을 따르는 난수 생성기를 통해 행렬을 구성하라.
// 행렬의 eigenvector 및 eigenvalue를 구하여라
#include <iostream>
#include <algorithm>
#include <time.h>
#include <vector>
#include "../recipes_cpp/gasdev.cpp"
#include "../recipes_cpp/ran1.cpp"
#include "../recipes_cpp/jacobi.cpp"
#include "../recipes_cpp/eigsrt.cpp"

using namespace std;

int main() {

    Mat_IO_DP a(11,11);
    Vec_O_DP d(11);
    Mat_O_DP v(11,11);
    int nrot;

    int x=clock()*clock()%10000007%100000;

    for(int i=0; i<11; i++) {
        for(int j=i; j<11; j++ )
            a[i][j] = a[j][i] = NR::gasdev(x);
    }
    cout <<"A" << endl;

    for(int i=0; i<11; i++) {
        for(int j=0; j<11; j++)
            cout << a[i][j] << " ";
        cout << endl;
    }
    cout << endl;

    NR::jacobi(a,d,v,nrot);
    NR::eigsrt(d,v);


    for(int i=0; i<11; i++) {
        cout << "고유치 : " << d[i] << endl;
        cout << "고유벡터 : (";
        for(int j=0; j<10; j++)
            cout << v[i][j] << ",";
        cout << v[i][10] << ")" << endl;
    }
    cout << endl;

    DP sum;
    for(int i=0; i<11; i++) {
        for(int j=i+1; j<11; j++) {
            sum = 0.0;
            for(int k=0; k<11; k++) {
                sum+=v[i][k]*v[j][k];
            }
            cout << i+1 << "번째 벡터 * " << j+1 << "번째 벡터 = " << sum << endl;
        }
    }

    return 0;
}