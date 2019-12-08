// 주어진 nonlinear model 에 data를 fitting한다.
// Levenberg-Marquardt Method를 사용하여 점진적으로 최적값에 도달하도록 구현한다.

#include <stdio.h>
#include <vector>

#include "../recipes_cpp/ludcmp.cpp"
#include "../recipes_cpp/lubksb.cpp"

using namespace std;

vector<DP> v[4];

void Input() {
    DP a, b, c, d;
    int i=0;

    FILE *fp=fopen("../data.dat", "r");

    while(fscanf(fp, "%lf %lf %lf %lf", &a, &b, &c, &d) != EOF)
    {
        v[0].push_back(a);
        v[1].push_back(b);
        v[2].push_back(c);
        v[3].push_back(d);
    }

    fclose(fp);
}

void NonlinearFitting() {

    DP a11, a12, a13, a31, a32;
    DP l, ll;

    // Guess a
    Mat_DP mat(5,5);
    Vec_INT indx(5);
    DP d;
    Vec_DP b(5);

    // data1 fitting

    // 초기화
    indx[0]=indx[1]=indx[2]=indx[3]=indx[4]=0;


    // Fitting data to a straight line
    for(int i=0; i<5; i++) {
        mat[i][0] = v[0][i];
        mat[i][1] = v[1][i];
        mat[i][2] = 1;
        mat[i][3] = -1 * v[0][i] * v[3][i];
        mat[i][4] = -1 * v[1][i] * v[3][i];

        b[i] = v[3][i];

    }

    NR::ludcmp(mat,indx,d);
    NR::lubksb(mat,indx,b);

    a11 = b[0]; a12 = b[1]; a13 = b[2]; a31 = b[3]; a32 = b[4];
    l = 0.0;
    for (int i = 0; i < v[0].size(); i++) {
        DP k = v[2][i] - (a11 * v[0][i] + a12 * v[1][i] + a13) / (a31 * v[0][i] + a32 * v[1][i] + 1);
        l += k * k;
    }

    // Gradient Update loop
    int cnt=0;
    double lamda = 0.001;
    while (cnt < 100) {

        // H' 초기화 (lamda * I)
        for (int i = 0; i < 5; i++) {
            b[i] = 0;
            indx[i] = 0;
            for (int j = i+1; j < 5; j++)
                mat[i][j] = mat[j][i] = 0;
            mat[i][i] = lamda;
        }
        // H' 구하기 / X2 Gradient 구하기
        for (int i = 0; i < v[0].size(); i++) {

            DP p, q;
            p = a11 * v[0][i] + a12 * v[1][i] + a13;
            q = a31 * v[0][i] + a32 * v[1][i] + 1;

            mat[0][3] = mat[3][0] += -2*v[0][i]*v[0][i]/(q*q);
            mat[0][4] = mat[4][0] += -2*v[0][i]*v[1][i]/(q*q);

            mat[1][3] = mat[3][1] += -2*v[0][i]*v[1][i]/(q*q);
            mat[1][4] = mat[4][1] += -2*v[1][i]*v[1][i]/(q*q);

            mat[2][3] = mat[3][2] += -2*v[0][i]/(q*q);
            mat[2][4] = mat[4][2] += -2*v[1][i]/(q*q);

            mat[3][4] = mat[4][3] += 4*v[0][i]*v[1][i]*p/(q*q*q);

            mat[3][3] += 4*v[0][i]*v[0][i]*p/(q*q*q);
            mat[4][4] += 4*v[1][i]*v[1][i]*p/(q*q*q);

            b[0] += 2*(v[2][i] - p/q) * v[0][i]/q;
            b[1] += 2*(v[2][i] - p/q) * v[1][i]/q;
            b[2] += 2*(v[2][i] - p/q) / q;
            b[3] += -2*(v[2][i] - p/q) * v[0][i] * p/q;
            b[4] += -2*(v[2][i] - p/q) * v[1][i] * p/q;
        }

        NR::ludcmp(mat,indx,d);
        NR::lubksb(mat,indx,b);

        ll = 0.0;
        for (int i = 0; i < v[0].size(); i++) {
            DP k = v[2][i] - ((a11+b[0]) * v[0][i] + (a12+b[1]) * v[1][i] + a13 + b[2]) / ((a31+b[3]) * v[0][i] + (a32+b[4]) * v[1][i] + 1);
            ll += k * k;
        }
        // 업데이트 여부 판단
        if(ll < l) {
            a11+=b[0];
            a12+=b[1];
            a13+=b[2];
            a31+=b[3];
            a32+=b[4];

            l=ll;
            lamda/=10;
        }
        else if(ll > l)
            lamda*=10;
        cnt++;

    }
    cout << "a11 : " << a11 << endl;
    cout << "a12 : " << a12 << endl;
    cout << "a13 : " << a13 << endl;
    cout << "a31 : " << a31 << endl;
    cout << "a32 : " << a32 << endl;

    Mat_DP mat2(3,3);
    Vec_INT indx2(3);
    DP d2;
    Vec_DP b2(3);

    for (int i = 0; i < 3; i++) {
        b2[i] = 0;
        indx2[i] = 0;
        for (int j = i+1; j < 3; j++)
            mat2[i][j] = mat2[j][i] = 0;
        mat2[i][i] = 0;
    }
    // Fitting data to a straight line
    for(int i=0; i<v[0].size(); i++) {
        mat2[0][0] += v[0][i];
        mat2[0][1] += v[1][i];
        mat2[0][2] += 1;

        mat2[1][0] += v[0][i] * v[0][i];
        mat2[1][1] += v[1][i] * v[0][i];
        mat2[1][2] += v[0][i];

        mat2[2][0] += v[0][i] * v[1][i];
        mat2[2][1] += v[1][i] * v[1][i];
        mat2[2][2] += v[1][i];

        DP q = a31 * v[0][i] + a32 * v[1][i] + 1;

        b2[0] += v[3][i] * q;
        b2[1] += v[3][i] * q * v[0][i];
        b2[2] += v[3][i] * q * v[1][i];
    }
    // LU 분해를 통한 해 구하기
    NR::ludcmp(mat2,indx2,d2);
    NR::lubksb(mat2,indx2,b2);

    cout << "a21 : " << b2[0] << endl;
    cout << "a22 : " << b2[1] << endl;
    cout << "a23 : " << b2[2] << endl;
}

int main() {

    Input();
    NonlinearFitting();

    return 0;
}