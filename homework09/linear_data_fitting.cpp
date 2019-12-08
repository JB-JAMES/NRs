// 주어진 데이터 쌍 간의 선형 변환식을 linear data fitting을 통해 구하여라.

#include <stdio.h>
#include <vector>

#include "../recipes_cpp/ludcmp.cpp"
#include "../recipes_cpp/lubksb.cpp"

using namespace std;

vector<DP> v[4];
// [0] : x, [1] : y, [2] : x', [3] : y'

void Input() {
    double a, b, c, d;

    FILE *fp=fopen("../fitdata/fitdata3.dat", "r");

    while(fscanf(fp, "%lf %lf %lf %lf", &a, &b, &c, &d) != EOF)
    {
        v[0].push_back(a);
        v[1].push_back(b);
        v[2].push_back(c);
        v[3].push_back(d);
    }

    fclose(fp);
}

void LinearFitting() {
    Mat_DP mat(3,3);
    Vec_INT indx(3);
    DP d;
    Vec_DP b(3);

    // data1 fitting

    // 초기화
    mat[0][0]=mat[0][1]=mat[0][2]=mat[1][0]=mat[1][1]=mat[1][2]=mat[2][0]=mat[2][1]=mat[2][2] = 0;
    indx[0]=indx[1]=indx[2]=0;
    b[0]=b[1]=b[2]=0;

    // Fitting data to a straight line
    for(int i=0; i<v[0].size(); i++) {
        mat[0][0] += v[0][i];
        mat[0][1] += v[1][i];
        mat[0][2] += 1;

        mat[1][0] += v[0][i] * v[0][i];
        mat[1][1] += v[1][i] * v[0][i];
        mat[1][2] += v[0][i];

        mat[2][0] += v[0][i] * v[1][i];
        mat[2][1] += v[1][i] * v[1][i];
        mat[2][2] += v[1][i];

        b[0] += v[2][i];
        b[1] += v[2][i] * v[0][i];
        b[2] += v[2][i] * v[1][i];
    }
    // LU 분해를 통한 해 구하기
    NR::ludcmp(mat,indx,d);
    NR::lubksb(mat,indx,b);

    cout << "a1 : " << b[0] << endl;
    cout << "a2 : " << b[1] << endl;
    cout << "a3 : " << b[2] << endl;

    // 초기화
    mat[0][0]=mat[0][1]=mat[0][2]=mat[1][0]=mat[1][1]=mat[1][2]=mat[2][0]=mat[2][1]=mat[2][2] = 0;
    indx[0]=indx[1]=indx[2]=0;
    b[0]=b[1]=b[2]=0;
    // Fitting data to a straight line
    for(int i=0; i<v[0].size(); i++) {
        mat[0][0] += v[0][i];
        mat[0][1] += v[1][i];
        mat[0][2] += 1;

        mat[1][0] += v[0][i] * v[0][i];
        mat[1][1] += v[1][i] * v[0][i];
        mat[1][2] += v[0][i];

        mat[2][0] += v[0][i] * v[1][i];
        mat[2][1] += v[1][i] * v[1][i];
        mat[2][2] += v[1][i];

        b[0] += v[3][i];
        b[1] += v[3][i] * v[0][i];
        b[2] += v[3][i] * v[1][i];
    }
    // LU 분해를 통한 해 구하기
    NR::ludcmp(mat,indx,d);
    NR::lubksb(mat,indx,b);

    cout << "a4 : " << b[0] << endl;
    cout << "a5 : " << b[1] << endl;
    cout << "a6 : " << b[2] << endl;

}

int main() {

    Input();
    LinearFitting();

    return 0;
}