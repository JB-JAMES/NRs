// double와 float의 Machine Accuracy를 직접 구해본다.

#include <iostream>

using namespace std;

template <typename DP>
int get_eps(DP &eps)
{
    int n;
    DP one, temp, two;
    n = 0;
    one = DP(1);
    two = DP(2);
    eps = DP(1);

    while(one + (eps/two) != one) {
        n++;
        eps /= two;
    }
    return n;
}


int main() {
    float eps_float;
    double eps_double;
    int n_float, n_double;

    n_float = get_eps(eps_float);
    n_double = get_eps(eps_double);

    cout << "Method 2" << endl;
    cout << "machine accuracy of “float” : " << eps_float << endl;        // 부동 소수점 정밀도
    cout << "minimum n : " << n_float << endl;
    cout << "machine accuracy of “double” : " << eps_double << endl;        // 부동 소수점 정밀도
    cout << "minimum n : " << n_double << endl;

    // cout << "<" << b << "," << temp << "," << itemp << ">" << endl;
    return 0;
}