// 다양한 root finding 함수를 사용하여 주어진 함수의 solution을 구해본다.

#include <stdio.h>
#include "../other_c/nrutil.c"

#include "../recipes_c/bessj0.c" // 베슬 함수 (J0)
#include "../recipes_c/bessj1.c" // 베슬 함수의 도함수 (J1)

// zbrak 함수는 어떠한 제약 조건 없이 실행 가능하다.
#include "../recipes_c/zbrak.c" // bracketing method (단위별로 자른뒤 하나씩 증가하면서 근이 존재하는 구간을 추적한다.


// rtbis / rtflsp 함수는 구간 내 근이 존재해야 실행가능한 함수들이다. 그렇기에 시작 시 중간값 정리를 통해 근의 존재여부 부터 판별한다.
#include "../recipes_c/rtbis.c" // 전체 구간을 절반씩 줄여나가면서 근이 존재하는 구간을 추적한다.
#include "../recipes_c/rtflsp.c" // 양 극점을 잇는 직선의 x절편을 x1이라고 하면, (x1,f(x1))을 새로운 극점으로 잡는다. 나머지 극점은 원래의 극점 중 새 극점과 부호가 다른 극점으로 한다.

#include "../recipes_c/rtnewt.c" // 뉴턴 랩슨 메소드, 접선의 x절편을 x1이라고 하면, (x1,f(x1))에서 새로운 접선을 긋는 작업을 반복한다.
#include "../recipes_c/rtsafe.c" // 뉴턴 렙슨 메소드와 비셉션 메소드를 섞어 놓았다.

#include "../recipes_c/rtsec.c" // P0, P1을 양 극점으로 초기값을 잡는다. 자연수 n에 대하여 pn-1, pn을 잇는 직선의 x절편을 x1이라고 하면, (x1,f(x1))을 새로운 pn+1로 잡고 이 과정을 반복한다.
                   // x절편이 존재하지 않거나 무수히 존재하는 경우(직선의 기울기가 0이면) 실행을 중단한다.

#include "muller.c" // Secant가 직선을 그린다면, muller's method는 세점을 지나는 포물선을 그린다.



// 4개의 함수들


float F1(float x) { // 1번 문제
    return 10.0*exp(-1.0*x)*sin(2*M_PI*x) - 2.0;
}
float F1_Prime(float x) {
    return -10.0*exp(-1.0*x)*sin(2*M_PI*x) + 20.0*M_PI*exp(-1.0*x)*cos(2*M_PI*x);
}

void F1_Pointer(float rtn, float *f, float *df) {
    *f = F1(rtn);
    *df = F1_Prime(rtn);
}


float F2(float x) { // 중근을 가지는 녀석이므로 정리해서 넣어준다.
    return x - exp(-1.0*x);
}

float F2_Prime(float x) {
    return 1.0 - exp(-1.0*x);
}

void F2_Pointer(float rtn, float *f, float *df) {
    *f = F2(rtn);
    *df = F2_Prime(rtn);
}

float F3(float x) { // 3번 문제
    return cos(x+sqrt(2.0)) + x*(x/2.0 + sqrt(2.0));
}

float F3_Prime(float x) {
    return -1.0*sin(x+sqrt(2.0)) + x +sqrt(2.0);
}

void F3_Pointer(float rtn, float *f, float *df) {
    *f = F3(rtn);
    *df = F3_Prime(rtn);
}

float F4(float x) { // 자작 문제
    return sin(x*x) - cos(x);
}

float F4_Prime(float x) {
    return 2*x*cos(x*x) + sin(x);
}

void F4_Pointer(float rtn, float *f, float *df) {
    *f = F4(rtn);
    *df = F4_Prime(rtn);
}


void Func(float rtn, float *f, float *df) {
    *f = bessj0(rtn);
    *df = (-1.0)*bessj1(rtn);
}

float xb1[20], xb2[20]; // 0.0으로 초기화

int main() {

    int n = 9000000, *nb;
    *nb = 20;

    zbrak(bessj0,1.0,10.0,n,xb1,xb2,nb);
    printf("zbrak  결과 : \n");
    for(int i=1; i<=*nb; i++) {
        printf("근%d : %f\n", i, xb2[i]);
    }


    printf("rtbis  결과 : %f\n", rtbis(bessj0, 1.0,10.0,0.000001));
    printf("rtflsp 결과 : %f\n", rtflsp(bessj0, 1.0,10.0,0.000001));

    void (*FuncPtr)(float, float *, float *) = Func; // 베슬 함수값 및 미분값을 넘겨주는 함수 Func의 포인터

    printf("rtnewt 결과 : %f\n", rtnewt(FuncPtr, 1.0,10.0,0.000001));
    printf("rtsafe 결과 : %f\n", rtsafe(FuncPtr, 1.0,10.0,0.000001));

    printf("rtsec  결과 : %f\n", rtsec(bessj0, 1.0,10.0,0.000001));

    printf("muller 결과 : %f\n", muller(bessj0, 1.0,10.0,0.000001));


    FuncPtr = F1_Pointer;
    printf("F1 rtsafe 결과 : %f\n", rtsafe(FuncPtr, 0.1,1.0,0.000001));
    FuncPtr = F2_Pointer;
    printf("F2 rtsafe 결과 : %f\n", rtsafe(FuncPtr, 0.0,1.0,0.000001));
    // F3은 양 극점에서 함수값이 같이 부호이므로 에러가 발생하였다.
//    FuncPtr = F3_Pointer;
//    printf("F3 rtsafe 결과 : %f\n", rtsafe(FuncPtr, -2.0,-1.0,0.000001));
    FuncPtr = F4_Pointer;
    printf("F4 rtsafe 결과 : %f\n", rtsafe(FuncPtr, 0.0,3.0,0.000001));


    return 0;
}