// 행렬을 이용한 linear equation의 solution 을 3가지의 method를 통해 시도해본다.

#include <stdio.h>
#include <stdlib.h>

#include "../recipes_c/gaussj.c"

#include "../recipes_c/ludcmp.c"
#include "../recipes_c/lubksb.c"

#include "../recipes_c/svdcmp.c"
#include "../recipes_c/svbksb.c"

#include "../recipes_c/mprove.c"

// Ax = c
int trash;
int n1, n2, n3; // 행렬 'A'의 크기
float **a1, **a2, **a3; // 행렬 'A'
float **c1, **c2, **c3; // Right Hand Vector 'c'

void Input() { // 파일 입력
    FILE *fp=fopen("../lineq1.dat", "r");

    fscanf(fp, "%d %d", &n1, &trash);

    a1 = (float**)malloc(sizeof(float*)*(n1+1));
    c1 = (float**)malloc(sizeof(float*)*(n1+1));

    for(int i=1; i<=n1; i++)
    {
        a1[i] = (float*)malloc(sizeof(float)*(n1+1));
        for (int j = 1; j <= n1; j++)
            fscanf(fp, "%f", &a1[i][j]);
    }

    for(int i=1; i<=n1; i++) {
        c1[i] = (float*)malloc(sizeof(float)*(2));
        fscanf(fp, "%f", &c1[i][1]);
    }

    fclose(fp);

    fp=fopen("../lineq2.dat", "r");
    fscanf(fp, "%d %d", &n2, &trash);
    a2 = (float**)malloc(sizeof(float*)*(n2+1));
    c2 = (float**)malloc(sizeof(float*)*(n2+1));

    for(int i=1; i<=n2; i++)
    {
        a2[i] = (float*)malloc(sizeof(float)*(n2+1));
        for (int j = 1; j <= n2; j++)
            fscanf(fp, "%f", &a2[i][j]);
    }

    for(int i=1; i<=n2; i++) {
        c2[i] = (float*)malloc(sizeof(float)*(2));
        fscanf(fp, "%f", &c2[i][1]);
    }

    fclose(fp);

    fp=fopen("../lineq3.dat", "r");
    fscanf(fp, "%d %d", &n3, &trash);
    a3 = (float**)malloc(sizeof(float*)*(n3+1));
    c3 = (float**)malloc(sizeof(float*)*(n3+1));

    for(int i=1; i<=n3; i++)
    {
        a3[i] = (float*)malloc(sizeof(float)*(n3+1));
        for (int j = 1; j <= n3; j++)
            fscanf(fp, "%f", &a3[i][j]);
    }

    for(int i=1; i<=n3; i++) {
        c3[i] = (float*)malloc(sizeof(float)*(2));
        fscanf(fp, "%f", &c3[i][1]);
    }

    fclose(fp);
}

void Print1(int n, float **a, float **c) { // 입력받은 값 A 및 c 그대로 출력
    printf("%d %d\n", n, n);
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=n; j++)
            printf("%f ", a[i][j]);
        printf("\n");
    }
    printf("\n");
    printf("RHV\n");
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=1; j++)
            printf("%f ", c[i][j]);
        printf("\n");
    }
    printf("\n");
}

float* Mprove(int n, float **aa, float **aalud, int *index, float **bb, float *xx) { // iterative improvement method

    int *indx = (int*)malloc(sizeof(int)*(n+1));
    float **mat, **matlud, *b, *x;

    mat = (float**)malloc(sizeof(float*)*(n+1));
    matlud = (float**)malloc(sizeof(float*)*(n+1));
    b = (float*)malloc(sizeof(float)*(n+1));
    x = (float*)malloc(sizeof(float)*(n+1));
    for(int i=1; i<=n; i++) {
        indx[i] = index[i];
        b[i]=bb[i][1];
        x[i]=xx[i];
        mat[i] = (float*)malloc(sizeof(float)*(n+1));
        matlud[i] = (float*)malloc(sizeof(float)*(n+1));
        for(int j=1; j<=n; j++) {
            mat[i][j] = aa[i][j];
            matlud[i][j] = aalud[i][j];
        }
    }

    mprove(mat,matlud,n,indx,b,x);

    printf("해\n");
    for(int i=1; i<=n; i++)
        printf("%f\n", x[i]);
    printf("\n");

    return x;

}

void Gaussj(int n, float **a, float **c) { // Gauss-Jordan Elimination
    printf("************GAUSSJ*************\n");
    Print1(n,a,c);

    float **mat, **b;

    mat = (float**)malloc(sizeof(float*)*(n+1));
    b = (float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++) {

        b[i] = (float*)malloc(sizeof(float)*(2));
        mat[i] = (float*)malloc(sizeof(float)*(n+1));
        for(int j=1; j<=n; j++)
            mat[i][j] = a[i][j];
        b[i][1] = c[i][1];
    }

    gaussj(mat,n,b,1);

    printf("해\n");
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=1; j++)
            printf("%f ", b[i][j]);
        printf("\n");
    }
    printf("\n");

    printf("검산결과\n");
    for(int i=1; i<=n; i++) {
        float sum=0.0;
        for(int j=1; j<=n; j++)
            sum += b[j][1]*a[i][j];
        printf("%f\n", sum);
    }

    printf("\n");
    printf("A의 역행렬\n");
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=n; j++)
            printf("%f ", mat[i][j]);
        printf("\n");
    }
    printf("\n");

}

void Ludcmp(int n, float **a, float **c) { // LU decomposition
    printf("************LUDCMP*************\n");
    Print1(n,a,c);

    int *indx = (int*)malloc(sizeof(int)*(n+1));
    float **mat, *b, d, **inverse;
    float determinant;
    float *xx;
    char ck;

    mat = (float**)malloc(sizeof(float*)*(n+1));
    b = (float*)malloc(sizeof(float)*(n+1));
    for(int i=1; i<=n; i++) {
        indx[i] = i;
        b[i]=c[i][1];
        mat[i] = (float*)malloc(sizeof(float)*(n+1));
        for(int j=1; j<=n; j++)
            mat[i][j] = a[i][j];
    }

    ludcmp(mat,n,indx,&d);

    determinant = d;
    for(int i=1; i<=n; i++) determinant*=mat[i][i];

    lubksb(mat,n,indx,b);

    printf("해\n");
    for(int i=1; i<=n; i++)
        printf("%f\n", b[i]);
    printf("\n");

    printf("검산결과\n");
    for(int i=1; i<=n; i++) {
        float sum=0.0;
        for(int j=1; j<=n; j++)
            sum += b[j]*a[i][j];
        printf("%f\n", sum);
    }
    printf("\n");

    printf("mprove method 사용을 통해 해 검산\n");
    printf("10^-4 이하의 오차는 무시한다\n");
    xx=Mprove(n,a,mat,indx,c,b);
    ck=0;
    for(int i=1; i<=n; i++) {
        if(fabs(b[i]-xx[i]) > 0.0001)
            ck=1;
    }
    if(ck==0) {
        printf("mprove 사용결과 동일!\n\n");
    }

    printf("A의 determinant\n%f\n\n",determinant);
    if(fabs(determinant -0.0) < 0.0001) {
        printf("A의 역행렬은 존재하지 않는다\n");
        return;
    }


    inverse = (float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++) {
        inverse[i]=(float*)malloc(sizeof(float)*(n+1));
        indx[i] = i;
        for(int j=1; j<=n; j++)
            mat[i][j] = a[i][j];
    }
    ludcmp(mat,n,indx,&d);
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=n; j++) {
            b[j]=0.0;
        }
        b[i]=1.0;
        lubksb(mat,n,indx,b);
        for(int j=1; j<=n; j++) inverse[j][i]=b[j];
    }
    printf("A의 역행렬\n");
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=n; j++)
            printf("%f ", inverse[i][j]);
        printf("\n");
    }
    printf("\n");
}

void Svdcmp(int n, float **a, float **c) { // Singular Value Decomposition

    printf("************SVDCMP*************\n");
    Print1(n,a,c);

    float *b = (float*)malloc(sizeof(float)*(n+1));
    float *x = (float*)malloc(sizeof(float)*(n+1));
    float *w = (float*)malloc(sizeof(float)*(n+1));
    float **v = (float**)malloc(sizeof(float*)*(n+1));
    float **mat = (float**)malloc(sizeof(float*)*(n+1));

    for(int i=1; i<=n; i++) {
        w[i]=b[i]=c[i][1];
        v[i] = (float*)malloc(sizeof(float)*(n+1));
        mat[i] = (float*)malloc(sizeof(float)*(n+1));
        for (int j=1;j<=n;j++) mat[i][j]=a[i][j];
    }

    svdcmp(mat,n,n,w,v);

    float wmax=0.0;
    for(int j=1;j<=n;j++) if (w[j] > wmax) wmax=w[j];
    float wmin=wmax*1.0e-6;;
    for(int j=1;j<=n;j++) if (w[j] < wmin) w[j]=0.0;
    svbksb(mat,w,v,n,n,b,x);

    printf("해\n");
    for(int i=1; i<=n; i++)
        printf("%f\n", x[i]);

    printf("\n");
    printf("검산결과\n");
    for(int i=1; i<=n; i++) {
        float sum=0.0;
        for(int j=1; j<=n; j++)
            sum += x[j]*a[i][j];
        printf("%f\n", sum);
    }
    printf("\n");

    for(int i=1; i<=n; i++) {
        if (fabs(w[i] - 0.0) < 0.0001) {
            printf("A의 역행렬은 존재하지 않는다\n");
            return;
        }
    }

    printf("A의 역행렬\n");
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=n; j++)
            v[i][j] /= w[j];
    }
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=n; j++) {
            float sum=0.0;
            for(int k=1; k<=n; k++)
                sum += v[i][k] * mat[j][k];
            printf("%f ", sum);
        }
        printf("\n");
    }


}

int main()
{
    Input();
    printf("*************************Problem 1\n");
//    Gaussj(n1,a1,c1); 1번 문제는 역행렬이 존재하지 않는 Singular Matrix이므로 작동X
    Ludcmp(n1,a1,c1);
    Svdcmp(n1,a1,c1);

    printf("*************************Problem 2\n");
    Gaussj(n2,a2,c2);
    Ludcmp(n2,a2,c2);
    Svdcmp(n2,a2,c2);

    printf("*************************Problem 3\n");
    Gaussj(n3,a3,c3);
    Ludcmp(n3,a3,c3);
    Svdcmp(n3,a3,c3);
    return 0;
}

