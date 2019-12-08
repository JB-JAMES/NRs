
#include <math.h>

#define MAXIT 100

float muller(func,x1,x2,xacc)
float (*func)(),x1,x2,xacc;
{
    void nrerror();
    int j;
    float fl,f,fa,dx,swap,xl,rts,xa;

    fl=(*func)(x1);
    f=(*func)(x2);

    if (fabs(fl) < fabs(f)) {
        rts=x1;
        xl=x2;
        swap=fl;
        fl=f;
        f=swap;
    } else {
        xl=x1;
        rts=x2;
    }
    if (x1 == xa || xa == x2)
        nrerror("Interval ERROR");

    xa = (xl+rts) * 0.5;
    fa=(*func)(xa);
    for (j=1;j<=MAXIT;j++) {
        float a, b, c; // 계산식
        c = f;
        b = ((xl-rts)*(xl-rts)*(fa-f) - (xa-rts)*(xa-rts)*(fl-f)) / ((xl-rts)*(xa-rts)*(xl-xa));
        a = ((xa-rts)*(fl-f) - (xl-rts)*(fa-f)) / ((xl-rts)*(xa-rts)*(xl-xa));

        if (b*b < 4*a*c)
            nrerror("ERROR");

        if(2*a + b > 0 && x1 < rts) // 위로 볼록하면 순서를 역순으로 한다.
        {
            swap=xl;
            xl=rts;
            rts=swap;
            swap=fl;
            fl=f;
            f=swap;
        }

        dx = -2.0*c/(b+((b>0) ? 1.0 : -1.0)*sqrt(b*b-4*a*c));
        xl=xa;
        fl=fa;
        xa=rts;
        fa=f;
        rts += dx;
        f=(*func)(rts);

        if (fabs(dx) < xacc || f == 0.0) {
            printf("\n돈 횟수: %d\n", j);
            return rts;
        }

    }
    nrerror("Maximum number of iterations exceeded in muller");
    return 0.0;
}
#undef MAXIT
