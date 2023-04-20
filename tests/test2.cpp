#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif 

#include "lsoda.h"
#include "common.h"

#ifdef __cplusplus
}
#endif

int fex(double t, double *y, double *ydot, void *data)
{
    // todo: observar que os n devem ser avaliados como inteiros
    double max[] = {2.96, 1.8768, 1.0653, 1.0101, 1.4608};
    double tau[] = {1.25, 4, 1.02, 1.57, 3.43};
    double n[] = {13, 4, 3, 4, 16};
    double k[] = {0.72, 0.50, 0.45, 0.51, 0.52};

    ydot[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[4]), n[0]) + pow(k[0], n[0]))) - (y[0] / max[0])) / tau[0];

    ydot[1] = (((pow((y[0] / max[0]), n[1])) / (pow((y[0] / max[0]), n[1]) + pow(k[1], n[1]))) - (y[1] / max[1])) /
              tau[1];

    ydot[2] = (((pow((y[1] / max[1]), n[2])) / (pow((y[1] / max[1]), n[2]) + pow(k[2], n[2]))) - (y[2] / max[2])) /
              tau[2];

    ydot[3] = (((pow((y[2] / max[2]), n[3])) / (pow((y[2] / max[2]), n[3]) + pow(k[3], n[3]))) - (y[3] / max[3])) /
              tau[3];

    ydot[4] = (((pow((y[3] / max[3]), n[4])) / (pow((y[3] / max[3]), n[4]) + pow(k[4], n[4]))) - (y[4] / max[4])) /
              tau[4];

    return 0;
}

int test(void)
{
    double tspan[] = {0.0, 72.0};
    int n = 49;
	double          atol[5], rtol[5], t, tout, y[5];
	int             neq = 5;
	int             iout;
    double dt = (tspan[1] - tspan[0]) / (double)(n);
	y[0] = 0.7095;
	y[1] = 0.1767;
	y[2] = 0.1522;
	y[3] = 0.0806;
	y[4] = 0.7095;
	t = 0.0E0;
	tout = dt;
	struct lsoda_opt_t opt = {0};
	opt.ixpr = 0;
	opt.rtol = rtol;
	opt.atol = atol;
	opt.itask = 1;

    rtol[0] = atol[0] = 1.49012e-8;
    rtol[1] = atol[1] = 1.49012e-8;
    rtol[2] = atol[2] = 1.49012e-8;
    rtol[3] = atol[3] = 1.49012e-8;
    rtol[4] = atol[4] = 1.49012e-8;


	struct lsoda_context_t ctx = {
		.function = fex,
		.data = NULL,
		.neq = neq,
		.state = 1,
	};

	lsoda_prepare(&ctx, &opt);

	for (iout = 1; iout <= n; iout++) {
		lsoda(&ctx, y, &t, tout);
		printf(" at t= %12.4e y= %14.6e %14.6e %14.6e %14.6e %14.6e\n", t, y[0], y[1], y[2], y[3], y[4]);
		if (ctx.state <= 0) {
			printf("error istate = %d\n", ctx.state);
			exit(0);
		}
		tout = tout + dt;
	}
	lsoda_free(&ctx);
	return(0);
}

int main(void) {
	int i;
	for(i = 0; i < 1; i++) {
		test();
	}
	return(0);
}

/*
 The correct answer (up to certain precision):

 at t=   4.0000e-01 y=   9.851712e-01   3.386380e-05   1.479493e-02
 at t=   4.0000e+00 y=   9.055333e-01   2.240655e-05   9.444430e-02
 at t=   4.0000e+01 y=   7.158403e-01   9.186334e-06   2.841505e-01
 at t=   4.0000e+02 y=   4.505250e-01   3.222964e-06   5.494717e-01
 at t=   4.0000e+03 y=   1.831976e-01   8.941773e-07   8.168015e-01
 at t=   4.0000e+04 y=   3.898729e-02   1.621940e-07   9.610125e-01
 at t=   4.0000e+05 y=   4.936362e-03   1.984221e-08   9.950636e-01
 at t=   4.0000e+06 y=   5.161833e-04   2.065787e-09   9.994838e-01
 at t=   4.0000e+07 y=   5.179804e-05   2.072027e-10   9.999482e-01
 at t=   4.0000e+08 y=   5.283675e-06   2.113481e-11   9.999947e-01
 at t=   4.0000e+09 y=   4.658667e-07   1.863468e-12   9.999995e-01
 at t=   4.0000e+10 y=   1.431100e-08   5.724404e-14   1.000000e+00
 */
