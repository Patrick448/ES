//
// Created by patri on 07/02/2023.
//

#ifndef ES_RK4_H
#define ES_RK4_H

#ifdef __cplusplus
extern "C" {
#endif

extern void rk4 ( int dydt (double t, double *y, double *ydot, void *_data), double tspan[2],
                  double y0[], int n, int m, double t[], double* yout, void *ctx );
//extern void timestamp();

#ifdef __cplusplus
}
#endif

#endif
