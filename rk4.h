//
// Created by patri on 07/02/2023.
//

#ifndef ES_RK4_H
#define ES_RK4_H

#ifdef __cplusplus
extern "C" {
#endif

extern void rk4 ( void dydt ( double t, double u[], double f[] ), double tspan[2], double y0[], int n, int m, double t[], double y[] );
//extern void timestamp();

#ifdef __cplusplus
}
#endif

#endif
