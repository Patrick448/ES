//
// Created by patrick on 04/01/24.
//

#ifndef ES_GRN5NEWMODEL_H
#define ES_GRN5NEWMODEL_H

#include "dependencies.h"
#include "GRNModel.h"

class GRN5NewModel:GRNModel {

public:
    static const int TAU_SIZE = 5;
    static const int N_SIZE = 8;
    static const int K_SIZE = 8;
    static const int MODEL_DIMENSION = TAU_SIZE + N_SIZE + K_SIZE;
    static constexpr double MIN_K = 0.1;
    static constexpr double MAX_K = 1;
    static constexpr double MIN_N = 1;
    static constexpr double MAX_N = 25;
    static constexpr double MIN_TAU = 0.1;
    static constexpr double MAX_TAU = 5;

    int modelFunction(double t, double *y, double *ydot, void *context) {
        appContext *ctx = (appContext *) context;
        double *individual = ctx->individual;
        double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();//maxA, maxB, maxC, maxD, maxE
        double *tau = &individual[0]; //tau[0], tau[1], tau[2], tau[3], tau[4]
        double *k = &individual[TAU_SIZE]; //k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]
        double *n = &individual[TAU_SIZE +K_SIZE]; //(int)n[0], (int)n[1], (int)n[2], (int)n[3], (int)n[4], (int)n[5], (int)n[6], (int)n[7]

        ydot[0] = ((1 - ((pow(y[4] / maxValues[4], (int) n[0])) /
                         (pow(y[4] / maxValues[4], (int) n[0]) + pow(k[0], (int) n[0])))) - (y[0] / maxValues[0])) /
                  tau[0];

        ydot[1] = ((1 - ((pow(y[4] / maxValues[4], (int) n[1])) /
                         (pow(y[4] / maxValues[4], (int) n[1]) + pow(k[1], (int) n[1])))) - (y[1] / maxValues[1])) /
                  tau[1];

        ydot[2] = ((((pow(y[1] / maxValues[1], (int) n[2])) /
                     (pow(y[1] / maxValues[1], (int) n[2]) + pow(k[2], (int) n[2]))) +
                    ((pow(y[2] / maxValues[2], (int) n[3])) /
                     (pow(y[2] / maxValues[2], (int) n[3]) + pow(k[3], (int) n[3])))) - (y[2] / maxValues[2])) / tau[2];

        ydot[3] = ((((pow(y[1] / maxValues[1], (int) n[4])) /
                     (pow(y[1] / maxValues[1], (int) n[4]) + pow(k[4], (int) n[4]))) +
                    ((pow(y[2] / maxValues[2], (int) n[5])) /
                     (pow(y[2] / maxValues[2], (int) n[5]) + pow(k[5], (int) n[5])))) - (y[3] / maxValues[3])) / tau[3];

        ydot[4] = (((pow(y[2] / maxValues[2], (int) n[6])) /
                    (pow(y[2] / maxValues[2], (int) n[6]) + pow(k[6], (int) n[6]))) *
                   (((pow(y[3] / maxValues[3], (int) n[7])) /
                     (pow(y[3] / maxValues[3], (int) n[7]) + pow(k[7], (int) n[7])))) - (y[4] / maxValues[4])) / tau[4];


        return 0;
    }

};

#endif //ES_GRN5NEWMODEL_H
