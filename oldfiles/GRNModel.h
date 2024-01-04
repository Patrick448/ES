//
// Created by patrick on 04/01/24.
//

#ifndef ES_GRNMODEL_H
#define ES_GRNMODEL_H
#include "dependencies.h"
#include "GRN5Model.h"

class GRNModel {

public:

    static const int TAU_SIZE;
    static const int N_SIZE;
    static const int K_SIZE;
    static const int MODEL_DIMENSION;
    static constexpr double MIN_K = 0.1;
    static constexpr double MAX_K = 1;
    static constexpr double MIN_N = 1;
    static constexpr double MAX_N = 25;
    static constexpr double MIN_TAU = 0.1;
    static constexpr double MAX_TAU = 5;
    virtual int modelFunction(double t, double *y, double *ydot, void *context);

};

#endif //ES_GRNMODEL_H
