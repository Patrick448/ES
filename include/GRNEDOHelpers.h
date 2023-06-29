//
// Created by patrick on 29/06/23.
//

#ifndef ES_GRNEDOHELPERS_H
#define ES_GRNEDOHELPERS_H
#include "appCtx.h"
#include "dependencies.h"

namespace GRNEDOHelpers{
    int twoBody5VarLSODA(double t, double *y, double *ydot, void *_data);
    int twoBody10VarLSODA(double t, double *y, double *ydot, void *_data);
    double difference(double *actual, double **expected, int numVariables, int numElements);
    double difference(double *actual, double **expected, int numVariables, int start, int end, int numSteps);
    void initializeGRN5Context(appContext* ctx, int mode, int granularity);
    void initializeGRN10Context(appContext* ctx, int mode, int granularity);
    void clearContext(appContext* ctx);
    double lsodaWrapper(int dydt(double t, double *y, double *ydot, void *data), appContext *appCtx, double *_yout);
    double getMaxValue(double *values, int numElements);
    void getMaxValues(double **data, double *outMaxValues, int numVariables, int numElements);
    double getMaxValue(double *values, int start, int end);
    void getMaxValues(double **data, double *outMaxValues, int numVariables, int start, int end);
    void printContext(appContext* ctx);
    void twoBody(double t, double y[], double max[], double tau[], double n[], double k[], double yp[]);
    string vectorToString(double *vec, int start, int end);
    void printGRNVector(double **vec, int rows, int cols);
    int twoBodyFixedLSODA(double t, double *y, double *ydot, void *data);
    void twoBodyFixed(double t, double y[], double *dim, double yp[]);
    double grn5EvaluationLSODA(void *ind, void* data);
    double grn5EvaluationRK4(void *ind, void* data);
    double grn10EvaluationLSODA(void *ind, void *data);
    double grn10EvaluationRK4(void *ind, void* data);




}

#endif //ES_GRNEDOHELPERS_H
