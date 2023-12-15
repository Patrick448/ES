//
// Created by patrick on 29/06/23.
//

#ifndef ES_GRNEDOHELPERS_H
#define ES_GRNEDOHELPERS_H
#include "appCtx.h"
#include "dependencies.h"


/// This namespace contains all the functions that are used to create the context of the problem
/// and to evaluate the individuals.
namespace GRNEDOHelpers{
    /// 5 variable GRN model.
    /// todo: rename to grn5Model
    int twoBody5VarLSODA(double t, double *y, double *ydot, void *_data);
    /// 10 variable GRN model.
    /// todo: rename to grn10Model
    int twoBody10VarLSODA(double t, double *y, double *ydot, void *_data);
    double difference(double *actual, double **expected, int numVariables, int numElements);
    double difference(double *actual, double **expected, int numVariables, int start, int end, int numSteps);

    /// initializes the context of the problem. Expects the GRN5.txt file to be in the same directory as the executable.
    /// @param ctx the context to be initialized.
    /// @param mode the mode of the problem, which can be training, validation, test
    /// or 'single set mode'. The latter means there are no separate sets.
    /// @param granularity how many steps to take between each time step:
    /// the data is divided into 49 time points, so if granularity is 1, there will be no steps between each time point.
    /// If granularity is 2, there will be 1 step between each time point, and so on.
    void initializeGRN5Context(appContext* ctx, int mode, int granularity);

    /// same as initializeGRN5Context, but for the 10 variable model.
    /// @see initializeGRN5Context
    void initializeGRN10Context(appContext* ctx, int mode, int granularity);
    void initializeGRNContext(appContext* ctx, int granularity, int numVariables, int numTau, int numN, int numK, int setStart, int setEnd, double** vectors, double * maxValues);

    /// clears the context of the problem, deallocating memory.
    /// should always be called at some point after initializing the context.
    void clearContext(appContext* ctx);
    void clearContext2Test(appContext* ctx);

    /// numerical integration of the given model using the LSODA algorithm.
    /// todo: rename to something more descriptive, like LSODAIntegration
    /// @param dydt the model to be integrated.
    /// @param appCtx the context of the problem.
    /// @param _yout the output of the integration (meaning the function values at each time step).
    double lsodaWrapper(int dydt(double t, double *y, double *ydot, void *data), appContext *appCtx, double *_yout);

    /// evaluates the given individual (5 variables) using the LSODA algorithm.
    double grn5EvaluationLSODA(void *ind, void* data);
    /// evaluates the given individual (5 variables) using the RK4 algorithm.
    double grn5EvaluationRK4(void *ind, void* data);

    /// evaluates the given individual (10 variables) using the LSODA algorithm.
    double grn10EvaluationLSODA(void *ind, void *data);

    /// evaluates the given individual (10 variables) using the RK4 algorithm.
    double grn10EvaluationRK4(void *ind, void* data);

    /// changes the mode of the problem.
    /// @see initializeGRN5Context for more details.
    void setMode(appContext* ctx, int mode);

    // other helpers
    string vectorToString(double *vec, int start, int end);
    void printGRNVector(double **vec, int rows, int cols);
    double getMaxValue(double *values, int numElements);
    void getMaxValues(double **data, double *outMaxValues, int numVariables, int numElements);
    double getMaxValue(double *values, int start, int end);
    void getMaxValues(double **data, double *outMaxValues, int numVariables, int start, int end);
    void printContext(appContext* ctx);

    // legacy functions
    void twoBody(double t, double y[], double max[], double tau[], double n[], double k[], double yp[]);
    int twoBodyFixedLSODA(double t, double *y, double *ydot, void *data);
    void twoBodyFixed(double t, double y[], double *dim, double yp[]);






}

#endif //ES_GRNEDOHELPERS_H
