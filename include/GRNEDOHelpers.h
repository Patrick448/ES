//
// Created by patrick on 29/06/23.
//

#ifndef ES_GRNEDOHELPERS_H
#define ES_GRNEDOHELPERS_H
#include "appCtx.h"
#include "dependencies.h"
#include "ProblemDescription.h"


/// This namespace contains all the functions that are used to create the context of the problem
/// and to evaluate the individuals.
namespace GRNEDOHelpers{

    static ProblemDescription grn5ProblemDescription = {
            .IND_SIZE = 19,        // Tamanho do indivíduo (quantidade de coeficientes)
            .MIN_K = 0.01, //0.1,       // Menor valor que K pode assumir
            .MAX_K = 1,       // Maior valor que K pode assumir
            .MIN_N = 1,        // Menor valor que N pode assumir
            .MAX_N = 30, //25       // Maior valor que N pode assumir
            .MIN_TAU = 0.1,      // Menor valor que TAU pode assumir
            .MAX_TAU = 6,  //5     // Maior valor que TAU pode assumir
            .TAU_SIZE = 5,
            .N_SIZE = 7,
            .K_SIZE = 7
    } ;
    /// 5 variable GRN model.
    /// todo: rename to grn5Model
    int twoBody5VarLSODA(double t, double *y, double *ydot, void *_data);
    int grn5Model(double t, double *y, double *ydot, void *data);
    int grn10Model(double t, double *y, double *ydot, void *data);

    /// 10 variable GRN model.
    /// todo: rename to grn10Model
    double difference(double *actual, double **expected, int numVariables, int numElements);
    double difference(double *actual, double **expected, int numVariables, int start, int end, int numSteps);
    double differenceTest(double *actual, double **expected, int numElements, int numVariables, int granularity);

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
    double lsodaWrapperTest(int dydt(double t, double *y, double *ydot, void *data), double* tspan, double* y_0, int totalSteps, int nVariables, double* times, double *_yout, void* context);


    double grnEvaluationLSODA(void* individual, void* context);
    double grnEvaluationRK4(void* individual, void* context);
    double grn10EvaluationLSODA(void* individual, void* context);
    double grn10EvaluationRK4(void* individual, void* context);
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
