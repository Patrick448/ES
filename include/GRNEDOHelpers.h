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
    int grn5Model(double t, double *y, double *ydot, void *data);
    int grn5NewModel(double t, double *y, double *ydot, void *data);
    int grn10Model(double t, double *y, double *ydot, void *data);
    int grn10NewModel(double t, double *y, double *ydot, void *data);
    int grn5NCYCModel(double t, double *y, double *ydot, void *context);
    int grn4NCYCModel(double t, double *y, double *ydot, void *context);
    double difference(double *actual, double **expected, int numElements, int numVariables, int granularity);

    /// todo: rename to something more descriptive, like LSODAIntegration
    /// @param dydt the model to be integrated.
    /// @param appCtx the context of the problem.
    /// @param _yout the output of the integration (meaning the function values at each time step).
    double lsodaWrapper(int dydt(double t, double *y, double *ydot, void *data), double* tspan, double* y_0, int totalSteps, int nVariables, double* times, double *_yout, void* context);
    double lsodaWrapperTest(int dydt(double t, double *y, double *ydot, void *data), double* tspan, double* y_0, int totalSteps, int nVariables, double* times, double *_yout, void* context);

    double grnEvaluationLSODA(void* individual, void* context);
    double grnEvaluationRK4(void* individual, void* context);
    void printODEIntSeries(void* individual, void* context, const string& outputFile, int solver);

    // other helpers
    string vectorToString(double *vec, int start, int end);
    void printGRNVector(double **vec, int rows, int cols);
    void printContext(appContext* ctx);







}

#endif //ES_GRNEDOHELPERS_H
