//
// Created by patrick on 05/06/23.
//

#ifndef ES_APPCTX_H
#define ES_APPCTX_H


struct appContext {

    int mode;
    int IND_SIZE;        // Tamanho do indivíduo (quantidade de coeficientes)
    double MIN_K;        // Menor valor que K pode assumir
    double MAX_K;        // Maior valor que K pode assumir
    double MIN_N;        // Menor valor que N pode assumir
    double MAX_N;        // Maior valor que N pode assumir
    double MIN_TAU;      // Menor valor que TAU pode assumir
    double MAX_TAU;      // Maior valor que TAU pode assumir
    double MIN_STRATEGY; // Menor valor que a estratégia pode assumir
    double MAX_STRATEGY; // Maior valor que a estratégia pode assumir
    int TAU_SIZE;
    int N_SIZE;
    int K_SIZE;
    double *maxValues;
    int granularity;
    int nVariables;
    int nSteps;
    int dataSetSize;
    int fullSetStart;
    int fullSetEnd;
    int totalSteps;
    int trainingSteps;
    int testSteps;
    int validationSteps;
    int trainingSetStart;
    int trainingSetEnd;
    int testSetStart;
    int testSetEnd;
    int validationSetStart;
    int validationSetEnd;
    int setStart;
    int setEnd;
    double tspan[2];
    double trainingTSpan[2];
    double *yout;
    double *y_0;
    double **vectors;
    double **expectedResult;
    double* individual;
    void *series;

} ;
#endif //ES_APPCTX_H
