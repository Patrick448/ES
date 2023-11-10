#include "GRNEDOHelpers.h"
#include "appCtx.h"
#include <cmath>
#include <cfloat>
#include "dependencies.h"
#ifdef __cplusplus
extern "C"
{
#endif

#include "lsoda.h"
#include "algModes.h"

#ifdef __cplusplus
}
#endif

using namespace std;
using namespace algModes;

int GRNEDOHelpers::twoBody5VarLSODA(double t, double *y, double *ydot, void *_data)
{

    appContext *ctx = (appContext *)_data;
    double* data = ctx->individual;
    double* maxValues = ctx->maxValues;
    double *tau = &data[0];
    double *k = &data[ctx->TAU_SIZE];
    double *n = &data[ctx->TAU_SIZE + ctx->N_SIZE];


    ydot[0] = ((1 - (pow((y[4] / maxValues[4]), (int)n[0])) /
                    (pow((y[4] / maxValues[4]), (int)n[0]) + pow(k[0], (int)n[0]))) -
               (y[0] / maxValues[0])) /
              tau[0];

    ydot[1] = (((pow((y[0] / maxValues[0]), (int)n[1])) /
                (pow((y[0] / maxValues[0]), (int)n[1]) + pow(k[1], (int)n[1]))) -
               (y[1] / maxValues[1])) /
              tau[1];

    ydot[2] = (((pow((y[1] / maxValues[1]), (int)n[2])) /
                (pow((y[1] / maxValues[1]), (int)n[2]) + pow(k[2], (int)n[2]))) -
               (y[2] / maxValues[2])) /
              tau[2];

    ydot[3] = (((pow((y[2] / maxValues[2]), (int)n[3])) /
                (pow((y[2] / maxValues[2]), (int)n[3]) + pow(k[3], (int)n[3]))) -
               (y[3] / maxValues[3])) /
              tau[3];

    ydot[4] = ((((pow(y[1] / maxValues[1], (int)n[4]) / (pow(y[1] / maxValues[1], (int)n[4]) + pow(k[4], (int)n[4]))) * (
            pow(y[3] / maxValues[3], (int)n[5]) / (pow(y[3] / maxValues[3], (int)n[5]) + pow(k[5], (int)n[5])))) + (
                        (pow(y[3] / maxValues[3], (int)n[5]) / (pow(y[3] / maxValues[3], (int)n[5]) + pow(k[5], (int)n[5]))) * (
                                pow(y[4] / maxValues[4], (int)n[6]) / (pow(y[4] / maxValues[4], (int)n[6]) + pow(k[6], (int)n[6]))))) - (
                       y[4] / maxValues[4])) / tau[4];

    return 0;
}


int GRNEDOHelpers::twoBody10VarLSODA(double t, double *y, double *ydot, void *_data)
{

    appContext *ctx = (appContext *)_data;
    double* data = ctx->individual;
    double* maxValues = ctx->maxValues;    double *tau = &data[0];

    double *k = &data[ctx->TAU_SIZE];
    double *n = &data[ctx->TAU_SIZE + ctx->N_SIZE];
    double maximo_A = maxValues[0];
    double maximo_B = maxValues[1];
    double maximo_C = maxValues[2];
    double maximo_D = maxValues[3];
    double maximo_E = maxValues[4];
    double maximo_F = maxValues[5];
    double maximo_G = maxValues[6];
    double maximo_H = maxValues[7];
    double maximo_I = maxValues[8];
    double maximo_J = maxValues[9];
    double tauA = data[0];
    double tauB = data[1];
    double tauC = data[2];
    double tauD = data[3];
    double tauE = data[4];
    double tauF = data[5];
    double tauG = data[6];
    double tauH = data[7];
    double tauI = data[8];
    double tauJ = data[9];
    double kAJ = data[10];
    double kBE = data[11];
    double kCB = data[12];
    double kCF = data[13];
    double kCA = data[14];
    double kDF = data[15];
    double kEJ = data[16];
    double kFA = data[17];
    double kGB = data[18];
    double kGF = data[19];
    double kGA = data[20];
    double kHF = data[21];
    double kIG = data[22];
    double kIH = data[23];
    double kJI = data[24];
    double nAJ = (int)data[25];
    double nBE = (int)data[26];
    double nCB = (int)data[27];
    double nCF = (int)data[28];
    double nCA = (int)data[29];
    double nDF = (int)data[30];
    double nEJ = (int)data[31];
    double nFA = (int)data[32];
    double nGB = (int)data[33];
    double nGF = (int)data[34];
    double nGA = (int)data[35];
    double nHF = (int)data[36];
    double nIG = (int)data[37];
    double nIH = (int)data[38];
    double nJI = (int)data[39];

    ydot[0] = ((1 - pow(y[9] / maximo_J, nAJ) / (pow(y[9] / maximo_J, nAJ) + pow(kAJ, nAJ))) - (y[0] / maximo_A)) / tauA;

    ydot[1] = (pow(y[4] / maximo_E, nBE) / (pow(y[4] / maximo_E, nBE) + pow(kBE, nBE)) - (y[1] / maximo_B)) / tauB;

    ydot[2] = (pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB)) *
               (1 - pow(y[5] / maximo_F, nCF) / (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF))) *
               (1 - pow(y[0] / maximo_A, nCA) / (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA))) +
               (1 - pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB))) * pow(y[5] / maximo_F, nCF) /
               (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF)) *
               (1 - pow(y[0] / maximo_A, nCA) / (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA))) +
               (1 - pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB))) *
               (1 - pow(y[5] / maximo_F, nCF) / (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF))) * pow(y[0] / maximo_A, nCA) /
               (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) +
               pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB)) *
               (1 - pow(y[5] / maximo_F, nCF) / (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF))) * pow(y[0] / maximo_A, nCA) /
               (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) +
               (1 - pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB))) * pow(y[5] / maximo_F, nCF) /
               (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF)) * pow(y[0] / maximo_A, nCA) /
               (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) +
               pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB)) * pow(y[5] / maximo_F, nCF) /
               (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF)) * pow(y[0] / maximo_A, nCA) /
               (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) -
               (y[2] / maximo_C)) /
              tauC;

    ydot[3] = (pow(y[5] / maximo_F, nDF) / (pow(y[4] / maximo_E, nDF) + pow(kDF, nDF)) - (y[3] / maximo_D)) / tauD;

    ydot[4] = (1 - pow(y[9] / maximo_J, nEJ) / (pow(y[9] / maximo_J, nEJ) + pow(kEJ, nEJ)) - (y[4] / maximo_E)) / tauE;

    ydot[5] = (pow(y[0] / maximo_A, nFA) / (pow(y[0] / maximo_A, nFA) + pow(kFA, nFA)) - (y[5] / maximo_F)) / tauF;

    ydot[6] = (pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB)) *
               (1 - pow(y[5] / maximo_F, nGF) / (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF))) *
               (1 - pow(y[0] / maximo_A, nGA) / (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA))) +
               (1 - pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB))) * pow(y[5] / maximo_F, nGF) /
               (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF)) *
               (1 - pow(y[0] / maximo_A, nGA) / (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA))) +
               (1 - pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB))) *
               (1 - pow(y[5] / maximo_F, nGF) / (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF))) * pow(y[0] / maximo_A, nGA) /
               (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) +
               pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB)) *
               (1 - pow(y[5] / maximo_F, nGF) / (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF))) * pow(y[0] / maximo_A, nGA) /
               (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) +
               (1 - pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB))) * pow(y[5] / maximo_F, nGF) /
               (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF)) * pow(y[0] / maximo_A, nGA) /
               (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) +
               pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB)) * pow(y[5] / maximo_F, nGF) /
               (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF)) * pow(y[0] / maximo_A, nGA) /
               (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) -
               (y[6] / maximo_G)) /
              tauG;

    ydot[7] = (pow(y[5] / maximo_F, nHF) / (pow(y[5] / maximo_F, nHF) + pow(kHF, nHF)) - (y[7] / maximo_H)) / tauH;

    ydot[8] = (pow(y[6] / maximo_G, nIG) / (pow(y[6] / maximo_G, nIG) + pow(kIG, nIG)) * pow(y[7] / maximo_H, nIH) /
               (pow(y[7] / maximo_H, nIH) + pow(kIH, nIH)) -
               (y[8] / maximo_I)) /
              tauI;

    ydot[9] = (pow(y[8] / maximo_I, nJI) / (pow(y[8] / maximo_I, nJI) + pow(kJI, nJI)) - (y[9] / maximo_J)) / tauJ;

    ydot[3] = (pow(y[5] / maximo_F, nDF) / (pow(y[4] / maximo_E, nDF) + pow(kDF, nDF)) - (y[3] / maximo_D)) / tauD;

    //    cout << "#########################" << endl;
    //    cout << pow(y[5] / maximo_F, nDF) << endl;
    //    cout << pow(y[4] / maximo_E, nDF) << endl;
    //    cout << pow(kDF, nDF) << endl;
    //    cout << (y[3] / maximo_D) << endl;
    //    cout << tauD << endl;
    //    cout << y[3] << endl;
    //    cout << maximo_D << endl;
    //    cout << "#########################" << endl;

    return 0;
}

double GRNEDOHelpers::difference(double *actual, double **expected, int numVariables, int numElements)
{
    // todo: testar função e remover essa linha
    // double *dif = new double [numVariables];
    double difTotal = 0.0;

    for (int i = 0; i < numVariables; i++)
    {
        // dif[i] = 0;
        for (int j = 0; j < numElements; j++)
        {
            //  dif[i] += fabs(actual[i][j] - expected[i][j]);
            // if(i==3)continue;
            difTotal += fabs(actual[j * numVariables + i] - expected[i][j]);
        }
    }

    if (isnan(difTotal))
    {
        return DBL_MAX;
    }

    return difTotal;
}
//todo: melhorar vetores, padronizar o formato e melhorar os acessos
double GRNEDOHelpers::difference(double *actual, double **expected, int numVariables, int start, int end, int numSteps)
{
    double difTotal = 0.0;
    int numElements = end - start + 1;
   // printGRNVector(expected, numVariables, 50);
  //  cout << "\n\n" << endl;
    int jump = numSteps / (numElements-1);

    for (int i = 0; i < numVariables; i++)
    {
        for (int j = start; j <= end; j++)
        {
            int index = jump*j* numVariables + i;
            difTotal += fabs(actual[index] - expected[i][j]);
            //cout << "(" << j << "," << i << ")"<< " - ";
            //cout <<  index << " - ";
           // cout << actual[index] << ",";

            //cout << actual[j*jump* numVariables + i]<< " ";
        }

      //  cout << endl;
    }
    //cout << endl;
    if (isnan(difTotal))
    {
        return DBL_MAX;
    }

    return difTotal;
}

void GRNEDOHelpers::setMode(appContext* ctx, int mode){
    if(mode == TRAINING_MODE){
        ctx->setStart = ctx->trainingSetStart;
        ctx->setEnd = ctx->trainingSetEnd;
        ctx->nSteps = ctx->trainingSteps;
    } else if(mode == VALIDATION_MODE){
        ctx->setStart = ctx->validationSetStart;
        ctx->setEnd = ctx->validationSetEnd;
        ctx->nSteps = ctx->validationSteps;
    }else if(mode == TEST_MODE){
        ctx->setStart = ctx->testSetStart;
        ctx->setEnd = ctx->testSetEnd;
        ctx->nSteps = ctx->testSteps;
    }else{
        ctx->setStart = 0;
        ctx->setEnd = ctx->dataSetSize - 1;
        ctx->nSteps = (ctx->dataSetSize - 1)*ctx->granularity;
    }
}

void GRNEDOHelpers::initializeGRN5Context(appContext* ctx, int mode, int granularity)
{
    //appContext *ctx = new appContext;
    /*ctx->TRAINING_MODE = 0;
    ctx->TEST_MODE = 2;
    ctx->VALIDATION_MODE = 1;
    ctx->SINGLE_SET_MODE = 3;*/
    ctx->IND_SIZE = 19;      // Tamanho do indivíduo (quantidade de coeficientes)
    ctx->MIN_K = 0.01; //0.1        // Menor valor que K pode assumir
    ctx->MAX_K = 1;          // Maior valor que K pode assumir
    ctx->MIN_N = 1;          // Menor valor que N pode assumir
    ctx->MAX_N = 30; //25        // Maior valor que N pode assumir
    ctx->MIN_TAU = 0.1;      // Menor valor que TAU pode assumir
    ctx->MAX_TAU = 6;//6        // Maior valor que TAU pode assumir
    ctx->MIN_STRATEGY = 0.1; // Menor valor que a estratégia pode assumir
    ctx->MAX_STRATEGY = 10;  // Maior valor que a estratégia pode assumir
    ctx->TAU_SIZE = 5;
    ctx->N_SIZE = 7;
    ctx->K_SIZE = 7;
    ctx->granularity = granularity;
    ctx->nVariables = 5;
    ctx->nSteps = 49;
    ctx->dataSetSize = 50;
    ctx->fullSetStart = 0;
    ctx->fullSetEnd = 49;
    ctx->totalSteps = granularity*49;
    ctx->trainingSetStart = 0;
    ctx->trainingSetEnd = 34;
    ctx->trainingSteps = granularity*35;
    ctx->validationSetStart = 30;
    ctx->validationSetEnd = 39;
    ctx->validationSteps = granularity*15;
    ctx->testSetStart = 35;
    ctx->testSetEnd = 49;
    ctx->testSteps = granularity*15;
    ctx->mode = mode;

    if(mode == TRAINING_MODE){
        ctx->setStart = ctx->trainingSetStart;
        ctx->setEnd = ctx->trainingSetEnd;
        ctx->nSteps = ctx->trainingSteps;
    } else if(mode == VALIDATION_MODE){
        ctx->setStart = ctx->validationSetStart;
        ctx->setEnd = ctx->validationSetEnd;
        ctx->nSteps = ctx->validationSteps;
    }else if(mode == TEST_MODE){
        ctx->setStart = ctx->testSetStart;
        ctx->setEnd = ctx->testSetEnd;
        ctx->nSteps = ctx->testSteps;
    }else{
        ctx->setStart = 0;
        ctx->setEnd = ctx->dataSetSize - 1;
        ctx->nSteps = (ctx->dataSetSize - 1)*ctx->granularity;
    }

    ctx->maxValues = new double[ctx->nVariables];
    ctx->yout = new double[(ctx->totalSteps + 1) * ctx->nVariables];
    ctx->vectors = new double *[ctx->nVariables + 1];
    for (int i = 0; i < ctx->nVariables + 1; i++)
    {
        ctx->vectors[i] = new double[50];
    }

    readGRNFileToVectors("GRN5.txt", ctx->nVariables + 1, ctx->vectors);
    getMaxValues(ctx->vectors, ctx->maxValues, ctx->nVariables, ctx->trainingSetStart, ctx->trainingSetEnd);



    ctx->y_0 = new double[ctx->nVariables];
    ctx->expectedResult = &ctx->vectors[1];
    //printGRNVector(ctx->expectedResult, ctx->nVariables, ctx->dataSetSize);
    //cout << "\n\n";
    //todo: verificar se o y_0 está correto para
    // o caso de inicio fora do zero
    for (int i = 0; i < ctx->nVariables; i++)
    {
        ctx->y_0[i] = ctx->vectors[i + 1][ctx->fullSetStart];
    }

    ctx->tspan[0] = ctx->vectors[0][ctx->fullSetStart];
    ctx->tspan[1] = ctx->vectors[0][ctx->fullSetEnd];


    /////////////////////////////

    /*ctx->nSteps = 49;
    ctx->tspan[0] = 0.0;
    ctx->tspan[1] = 72.0;*/

    //ctx->trainingTSpan[0] = 0.0;
    //ctx->trainingTSpan[1] = (ctx->tspan[1]/ctx->nSteps)*ctx->trainingSteps;
}

void GRNEDOHelpers::initializeGRN10Context(appContext* ctx, int mode, int granularity)
{

    /*ctx->TRAINING_MODE = 0;
    ctx->TEST_MODE = 2;
    ctx->VALIDATION_MODE = 1;
    ctx->SINGLE_SET_MODE = 3;*/
    ctx->IND_SIZE = 40;      // Tamanho do indivíduo (quantidade de coeficientes)
    ctx->MIN_K = 0.1;        // Menor valor que K pode assumir
    ctx->MAX_K = 1;          // Maior valor que K pode assumir
    ctx->MIN_N = 1;          // Menor valor que N pode assumir
    ctx->MAX_N = 25;         // Maior valor que N pode assumir
    ctx->MIN_TAU = 0.1;      // Menor valor que TAU pode assumir
    ctx->MAX_TAU = 5;        // Maior valor que TAU pode assumir
    ctx->MIN_STRATEGY = 0.1; // Menor valor que a estratégia pode assumir
    ctx->MAX_STRATEGY = 10;  // Maior valor que a estratégia pode assumir
    ctx->TAU_SIZE = 10;
    ctx->N_SIZE = 15;
    ctx->K_SIZE = 15;
    ctx->granularity = granularity;
    ctx->nVariables = 10;
    ctx->nSteps = 49;
    ctx->dataSetSize = 50;
    ctx->fullSetStart = 0;
    ctx->fullSetEnd = 49;
    ctx->totalSteps = 49*granularity;
    ctx->trainingSetStart = 0;
    ctx->trainingSetEnd = 49;
    ctx->trainingSteps = granularity*49;
    ctx->validationSetStart = 20;
    ctx->validationSetEnd = 34;
    ctx->validationSteps = granularity*15;
    ctx->testSetStart = 35;
    ctx->testSetEnd = 49;
    ctx->testSteps = granularity*15;
    ctx->mode = mode;

    if(mode == TRAINING_MODE){
        ctx->setStart = ctx->trainingSetStart;
        ctx->setEnd = ctx->trainingSetEnd;
        ctx->nSteps = ctx->trainingSteps;
    } else if(mode == VALIDATION_MODE){
        ctx->setStart = ctx->validationSetStart;
        ctx->setEnd = ctx->validationSetEnd;
        ctx->nSteps = ctx->validationSteps;
    }else if(mode == TEST_MODE){
        ctx->setStart = ctx->testSetStart;
        ctx->setEnd = ctx->testSetEnd;
        ctx->nSteps = ctx->testSteps;
    }else{
        ctx->setStart = 0;
        ctx->setEnd = ctx->dataSetSize - 1;
        ctx->nSteps = (ctx->dataSetSize - 1)*ctx->granularity;
    }

    ctx->maxValues = new double[ctx->nVariables];
    ctx->yout = new double[(ctx->totalSteps + 1) * ctx->nVariables];
    ctx->vectors = new double *[ctx->nVariables + 1];
    for (int i = 0; i < ctx->nVariables + 1; i++)
    {
        ctx->vectors[i] = new double[50];
    }

    readGRNFileToVectors("GRN10.txt", ctx->nVariables + 1, ctx->vectors);
    getMaxValues(ctx->vectors, ctx->maxValues, ctx->nVariables, ctx->trainingSetStart, ctx->trainingSetEnd);


    ctx->y_0 = new double[ctx->nVariables];
    ctx->expectedResult = &ctx->vectors[1];
    //todo: verificar se o y_0 está correto para
    // o caso de inicio fora do zero
    for (int i = 0; i < ctx->nVariables; i++)
    {
        ctx->y_0[i] = ctx->vectors[i + 1][ctx->fullSetStart];
    }

    ctx->tspan[0] = ctx->vectors[0][ctx->fullSetStart];
    ctx->tspan[1] = ctx->vectors[0][ctx->fullSetEnd];

}

void GRNEDOHelpers::clearContext(appContext* ctx)
{
    // free ( t );
    delete[] ctx->yout;
    delete[] ctx->y_0;
    // todo: algum problema na desalocação, investigar
    delete[] ctx->maxValues;
    for (int i = 0; i < ctx->nVariables+1; i++)
    {
        delete[] ctx->vectors[i];
    }

    delete [] ctx->vectors;
}

double GRNEDOHelpers::lsodaWrapper(int dydt(double t, double *y, double *ydot, void *data), appContext *appCtx, double *_yout)
{

    // todo: tentar colocar essa alocação fora da função
    //  esses vetores serão alocados toda vez, e essa função será chamada
    //  a cada avaliação de indivíduo

    double *atol = new double[appCtx->nVariables];
    double *rtol = new double[appCtx->nVariables];
    double t, tout, dt;

    double *y = new double[appCtx->nVariables];
    int iout;

    for (int i = 0; i < appCtx->nVariables; i++)
    {
        y[i] = appCtx->y_0[i];
        rtol[i] = atol[i] = 1.49012e-4;
        _yout[i] = appCtx->y_0[i];
        //cout << appCtx->y_0[i] <<endl;
    }

    t = 0.0E0;
    dt = (appCtx->tspan[1] - appCtx->tspan[0]) / (double)(appCtx->totalSteps);
    tout = dt;

    struct lsoda_opt_t opt = {0};
    opt.ixpr = 0;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.itask = 1;
    opt.mxstep = 500;


    struct lsoda_context_t ctx = {
            .function = dydt,
            .data = appCtx,
            .neq = appCtx->nVariables,
            .state = 1,
    };
    //ctx.data = coefficients;

    lsoda_prepare(&ctx, &opt);

    for (iout =1; iout <= appCtx->totalSteps; iout++)
    {
        lsoda(&ctx, y, &t, tout);
        //printf(" at t= %12.4e y= %14.6e %14.6e %14.6e %14.6e %14.6e\n", t, y[0], y[1], y[2], y[3], y[4]);

        for(int i=0; i<appCtx->nVariables; i++) {
            int outIndex = appCtx->nVariables * iout + i;
            _yout[outIndex] = y[i];

        }

        //cout << "iout: " + to_string(iout) << endl;
        //cout << "t: " + to_string(t) << endl;
        //cout << "tout: " + to_string(tout) << endl;
        //cout << "t=" + to_string(t) +" y="+ vectorToString(y, 0, 4) <<endl;
        //cout << "" + to_string(t) +" "+ vectorToString(y, 0, 4) <<endl;

        if (ctx.state <= 0)
        { //todo: ver se devo abortar ou não.
            //todo: entender esse limite de passos e pq está sendo atingido mesmo com tamanho de passo pequeno
            // outputToFile("problematicInds.txt", vectorToString(appCtx->individual, 0, appCtx->IND_SIZE-1) + "\n", true);
            //cout << vectorToString(appCtx->individual, 0, appCtx->IND_SIZE-1)<<endl;
            printf("error istate = %d\n", ctx.state);
            for(int i=0; i<appCtx->nVariables; i++) {
                int outIndex = appCtx->nVariables * iout + i;
                _yout[outIndex] = INFINITY;
            }
            break;

            //cout << y[0] << endl;
            //cout << y[1] << endl;
            //cout << y[2] << endl;
            //cout << y[3] << endl;
            //cout << y[4] << endl;
            //break;
            //exit(0);
            //cout << "_----------------------------" << endl;
        }
        tout = tout + dt;
    }

    delete[] rtol;
    delete[] atol;
    delete[] y;
    lsoda_free(&ctx);

    return 0;
}

double GRNEDOHelpers::getMaxValue(double *values, int numElements)
{
    double maxValue = 0;
    for (int i = 0; i < numElements; i++)
    {
        if (values[i] > maxValue)
        {
            maxValue = values[i];
        }
    }

    return maxValue;
}

void GRNEDOHelpers::getMaxValues(double **data, double *outMaxValues, int numVariables, int numElements)
{

    for (int i = 1; i < numVariables + 1; i++)
    {
        outMaxValues[i - 1] = getMaxValue(data[i], numElements);
    }
}

double GRNEDOHelpers::getMaxValue(double *values, int start, int end)
{
    double maxValue = 0;
    for (int i = start; i <= end; i++)
    {
        if (values[i] > maxValue)
        {
            maxValue = values[i];
        }
    }

    return maxValue;
}

void GRNEDOHelpers::getMaxValues(double **data, double *outMaxValues, int numVariables, int start, int end)
{

    for (int i = 1; i < numVariables + 1; i++)
    {
        outMaxValues[i - 1] = getMaxValue(data[i], start, end);
    }
}

void GRNEDOHelpers::printContext(appContext* ctx){
    printf("IND_SIZE: %d\n", ctx->IND_SIZE);
    printf("MIN_K: %f\n", ctx->MIN_K);
    printf("MAX_K: %f\n", ctx->MAX_K);
    printf("MIN_N: %f\n", ctx->MIN_N);
    printf("MAX_N: %f\n", ctx->MAX_N);
    printf("MIN_TAU: %f\n", ctx->MIN_TAU);
    printf("MAX_TAU: %f\n", ctx->MAX_TAU);
    printf("MIN_STRATEGY: %f\n", ctx->MIN_STRATEGY);
    printf("MAX_STRATEGY: %f\n", ctx->MAX_STRATEGY);
    printf("TAU_SIZE: %d\n", ctx->TAU_SIZE);
    printf("N_SIZE: %d\n", ctx->N_SIZE);
    printf("K_SIZE: %d\n", ctx->K_SIZE);
    printf("nVariables: %d\n", ctx->nVariables);
    printf("nSteps: %d\n", ctx->nSteps);
    printf("setStart: %d\n", ctx->setStart);
    printf("setEnd: %d\n", ctx->setEnd);
    printf("dataSetSize: %d\n", ctx->dataSetSize);
    printf("trainingSetStart: %d\n", ctx->trainingSetStart);
    printf("trainingSetEnd: %d\n", ctx->trainingSetEnd);
    printf("trainingSteps: %d\n", ctx->trainingSteps);
    printf("validationSetStart: %d\n", ctx->validationSetStart);
    printf("validationSetEnd: %d\n", ctx->validationSetEnd);
    printf("validationSteps: %d\n", ctx->validationSteps);
    printf("testSetStart: %d\n", ctx->testSetStart);
    printf("testSetEnd: %d\n", ctx->testSetEnd);
    printf("testSteps: %d\n", ctx->testSteps);
    printf("mode: %d\n", ctx->mode);
    printf("tspan[0]: %f\n", ctx->tspan[0]);
    printf("tspan[1]: %f\n", ctx->tspan[1]);
    printf("maxValues: %s\n",  vectorToString(ctx->maxValues, 0,  ctx->nVariables-1).c_str());
    printf("y_0: %s\n",  vectorToString(ctx->y_0, 0,  ctx->nVariables-1).c_str());


}

void GRNEDOHelpers::twoBody(double t, double y[], double max[], double tau[], double n[], double k[], double yp[])
{

    yp[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[0]), n[0]) + pow(k[0], n[0]))) - (y[0] / max[0])) / tau[0];

    yp[1] = (((pow((y[0] / max[0]), (int)n[1])) / (pow((y[0] / max[0]), (int)n[1]) + pow(k[1], (int)n[1]))) -
             (y[1] / max[1])) /
            tau[1];

    yp[2] = (((pow((y[1] / max[1]), n[2])) / (pow((y[1] / max[1]), (int)n[2]) + pow(k[2], (int)n[2]))) -
             (y[2] / max[2])) /
            tau[2];

    yp[3] = (((pow((y[2] / max[2]), (int)n[3])) / (pow((y[2] / max[2]), (int)n[3]) + pow(k[3], (int)n[3]))) -
             (y[3] / max[3])) /
            tau[3];

    yp[4] = (((pow((y[3] / max[3]), (int)n[4])) / (pow((y[3] / max[3]), (int)n[4]) + pow(k[4], (int)n[4]))) -
             (y[4] / max[4])) /
            tau[4];
}

string GRNEDOHelpers::vectorToString(double *vec, int start, int end)
{
    string s = "";
    for (int i = start; i <= end; i++)
    {
        s.append(to_string(vec[i]) + " " );
    }

    return s;
}

void GRNEDOHelpers::printGRNVector(double **vec, int rows, int cols)
{
    //todo: parece que rows e cols estão invertidos, ver se renomeio as variáveis na assinatura
    for(int j = 0; j < cols; j++){
        for(int i = 0; i < rows; i++){
            cout << vec[i][j] << ",";
        }
        cout << "\n";
    }

}

int GRNEDOHelpers::twoBodyFixedLSODA(double t, double *y, double *ydot, void *data)
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

void GRNEDOHelpers::twoBodyFixed(double t, double y[], double *dim, double yp[])
{
    // todo: observar que os n devem ser avaliados como inteiros
    double max[] = {2.96, 1.8768, 1.0653, 1.0101, 1.4608};
    double tau[] = {1.25, 4, 1.02, 1.57, 3.43};
    double n[] = {13, 4, 3, 4, 16};
    double k[] = {0.72, 0.50, 0.45, 0.51, 0.52};

    yp[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[4]), n[0]) + pow(k[0], n[0]))) - (y[0] / max[0])) / tau[0];

    yp[1] = (((pow((y[0] / max[0]), n[1])) / (pow((y[0] / max[0]), n[1]) + pow(k[1], n[1]))) - (y[1] / max[1])) /
            tau[1];

    yp[2] = (((pow((y[1] / max[1]), n[2])) / (pow((y[1] / max[1]), n[2]) + pow(k[2], n[2]))) - (y[2] / max[2])) /
            tau[2];

    yp[3] = (((pow((y[2] / max[2]), n[3])) / (pow((y[2] / max[2]), n[3]) + pow(k[3], n[3]))) - (y[3] / max[3])) /
            tau[3];

    yp[4] = (((pow((y[3] / max[3]), n[4])) / (pow((y[3] / max[3]), n[4]) + pow(k[4], n[4]))) - (y[4] / max[4])) /
            tau[4];
}

//todo: renomear quando concluir
double GRNEDOHelpers::grn5EvaluationLSODA(void *ind, void* data)
{
    appContext* ctx = (appContext*)(data);
    double* _ind = (double *)ind;
    ctx->individual = _ind;
    lsodaWrapper(twoBody5VarLSODA, ctx, ctx->yout);

    double eval = difference(ctx->yout, ctx->expectedResult, ctx->nVariables, ctx->setStart, ctx->setEnd, ctx->nSteps);

    return eval;

}

double GRNEDOHelpers::grn5EvaluationRK4(void *ind, void* data)

{   appContext* ctx = (appContext*)(data);
    double* _ind = (double *)ind;
    ctx->individual = _ind;
    double *t = new double [ctx->totalSteps+1];
    rk4(twoBody5VarLSODA, ctx->tspan, ctx->y_0, ctx->totalSteps, ctx->nVariables, t, _ind, ctx);
    delete [] t;
    return difference(ctx->yout, ctx->expectedResult, ctx->nVariables, ctx->setStart, ctx->setEnd, ctx->nSteps);

}

double GRNEDOHelpers::grn10EvaluationLSODA(void *ind, void *data)
{
    appContext* ctx = (appContext*)(data);
    double* _ind = (double *)ind;
    ctx->individual = _ind;
    lsodaWrapper(twoBody10VarLSODA, ctx, ctx->yout);

    return difference(ctx->yout, ctx->expectedResult, ctx->nVariables, ctx->setStart, ctx->setEnd, ctx->nSteps);

}

double GRNEDOHelpers::grn10EvaluationRK4(void *ind, void* data)

{   appContext* ctx = (appContext*)(data);
    double* _ind = (double *)ind;
    ctx->individual = _ind;
    double *t = new double [ctx->totalSteps+1];
    rk4(twoBody10VarLSODA, ctx->tspan, ctx->y_0, ctx->totalSteps, ctx->nVariables, t, _ind, ctx);
    delete [] t;
    return difference(ctx->yout, ctx->expectedResult, ctx->nVariables, ctx->setStart, ctx->setEnd, ctx->nSteps);

}

