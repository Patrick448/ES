#include "GRNEDOHelpers.h"
#include "GRNSeries.h"
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

//Modelo do DadosZeEduardo.txt
int GRNEDOHelpers::grn5NCYCModel(double t, double *y, double *ydot, void *context) {
    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    //double x[] = {5.000000,5.000000,5.000000,0.100000,0.803348,0.100000,0.100000,1.000000,0.100000,0.100000,0.499194,0.100000,1.000000,1.000000,0.828747,0.181154,0.455022,0.100000,9.677138,4.894019,5.414821,20.363100,24.478013,18.839608,6.112694,8.810225,15.297329,16.423215,1.000000,13.097807,7.390248};
    //individual = x;
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];

    ydot[0] = (((((pow(y[2] / maxValues[2], (int) n[0])) /
                  (pow(y[2] / maxValues[2], (int) n[0]) + pow(k[0], (int) n[0]))) *
                 ((pow(y[4] / maxValues[4], (int) n[1])) /
                  (pow(y[4] / maxValues[4], (int) n[1]) + pow(k[1], (int) n[1]))))
                +
                (((pow(y[1] / maxValues[1], (int) n[2])) /
                  (pow(y[1] / maxValues[1], (int) n[2]) + pow(k[2], (int) n[2])))
                 *
                 ((pow(y[4] / maxValues[4], (int) n[1])) /
                  (pow(y[4] / maxValues[4], (int) n[1]) + pow(k[1], (int) n[1]))))
                +
                ((1 - ((pow(y[1] / maxValues[1], (int) n[2])) /
                       (pow(y[1] / maxValues[1], (int) n[2]) + pow(k[2], (int) n[2]))))
                 *
                 (1 - ((pow(y[2] / maxValues[2], (int) n[0])) /
                       (pow(y[2] / maxValues[2], (int) n[0]) + pow(k[0], (int) n[0]))))
                 *
                 (1 - ((pow(y[4] / maxValues[4], (int) n[1])) /
                       (pow(y[4] / maxValues[4], (int) n[1]) + pow(k[1], (int) n[1])))))) - (y[0] / maxValues[0])) /
              tau[0];


    ydot[1] = ((((pow(y[2] / maxValues[2], (int) n[3])) /
                 (pow(y[2] / maxValues[2], (int) n[3]) + pow(k[3], (int) n[3]))) +
                ((pow(y[4] / maxValues[4], (int) n[4])) /
                 (pow(y[4] / maxValues[4], (int) n[4]) + pow(k[4], (int) n[4])))) - (y[1] / maxValues[1])) / tau[1];


    ydot[2] = (((((pow(y[2] / maxValues[2], (int) n[6])) /
                  (pow(y[2] / maxValues[2], (int) n[6]) + pow(k[6], (int) n[6]))) *
                 ((pow(y[4] / maxValues[4], (int) n[7])) /
                  (pow(y[4] / maxValues[4], (int) n[7]) + pow(k[7], (int) n[7]))))
                +
                (((pow(y[1] / maxValues[1], (int) n[5])) /
                  (pow(y[1] / maxValues[1], (int) n[5]) + pow(k[5], (int) n[5])))
                 *
                 ((pow(y[4] / maxValues[4], (int) n[7])) /
                  (pow(y[4] / maxValues[4], (int) n[7]) + pow(k[7], (int) n[7]))))
                +
                ((1 - ((pow(y[1] / maxValues[1], (int) n[5])) /
                       (pow(y[1] / maxValues[1], (int) n[5]) + pow(k[5], (int) n[5]))))
                 *
                 (1 - ((pow(y[2] / maxValues[2], (int) n[6])) /
                       (pow(y[2] / maxValues[2], (int) n[6]) + pow(k[6], (int) n[6]))))
                 *
                 (1 - ((pow(y[4] / maxValues[4], (int) n[7])) /
                       (pow(y[4] / maxValues[4], (int) n[7]) + pow(k[7], (int) n[7])))))) - (y[2] / maxValues[2])) /
              tau[2];


    ydot[3] = ((((pow(y[2] / maxValues[2], (int) n[8])) /
                 (pow(y[2] / maxValues[2], (int) n[8]) + pow(k[8], (int) n[8]))) +
                ((pow(y[4] / maxValues[4], (int) n[9])) /
                 (pow(y[4] / maxValues[4], (int) n[9]) + pow(k[9], (int) n[9])))) - (y[3] / maxValues[3])) / tau[3];


    ydot[4] = (((((pow(y[2] / maxValues[2], (int) n[11])) /
                  (pow(y[2] / maxValues[2], (int) n[11]) + pow(k[11], (int) n[11]))) *
                 ((pow(y[4] / maxValues[4], (int) n[12])) /
                  (pow(y[4] / maxValues[4], (int) n[12]) + pow(k[12], (int) n[12]))))
                +
                (((pow(y[1] / maxValues[1], (int) n[10])) /
                  (pow(y[1] / maxValues[1], (int) n[10]) + pow(k[10], (int) n[10])))
                 *
                 ((pow(y[4] / maxValues[4], (int) n[12])) /
                  (pow(y[4] / maxValues[4], (int) n[12]) + pow(k[12], (int) n[12]))))
                +
                ((1 - ((pow(y[1] / maxValues[1], (int) n[10])) /
                       (pow(y[1] / maxValues[1], (int) n[10]) + pow(k[10], (int) n[10]))))
                 *
                 (1 - ((pow(y[2] / maxValues[2], (int) n[11])) /
                       (pow(y[2] / maxValues[2], (int) n[11]) + pow(k[11], (int) n[11]))))
                 *
                 (1 - ((pow(y[4] / maxValues[4], (int) n[12])) /
                       (pow(y[4] / maxValues[4], (int) n[12]) + pow(k[12], (int) n[12])))))) - (y[4] / maxValues[4])) /
              tau[4];

    return 0;
}

int GRNEDOHelpers::grn5Model(double t, double *y, double *ydot, void *context) {
    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];


    ydot[0] = ((1 - (pow((y[4] / maxValues[4]), (int) n[0])) /
                    (pow((y[4] / maxValues[4]), (int) n[0]) + pow(k[0], (int) n[0]))) -
               (y[0] / maxValues[0])) /
              tau[0];

    ydot[1] = (((pow((y[0] / maxValues[0]), (int) n[1])) /
                (pow((y[0] / maxValues[0]), (int) n[1]) + pow(k[1], (int) n[1]))) -
               (y[1] / maxValues[1])) /
              tau[1];

    ydot[2] = (((pow((y[1] / maxValues[1]), (int) n[2])) /
                (pow((y[1] / maxValues[1]), (int) n[2]) + pow(k[2], (int) n[2]))) -
               (y[2] / maxValues[2])) /
              tau[2];

    ydot[3] = (((pow((y[2] / maxValues[2]), (int) n[3])) /
                (pow((y[2] / maxValues[2]), (int) n[3]) + pow(k[3], (int) n[3]))) -
               (y[3] / maxValues[3])) /
              tau[3];

    ydot[4] =
            ((((pow(y[1] / maxValues[1], (int) n[4]) / (pow(y[1] / maxValues[1], (int) n[4]) + pow(k[4], (int) n[4]))) *
               (
                       pow(y[3] / maxValues[3], (int) n[5]) /
                       (pow(y[3] / maxValues[3], (int) n[5]) + pow(k[5], (int) n[5])))) + (
                      (pow(y[3] / maxValues[3], (int) n[5]) /
                       (pow(y[3] / maxValues[3], (int) n[5]) + pow(k[5], (int) n[5]))) * (
                              pow(y[4] / maxValues[4], (int) n[6]) /
                              (pow(y[4] / maxValues[4], (int) n[6]) + pow(k[6], (int) n[6]))))) - (
                     y[4] / maxValues[4])) / tau[4];

    return 0;
}

int GRNEDOHelpers::grn5NewModel(double t, double *y, double *ydot, void *context) {
    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];

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

int GRNEDOHelpers::grn10Model(double t, double *y, double *ydot, void *context) {

    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];

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
    double tauA = individual[0];
    double tauB = individual[1];
    double tauC = individual[2];
    double tauD = individual[3];
    double tauE = individual[4];
    double tauF = individual[5];
    double tauG = individual[6];
    double tauH = individual[7];
    double tauI = individual[8];
    double tauJ = individual[9];
    double kAJ = individual[10];
    double kBE = individual[11];
    double kCB = individual[12];
    double kCF = individual[13];
    double kCA = individual[14];
    double kDF = individual[15];
    double kEJ = individual[16];
    double kFA = individual[17];
    double kGB = individual[18];
    double kGF = individual[19];
    double kGA = individual[20];
    double kHF = individual[21];
    double kIG = individual[22];
    double kIH = individual[23];
    double kJI = individual[24];
    double nAJ = (int) individual[25];
    double nBE = (int) individual[26];
    double nCB = (int) individual[27];
    double nCF = (int) individual[28];
    double nCA = (int) individual[29];
    double nDF = (int) individual[30];
    double nEJ = (int) individual[31];
    double nFA = (int) individual[32];
    double nGB = (int) individual[33];
    double nGF = (int) individual[34];
    double nGA = (int) individual[35];
    double nHF = (int) individual[36];
    double nIG = (int) individual[37];
    double nIH = (int) individual[38];
    double nJI = (int) individual[39];

    ydot[0] =
            ((1 - pow(y[9] / maximo_J, nAJ) / (pow(y[9] / maximo_J, nAJ) + pow(kAJ, nAJ))) - (y[0] / maximo_A)) / tauA;

    ydot[1] = (pow(y[4] / maximo_E, nBE) / (pow(y[4] / maximo_E, nBE) + pow(kBE, nBE)) - (y[1] / maximo_B)) / tauB;

    ydot[2] = (pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB)) *
               (1 - pow(y[5] / maximo_F, nCF) / (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF))) *
               (1 - pow(y[0] / maximo_A, nCA) / (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA))) +
               (1 - pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB))) *
               pow(y[5] / maximo_F, nCF) /
               (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF)) *
               (1 - pow(y[0] / maximo_A, nCA) / (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA))) +
               (1 - pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB))) *
               (1 - pow(y[5] / maximo_F, nCF) / (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF))) *
               pow(y[0] / maximo_A, nCA) /
               (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) +
               pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB)) *
               (1 - pow(y[5] / maximo_F, nCF) / (pow(y[5] / maximo_F, nCF) + pow(kCF, nCF))) *
               pow(y[0] / maximo_A, nCA) /
               (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) +
               (1 - pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB))) *
               pow(y[5] / maximo_F, nCF) /
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
               (1 - pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB))) *
               pow(y[5] / maximo_F, nGF) /
               (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF)) *
               (1 - pow(y[0] / maximo_A, nGA) / (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA))) +
               (1 - pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB))) *
               (1 - pow(y[5] / maximo_F, nGF) / (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF))) *
               pow(y[0] / maximo_A, nGA) /
               (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) +
               pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB)) *
               (1 - pow(y[5] / maximo_F, nGF) / (pow(y[5] / maximo_F, nGF) + pow(kGF, nGF))) *
               pow(y[0] / maximo_A, nGA) /
               (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) +
               (1 - pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB))) *
               pow(y[5] / maximo_F, nGF) /
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

//modelo do arquivo genesABCD
int GRNEDOHelpers::grn4NCYCModel(double t, double *y, double *ydot, void *context) {
    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];


    ydot[0] = (((1 - ((pow(y[1] / maxValues[1], (int) n[0])) /
                      (pow(y[1] / maxValues[1], (int) n[0]) + pow(k[0], (int) n[0]))))
                *
                (((pow(y[2] / maxValues[2], (int) n[1])) /
                  (pow(y[2] / maxValues[2], (int) n[1]) + pow(k[1], (int) n[1])))))
               - (y[0] / maxValues[0])) / tau[0];


    ydot[1] = (((((pow(y[0] / maxValues[0], (int) n[2])) /
                  (pow(y[0] / maxValues[0], (int) n[2]) + pow(k[2], (int) n[2]))))
                +
                (((pow(y[1] / maxValues[1], (int) n[3])) /
                  (pow(y[1] / maxValues[1], (int) n[3]) + pow(k[3], (int) n[3]))))
                +
                (1 - ((pow(y[2] / maxValues[2], (int) n[4])) /
                      (pow(y[2] / maxValues[2], (int) n[4]) + pow(k[4], (int) n[4]))))) - (y[1] / maxValues[1])) /
              tau[1];


    ydot[2] = (((((pow(y[0] / maxValues[0], (int) n[5])) /
                  (pow(y[0] / maxValues[0], (int) n[5]) + pow(k[5], (int) n[5]))))
                +
                (((pow(y[2] / maxValues[2], (int) n[6])) /
                  (pow(y[2] / maxValues[2], (int) n[6]) + pow(k[6], (int) n[6]))))) - (y[2] / maxValues[2])) / tau[2];


    ydot[3] = ((((1 - ((pow(y[0] / maxValues[0], (int) n[7])) /
                       (pow(y[0] / maxValues[0], (int) n[7]) + pow(k[7], (int) n[7]))))
                 *
                 (((pow(y[1] / maxValues[1], (int) n[8])) /
                   (pow(y[1] / maxValues[1], (int) n[8]) + pow(k[8], (int) n[8]))))
                 *
                 (1 - ((pow(y[3] / maxValues[3], (int) n[10])) /
                       (pow(y[3] / maxValues[3], (int) n[10]) + pow(k[10], (int) n[10])))))
                +
                (
                        (1 - ((pow(y[0] / maxValues[0], (int) n[7])) /
                              (pow(y[0] / maxValues[0], (int) n[7]) + pow(k[7], (int) n[7]))))
                        *
                        (1 - ((pow(y[1] / maxValues[1], (int) n[8])) /
                              (pow(y[1] / maxValues[1], (int) n[8]) + pow(k[8], (int) n[8]))))
                        *
                        (1 - ((pow(y[2] / maxValues[2], (int) n[9])) /
                              (pow(y[2] / maxValues[2], (int) n[9]) + pow(k[9], (int) n[9]))))
                        *
                        (((pow(y[3] / maxValues[3], (int) n[10])) /
                          (pow(y[3] / maxValues[3], (int) n[10]) + pow(k[10], (int) n[10]))))
                )) - (y[3] / maxValues[3])) / tau[3];

    return 0;
}

int GRNEDOHelpers::grn10NewModel(double t, double *y, double *ydot, void *context) {

    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];

    ydot[0] = ((1 - ((pow(y[4] / maxValues[4], (int) n[0])) /
                     (pow(y[4] / maxValues[4], (int) n[0]) + pow(k[0], (int) n[0])))) - (y[0] / maxValues[0])) / tau[0];
    ydot[1] = ((1 - ((pow(y[4] / maxValues[4], (int) n[1])) /
                     (pow(y[4] / maxValues[4], (int) n[1]) + pow(k[1], (int) n[1])))) - (y[1] / maxValues[1])) / tau[1];
    ydot[2] = ((((pow(y[0] / maxValues[0], (int) n[2])) /
                 (pow(y[0] / maxValues[0], (int) n[2]) + pow(k[2], (int) n[2]))) *
                ((pow(y[6] / maxValues[6], (int) n[3])) /
                 (pow(y[6] / maxValues[6], (int) n[3]) + pow(k[3], (int) n[3])))) - (y[2] / maxValues[2])) / tau[2];
    ydot[3] = ((((pow(y[1] / maxValues[1], (int) n[4])) /
                 (pow(y[1] / maxValues[1], (int) n[4]) + pow(k[4], (int) n[4]))) +
                ((pow(y[3] / maxValues[3], (int) n[5])) /
                 (pow(y[3] / maxValues[3], (int) n[5]) + pow(k[5], (int) n[5])))) - (y[3] / maxValues[3])) / tau[3];
    ydot[4] = ((((pow(y[0] / maxValues[0], (int) n[6])) /
                 (pow(y[0] / maxValues[0], (int) n[6]) + pow(k[6], (int) n[6]))) *
                ((pow(y[6] / maxValues[6], (int) n[7])) /
                 (pow(y[6] / maxValues[6], (int) n[7]) + pow(k[7], (int) n[7]))) +
                ((pow(y[3] / maxValues[3], (int) n[8])) /
                 (pow(y[3] / maxValues[3], (int) n[8]) + pow(k[8], (int) n[8])))) - (y[4] / maxValues[4])) / tau[4];
    ydot[5] = ((1 - ((pow(y[4] / maxValues[4], (int) n[9])) /
                     (pow(y[4] / maxValues[4], (int) n[9]) + pow(k[9], (int) n[9])))) - (y[5] / maxValues[5])) / tau[5];
    ydot[6] = ((1 - ((pow(y[4] / maxValues[4], (int) n[10])) /
                     (pow(y[4] / maxValues[4], (int) n[10]) + pow(k[10], (int) n[10])))) - (y[6] / maxValues[6])) /
              tau[6];
    ydot[7] = ((((pow(y[4] / maxValues[4], (int) n[11])) /
                 (pow(y[4] / maxValues[4], (int) n[11]) + pow(k[11], (int) n[11]))) +
                ((pow(y[5] / maxValues[5], (int) n[12])) /
                 (pow(y[5] / maxValues[5], (int) n[12]) + pow(k[12], (int) n[12])))) - (y[7] / maxValues[7])) / tau[7];
    ydot[8] = ((((pow(y[0] / maxValues[0], (int) n[13])) /
                 (pow(y[0] / maxValues[0], (int) n[13]) + pow(k[13], (int) n[13]))) *
                ((pow(y[6] / maxValues[6], (int) n[14])) /
                 (pow(y[6] / maxValues[6], (int) n[14]) + pow(k[14], (int) n[14]))) +
                ((pow(y[3] / maxValues[3], (int) n[15])) /
                 (pow(y[3] / maxValues[3], (int) n[15]) + pow(k[15], (int) n[15])))) - (y[8] / maxValues[8])) / tau[8];
    ydot[9] = ((((pow(y[0] / maxValues[0], (int) n[16])) /
                 (pow(y[0] / maxValues[0], (int) n[16]) + pow(k[16], (int) n[16]))) *
                ((pow(y[6] / maxValues[6], (int) n[17])) /
                 (pow(y[6] / maxValues[6], (int) n[17]) + pow(k[17], (int) n[17]))) +
                ((pow(y[3] / maxValues[3], (int) n[18])) /
                 (pow(y[3] / maxValues[3], (int) n[18]) + pow(k[18], (int) n[18])))) - (y[9] / maxValues[9])) / tau[9];


    return 0;
}

int GRNEDOHelpers::grn10New2Model(double t, double *y, double *ydot, void *context) {

    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];

    //tau[0], k[0], (int)n[0], k[1], (int)n[1]

    ydot[0] = ((((pow(y[7] / maxValues[7], (int) n[0])) /
                 (pow(y[7] / maxValues[7], (int) n[0]) + pow(k[0], (int) n[0]))) +
                ((pow(y[1] / maxValues[1], (int) n[1])) /
                 (pow(y[1] / maxValues[1], (int) n[1]) + pow(k[1], (int) n[1])))) - (y[0] / maxValues[0])) / tau[0];

//tau[1], k[2], (int)n[2], k[3], (int)n[3]

    ydot[1] = ((((1 - ((pow(y[1] / maxValues[1], (int) n[2])) /
                       (pow(y[1] / maxValues[1], (int) n[2]) + pow(k[2], (int) n[2])))) *
                 ((pow(y[5] / maxValues[5], (int) n[3])) /
                  (pow(y[5] / maxValues[5], (int) n[3]) + pow(k[3], (int) n[3])))) + ((1 - ((pow(y[5] / maxValues[5],
                                                                                                 (int) n[3])) /
                                                                                            (pow(y[5] / maxValues[5],
                                                                                                 (int) n[3]) +
                                                                                             pow(k[3], (int) n[3])))) *
                                                                                      ((pow(y[1] / maxValues[1],
                                                                                            (int) n[2])) /
                                                                                       (pow(y[1] / maxValues[1],
                                                                                            (int) n[2]) +
                                                                                        pow(k[2], (int) n[2]))))) -
               (y[1] / maxValues[1])) / tau[1];

//tau[2], k[4], k[5], (int)n[4], (int)n[5]

    ydot[2] = ((((pow(y[0] / maxValues[0], (int) n[4])) /
                 (pow(y[0] / maxValues[0], (int) n[4]) + pow(k[4], (int) n[4]))) +
                ((pow(y[2] / maxValues[2], (int) n[5])) /
                 (pow(y[2] / maxValues[2], (int) n[5]) + pow(k[5], (int) n[5])))) - (y[2] / maxValues[2])) / tau[2];

//tau[3], k[6], k[7], k[8], k[9], (int)n[6], (int)n[7], (int)n[8], (int)n[9]

    ydot[3] =
            (((((pow(y[1] / maxValues[1], (int) n[6])) / (pow(y[1] / maxValues[1], (int) n[6]) + pow(k[6], (int) n[6])))
               *
               (1 - ((pow(y[4] / maxValues[4], (int) n[8])) /
                     (pow(y[4] / maxValues[4], (int) n[8]) + pow(k[8], (int) n[8]))))
               *
               ((pow(y[7] / maxValues[7], (int) n[9])) /
                (pow(y[7] / maxValues[7], (int) n[9]) + pow(k[9], (int) n[9]))))

              +

              (((pow(y[1] / maxValues[1], (int) n[6])) / (pow(y[1] / maxValues[1], (int) n[6]) + pow(k[6], (int) n[6])))
               *
               (1 - ((pow(y[3] / maxValues[3], (int) n[7])) /
                     (pow(y[3] / maxValues[3], (int) n[7]) + pow(k[7], (int) n[7]))))
               *
               ((pow(y[7] / maxValues[7], (int) n[9])) /
                (pow(y[7] / maxValues[7], (int) n[9]) + pow(k[9], (int) n[9]))))

              +

              (((pow(y[1] / maxValues[1], (int) n[6])) / (pow(y[1] / maxValues[1], (int) n[6]) + pow(k[6], (int) n[6])))
               *
               ((pow(y[3] / maxValues[3], (int) n[7])) / (pow(y[3] / maxValues[3], (int) n[7]) + pow(k[7], (int) n[7])))
               *
               ((pow(y[4] / maxValues[4], (int) n[8])) / (pow(y[4] / maxValues[4], (int) n[8]) + pow(k[8], (int) n[8])))
               *
               (1 - ((pow(y[7] / maxValues[7], (int) n[9])) /
                     (pow(y[7] / maxValues[7], (int) n[9]) + pow(k[9], (int) n[9])))))) - (y[3] / maxValues[3])) /
            tau[3];

//tau[4], k[10], (int)n[10]

    ydot[4] = (((pow(y[7] / maxValues[7], (int) n[10])) /
                (pow(y[7] / maxValues[7], (int) n[10]) + pow(k[10], (int) n[10]))) - (y[4] / maxValues[4])) / tau[4];

//tau[5], k[11], k[12], k[13], (int)n[11], (int)n[12], (int)n[13]

    ydot[5] = ((((1 - ((pow(y[0] / maxValues[0], (int) n[11])) /
                       (pow(y[0] / maxValues[0], (int) n[11]) + pow(k[11], (int) n[11]))))
                 *
                 (1 - ((pow(y[2] / maxValues[2], (int) n[12])) /
                       (pow(y[2] / maxValues[2], (int) n[12]) + pow(k[12], (int) n[12]))))
                 *
                 ((pow(y[8] / maxValues[8], (int) n[13])) /
                  (pow(y[8] / maxValues[8], (int) n[13]) + pow(k[13], (int) n[13]))))

                +

                (((pow(y[0] / maxValues[0], (int) n[11])) /
                  (pow(y[0] / maxValues[0], (int) n[11]) + pow(k[11], (int) n[11])))
                 *
                 (1 - ((pow(y[8] / maxValues[8], (int) n[13])) /
                       (pow(y[8] / maxValues[8], (int) n[13]) + pow(k[13], (int) n[13])))))

                +

                (((pow(y[2] / maxValues[2], (int) n[12])) /
                  (pow(y[2] / maxValues[2], (int) n[12]) + pow(k[12], (int) n[12])))
                 *
                 (1 - ((pow(y[8] / maxValues[8], (int) n[13])) /
                       (pow(y[8] / maxValues[8], (int) n[13]) + pow(k[13], (int) n[13])))))) - (y[5] / maxValues[5])) /
              tau[5];

//tau[6], k[14], k[15], k[16], (int)n[14], (int)n[15], (int)n[16]

    ydot[6] = ((((1 - ((pow(y[4] / maxValues[4], (int) n[15])) /
                       (pow(y[4] / maxValues[4], (int) n[15]) + pow(k[15], (int) n[15]))))
                 *
                 ((pow(y[7] / maxValues[7], (int) n[16])) /
                  (pow(y[7] / maxValues[7], (int) n[16]) + pow(k[16], (int) n[16]))))

                +

                ((1 - ((pow(y[3] / maxValues[3], (int) n[14])) /
                       (pow(y[3] / maxValues[3], (int) n[14]) + pow(k[14], (int) n[14]))))
                 *
                 ((pow(y[7] / maxValues[7], (int) n[16])) /
                  (pow(y[7] / maxValues[7], (int) n[16]) + pow(k[16], (int) n[16]))))

                +

                (((pow(y[3] / maxValues[3], (int) n[14])) /
                  (pow(y[3] / maxValues[3], (int) n[14]) + pow(k[14], (int) n[14])))
                 *
                 ((pow(y[4] / maxValues[4], (int) n[15])) /
                  (pow(y[4] / maxValues[4], (int) n[15]) + pow(k[15], (int) n[15])))
                 *
                 ((pow(y[7] / maxValues[7], (int) n[16])) /
                  (pow(y[7] / maxValues[7], (int) n[16]) + pow(k[16], (int) n[16]))))) - (y[6] / maxValues[6])) /
              tau[6];

//tau[7], k[17], (int)n[17], k[18], (int)n[18]

    ydot[7] = ((((pow(y[7] / maxValues[7], (int) n[17])) /
                 (pow(y[7] / maxValues[7], (int) n[17]) + pow(k[17], (int) n[17]))) +
                ((pow(y[1] / maxValues[1], (int) n[18])) /
                 (pow(y[1] / maxValues[1], (int) n[18]) + pow(k[18], (int) n[18])))) - (y[7] / maxValues[7])) / tau[7];

//tau[8], k[19], k[20], k[21], k[22], (int)n[19], (int)n[20], (int)n[21], (int)n[22]

    ydot[8] = (((((pow(y[1] / maxValues[1], (int) n[19])) /
                  (pow(y[1] / maxValues[1], (int) n[19]) + pow(k[19], (int) n[19])))
                 *
                 (1 - ((pow(y[4] / maxValues[4], (int) n[21])) /
                       (pow(y[4] / maxValues[4], (int) n[21]) + pow(k[21], (int) n[21]))))
                 *
                 ((pow(y[7] / maxValues[7], (int) n[22])) /
                  (pow(y[7] / maxValues[7], (int) n[22]) + pow(k[22], (int) n[22]))))

                +

                (((pow(y[1] / maxValues[1], (int) n[19])) /
                  (pow(y[1] / maxValues[1], (int) n[19]) + pow(k[19], (int) n[19])))
                 *
                 (1 - ((pow(y[3] / maxValues[3], (int) n[20])) /
                       (pow(y[3] / maxValues[3], (int) n[20]) + pow(k[20], (int) n[20]))))
                 *
                 ((pow(y[7] / maxValues[7], (int) n[22])) /
                  (pow(y[7] / maxValues[7], (int) n[22]) + pow(k[22], (int) n[22]))))

                +

                (((pow(y[1] / maxValues[1], (int) n[19])) /
                  (pow(y[1] / maxValues[1], (int) n[19]) + pow(k[19], (int) n[19])))
                 *
                 ((pow(y[3] / maxValues[3], (int) n[20])) /
                  (pow(y[3] / maxValues[3], (int) n[20]) + pow(k[20], (int) n[20])))
                 *
                 ((pow(y[4] / maxValues[4], (int) n[21])) /
                  (pow(y[4] / maxValues[4], (int) n[21]) + pow(k[21], (int) n[21])))
                 *
                 (1 - ((pow(y[7] / maxValues[7], (int) n[22])) /
                       (pow(y[7] / maxValues[7], (int) n[22]) + pow(k[22], (int) n[22])))))) - (y[8] / maxValues[8])) /
              tau[8];

//tau[9], k[23], (int)n[23], k[24], (int)n[24]

    ydot[9] = ((((1 - ((pow(y[1] / maxValues[1], (int) n[23])) /
                       (pow(y[1] / maxValues[1], (int) n[23]) + pow(k[23], (int) n[23])))) *
                 ((pow(y[5] / maxValues[5], (int) n[24])) /
                  (pow(y[5] / maxValues[5], (int) n[24]) + pow(k[24], (int) n[24])))) + ((1 - ((pow(y[5] / maxValues[5],
                                                                                                    (int) n[24])) /
                                                                                               (pow(y[5] / maxValues[5],
                                                                                                    (int) n[24]) +
                                                                                                pow(k[24],
                                                                                                    (int) n[24])))) *
                                                                                         ((pow(y[1] / maxValues[1],
                                                                                               (int) n[23])) /
                                                                                          (pow(y[1] / maxValues[1],
                                                                                               (int) n[23]) +
                                                                                           pow(k[23], (int) n[23]))))) -
               (y[9] / maxValues[9])) / tau[9];


    return 0;
}

//todo: melhorar vetores, padronizar o formato e melhorar os acessos
double GRNEDOHelpers::difference(double *actual, double **expected, int numElements, int numVariables, int granularity) {
    double difTotal = 0.0;
    for (int i = 0; i < numVariables; i++) {
        for (int j = 0; j <= numElements; j++) {
            int index = granularity * j * numVariables + i;
            difTotal += fabs(actual[index] - expected[i][j]);
            // cout << actual[index] << " ";
        }
        //cout << endl;
    }
    if (isnan(difTotal)) {
        return DBL_MAX;
    }
    return difTotal;
}


double GRNEDOHelpers::differenceNormalized(double *actual, double **expected, int numElements, int numVariables, int granularity, double* maxValues, double* minValues) {
    double difTotal = 0.0;
    for (int i = 0; i < numVariables; i++) {
        for (int j = 0; j <= numElements; j++) {
            int index = granularity * j * numVariables + i;
            double dif = actual[index] - expected[i][j];
            double normalizedDif = fabs(dif / (maxValues[i] - minValues[i]));
            difTotal += normalizedDif;
        }
    }
    if (isnan(difTotal)) {
        return DBL_MAX;
    }
    return difTotal;
}

double GRNEDOHelpers::differenceNormalized2(double *actual, double **expected, int numElements, int numVariables, int granularity, double* maxValues, double* minValues) {
    double difTotal = 0.0;
    for (int i = 0; i < numVariables; i++) {
        for (int j = 0; j <= numElements; j++) {
            int index = granularity * j * numVariables + i;
            double dif = actual[index] - expected[i][j];
            double normalizedDif = fabs(dif / maxValues[i]);
            difTotal += normalizedDif;
        }
    }
    if (isnan(difTotal)) {
        return DBL_MAX;
    }
    return difTotal;
}

double GRNEDOHelpers::lsodaWrapper(int dydt(double t, double *y, double *ydot, void *data), double *tspan, double *y_0,
                                   int totalSteps, int nVariables, double *times, double *_yout, void *context) {

    // todo: tentar colocar essa alocação fora da função
    //  esses vetores serão alocados toda vez, e essa função será chamada
    //  a cada avaliação de indivíduo

    double *atol = new double[nVariables];
    double *rtol = new double[nVariables];
    double t, tout, dt;

    double *y = new double[nVariables];
    int iout;

    for (int i = 0; i < nVariables; i++) {
        y[i] = y_0[i];
        rtol[i] = atol[i] = 1.49012e-4;
        _yout[i] = y_0[i];
    }

    t = tspan[0];
    dt = (tspan[1] - tspan[0]) / (double) (totalSteps);
    tout = t + dt;

    struct lsoda_opt_t opt = {0};
    opt.ixpr = 0;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.itask = 1;
    opt.mxstep = 500;


    struct lsoda_context_t ctx = {
            .function = dydt,
            .data = context,
            .neq = nVariables,
            .state = 1,
    };
    //ctx.data = coefficients;

    lsoda_prepare(&ctx, &opt);

    for (iout = 1; iout <= totalSteps; iout++) {
        lsoda(&ctx, y, &t, tout);
        //printf(" at t= %12.4e y= %14.6e %14.6e %14.6e %14.6e %14.6e\n", t, y[0], y[1], y[2], y[3], y[4]);

        for (int i = 0; i < nVariables; i++) {
            int outIndex = nVariables * iout + i;
            _yout[outIndex] = y[i];

        }

        if (ctx.state <= 0) { //todo: ver se devo abortar ou não.
            //todo: entender esse limite de passos e pq está sendo atingido mesmo com tamanho de passo pequeno
            // outputToFile("problematicInds.txt", vectorToString(appCtx->individual, 0, appCtx->IND_SIZE-1) + "\n", true);
            //cout << vectorToString(appCtx->individual, 0, appCtx->IND_SIZE-1)<<endl;
            printf("error istate = %d\n", ctx.state);
            for (int i = 0; i < nVariables; i++) {
                int outIndex = nVariables * iout + i;
                _yout[outIndex] = INFINITY;
            }
            break;
        }
        tout = tout + dt;
    }

    delete[] rtol;
    delete[] atol;
    delete[] y;
    lsoda_free(&ctx);

    return 0;
}

string GRNEDOHelpers::vectorToString(double *vec, int start, int end) {
    string s = "";
    for (int i = start; i <= end; i++) {
        s.append(to_string(vec[i]) + " ");
    }
    s.pop_back();

    return s;
}

double GRNEDOHelpers::lsodaWrapperTest(int dydt(double t, double *y, double *ydot, void *data), double *tspan, double *y_0,
                                int totalSteps, int nVariables, double *times, double *_yout, void *context) {

    // todo: tentar colocar essa alocação fora da função
    //  esses vetores serão alocados toda vez, e essa função será chamada
    //  a cada avaliação de indivíduo

    double *atol = new double[nVariables];
    double *rtol = new double[nVariables];
    double t, tout, dt;

    double *y = new double[nVariables];
    int iout;

    for (int i = 0; i < nVariables; i++) {
        y[i] = y_0[i];
        rtol[i] = atol[i] = 1.49012e-4;
        _yout[i] = y_0[i];
    }

    t = times[0];
    //dt = (tspan[1] - tspan[0]) / (double)(totalSteps);
    tout = times[1];

    struct lsoda_opt_t opt = {0};
    opt.ixpr = 0;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.itask = 1;
    opt.mxstep = 500;

    struct lsoda_context_t ctx = {.function = dydt, .data = context, .neq = nVariables, .state = 1};
    lsoda_prepare(&ctx, &opt);

    for (iout = 1; iout <= totalSteps; iout++) {
        lsoda(&ctx, y, &t, tout);
        //printf(" at t= %12.4e y= %14.6e %14.6e %14.6e %14.6e %14.6e\n", t, y[0], y[1], y[2], y[3], y[4]);

        for (int i = 0; i < nVariables; i++) {
            int outIndex = nVariables * iout + i;
            _yout[outIndex] = y[i];
        }

        if (ctx.state <= 0) { //todo: ver se devo abortar ou não.
            //todo: entender esse limite de passos e pq está sendo atingido mesmo com tamanho de passo pequeno
            // outputToFile("problematicInds.txt", vectorToString(appCtx->individual, 0, appCtx->IND_SIZE-1) + "\n", true);
            //cout << vectorToString(appCtx->individual, 0, appCtx->IND_SIZE-1)<<endl;
            printf("error istate = %d\n", ctx.state);
            appContext *appCtx = (appContext *) context;
            printf("coefs: %s\n", vectorToString(appCtx->individual,0, appCtx->description->IND_SIZE-1).c_str());
            for (int i = 0; i < nVariables; i++) {
                int outIndex = nVariables * iout + i;
                _yout[outIndex] = INFINITY;
            }
            break;
        }
        if (iout < totalSteps) {
            tout = times[iout + 1];
        }

    }

    delete[] rtol;
    delete[] atol;
    delete[] y;
    lsoda_free(&ctx);

    return 0;
}


void GRNEDOHelpers::printGRNVector(double **vec, int rows, int cols) {
    //todo: parece que rows e cols estão invertidos, ver se renomeio as variáveis na assinatura
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            cout << vec[i][j] << ",";
        }
        cout << "\n";
    }

}

//todo: generalizar função para qualquer modelo
double GRNEDOHelpers::grnEvaluationLSODA(void *individual, void *context) {
    appContext *ctx = (appContext *) (context);
    GRNSeries *evalSeries = (GRNSeries *) (ctx->series);
    double *_ind = (double *) individual;
    ctx->individual = _ind;

    //todo: ver se volto granularity para o context
    int granularity = 1;
    int totalSteps = (evalSeries->getNumTimeSteps() - 1) * granularity;
    //double* yout = new double[(evalSeries->getNumTimeSteps()) * evalSeries->getNumVariables()];
    double *yout = new double[(totalSteps + 1) * evalSeries->getNumVariables()];

    double tspan[]{evalSeries->getStartTime(), evalSeries->getEndTime()};
    double **expectedResult = &evalSeries->getVectors()[1];
    double *t = &evalSeries->getVectors()[0][0];
    int nVariables = evalSeries->getNumVariables();
    double *y_0 = evalSeries->getInitialValues();

    lsodaWrapperTest(ctx->description->modelFunction, tspan, y_0, totalSteps, nVariables, t, yout, ctx);

    double eval = difference(yout, expectedResult, evalSeries->getNumTimeSteps() - 1, nVariables, granularity/*,
                                       evalSeries->getMaxValues(), evalSeries->getMinValues()*/);

    delete[] yout;
    return eval;

}

//todo: generalizar função para qualquer modelo
double GRNEDOHelpers::grnEvaluationRK4(void *individual, void *context) {
    appContext *ctx = (appContext *) (context);
    GRNSeries *evalSeries = (GRNSeries *) (ctx->series);
    double *_ind = (double *) individual;
    ctx->individual = _ind;

    //todo: ver se volto granularity para o context
    int granularity = 20;
    int totalSteps = (evalSeries->getNumTimeSteps() - 1) * granularity;
    double *yout = new double[(totalSteps + 1) * evalSeries->getNumVariables()];
    double tspan[]{evalSeries->getStartTime(), evalSeries->getEndTime()};
    double **expectedResult = &evalSeries->getVectors()[1];
    int nVariables = evalSeries->getNumVariables();
    double *y_0 = evalSeries->getInitialValues();
    double *t = new double[totalSteps + 1];

    rk4(ctx->description->modelFunction, tspan, y_0, totalSteps, nVariables, t, yout, ctx);

    double eval = difference(yout, expectedResult, evalSeries->getNumTimeSteps() - 1, nVariables, granularity);

    delete[] yout;
    delete[] t;

    return eval;

}

/// helper for outputting text to file
void outputToFile(string path, string text, bool append) {
    ofstream outputf;

    if (append) {
        outputf.open(path, std::ios_base::app);
    } else {
        outputf.open(path);
    }

    outputf << text;
    outputf.close();
}

//todo: mudar de void* para double* (individual)
void GRNEDOHelpers::printODEIntSeries(void *individual, void *context, const string &outputFile, int solver) {
    appContext *ctx = (appContext *) (context);
    GRNSeries *evalSeries = (GRNSeries *) (ctx->series);
    double *_ind = (double *) individual;
    ctx->individual = _ind;

    //todo: ver se volto granularity para o context
    int granularity = 1;
    int totalSteps = (evalSeries->getNumTimeSteps() - 1) * granularity;
    double *yout = new double[(totalSteps + 1) * evalSeries->getNumVariables()];
    double *t = new double[totalSteps + 1];

    double tspan[]{evalSeries->getStartTime(), evalSeries->getEndTime()};
    double **expectedResult = &evalSeries->getVectors()[1];
    int nVariables = evalSeries->getNumVariables();
    double *y_0 = evalSeries->getInitialValues();
    double *times = &evalSeries->getVectors()[0][0];

    if (solver == 1) {
        rk4(ctx->description->modelFunction, tspan, y_0, totalSteps, nVariables, t, yout, ctx);
    } else if (solver == 0) {
        lsodaWrapperTest(ctx->description->modelFunction, tspan, y_0, totalSteps, nVariables, times, yout, ctx);
    }

    GRNSeries resultSeries = GRNSeries(evalSeries->getNumTimeSteps(), evalSeries->getNumVariables(), yout, times,
                                       granularity);
    if (!outputFile.empty()) {
        outputToFile(outputFile, resultSeries.toString(), false);
    } else {
        cout << resultSeries.toString() << endl;
    }
    delete[] yout;
    delete[] t;

}