//
// Created by patrick on 04/01/24.
//

#ifndef ES_GRN10MODEL_H
#define ES_GRN10MODEL_H
#include "dependencies.h"

class GRN10Model {

public:
    static const int TAU_SIZE = 10;
    static const int N_SIZE = 15;
    static const int K_SIZE = 15;
    static const int MODEL_DIMENSION = TAU_SIZE + N_SIZE + K_SIZE;
    static constexpr double MIN_K = 0.1;
    static constexpr double MAX_K = 1;
    static constexpr double MIN_N = 1;
    static constexpr double MAX_N = 25;
    static constexpr double MIN_TAU = 0.1;
    static constexpr double MAX_TAU = 5;
    static int modelFunction(double t, double *y, double *ydot, void *context) {
        appContext *ctx = (appContext *)context;
        double* individual = ctx->individual;
        double* maxValues = ((GRNSeries*)ctx->series)->getMaxValues();
        double *tau = &individual[0];
        double *k = &individual[10];
        double *n = &individual[25];

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
        double nAJ = (int)individual[25];
        double nBE = (int)individual[26];
        double nCB = (int)individual[27];
        double nCF = (int)individual[28];
        double nCA = (int)individual[29];
        double nDF = (int)individual[30];
        double nEJ = (int)individual[31];
        double nFA = (int)individual[32];
        double nGB = (int)individual[33];
        double nGF = (int)individual[34];
        double nGA = (int)individual[35];
        double nHF = (int)individual[36];
        double nIG = (int)individual[37];
        double nIH = (int)individual[38];
        double nJI = (int)individual[39];

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


        return 0;
    }

};

#endif //ES_GRN10MODEL_H
