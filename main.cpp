//
// Created by patri on 27/01/2023.
//

#include <cstring>
#include "dependencies.h"

#include "lsoda.h"


using namespace std;
int IND_SIZE;  // Tamanho do indivíduo (quantidade de coeficientes)
double MIN_K;  // Menor valor que K pode assumir
double MAX_K;  // Maior valor que K pode assumir
double MIN_N; // Menor valor que N pode assumir
double MAX_N; // Maior valor que N pode assumir
double MIN_TAU; // Menor valor que TAU pode assumir
double MAX_TAU; // Maior valor que TAU pode assumir
double MIN_STRATEGY; // Menor valor que a estratégia pode assumir
double MAX_STRATEGY; // Maior valor que a estratégia pode assumir
int TAU_SIZE;
int N_SIZE;
int K_SIZE;
double *maxValues;

int nVariables;
int nSteps;
double tspan[2];
double *yout;
double *y_0;
double **vectors;
double **expectedResult;

void twoBody(double t, double y[], double max[], double tau[], double n[], double k[], double yp[]) {

    yp[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[0]), n[0]) + pow(k[0], n[0]))) - (
            y[0] / max[0])) / tau[0];

    yp[1] = (((pow((y[0] / max[0]), (int) n[1])) / (pow((y[0] / max[0]), (int) n[1]) + pow(k[1], (int) n[1]))) -
             (y[1] / max[1])) / tau[1];

    yp[2] = (((pow((y[1] / max[1]), n[2])) / (pow((y[1] / max[1]), (int) n[2]) + pow(k[2], (int) n[2]))) -
             (y[2] / max[2])) / tau[2];

    yp[3] = (((pow((y[2] / max[2]), (int) n[3])) / (pow((y[2] / max[2]), (int) n[3]) + pow(k[3], (int) n[3]))) -
             (y[3] / max[3])) / tau[3];

    yp[4] = (((pow((y[3] / max[3]), (int) n[4])) / (pow((y[3] / max[3]), (int) n[4]) + pow(k[4], (int) n[4]))) -
             (y[4] / max[4])) / tau[4];

}

void printVector(double *vec, int size) {
    for (int i = 0; i < size; i++) {
        cout << vec[i] << " ";
    }
    cout << "\n";
}

void twoBody5Var(double t, double y[], double *dim, double yp[]) {

    double *tau = &dim[0];
    double *k = &dim[TAU_SIZE];
    double *n = &dim[TAU_SIZE + N_SIZE];

    yp[0] = ((1 - (pow((y[4] / maxValues[4]), (int) n[0])) /
                  (pow((y[4] / maxValues[4]), (int) n[0]) + pow(k[0], (int) n[0]))) - (y[0] / maxValues[0])) / tau[0];

    yp[1] = (((pow((y[0] / maxValues[0]), (int) n[1])) /
              (pow((y[0] / maxValues[0]), (int) n[1]) + pow(k[1], (int) n[1]))) - (y[1] / maxValues[1])) / tau[1];

    yp[2] = (((pow((y[1] / maxValues[1]), (int) n[2])) /
              (pow((y[1] / maxValues[1]), (int) n[2]) + pow(k[2], (int) n[2]))) - (y[2] / maxValues[2])) / tau[2];

    yp[3] = (((pow((y[2] / maxValues[2]), (int) n[3])) /
              (pow((y[2] / maxValues[2]), (int) n[3]) + pow(k[3], (int) n[3]))) - (y[3] / maxValues[3])) / tau[3];

    yp[4] = (((pow((y[3] / maxValues[3]), (int) n[4])) /
              (pow((y[3] / maxValues[3]), (int) n[4]) + pow(k[4], (int) n[4]))) - (y[4] / maxValues[4])) / tau[4];


    return;
}

/*void twoBodyLSODA(double t, double *y, double *ydot, void *dim) {

    double* tau = &dim[0];
    double* k = &dim[TAU_SIZE];
    double* n = &dim[TAU_SIZE+N_SIZE];

    ydot[0] = ((1 - (pow((y[4] / maxValues[4]), (int)n[0])) / (pow((y[4] / maxValues[4]),  (int)n[0]) + pow(k[0],  (int)n[0]))) - (y[0] / maxValues[0])) / tau[0];

    ydot[1] = (((pow((y[0] / maxValues[0]), (int)n[1])) / (pow((y[0] / maxValues[0]),  (int)n[1]) + pow(k[1],  (int)n[1]))) - (y[1] / maxValues[1])) / tau[1];

    ydot[2] = (((pow((y[1] / maxValues[1]),  (int)n[2])) / (pow((y[1] / maxValues[1]),  (int)n[2]) + pow(k[2],  (int)n[2]))) - (y[2] / maxValues[2])) / tau[2];

    ydot[3] = (((pow((y[2] / maxValues[2]),  (int)n[3])) / (pow((y[2] / maxValues[2]),  (int)n[3]) + pow(k[3],  (int)n[3]))) - (y[3] / maxValues[3])) / tau[3];

    ydot[4] = (((pow((y[3] / maxValues[3]),  (int)n[4])) / (pow((y[3] / maxValues[3]),  (int)n[4]) + pow(k[4],  (int)n[4]))) - (y[4] / maxValues[4])) / tau[4];


    return;
}*/

int twoBodyFixedLSODA(double t, double *y, double *ydot, double *data) {
    //todo: observar que os n devem ser avaliados como inteiros
    double max[] = {2.96, 1.8768, 1.0653, 1.0101, 1.4608};
    double tau[] = {1.25, 4, 1.02, 1.57, 3.43};
    double n[] = {13, 4, 3, 4, 16};
    double k[] = {0.72, 0.50, 0.45, 0.51, 0.52};

    ydot[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[4]), n[0]) + pow(k[0], n[0]))) - (
            y[0] / max[0])) / tau[0];

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

int twoBody5VarLSODA(double t, double *y, double *ydot, void *_data) {

    double *data = (double *) _data;
    double *tau = &data[0];
    double *k = &data[TAU_SIZE];
    double *n = &data[TAU_SIZE + N_SIZE];
    //double max[] = {2.96, 1.8768, 1.0653, 1.0101, 1.4608};

    //maxValues = max;

    ydot[0] = ((1 - (pow((y[4] / maxValues[4]), (int) n[0])) /
                  (pow((y[4] / maxValues[4]), (int) n[0]) + pow(k[0], (int) n[0]))) - (y[0] / maxValues[0])) / tau[0];

    ydot[1] = (((pow((y[0] / maxValues[0]), (int) n[1])) /
              (pow((y[0] / maxValues[0]), (int) n[1]) + pow(k[1], (int) n[1]))) - (y[1] / maxValues[1])) / tau[1];

    ydot[2] = (((pow((y[1] / maxValues[1]), (int) n[2])) /
              (pow((y[1] / maxValues[1]), (int) n[2]) + pow(k[2], (int) n[2]))) - (y[2] / maxValues[2])) / tau[2];

    ydot[3] = (((pow((y[2] / maxValues[2]), (int) n[3])) /
              (pow((y[2] / maxValues[2]), (int) n[3]) + pow(k[3], (int) n[3]))) - (y[3] / maxValues[3])) / tau[3];

    ydot[4] = (((pow((y[3] / maxValues[3]), (int) n[4])) /
              (pow((y[3] / maxValues[3]), (int) n[4]) + pow(k[4], (int) n[4]))) - (y[4] / maxValues[4])) / tau[4];


    return 0;
}

int twoBody10VarLSODA(double t, double *y, double *ydot, void *_data) {
    double *data = (double *) _data;
    double *tau = &data[0];
    double *k = &data[TAU_SIZE];
    double *n = &data[TAU_SIZE + N_SIZE];
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
    double nAJ = (int) data[25];
    double nBE = (int) data[26];
    double nCB = (int) data[27];
    double nCF = (int) data[28];
    double nCA = (int) data[29];
    double nDF = (int) data[30];
    double nEJ = (int) data[31];
    double nFA = (int) data[32];
    double nGB = (int) data[33];
    double nGF = (int) data[34];
    double nGA = (int) data[35];
    double nHF = (int) data[36];
    double nIG = (int) data[37];
    double nIH = (int) data[38];
    double nJI = (int) data[39];

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
             (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) - (y[2] / maximo_C)) / tauC;

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
             (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) - (y[6] / maximo_G)) / tauG;

    ydot[7] = (pow(y[5] / maximo_F, nHF) / (pow(y[5] / maximo_F, nHF) + pow(kHF, nHF)) - (y[7] / maximo_H)) / tauH;

    ydot[8] = (pow(y[6] / maximo_G, nIG) / (pow(y[6] / maximo_G, nIG) + pow(kIG, nIG)) * pow(y[7] / maximo_H, nIH) /
             (pow(y[7] / maximo_H, nIH) + pow(kIH, nIH)) - (y[8] / maximo_I)) / tauI;

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


void twoBody10Var(double t, double y[], double *dim, double yp[]) {

    double *tau = &dim[0];
    double *k = &dim[TAU_SIZE];
    double *n = &dim[TAU_SIZE + N_SIZE];
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
    double tauA = dim[0];
    double tauB = dim[1];
    double tauC = dim[2];
    double tauD = dim[3];
    double tauE = dim[4];
    double tauF = dim[5];
    double tauG = dim[6];
    double tauH = dim[7];
    double tauI = dim[8];
    double tauJ = dim[9];
    double kAJ = dim[10];
    double kBE = dim[11];
    double kCB = dim[12];
    double kCF = dim[13];
    double kCA = dim[14];
    double kDF = dim[15];
    double kEJ = dim[16];
    double kFA = dim[17];
    double kGB = dim[18];
    double kGF = dim[19];
    double kGA = dim[20];
    double kHF = dim[21];
    double kIG = dim[22];
    double kIH = dim[23];
    double kJI = dim[24];
    double nAJ = (int) dim[25];
    double nBE = (int) dim[26];
    double nCB = (int) dim[27];
    double nCF = (int) dim[28];
    double nCA = (int) dim[29];
    double nDF = (int) dim[30];
    double nEJ = (int) dim[31];
    double nFA = (int) dim[32];
    double nGB = (int) dim[33];
    double nGF = (int) dim[34];
    double nGA = (int) dim[35];
    double nHF = (int) dim[36];
    double nIG = (int) dim[37];
    double nIH = (int) dim[38];
    double nJI = (int) dim[39];

    yp[0] = ((1 - pow(y[9] / maximo_J, nAJ) / (pow(y[9] / maximo_J, nAJ) + pow(kAJ, nAJ))) - (y[0] / maximo_A)) / tauA;

    yp[1] = (pow(y[4] / maximo_E, nBE) / (pow(y[4] / maximo_E, nBE) + pow(kBE, nBE)) - (y[1] / maximo_B)) / tauB;

    yp[2] = (pow(y[1] / maximo_B, nCB) / (pow(y[1] / maximo_B, nCB) + pow(kCB, nCB)) *
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
             (pow(y[0] / maximo_A, nCA) + pow(kCA, nCA)) - (y[2] / maximo_C)) / tauC;

    yp[3] = (pow(y[5] / maximo_F, nDF) / (pow(y[4] / maximo_E, nDF) + pow(kDF, nDF)) - (y[3] / maximo_D)) / tauD;

    yp[4] = (1 - pow(y[9] / maximo_J, nEJ) / (pow(y[9] / maximo_J, nEJ) + pow(kEJ, nEJ)) - (y[4] / maximo_E)) / tauE;

    yp[5] = (pow(y[0] / maximo_A, nFA) / (pow(y[0] / maximo_A, nFA) + pow(kFA, nFA)) - (y[5] / maximo_F)) / tauF;

    yp[6] = (pow(y[1] / maximo_B, nGB) / (pow(y[1] / maximo_B, nGB) + pow(kGB, nGB)) *
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
             (pow(y[0] / maximo_A, nGA) + pow(kGA, nGA)) - (y[6] / maximo_G)) / tauG;

    yp[7] = (pow(y[5] / maximo_F, nHF) / (pow(y[5] / maximo_F, nHF) + pow(kHF, nHF)) - (y[7] / maximo_H)) / tauH;

    yp[8] = (pow(y[6] / maximo_G, nIG) / (pow(y[6] / maximo_G, nIG) + pow(kIG, nIG)) * pow(y[7] / maximo_H, nIH) /
             (pow(y[7] / maximo_H, nIH) + pow(kIH, nIH)) - (y[8] / maximo_I)) / tauI;

    yp[9] = (pow(y[8] / maximo_I, nJI) / (pow(y[8] / maximo_I, nJI) + pow(kJI, nJI)) - (y[9] / maximo_J)) / tauJ;


    yp[3] = (pow(y[5] / maximo_F, nDF) / (pow(y[4] / maximo_E, nDF) + pow(kDF, nDF)) - (y[3] / maximo_D)) / tauD;

//    cout << "#########################" << endl;
//    cout << pow(y[5] / maximo_F, nDF) << endl;
//    cout << pow(y[4] / maximo_E, nDF) << endl;
//    cout << pow(kDF, nDF) << endl;
//    cout << (y[3] / maximo_D) << endl;
//    cout << tauD << endl;
//    cout << y[3] << endl;
//    cout << maximo_D << endl;
//    cout << "#########################" << endl;


}

void twoBodyFixed(double t, double y[], double *dim, double yp[]) {
    //todo: observar que os n devem ser avaliados como inteiros
    double max[] = {2.96, 1.8768, 1.0653, 1.0101, 1.4608};
    double tau[] = {1.25, 4, 1.02, 1.57, 3.43};
    double n[] = {13, 4, 3, 4, 16};
    double k[] = {0.72, 0.50, 0.45, 0.51, 0.52};

    yp[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[4]), n[0]) + pow(k[0], n[0]))) - (
            y[0] / max[0])) / tau[0];

    yp[1] = (((pow((y[0] / max[0]), n[1])) / (pow((y[0] / max[0]), n[1]) + pow(k[1], n[1]))) - (y[1] / max[1])) /
            tau[1];

    yp[2] = (((pow((y[1] / max[1]), n[2])) / (pow((y[1] / max[1]), n[2]) + pow(k[2], n[2]))) - (y[2] / max[2])) /
            tau[2];

    yp[3] = (((pow((y[2] / max[2]), n[3])) / (pow((y[2] / max[2]), n[3]) + pow(k[3], n[3]))) - (y[3] / max[3])) /
            tau[3];

    yp[4] = (((pow((y[3] / max[3]), n[4])) / (pow((y[3] / max[3]), n[4]) + pow(k[4], n[4]))) - (y[4] / max[4])) /
            tau[4];

}

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


void phase_plot(int n, int m, double t[], double y[])
/*
  Purpose:

    predator_phase_plot makes a phase plot of the results.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2020

  Author:

    John Burkardt

  Input:

    int n: the number of steps to take.

    int m: the number of variables.

    double t[n+1], y[(n+1)*m]: the times and solution values.
*/
{
    char command_filename[80];
    FILE *command;
    char data_filename[80];
    FILE *data;
    char header[] = "predator";
    int j;

    printf("\n");
    printf("phase_plot:\n");
    printf("  Write command and data files that can be used\n");
    printf("  by gnuplot for a predator-prey phase plot.\n");
/*
  Create the data file.
*/
    strcpy(data_filename, header);
    strcat(data_filename, "_data.txt");

    data = fopen(data_filename, "wt");

    for (j = 0; j <= n; j++) {
        fprintf(data, "%g", t[j]);

        for (int i = 0; i < m; i++) {
            fprintf(data, " %g", y[i + j * m]);
        }

        fprintf(data, "\n");

        //fprintf ( data, "  %g  %g  %g\n", t[j], y[0+j*m], y[1+j*m] );
    }

    fclose(data);

    printf("\n");
    printf("  phase_plot: data stored in \"%s\".\n", data_filename);
/*
  Create the command file.
*/
    strcpy(command_filename, header);
    strcat(command_filename, "_commands.txt");

    command = fopen(command_filename, "wt");

    fprintf(command, "# %s\n", command_filename);
    fprintf(command, "#\n");
    fprintf(command, "# Usage:\n");
    fprintf(command, "#  gnuplot < %s\n", command_filename);
    fprintf(command, "#\n");
    fprintf(command, "set term png\n");
    fprintf(command, "set output '%s.png'\n", header);
    fprintf(command, "set xlabel '<-- Prey -->'\n");
    fprintf(command, "set ylabel '<-- Predator -->'\n");
    fprintf(command, "set title 'Predator-prey solution by rk4'\n");
    fprintf(command, "set grid\n");
    fprintf(command, "set style data lines\n");
    fprintf(command, "plot '%s' using 2:3 with lines\n", data_filename);
    fprintf(command, "quit\n");

    fclose(command);

    printf("  phase_plot: plot commands stored in \"%s\".\n",
           command_filename);

    return;
}


double difference(double *actual, double **expected, int numVariables, int numElements) {
    //todo: testar função e remover essa linha
    //double *dif = new double [numVariables];
    double difTotal = 0.0;

    for (int i = 0; i < numVariables; i++) {
        // dif[i] = 0;
        for (int j = 0; j < numElements; j++) {
            //  dif[i] += fabs(actual[i][j] - expected[i][j]);
            //if(i==3)continue;
            difTotal += fabs(actual[j * numVariables + i] - expected[i][j]);
        }
    }

    if (isnan(difTotal)) {
        return DBL_MAX;
    }

    return difTotal;
}

double grn5Evaluation(double *dim) {
    rk4(twoBody5Var, tspan, y_0, nSteps, nVariables, vectors[0], dim, yout);
    return difference(yout, expectedResult, nVariables, 50);
}

double grn10Evaluation(double *dim) {
    rk4(twoBody10Var, tspan, y_0, nSteps, nVariables, vectors[0], dim, yout);
    return difference(yout, expectedResult, nVariables, 50);
}

void grn_test() {
    //int m=5;
    //int n = 49;
    double *t;
    //double tspan[2];
    // double *y;
    double *y0;
    //double **vectors = new double*[m+1];
    yout = new double[(nSteps + 1) * nVariables];
    vectors = new double *[nVariables + 1];

    for (int i = 0; i < nVariables + 1; i++) {
        vectors[i] = new double[50];
    }

    y0 = new double[nVariables];
    readFileToVectors("../GRN5.txt", nVariables + 1, vectors);

    for (int i = 0; i < nVariables; i++) {
        y0[i] = vectors[i + 1][0];
    }


    // t = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    //yout = (double *) malloc((nSteps + 1) * nVariables * sizeof(double));
    //y0 = ( double * ) malloc ( m * sizeof ( double ) );

    printf("\n");
    printf("twoBodyFixed\n");
    printf("  Use rk4() to solve the twoBody5Var ODE.\n");

    tspan[0] = 0.0;
    tspan[1] = 72.0;
    //y0[0] = 5000.0;
    //y0[1] = 100.0;
    double *coefficients;

    rk4(twoBodyFixed, tspan, y0, nSteps, nVariables, vectors[0], coefficients, yout);

    phase_plot(nSteps, nVariables, vectors[0], yout);
/*
  Free memory.
*/
    //free ( t );
    //free ( y );
    delete[] yout;
    delete[] y0;
    for (int i = 0; i < nVariables; i++) {
        delete[] vectors[i];
    }

    return;
}

double testFunc(Individual *ind) {
    double eval = 0;
    double x = ind->getDimension(0);
    double y = ind->getDimension(1);
    eval = pow(x + 2 * y - 7, 2) + pow(2 * x + y - 5, 2);

    return eval;
}

double getMaxValue(double *values, int numElements) {
    double maxValue = 0;
    for (int i = 0; i < numElements; i++) {
        if (values[i] > maxValue) {
            maxValue = values[i];
        }
    }

    return maxValue;
}

void getMaxValues(double **data, double *outMaxValues, int numVariables, int numElements) {

    for (int i = 1; i < numVariables + 1; i++) {
        outMaxValues[i - 1] = getMaxValue(data[i], numElements);
    }
}

void initializeGRN5() {
    IND_SIZE = 15;  // Tamanho do indivíduo (quantidade de coeficientes)
    MIN_K = 0.1;  // Menor valor que K pode assumir
    MAX_K = 1;  // Maior valor que K pode assumir
    MIN_N = 1; // Menor valor que N pode assumir
    MAX_N = 25; // Maior valor que N pode assumir
    MIN_TAU = 0.1; // Menor valor que TAU pode assumir
    MAX_TAU = 5; // Maior valor que TAU pode assumir
    MIN_STRATEGY = 0.1; // Menor valor que a estratégia pode assumir
    MAX_STRATEGY = 10; // Maior valor que a estratégia pode assumir
    TAU_SIZE = 5;
    N_SIZE = 5;
    K_SIZE = 5;
    nVariables = 5;
    nSteps = 49;

    maxValues = new double[nVariables];

    yout = new double[(nSteps + 1) * nVariables];

    vectors = new double *[nVariables + 1];
    for (int i = 0; i < nVariables + 1; i++) {
        vectors[i] = new double[50];
    }

    readFileToVectors("GRN5.txt", nVariables + 1, vectors);
    getMaxValues(vectors, maxValues, nVariables, nSteps + 1);

    y_0 = new double[nVariables];
    expectedResult = &vectors[1];
    for (int i = 0; i < nVariables; i++) {
        y_0[i] = vectors[i + 1][0];
    }

    tspan[0] = 0.0;
    tspan[1] = 72.0;

}

void initializeGRN10() {
    IND_SIZE = 40;  // Tamanho do indivíduo (quantidade de coeficientes)
    MIN_K = 0.1;  // Menor valor que K pode assumir
    MAX_K = 1;  // Maior valor que K pode assumir
    MIN_N = 1; // Menor valor que N pode assumir
    MAX_N = 25; // Maior valor que N pode assumir
    MIN_TAU = 0.1; // Menor valor que TAU pode assumir
    MAX_TAU = 5; // Maior valor que TAU pode assumir
    MIN_STRATEGY = 0.1; // Menor valor que a estratégia pode assumir
    MAX_STRATEGY = 10; // Maior valor que a estratégia pode assumir
    TAU_SIZE = 10;
    N_SIZE = 15;
    K_SIZE = 15;
    nVariables = 10;
    nSteps = 49;
    maxValues = new double[nVariables];

    yout = new double[(nSteps + 1) * nVariables];

    vectors = new double *[nVariables + 1];
    for (int i = 0; i < nVariables + 1; i++) {
        vectors[i] = new double[50];
    }

    readFileToVectors("GRN10.txt", nVariables + 1, vectors);
    getMaxValues(vectors, maxValues, nVariables, nSteps + 1);

    y_0 = new double[nVariables];
    expectedResult = &vectors[1];
    for (int i = 0; i < nVariables; i++) {
        y_0[i] = vectors[i + 1][0];
    }

    tspan[0] = 0.0;
    tspan[1] = 72.0;

}


void clearGRNData() {
    //free ( t );
    delete[] yout;
    delete[] y_0;
    //todo: algum problema na desalocação, investigar
    delete[] maxValues;
    for (int i = 0; i < nVariables; i++) {
        delete[] vectors[i];
    }
}


double lsodaWrapper(int dydt(double t, double *y, double *ydot, void *data), double tspan[2],
                    double y0[], int n, int m, double tVec[], double *coefficients, double _yout[]) {

    //todo: tentar colocar essa alocação fora da função
    // esses vetores serão alocados toda vez, e essa função será chamada
    // a cada avaliação de indivíduo

    double* atol = new double[m];
    double* rtol  = new double[m];
    double t, tout, dt;

    double*  iterYOut = new double[m];
    int iout;

    for(int i=0; i<m; i++){
        iterYOut[i] = y0[i];
        rtol[i] = atol[i] = 1.49012e-8;
        _yout[i] = y0[i];
    }

    t = 0.0E0;
    dt = ( tspan[1] - tspan[0] ) / ( double ) ( n );
    tout = dt;

    struct lsoda_opt_t opt = {0};
    opt.ixpr = 0;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.itask = 1;

    struct lsoda_context_t ctx = {
            .function = dydt,
            .data = NULL,
            .neq = m,
            .state = 1,
    };
    ctx.data = coefficients;

    lsoda_prepare(&ctx, &opt);

    for (iout = 1; iout <= n; iout++) {
        lsoda(&ctx, iterYOut, &t, tout);
        for(int i=0; i<m; i++){
           _yout[iout*m + i] = iterYOut[i];
        }

        if (ctx.state <= 0) {
            printf("error istate = %d\n", ctx.state);
            exit(0);
        }
        tout = tout + dt ;
    }

    delete [] rtol;
    delete [] atol;
    delete [] iterYOut;
    lsoda_free(&ctx);

    return 0;
}

double grn5EvaluationLSODA(double *dim) {
    lsodaWrapper(twoBody5VarLSODA, tspan, y_0, nSteps, nVariables, vectors[0], dim, yout);
    return difference(yout, expectedResult, nVariables, 50);
}

double grn10EvaluationLSODA(double *dim) {
    lsodaWrapper(twoBody10VarLSODA, tspan, y_0, nSteps, nVariables, vectors[0], dim, yout);
    return difference(yout, expectedResult, nVariables, 50);
}

double grn5EvaluationLSODA() {
    double atol[5], rtol[5], t, tout, y[5];
    int neq = 5;
    int iout;
    initializeGRN5();

    //esse y será y_0
    y[0] = y_0[0];
    y[1] = y_0[1];
    y[2] = y_0[2];
    y[3] = y_0[3];
    y[4] = y_0[4];
    t = 0.0E0;
    tout = 72.0 / 49.0;
    struct lsoda_opt_t opt = {0};
    opt.ixpr = 0;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.itask = 1;

    rtol[0] = atol[0] = 1.49012e-8;
    rtol[1] = atol[1] = 1.49012e-8;
    rtol[2] = atol[2] = 1.49012e-8;
    rtol[3] = atol[3] = 1.49012e-8;
    rtol[4] = atol[4] = 1.49012e-8;



    struct lsoda_context_t ctx = {
            .function = twoBodyFixedLSODA,
            .data = NULL,
            .neq = neq,
            .state = 1,
    };


    lsoda_prepare(&ctx, &opt);
    double *dummyDBL;
    double dim[] = {1.25, 4, 1.02, 1.57, 3.43, 0.72, 0.5, 0.45, 0.51, 0.52, 13, 4, 3, 4, 16};

    for (iout = 1; iout <= 49; iout++) {
        ctx.data = dim;
        lsoda(&ctx, y, &t, tout);
        printf(" at t= %f y= %14.6e %14.6e %14.6e %14.6e %14.6e\n", t, y[0], y[1], y[2], y[3], y[4]);
        if (ctx.state <= 0) {
            printf("error istate = %d\n", ctx.state);
            exit(0);
        }
        tout = tout + (72.0 / 49.0);
    }
    lsoda_free(&ctx);
    clearGRNData();

    //rk4 (twoBody, tspan, y_0, nSteps, nVariables, vectors[0], dim, y);
    return 0;//difference(y, expectedResult, nVariables, 50);
}

double grn5EvaluatioTest() {
    initializeGRN5();
    double dim[] = {1.25, 4, 1.02, 1.57, 3.43, 0.72, 0.5, 0.45, 0.51, 0.52, 13, 4, 3, 4, 16};
    //double dim[] = {5, 5, 0.1, 2.57215, 0.1, 0.1, 1, 0.1, 1, 1, 20.4753, 14.8202, 23.1872, 9.29585, 8.17558};
    //double dim[] = {5, 5, 0.3, 2.57215, 0.5, 0.1, 1, 0.1, 1, 1, 20.4753, 14.8202, 23.1872, 9.29585, 8.17558};
    //double dim[] = {0.1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    //double dim[] = {0.1, 5, 0.3, 2.57215, 0.5, 0.1, 1, 0.1, 1, 1, 20.4753, 14.8202, 23.1872, 9.29585, 8.17558};


    rk4(twoBody5Var, tspan, y_0, nSteps, nVariables, vectors[0], dim, yout);
    for (int i = 0; i < nVariables; i++) {
        for (int j = 0; j < nSteps; j++) {
            cout << expectedResult[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n\n";
    for (int i = 0; i < nVariables; i++) {
        for (int j = 0; j < nSteps; j++) {
            cout << yout[j * nVariables + i] << " ";
        }
        cout << "\n";
    }


    cout << "Diference: " << difference(yout, expectedResult, nVariables, 50) << endl;
    clearGRNData();

    return 0;
}


double grn10EvaluatioTest() {
    initializeGRN10();
    double dim[] = {1.73, 2, 0.81, 0.11, 1.23, 1.78, 1.14, 1.04, 3.47, 3.21,
                    0.45, 0.56, 0.99, 0.77, 0.71, 0.66, 0.46, 0.48, 0.66, 0.99, 0.85, 0.61, 0.55, 0.46, 0.17,
                    20, 9, 24, 12, 2, 2, 6, 4, 7, 24, 2, 7, 21, 20, 3};


    rk4(twoBody10Var, tspan, y_0, nSteps, nVariables, vectors[0], dim, yout);

//    for(int i=0; i<50; i++){
//        y[i*10 + 3] = 0;
//   }

    for (int i = 0; i < nVariables; i++) {
        for (int j = 0; j < nSteps + 1; j++) {
            cout << expectedResult[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n\n";
    for (int i = 0; i < nVariables; i++) {
        for (int j = 0; j < nSteps + 1; j++) {
            cout << yout[j * nVariables + i] << " ";
        }
        cout << "\n";
    }


    cout << "Diference: " << difference(yout, expectedResult, nVariables, 50) << "\n\n";
    clearGRNData();

    return 0;
}

void runGRN5ESComparisonExperiment() {
    initializeGRN5();

    ESAlgorithm esAlgorithm = ESAlgorithm(IND_SIZE);
    esAlgorithm.setEvaluationFunction(grn5EvaluationLSODA);
    esAlgorithm.setSigmaBounds(MIN_STRATEGY, MAX_STRATEGY);

    int cont = 0;
    for (int i = 0; i < TAU_SIZE; i++) {
        esAlgorithm.setBounds(i, MIN_TAU, MAX_TAU, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < TAU_SIZE + K_SIZE; i++) {
        esAlgorithm.setBounds(i, MIN_K, MAX_K, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < TAU_SIZE + K_SIZE + N_SIZE; i++) {
        esAlgorithm.setBounds(i, MIN_N, MAX_N, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
        cont = i;
    }

    int numRuns = 1;

    /*todo: usar representação contígua para essa matriz
     * e, para calcular posição, usar
     * #define R(i,j) results[i*5 + j]
    */
    double** results = new double*[5];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];
    results[3] = new double[numRuns];
    results[4] = new double[numRuns];

    vector<string> bestInds(5);
    vector<double> bestIndsEval(5, DBL_MAX);

    for (int i = 0; i < numRuns; i++) {
        cout << "Run " << to_string(i) << "\n";

        cout<< "1" <<"\n";
        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, 200);
        results[0][i] = esAlgorithm.getPopulation().back()->getEvaluation();
        if (results[0][i] < bestIndsEval[0]) {
            bestIndsEval[0] = results[0][i];
            bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }

        cout<< "2" <<"\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 10, 10, 20);
        results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[1][i] < bestIndsEval[1]) {
            bestIndsEval[1] = results[1][i];
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout<< "3" <<"\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 10, 10, 20);
        results[2][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[2][i] < bestIndsEval[2]) {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout<< "4" <<"\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 10, 5, 20);
        results[3][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[3][i] < bestIndsEval[3]) {
            bestIndsEval[3] = results[3][i];
            bestInds[3] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout<< "5" <<"\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 10, 5, 20);
        results[4][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[4][i] < bestIndsEval[4]) {
            bestIndsEval[4] = results[4][i];
            bestInds[4] = esAlgorithm.getPopulation()[0]->toCSVString();
        }
    }

    string csvOutput = "1+1,10+20-i,10+20-ni,5+10-i,5+10-ni\n";
    for (int j = 0; j < numRuns; j++) {
        csvOutput += to_string(results[0][j]) + ","
                     + to_string(results[1][j]) + ","
                     + to_string(results[2][j]) + ","
                     + to_string(results[3][j]) + ","
                     + to_string(results[4][j]) + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2] + "\n" + bestInds[3] + "\n" + bestInds[4];

    outputToFile("../comparison-30runs-200000it_LSODA.csv", csvOutput, false);
    outputToFile("../best-individuals_LSODA.txt", bestIndividuals, false);


    delete [] results[0];
    delete [] results[1];
    delete [] results[2];
    delete [] results[3];
    delete [] results[4];
    delete [] results;
    clearGRNData();
}

void runGRN10ESComparisonExperiment() {
    initializeGRN10();

    ESAlgorithm esAlgorithm = ESAlgorithm(IND_SIZE);
    esAlgorithm.setEvaluationFunction(grn10EvaluationLSODA);
    esAlgorithm.setSigmaBounds(MIN_STRATEGY, MAX_STRATEGY);

    int cont = 0;
    for (int i = 0; i < TAU_SIZE; i++) {
        esAlgorithm.setBounds(i, MIN_TAU, MAX_TAU, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < TAU_SIZE + K_SIZE; i++) {
        esAlgorithm.setBounds(i, MIN_K, MAX_K, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < TAU_SIZE + K_SIZE + N_SIZE; i++) {
        esAlgorithm.setBounds(i, MIN_N, MAX_N, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
        cont = i;
    }

    int numRuns = 10;
    /*vector<vector<double>> results(5);
    
    results[0].resize(numRuns);
    results[1].resize(numRuns);
    results[2].resize(numRuns);
    results[3].resize(numRuns);
    results[4].resize(numRuns);
    
    */
    double** results = new double*[5];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];
    results[3] = new double[numRuns];
    results[4] = new double[numRuns];

    vector<string> bestInds(5);
    vector<double> bestIndsEval(5, DBL_MAX);

    for (int i = 0; i < numRuns; i++) {
        cout << "Run " << to_string(i) << endl;

        cout << "0" << "\n";
        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, 200);
        results[0][i] = esAlgorithm.getPopulation().back()->getEvaluation();
        if (results[0][i] < bestIndsEval[0]) {
            bestIndsEval[0] = results[0][i];
            bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }
        cout << "1" << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 10, 10, 20);
        results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[1][i] < bestIndsEval[1]) {
            bestIndsEval[1] = results[1][i];
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "2" << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 10, 10, 20);
        results[2][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[2][i] < bestIndsEval[2]) {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "3" << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 10, 10, 20);
        results[3][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[3][i] < bestIndsEval[3]) {
            bestIndsEval[3] = results[3][i];
            bestInds[3] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "4" << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 10, 10, 20);
        results[4][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[4][i] < bestIndsEval[4]) {
            bestIndsEval[4] = results[4][i];
            bestInds[4] = esAlgorithm.getPopulation()[0]->toCSVString();
        }
    }

    string csvOutput = "1+1,10+20-i,10+20-ni,5+10-i,5+10-ni\n";
    for (int j = 0; j < numRuns; j++) {
        csvOutput += to_string(results[0][j]) + ","
                     + to_string(results[1][j]) + ","
                     + to_string(results[2][j]) + ","
                     + to_string(results[3][j]) + ","
                     + to_string(results[4][j]) + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2] + "\n" + bestInds[3] + "\n" + bestInds[4];

    outputToFile("../lsoda-comparison-GRN10-30runs-200000it.csv", csvOutput, false);
    outputToFile("../lsoda-best-individuals-GRN10.txt", bestIndividuals, false);


    delete [] results[0];
    delete [] results[1];
    delete [] results[2];
    delete [] results[3];
    delete [] results[4];
    delete [] results;


    clearGRNData();
}

int fex(double t, double *y, double *ydot, void *data) {
    ydot[0] = 1.0E4 * y[1] * y[2] - .04E0 * y[0];
    ydot[2] = 3.0E7 * y[1] * y[1];
    ydot[1] = -1.0 * (ydot[0] + ydot[2]);
    return (0);
}


int main() {

    //testa a saida da função difference com coeficientes pré-definidos
    //resultado: ~108.0
    initializeGRN5();
    double dim5[] = {1.25, 4, 1.02, 1.57, 3.43, 0.72, 0.5, 0.45, 0.51, 0.52, 13, 4, 3, 4, 16};
    cout << grn5EvaluationLSODA(dim5) << "\n";
    cout << "\n\n\n";
    clearGRNData();

    //resultado: ~56.0
    initializeGRN10();
    double dim10[] = {1.73,2,0.81,0.11, 1.23, 1.78, 1.14, 1.04, 3.47, 3.21,
                      0.45, 0.56, 0.99, 0.77, 0.71, 0.66, 0.46, 0.48, 0.66, 0.99, 0.85, 0.61, 0.55, 0.46, 0.17,
                      20, 9, 24, 12, 2, 2, 6, 4, 7, 24, 2, 7, 21, 20, 3};
    cout << grn10EvaluationLSODA(dim10) << "\n";
    cout << "\n\n\n";
    clearGRNData();

    return 0;


    //runGRN10ESComparisonExperiment();
    //runGRN5ESComparisonExperiment();
    //return 0;


}



