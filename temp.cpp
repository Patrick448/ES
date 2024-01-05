//
// Created by patrick on 05/01/24.
//

int grn4NCYCModel(double t, double *y, double *ydot, void *context)
{
    appContext *ctx = (appContext *)context;
    ProblemDescription* desc = ctx->description;
    double* individual = ctx->individual;
    double* maxValues = ((GRNSeries*)ctx->series)->getMaxValues();
    double *tau = &individual[0];//tau[0], tau[1], tau[2], tau[3]
    double *k = &individual[desc->TAU_SIZE];             //k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[10]
                                                        //nAB,nAC,nBA,nBB,nBC,nCA,nCC,nDA,nDB, (int)n[9],nDD
    double *n = &individual[desc->TAU_SIZE+desc->K_SIZE];//(int)n[0],(int)n[1],(int)n[2],(int)n[3],(int)n[4],(int)n[5],(int)n[6],(int)n[7],(int)n[8],(int)n[10]


    ydot[0] = (((1 - ((pow(y[1]/maxValues[1], (int)n[0])) / (pow(y[1]/maxValues[1], (int)n[0]) + pow(k[0], (int)n[0]))))
                *
                (((pow(y[2]/maxValues[2], (int)n[1])) / (pow(y[2]/maxValues[2], (int)n[1]) + pow(k[1], (int)n[1])))))
               - (y[0]/maxValues[0])) / tau[0];


    ydot[1] = (((((pow(y[0]/maxValues[0], (int)n[2])) / (pow(y[0]/maxValues[0], (int)n[2]) + pow(k[2], (int)n[2]))))
                +
                (((pow(y[1]/maxValues[1], (int)n[3])) / (pow(y[1]/maxValues[1], (int)n[3]) + pow(k[3], (int)n[3]))))
                +
                (1-((pow(y[2]/maxValues[2], (int)n[4])) / (pow(y[2]/maxValues[2], (int)n[4]) + pow(k[4], (int)n[4]))))) - (y[1]/maxValues[1])) / tau[1];

    ydot[2] = (((((pow(y[0]/maxValues[0], (int)n[5])) / (pow(y[0]/maxValues[0], (int)n[5]) + pow(k[5], (int)n[5]))))
                +
                (((pow(y[2]/maxValues[2], (int)n[6])) / (pow(y[2]/maxValues[2], (int)n[6]) + pow(k[6], (int)n[6]))))) - (y[2]/maxValues[2])) / tau[2];


    ydot[3] =  ((((1-((pow(y[0]/maxValues[0], (int)n[7])) / (pow(y[0]/maxValues[0], (int)n[7]) + pow(k[7], (int)n[7]))))
                  *
                  (((pow(y[1]/maxValues[1], (int)n[8])) / (pow(y[1]/maxValues[1], (int)n[8]) + pow(k[8], (int)n[8]))))
                  *
                  (1-((pow(y[3]/maxValues[3], (int)n[10])) / (pow(y[3]/maxValues[3], (int)n[10]) + pow(k[10], (int)n[10])))))
                 +
                 (
                         (1-((pow(y[0]/maxValues[0], (int)n[7])) / (pow(y[0]/maxValues[0], (int)n[7]) + pow(k[7], (int)n[7]))))
                         *
                         (1-((pow(y[1]/maxValues[1], (int)n[8])) / (pow(y[1]/maxValues[1], (int)n[8]) + pow(k[8], (int)n[8]))))
                         *
                         (1-((pow(y[2]/maxValues[2], (int)n[9])) / (pow(y[2]/maxValues[2], (int)n[9]) + pow(k[9], (int)n[9]))))
                         *
                         (((pow(y[3]/maxValues[3], (int)n[10])) / (pow(y[3]/maxValues[3], (int)n[10]) + pow(k[10], (int)n[10]))))
                 )) - (y[3]/maxValues[3])) / tau[3];

    return 0;
}