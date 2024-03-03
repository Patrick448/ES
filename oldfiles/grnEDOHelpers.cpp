
int GRNEDOHelpers::trpEcoliModel(double t, double *y, double *ydot, void *context) {
    appContext *ctx = (appContext *) context;
    ProblemDescription *desc = ctx->description;
    double *individual = ctx->individual;
    double *maxValues = ((GRNSeries *) ctx->series)->getMaxValues();
    double *tau = &individual[0];
    double *k = &individual[desc->TAU_SIZE];
    double *n = &individual[desc->TAU_SIZE + desc->K_SIZE];

    ydot[0] = ((((1 - ((pow(y[0] / maximo_A, nAA)) / (pow(y[0] / maximo_A, nAA) + pow(kAA, nAA))))
                 *
                 (1 - ((pow(y[3] / maximo_D, nAD)) / (pow(y[3] / maximo_D, nAD) + pow(kAD, nAD))))
                 *
                 (1 - ((pow(y[4] / maximo_E, nAE)) / (pow(y[4] / maximo_E, nAE) + pow(kAE, nAE)))))
                +
                ((1 - ((pow(y[0] / maximo_A, nAA)) / (pow(y[0] / maximo_A, nAA) + pow(kAA, nAA))))
                 *
                 (1 - ((pow(y[2] / maximo_C, nAC)) / (pow(y[2] / maximo_C, nAC) + pow(kAC, nAC))))
                 *
                 (((pow(y[4] / maximo_E, nAE)) / (pow(y[4] / maximo_E, nAE) + pow(kAE, nAE)))))
                +
                ((((pow(y[0] / maximo_A, nAA)) / (pow(y[0] / maximo_A, nAA) + pow(kAA, nAA))))
                 *
                 (((pow(y[3] / maximo_D, nAD)) / (pow(y[3] / maximo_D, nAD) + pow(kAD, nAD))))
                 *
                 (1 - ((pow(y[4] / maximo_E, nAE)) / (pow(y[4] / maximo_E, nAE) + pow(kAE, nAE)))))) -
               (y[0] / maximo_A)) / tauA

//tauB, kBA, kBC, kBD, kBE, nBA, nBC, nBD, nBE

    ydot[1] = ((((1 - ((pow(y[0] / maximo_A, nBA)) / (pow(y[0] / maximo_A, nBA) + pow(kBA, nBA))))
                 *
                 (1 - ((pow(y[3] / maximo_D, nBD)) / (pow(y[3] / maximo_D, nBD) + pow(kBD, nBD))))
                 *
                 (1 - ((pow(y[4] / maximo_E, nBE)) / (pow(y[4] / maximo_E, nBE) + pow(kBE, nBE)))))
                +
                ((1 - ((pow(y[0] / maximo_A, nBA)) / (pow(y[0] / maximo_A, nBA) + pow(kBA, nBA))))
                 *
                 (1 - ((pow(y[2] / maximo_C, nBC)) / (pow(y[2] / maximo_C, nBC) + pow(kBC, nBC))))
                 *
                 (((pow(y[4] / maximo_E, nBE)) / (pow(y[4] / maximo_E, nBE) + pow(kBE, nBE)))))
                +
                ((((pow(y[0] / maximo_A, nBA)) / (pow(y[0] / maximo_A, nBA) + pow(kBA, nBA))))
                 *
                 (((pow(y[3] / maximo_D, nBD)) / (pow(y[3] / maximo_D, nBD) + pow(kBD, nBD))))
                 *
                 (1 - ((pow(y[4] / maximo_E, nBE)) / (pow(y[4] / maximo_E, nBE) + pow(kBE, nBE)))))) -
               (y[0] / maximo_B)) / tauB


//tauC, kCD, kCE, nCD, nCE

    ydot[2] = (((((pow(y[3] / maximo_D, nCD)) / (pow(y[3] / maximo_D, nCD) + pow(kCD, nCD)))) *
                (((pow(y[4] / maximo_E, nCE)) / (pow(y[4] / maximo_E, nCE) + pow(kCE, nCE))))) - (y[2] / maximo_C)) /
              tauC


//tauD, kDC, nDC

    ydot[3] = ((((pow(y[2] / maximo_C, nDC)) / (pow(y[2] / maximo_C, nDC) + pow(kDC, nDC)))) - (y[3] / maximo_D)) / tauD

//tauE, kAE, kBE, kCE, kEE, nAE, nBE, nCE, nEE

    ydot[4] = (((((pow(y[0] / maximo_A, nEA)) / (pow(y[0] / maximo_A, nEA) + pow(kEA, nEA)))
                 *
                 ((pow(y[1] / maximo_B, nEB)) / (pow(y[1] / maximo_B, nEB) + pow(kEB, nEB)))
                 *
                 ((pow(y[2] / maximo_C, nEC)) / (pow(y[2] / maximo_C, nEC) + pow(kEC, nEC))))
                +
                ((pow(y[4] / maximo_E, nEE)) / (pow(y[4] / maximo_E, nEE) + pow(kEE, nEE)))) - (y[4] / maximo_E)) / tauE
}