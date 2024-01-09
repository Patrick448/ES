//
// Created by patrick on 05/01/24.
//


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

void GRNEDOHelpers::initializeGRNContext(appContext* ctx, int granularity, int numVariables, int numTau, int numN, int numK, int setStart, int setEnd, double** vectors, double * maxValues)
{
    ctx->IND_SIZE = numVariables;      // Tamanho do indivíduo (quantidade de coeficientes)
    ctx->MIN_K = 0.01; //0.1        // Menor valor que K pode assumir
    ctx->MAX_K = 1;          // Maior valor que K pode assumir
    ctx->MIN_N = 1;          // Menor valor que N pode assumir
    ctx->MAX_N = 30; //25        // Maior valor que N pode assumir
    ctx->MIN_TAU = 0.1;      // Menor valor que TAU pode assumir
    ctx->MAX_TAU = 6;//5       // Maior valor que TAU pode assumir
    ctx->MIN_STRATEGY = 0.1; // Menor valor que a estratégia pode assumir
    ctx->MAX_STRATEGY = 10;  // Maior valor que a estratégia pode assumir
    ctx->TAU_SIZE = numTau;
    ctx->N_SIZE = numN;
    ctx->K_SIZE = numK;
    ctx->granularity = granularity;
    ctx->nVariables = numVariables;
    ctx->totalSteps = granularity*49;
    ctx->setStart = setStart;
    ctx->setEnd = setEnd;
    ctx->nSteps = setEnd - setStart + 1;
    ctx->maxValues = maxValues;
    ctx->yout = new double[(ctx->totalSteps + 1) * ctx->nVariables];
    ctx->vectors = vectors;

    ctx->y_0 = new double[ctx->nVariables];
    ctx->expectedResult = &ctx->vectors[1];

    //todo: verificar se o y_0 está correto para
    // o caso de inicio fora do zero
    for (int i = 0; i < ctx->nVariables; i++)
    {
        ctx->y_0[i] = ctx->vectors[i + 1][ctx->setStart];
    }

    ctx->tspan[0] = ctx->vectors[0][ctx->setStart];
    ctx->tspan[1] = ctx->vectors[0][ctx->setEnd];

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
    ctx->MAX_TAU = 6;//5       // Maior valor que TAU pode assumir
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

}

void GRNEDOHelpers::initializeGRN10Context(appContext* ctx, int mode, int granularity)
{

    /*ctx->TRAINING_MODE = 0;
    ctx->TEST_MODE = 2;
    ctx->VALIDATION_MODE = 1;
    ctx->SINGLE_SET_MODE = 3;*/
    ctx->IND_SIZE = 40;      // Tamanho do indivíduo (quantidade de coeficientes)
    ctx->MIN_K = 0.01; //0.1        // Menor valor que K pode assumir
    ctx->MAX_K = 1;          // Maior valor que K pode assumir
    ctx->MIN_N = 1;          // Menor valor que N pode assumir
    ctx->MAX_N = 30; //25        // Maior valor que N pode assumir
    ctx->MIN_TAU = 0.1;      // Menor valor que TAU pode assumir
    ctx->MAX_TAU = 6;//5        // Maior valor que TAU pode assumir
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

void GRNEDOHelpers::clearContext2Test(appContext* ctx)
{
    delete[] ctx->yout;
    delete[] ctx->y_0;

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
        }
        tout = tout + dt;
    }

    delete[] rtol;
    delete[] atol;
    delete[] y;
    lsoda_free(&ctx);

    return 0;
}

double GRNEDOHelpers::difference(double *actual, double **expected, int numVariables, int start, int end, int numSteps)
{
    double difTotal = 0.0;
    int numElements = end - start + 1;
    int jump = numSteps / (numElements-1);

    for (int i = 0; i < numVariables; i++)
    {
        for (int j = start; j <= end; j++)
        {
            int index = jump*j* numVariables + i;
            difTotal += fabs(actual[index] - expected[i][j]);
        }
    }
    if (isnan(difTotal))
    {
        return DBL_MAX;
    }
    return difTotal;
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


void func(){
    appContext *ctx = (appContext *)context;
    ProblemDescription* desc = ctx->description;
    double* individual = ctx->individual;
    double* maxValues = ((GRNSeries*)ctx->series)->getMaxValues();
    double *tau = &individual[0];//tau[0], tau[1], tau[2], tau[3], tau[4], tau[5], tau[6], tau[7], tau[8], tau[9]
    double *k = &individual[desc->TAU_SIZE];              //k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[9],
    // k[10],k[11],k[12],k[13],k[14],k[15],k[16],k[17],k[18],k[19],k[20],k[21],k[22],k[23],k[24]
    double *n = &individual[desc->TAU_SIZE+desc->K_SIZE];//(int)n[0],(int)n[1],(int)n[2],(int)n[3],(int)n[4],(int)n[5],(int)n[6],(int)n[7],
    // (int)n[8],(int)n[9],(int)n[10],(int)n[11],(int)n[12],(int)n[13],(int)n[14],(int)n[15],(int)n[16],
    // (int)n[17],(int)n[18],(int)n[19],(int)n[20],(int)n[21],(int)n[22],(int)n[23],(int)n[24]

//tau[0], k[0], (int)n[0], k[1], (int)n[1]

    ydot[0] = ((((pow(y[7]/maxValues[7], (int)n[0])) / (pow(y[7]/maxValues[7], (int)n[0]) + pow(k[0], (int)n[0]))) + ((pow(y[1]/maxValues[1], (int)n[1])) / (pow(y[1]/maxValues[1], (int)n[1]) + pow(k[1], (int)n[1])))) - (y[0]/maxValues[0])) / tau[0];

//tau[1], k[2], (int)n[2], k[3], (int)n[3]

    ydot[1] = ((((1-((pow(y[1]/maxValues[1], (int)n[2])) / (pow(y[1]/maxValues[1], (int)n[2]) + pow(k[2], (int)n[2])))) * ((pow(y[5]/maxValues[5], (int)n[3])) / (pow(y[5]/maxValues[5], (int)n[3]) + pow(k[3], (int)n[3])))) + ((1-((pow(y[5]/maxValues[5], (int)n[3])) / (pow(y[5]/maxValues[5], (int)n[3]) + pow(k[3], (int)n[3])))) * ((pow(y[1]/maxValues[1], (int)n[2])) / (pow(y[1]/maxValues[1], (int)n[2]) + pow(k[2], (int)n[2]))))) - (y[1]/maxValues[1])) / tau[1];

//tau[2], k[4], k[5], (int)n[4], (int)n[5]

    ydot[2] = ((((pow(y[0]/maxValues[0], (int)n[4])) / (pow(y[0]/maxValues[0], (int)n[4]) + pow(k[4], (int)n[4]))) + ((pow(y[2]/maxValues[2], (int)n[5])) / (pow(y[2]/maxValues[2], (int)n[5]) + pow(k[5], (int)n[5])))) - (y[2]/maxValues[2])) / tau[2];

//tau[3], k[6], k[7], k[8], k[9], (int)n[6], (int)n[7], (int)n[8], (int)n[9]

    ydot[3] = (((((pow(y[1]/maxValues[1], (int)n[6])) / (pow(y[1]/maxValues[1], (int)n[6]) + pow(k[6], (int)n[6])))
                 *
                 (1-((pow(y[4]/maxValues[4], (int)n[8])) / (pow(y[4]/maxValues[4], (int)n[8]) + pow(k[8], (int)n[8]))))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[9])) / (pow(y[7]/maxValues[7], (int)n[9]) + pow(k[9], (int)n[9]))))

                +

                (((pow(y[1]/maxValues[1], (int)n[6])) / (pow(y[1]/maxValues[1], (int)n[6]) + pow(k[6], (int)n[6])))
                 *
                 (1-((pow(y[3]/maxValues[3], (int)n[7])) / (pow(y[3]/maxValues[3], (int)n[7]) + pow(k[7], (int)n[7]))))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[9])) / (pow(y[7]/maxValues[7], (int)n[9]) + pow(k[9], (int)n[9]))))

                +

                (((pow(y[1]/maxValues[1], (int)n[6])) / (pow(y[1]/maxValues[1], (int)n[6]) + pow(k[6], (int)n[6])))
                 *
                 ((pow(y[3]/maxValues[3], (int)n[7])) / (pow(y[3]/maxValues[3], (int)n[7]) + pow(k[7], (int)n[7])))
                 *
                 ((pow(y[4]/maxValues[4], (int)n[8])) / (pow(y[4]/maxValues[4], (int)n[8]) + pow(k[8], (int)n[8])))
                 *
                 (1-((pow(y[7]/maxValues[7], (int)n[9])) / (pow(y[7]/maxValues[7], (int)n[9]) + pow(k[9], (int)n[9])))))) - (y[3]/maxValues[3])) / tau[3];

//tau[4], k[10], (int)n[10]

    ydot[4] = (((pow(y[7]/maxValues[7], (int)n[10])) / (pow(y[7]/maxValues[7], (int)n[10]) + pow(k[10], (int)n[10]))) - (y[4]/maxValues[4])) / tau[4];

//tau[5], k[11], k[12], k[13], (int)n[11], (int)n[12], (int)n[13]

    ydot[5] = ((((1-((pow(y[0]/maxValues[0], (int)n[11])) / (pow(y[0]/maxValues[0], (int)n[11]) + pow(k[11], (int)n[11]))))
                 *
                 (1-((pow(y[2]/maxValues[2], (int)n[12])) / (pow(y[2]/maxValues[2], (int)n[12]) + pow(k[12], (int)n[12]))))
                 *
                 ((pow(y[8]/maxValues[8], (int)n[13])) / (pow(y[8]/maxValues[8], (int)n[13]) + pow(k[13], (int)n[13]))))

                +

                (((pow(y[0]/maxValues[0], (int)n[11])) / (pow(y[0]/maxValues[0], (int)n[11]) + pow(k[11], (int)n[11])))
                 *
                 (1-((pow(y[8]/maxValues[8], (int)n[13])) / (pow(y[8]/maxValues[8], (int)n[13]) + pow(k[13], (int)n[13])))))

                +

                (((pow(y[2]/maxValues[2], (int)n[12])) / (pow(y[2]/maxValues[2], (int)n[12]) + pow(k[12], (int)n[12])))
                 *
                 (1-((pow(y[8]/maxValues[8], (int)n[13])) / (pow(y[8]/maxValues[8], (int)n[13]) + pow(k[13], (int)n[13])))))) - (y[5]/maxValues[5])) / tau[5];

//tau[6], k[14], k[15], k[16], (int)n[14], (int)n[15], (int)n[16]

    ydot[6] = ((((1-((pow(y[4]/maxValues[4], (int)n[15])) / (pow(y[4]/maxValues[4], (int)n[15]) + pow(k[15], (int)n[15]))))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[16])) / (pow(y[7]/maxValues[7], (int)n[16]) + pow(k[16], (int)n[16]))))

                +

                ((1-((pow(y[3]/maxValues[3], (int)n[14])) / (pow(y[3]/maxValues[3], (int)n[14]) + pow(k[14], (int)n[14]))))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[16])) / (pow(y[7]/maxValues[7], (int)n[16]) + pow(k[16], (int)n[16]))))

                +

                (((pow(y[3]/maxValues[3], (int)n[14])) / (pow(y[3]/maxValues[3], (int)n[14]) + pow(k[14], (int)n[14])))
                 *
                 ((pow(y[4]/maxValues[4], (int)n[15])) / (pow(y[4]/maxValues[4], (int)n[15]) + pow(k[15], (int)n[15])))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[16])) / (pow(y[7]/maxValues[7], (int)n[16]) + pow(k[16], (int)n[16]))))) - (y[6]/maxValues[6])) / tau[6];

//tau[7], k[17], (int)n[17], k[18], (int)n[18]

    ydot[7] = ((((pow(y[7]/maxValues[7], (int)n[17])) / (pow(y[7]/maxValues[7], (int)n[17]) + pow(k[17], (int)n[17]))) + ((pow(y[1]/maxValues[1], (int)n[18])) / (pow(y[1]/maxValues[1], (int)n[18]) + pow(k[18], (int)n[18])))) - (y[7]/maxValues[7])) / tau[7];

//tau[8], k[19], k[20], k[21], k[22], (int)n[19], (int)n[20], (int)n[21], (int)n[22]

    ydot[8] = (((((pow(y[1]/maxValues[1], (int)n[19])) / (pow(y[1]/maxValues[1], (int)n[19]) + pow(k[19], (int)n[19])))
                 *
                 (1-((pow(y[4]/maxValues[4], (int)n[21])) / (pow(y[4]/maxValues[4], (int)n[21]) + pow(k[21], (int)n[21]))))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[22])) / (pow(y[7]/maxValues[7], (int)n[22]) + pow(k[22], (int)n[22]))))

                +

                (((pow(y[1]/maxValues[1], (int)n[19])) / (pow(y[1]/maxValues[1], (int)n[19]) + pow(k[19], (int)n[19])))
                 *
                 (1-((pow(y[3]/maxValues[3], (int)n[20])) / (pow(y[3]/maxValues[3], (int)n[20]) + pow(k[20], (int)n[20]))))
                 *
                 ((pow(y[7]/maxValues[7], (int)n[22])) / (pow(y[7]/maxValues[7], (int)n[22]) + pow(k[22], (int)n[22]))))

                +

                (((pow(y[1]/maxValues[1], (int)n[19])) / (pow(y[1]/maxValues[1], (int)n[19]) + pow(k[19], (int)n[19])))
                 *
                 ((pow(y[3]/maxValues[3], (int)n[20])) / (pow(y[3]/maxValues[3], (int)n[20]) + pow(k[20], (int)n[20])))
                 *
                 ((pow(y[4]/maxValues[4], (int)n[21])) / (pow(y[4]/maxValues[4], (int)n[21]) + pow(k[21], (int)n[21])))
                 *
                 (1-((pow(y[7]/maxValues[7], (int)n[22])) / (pow(y[7]/maxValues[7], (int)n[22]) + pow(k[22], (int)n[22])))))) - (y[8]/maxValues[8])) / tau[8];

//tau[9], k[23], (int)n[23], k[24], (int)n[24]

    ydot[9] = ((((1-((pow(y[1]/maxValues[1], (int)n[23])) / (pow(y[1]/maxValues[1], (int)n[23]) + pow(k[23], (int)n[23])))) * ((pow(y[5]/maxValues[5], (int)n[24])) / (pow(y[5]/maxValues[5], (int)n[24]) + pow(k[24], (int)n[24])))) + ((1-((pow(y[5]/maxValues[5], (int)n[24])) / (pow(y[5]/maxValues[5], (int)n[24]) + pow(k[24], (int)n[24])))) * ((pow(y[1]/maxValues[1], (int)n[23])) / (pow(y[1]/maxValues[1], (int)n[23]) + pow(k[23], (int)n[23]))))) - (y[9]/maxValues[9])) / tau[9];
}