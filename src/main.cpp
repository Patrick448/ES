//
// Created by patri on 27/01/2023.
//

#include <cstring>
#include <chrono>
#include "dependencies.h"
#include "appCtx.h"
#include "GRNEDOHelpers.h"
#include "algModes.h"
#include <pagmo/algorithms/cmaes.hpp>
#include <pagmo/population.hpp>
#include <pagmo/utils/constrained.hpp>
#include "GRNCoefProblem.h"
#include "GRNSeries.h"
#include "ProblemDescription.h"

#ifdef __cplusplus
extern "C"
{
#endif

#include "lsoda.h"

#ifdef __cplusplus
}
#endif

using namespace std;
using namespace GRNEDOHelpers;
using namespace algModes;
using namespace pagmo;
/*
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

int nVariables;
int nSteps;
double tspan[2];
double *yout;
double *y_0;
double **vectors;
double **expectedResult;
*/

/*struct appContext {
    int TRAINING_MODE = 0;
    int TEST_MODE = 2;
    int VALIDATION_MODE = 1;
    int mode=0;
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

    int nVariables;
    int nSteps;
    int dataSetSize;
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

};*/

/// helper for outputting text to file
void outputToFile(string path, string text, bool append)
{
    ofstream outputf;

    if (append)
    {
        outputf.open(path, std::ios_base::app);
    }
    else
    {
        outputf.open(path);
    }

    outputf << text;
    outputf.close();
}

/// Runs an experiment.
/// todo: parametrize the number of runs and the number of evaluations
/// todo: remove unnecessary experiments
/// @param grnMode: grn5 or grn10
/// @param evalMode: lsoda or rk4
/// @param expName: name of the experiment (names the folder where the results will be stored)
/// @param outputDir: path to the output directory where the the expName folder will be created
void runESComparisonExperiment(string grnMode, string evalMode, string expName, string outputDir)
{
    appContext ctx{};
    double (*func)(void*,void*);

    //one plus one
    int maxEvals = 200000;

    //firt part populational algorithms
    int numParents = 10;
    int numOffspring = 20;
    int numRuns = 30;

    if(grnMode == "grn5"){
        maxEvals = 200000;
        numParents = 10;
        numOffspring = 20;

        if(evalMode == "lsoda"){
            func = &grn5EvaluationLSODA;
            initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
        }
        else {
            func = &grn5EvaluationRK4;
            initializeGRN5Context(&ctx, SINGLE_SET_MODE, 10);
        }
    }
    else {
        maxEvals = 200000;
        numParents = 20;
        numOffspring = 40;

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODA;
            initializeGRN10Context(&ctx, SINGLE_SET_MODE, 1);
        }
        else {
            func = &grn10EvaluationRK4;
            initializeGRN10Context(&ctx, SINGLE_SET_MODE, 10);
        }
    }

    int maxGenerations = maxEvals / numOffspring;
    string experimentId = grnMode + "-" + evalMode + "-" + to_string(numRuns) + "runs-"+ to_string(maxEvals) + "evals";
    string experimentGroup = expName;

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(func);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);

    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }



    /*todo: usar representação contígua para essa matriz
     * e, para calcular posição, usar
     * #define R(i,j) results[i*5 + j]
     */

    double **results = new double *[5];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];
    results[3] = new double[numRuns];
    results[4] = new double[numRuns];

    vector<string> bestInds(5);
    vector<double> bestIndsEval(5, DBL_MAX);


    for (int i = 0; i < numRuns; i++)
    {
        auto beg = chrono::high_resolution_clock::now();

        cout << "Run " << to_string(i) << "\n";

        cout << "1"
             << "\n";
        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, maxEvals);
        results[0][i] = esAlgorithm.getPopulation().back()->getEvaluation();
        if (results[0][i] < bestIndsEval[0])
        {
            bestIndsEval[0] = results[0][i];
            bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }

        cout << "2"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[1][i] < bestIndsEval[1])
        {
            bestIndsEval[1] = results[1][i];
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "3"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        results[2][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[2][i] < bestIndsEval[2])
        {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "4"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        results[3][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[3][i] < bestIndsEval[3])
        {
            bestIndsEval[3] = results[3][i];
            bestInds[3] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "5"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        results[4][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[4][i] < bestIndsEval[4])
        {
            bestIndsEval[4] = results[4][i];
            bestInds[4] = esAlgorithm.getPopulation()[0]->toCSVString();
        }


        //temporização
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
        outputToFile(outputDir + experimentGroup + "/" + experimentId+" -time.csv", to_string(duration.count()) + ",", true);
        cout << "Elapsed Time: " << duration.count() << "\n";
    }

    //todo: modificar para refletir os parâmetros da execução
    string csvOutput = "1+1,10+20-i,10+20-ni,5+10-i,5+10-ni\n";

    for (int j = 0; j < numRuns; j++)
    {
        csvOutput += to_string(results[0][j]) + "," + to_string(results[1][j]) + "," + to_string(results[2][j]) + "," + to_string(results[3][j]) + "," + to_string(results[4][j]) + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2] + "\n" + bestInds[3] + "\n" + bestInds[4];

    outputToFile(outputDir + experimentGroup + "/"+experimentId+".csv", csvOutput, false);
    outputToFile(outputDir + experimentGroup + "/"+experimentId+"-best.csv", bestIndividuals, false);

    delete[] results[0];
    delete[] results[1];
    delete[] results[2];
    delete[] results[3];
    delete[] results[4];
    delete[] results;
    clearContext(&ctx);
}

void runCECComparisonExperiment(string grnMode, string evalMode)
{
    appContext ctx{};
    appContext validationContext{};
    double (*func)(void*,void*);

    //one plus one
    int maxEvals = 15 * 105*10000;

    //firt part populational algorithms
    int numParents = 15;
    int numOffspring = 105;
    int numRuns = 3;

    if(grnMode == "grn5"){
        maxEvals = 105*100;
        numParents = 15;
        numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn5EvaluationLSODA;
            initializeGRN5Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn5EvaluationRK4;
            initializeGRN5Context(&ctx, TRAINING_MODE, 20);
        }
    }
    else {
        maxEvals =  105*100;
        numParents = 15;
        numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODA;
            initializeGRN10Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn10EvaluationRK4;
            initializeGRN10Context(&ctx, TRAINING_MODE, 20);
        }
    }

    int maxGenerations = maxEvals / numOffspring;
    string experimentId = grnMode + "-" + evalMode + "-" + to_string(numRuns) + "runs-"+ to_string(maxEvals) + "evals";
    string experimentGroup = "exp8";

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(func);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);


    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }



    /*todo: usar representação contígua para essa matriz
     * e, para calcular posição, usar
     * #define R(i,j) results[i*5 + j]
     */

    double **results = new double *[3];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];

    vector<string> bestInds(3);
    vector<double> bestIndsEvalValidation(3, DBL_MAX);
    vector<double> bestIndsEvalTest(3, DBL_MAX);
    double bestIndEval;

    for (int i = 0; i < numRuns; i++)
    {
        auto beg = chrono::high_resolution_clock::now();

        cout << "Run " << to_string(i) << "\n";

        cout << "1"
             << "\n";

        setMode(&ctx, TRAINING_MODE);
        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, maxEvals);
        setMode(&ctx, VALIDATION_MODE);
        esAlgorithm.reevaluateAllNoCounter();
        bestIndEval = esAlgorithm.getPopulation()[0]->getEvaluation();
        setMode(&ctx, TEST_MODE);
        results[0][i] = esAlgorithm.getReevaluationByIndexNoCounter(0);

        if (bestIndEval < bestIndsEvalValidation[0])
        {
            bestIndsEvalValidation[0] = bestIndEval;
            bestInds[0] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "2"
             << "\n";
        setMode(&ctx, TRAINING_MODE);
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        setMode(&ctx, VALIDATION_MODE);
        esAlgorithm.reevaluateAllNoCounter();
        bestIndEval = esAlgorithm.getPopulation()[0]->getEvaluation();
        setMode(&ctx, TEST_MODE);
        results[1][i] = esAlgorithm.getReevaluationByIndexNoCounter(0);

        if (bestIndEval < bestIndsEvalValidation[1])
        {
            bestIndsEvalValidation[1] = bestIndEval;
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "3"
             << "\n";
        setMode(&ctx, TRAINING_MODE);
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        setMode(&ctx, VALIDATION_MODE);
        esAlgorithm.reevaluateAllNoCounter();
        bestIndEval = esAlgorithm.getPopulation()[0]->getEvaluation();
        setMode(&ctx, TEST_MODE);
        results[2][i] = esAlgorithm.getReevaluationByIndexNoCounter(0);

        if (bestIndEval < bestIndsEvalValidation[2])
        {
            bestIndsEvalValidation[2] = bestIndEval;
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        //temporização
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
        outputToFile("../results/" + experimentGroup + "/" + experimentId+" -time.csv", to_string(duration.count()) + ",", true);
        cout << "Elapsed Time: " << duration.count() << "\n";
    }

    //todo: modificar para refletir os parâmetros da execução
    string csvOutput = "1+1,15+105-i,15+105-ni\n";

    for (int j = 0; j < numRuns; j++)
    {
        csvOutput += to_string(results[0][j]) + "," + to_string(results[1][j]) + "," + to_string(results[2][j])  + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2];

    outputToFile("../results/" + experimentGroup + "/"+experimentId+".csv", csvOutput, false);
    outputToFile("../results/" + experimentGroup + "/"+experimentId+"-best.csv", bestIndividuals, false);

    delete[] results[0];
    delete[] results[1];
    delete[] results[2];
    delete[] results;
    clearContext(&ctx);
}

void runCECComparisonExperiment2(string grnMode, string evalMode, string expName, string outputDir)
{
    appContext ctx{};
    double (*func)(void*,void*);

    int maxEvals =  105*10000;

    //firt part populational algorithms
    int numParents = 15;
    int numOffspring = 105;
    int numRuns = 30;

    if(grnMode == "grn5"){
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn5EvaluationLSODA;
            initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
        }
        else {
            func = &grn5EvaluationRK4;
            initializeGRN5Context(&ctx, SINGLE_SET_MODE, 20);
        }
    }
    else {
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODA;
            initializeGRN10Context(&ctx, SINGLE_SET_MODE, 1);
        }
        else {
            func = &grn10EvaluationRK4;
            initializeGRN10Context(&ctx, SINGLE_SET_MODE, 20);
        }
    }

    int maxGenerations = maxEvals / numOffspring;
    string experimentId = grnMode + "-" + evalMode + "-" + to_string(numRuns) + "runs-"+ to_string(maxEvals) + "evals";
    string experimentGroup = expName;

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(func);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);

    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }



    /*todo: usar representação contígua para essa matriz
     * e, para calcular posição, usar
     * #define R(i,j) results[i*5 + j]
     */

    double **results = new double *[3];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];

    vector<string> bestInds(3);
    vector<double> bestIndsEval(3, DBL_MAX);

    for (int i = 0; i < numRuns; i++)
    {
        auto beg = chrono::high_resolution_clock::now();

        cout << "Run " << to_string(i) << "\n";
        outputToFile("output.txt", "Run " + to_string(i) + "\n", true);

        cout << "1"
             << "\n";

        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, maxEvals);
        results[0][i] = esAlgorithm.getPopulation().back()->getEvaluation();

        if (results[0][i] < bestIndsEval[0])
        {
            bestIndsEval[0] = results[0][i];
            bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }
        cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

        cout << "2"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[1][i] < bestIndsEval[1])
        {
            bestIndsEval[1] = results[1][i];
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }
        cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

        cout << "3"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        results[2][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[2][i] < bestIndsEval[2])
        {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }
        cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

        //temporização
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
        outputToFile(outputDir  + experimentGroup + "/" + experimentId+" -time.csv", to_string(duration.count()) + ",", true);
        cout << "Elapsed Time: " << duration.count() << "\n";
    }

    //todo: modificar para refletir os parâmetros da execução
    string csvOutput = "1+1,15+105-i,15+105-ni\n";

    for (int j = 0; j < numRuns; j++)
    {
        csvOutput += to_string(results[0][j]) + "," + to_string(results[1][j]) + "," + to_string(results[2][j])  + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2];

    outputToFile(outputDir +  experimentGroup + "/"+experimentId+".csv", csvOutput, false);
    outputToFile(outputDir +  experimentGroup + "/"+experimentId+"-best.csv", bestIndividuals, false);

    delete[] results[0];
    delete[] results[1];
    delete[] results[2];
    delete[] results;
    clearContext(&ctx);
}


void runCMAESComparisonExperimentTrainingTest(string grnMode, string evalMode, string expName, string outputDir)
{
    appContext ctx{};
    double (*func)(void*,void*);

    int maxEvals =  105*10000;

    //firt part populational algorithms
    int numParents = 15;
    int numOffspring = 105;
    int numRuns = 30;

    if(grnMode == "grn5"){
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn5EvaluationLSODA;
            initializeGRN5Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn5EvaluationRK4;
            initializeGRN5Context(&ctx, TRAINING_MODE, 20);
        }
    }
    else {
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODA;
            initializeGRN10Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn10EvaluationRK4;
            initializeGRN10Context(&ctx, TRAINING_MODE, 20);
        }
    }

    int maxGenerations = maxEvals / numOffspring;
    string experimentId = grnMode + "-" + evalMode + "-" + to_string(numRuns) + "runs-"+ to_string(maxEvals) + "evals";
    string experimentGroup = expName;

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(func);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);

    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }



    /*todo: usar representação contígua para essa matriz
     * e, para calcular posição, usar
     * #define R(i,j) results[i*5 + j]
     */

    double **results = new double *[3];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];

    vector<string> bestInds(3);
    vector<double> bestIndsEval(3, DBL_MAX);

    for (int i = 0; i < numRuns; i++)
    {
        auto beg = chrono::high_resolution_clock::now();

        cout << "Run " << to_string(i) << "\n";
        outputToFile("output.txt", "Run " + to_string(i) + "\n", true);

        cout << "1"
             << "\n";

       // esAlgorithm.runPopulationalIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
       // GRNEDOHelpers::setMode(&ctx, TEST_MODE);
       // esAlgorithm.evaluate(esAlgorithm.getPopulation()[0]);
        //results[0][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        results[0][i] = 0.0;
        //GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);

        if (results[0][i] < bestIndsEval[0])
        {
           // bestIndsEval[0] = results[0][i];
           // bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }
        //cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

        cout << "2"
             << "\n";
        //esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, maxEvals, numParents, numOffspring);
        //GRNEDOHelpers::setMode(&ctx, TEST_MODE);
       // esAlgorithm.evaluate(esAlgorithm.getPopulation()[0]);
        //results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        results[1][i] = 0.0;
       // GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);

        if (results[1][i] < bestIndsEval[1])
        {
            //bestIndsEval[1] = results[1][i];
            //bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }
        //cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

        cout << "3"
             << "\n";
        esAlgorithm.runCMAES(i, maxEvals, 40);
        GRNEDOHelpers::setMode(&ctx, TEST_MODE);
        esAlgorithm.evaluate(esAlgorithm.getBestIndividual());
        results[2][i] = esAlgorithm.getBestIndividual()->getEvaluation();
        GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);

        if (results[2][i] < bestIndsEval[2])
        {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getBestIndividual()->toCSVString();
        }
        cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

        //temporização
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
        outputToFile(outputDir  + experimentGroup + "/" + experimentId+"-train-test-time.csv", to_string(duration.count()) + ",", true);
        cout << "Elapsed Time: " << duration.count() << "\n";
    }

    //todo: modificar para refletir os parâmetros da execução
    string csvOutput = "15+105-i,15+105-ni,CMAES\n";

    for (int j = 0; j < numRuns; j++)
    {
        csvOutput += to_string(results[0][j]) + "," + to_string(results[1][j]) + "," + to_string(results[2][j])  + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2];

    outputToFile(outputDir +  experimentGroup + "/"+experimentId+"-train-test.csv", csvOutput, false);
    outputToFile(outputDir +  experimentGroup + "/"+experimentId+"-train-test-best.csv", bestIndividuals, false);

    delete[] results[0];
    delete[] results[1];
    delete[] results[2];
    delete[] results;
    clearContext(&ctx);
}

void runExperimentRound(string grnMode, string evalMode, string algName, int maxEvals, int seed)
{
    appContext ctx{};
    double (*func)(void*,void*);

    //int maxEvals =  105*10000;

    //firt part populational algorithms
    int numParents = 15;
    int numOffspring = 105;

    if(grnMode == "grn5"){
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn5EvaluationLSODA;
            initializeGRN5Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn5EvaluationRK4;
            initializeGRN5Context(&ctx, TRAINING_MODE, 20);
        }
    }
    else {
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODA;
            initializeGRN10Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn10EvaluationRK4;
            initializeGRN10Context(&ctx, TRAINING_MODE, 20);
        }
    }

    int maxGenerations = maxEvals / numOffspring;

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(func);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);

    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }


    //string resultCsv = "seed,eval,time,numEvals,ind\n";
    string resultCsv = "";
    auto beg = chrono::high_resolution_clock::now();
    Individual *bestInd = nullptr;
    double  bestEval = 0;

    //cout << "Running " << "\n";

    if(algName=="cmaes"){
        esAlgorithm.runCMAES(seed, maxEvals, 40);
    }else if(algName=="es-i"){
        esAlgorithm.runPopulationalIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }else if(algName=="es-ni"){
        esAlgorithm.runPopulationalNonIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }
    else if(algName=="1+1"){
        esAlgorithm.run1Plus1ES(seed, 0.5, 0.817, 10, maxEvals);
    }else if(algName=="de") {
        esAlgorithm.runDE(seed, maxEvals, 40);
    }else if(algName=="sade") {
        esAlgorithm.runSADE(seed, maxEvals, 40);
    }
    GRNEDOHelpers::setMode(&ctx, TEST_MODE);
    esAlgorithm.evaluate(esAlgorithm.getBestIndividual());
    bestInd = esAlgorithm.getBestIndividual();
    GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);

    //cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

    //temporização
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
    //cout << "Elapsed Time: " << duration.count() << "\n";
    resultCsv += to_string(seed) +","
            + to_string(bestInd->getEvaluation()) + ","
            + to_string(duration.count()) + ","
            + to_string(esAlgorithm.getEvaluations()) + ",["
            + bestInd->toCSVString()+ "]";

    cout << resultCsv << endl;

    clearContext(&ctx);
}


void runExperimentRoundTest(Algorithm& esAlgorithm, string algName, int maxEvals, int seed)
{
    //string resultCsv = "seed,eval,time,numEvals,ind\n";
    string resultCsv = "";
    auto beg = chrono::high_resolution_clock::now();
    Individual *bestInd = nullptr;
    double  bestEval = 0;

    if(algName=="cmaes"){
        esAlgorithm.runCMAES(seed, maxEvals, 40);
    }else if(algName=="es-i"){
        esAlgorithm.runPopulationalIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }else if(algName=="es-ni"){
        esAlgorithm.runPopulationalNonIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }
    else if(algName=="1+1"){
        esAlgorithm.run1Plus1ES(seed, 0.5, 0.817, 10, maxEvals);
    }else if(algName=="de") {
        esAlgorithm.runDE(seed, maxEvals, 40);
    }else if(algName=="sade") {
        esAlgorithm.runSADE(seed, maxEvals, 40);
    }

    //todo: uma função que reavalia população segundo conjunto de teste
    //  esAlgorithm.evaluate(esAlgorithm.getBestIndividual(), test);

    //GRNEDOHelpers::setMode(&ctx, TEST_MODE);
    //esAlgorithm.evaluate(esAlgorithm.getBestIndividual());
    esAlgorithm.reevaluateBestIndividualUsingTestSet();
    bestInd = esAlgorithm.getBestIndividual();
    //GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);


    //temporização
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
    resultCsv += to_string(seed) +","
                 + to_string(bestInd->getEvaluation()) + ","
                 + to_string(duration.count()) + ","
                 + to_string(esAlgorithm.getEvaluations()) + ",["
                 + bestInd->toCSVString()+ "]";

    cout << resultCsv << endl;

    //clearContext2Test(&ctx);
}

void runBoundlessExperimentRound(string grnMode, string evalMode, string algName, int maxEvals, int seed)
{
    appContext ctx{};
    double (*func)(void*,void*);

    //int maxEvals =  105*10000;

    //firt part populational algorithms
    int numParents = 15;
    int numOffspring = 105;

    if(grnMode == "grn5"){
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn5EvaluationLSODA;
            initializeGRN5Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn5EvaluationRK4;
            initializeGRN5Context(&ctx, TRAINING_MODE, 20);
        }
    }
    else {
        //maxEvals =  105*10000;
        //numParents = 15;
        //numOffspring = 105;

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODA;
            initializeGRN10Context(&ctx, TRAINING_MODE, 1);
        }
        else {
            func = &grn10EvaluationRK4;
            initializeGRN10Context(&ctx, TRAINING_MODE, 20);
        }
    }


    ctx.MIN_K = 0.0;


    int maxGenerations = maxEvals / numOffspring;

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(func);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);


    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }


    //string resultCsv = "seed,eval,time,numEvals,ind\n";
    string resultCsv = "";
    auto beg = chrono::high_resolution_clock::now();
    Individual *bestInd = nullptr;
    double  bestEval = 0;

    //cout << "Running " << "\n";

    if(algName=="cmaes"){
        esAlgorithm.runCMAES(seed, maxEvals, 40);
    }else if(algName=="es-i"){
        esAlgorithm.runPopulationalIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }else if(algName=="es-ni"){
        esAlgorithm.runPopulationalNonIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }
    else if(algName=="1+1"){
        esAlgorithm.run1Plus1ES(seed, 0.5, 0.817, 10, maxEvals);
    }

    GRNEDOHelpers::setMode(&ctx, TEST_MODE);
    esAlgorithm.evaluate(esAlgorithm.getBestIndividual());
    bestInd = esAlgorithm.getBestIndividual();
    GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);

    //cout << "Evals: " << esAlgorithm.getEvaluations() << "\n";

    //temporização
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
    //cout << "Elapsed Time: " << duration.count() << "\n";
    resultCsv += to_string(seed) +","
                 + to_string(bestInd->getEvaluation()) + ","
                 + to_string(duration.count()) + ","
                 + to_string(esAlgorithm.getEvaluations()) + ",["
                 + bestInd->toCSVString()+ "]";

    cout << resultCsv << endl;

    clearContext(&ctx);
}

void runGRN5ESComparisonExperiment(string evalMode, string experimentId)
{
    appContext ctx{};
    initializeGRN5Context(&ctx, TRAINING_MODE, 1);

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(grn5EvaluationLSODA);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);
    int numRuns = 1;


    string experimentPrefix = "grn5-" + evalMode + "-" + to_string(numRuns) + "runs-"+ to_string(2000) + "evals";

    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }


    /*todo: usar representação contígua para essa matriz
     * e, para calcular posição, usar
     * #define R(i,j) results[i*5 + j]
     */
    double **results = new double *[5];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];
    results[3] = new double[numRuns];
    results[4] = new double[numRuns];

    vector<string> bestInds(5);
    vector<double> bestIndsEval(5, DBL_MAX);


    for (int i = 0; i < numRuns; i++)
    {
        auto beg = chrono::high_resolution_clock::now();

        cout << "Run " << to_string(i) << "\n";

        cout << "1"
             << "\n";
        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, 2000);
        results[0][i] = esAlgorithm.getPopulation().back()->getEvaluation();
        if (results[0][i] < bestIndsEval[0])
        {
            bestIndsEval[0] = results[0][i];
            bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }

        cout << "2"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 100, 10, 20);
        results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[1][i] < bestIndsEval[1])
        {
            bestIndsEval[1] = results[1][i];
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "3"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 100, 10, 20);
        results[2][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[2][i] < bestIndsEval[2])
        {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "4"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 200, 5, 10);
        results[3][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[3][i] < bestIndsEval[3])
        {
            bestIndsEval[3] = results[3][i];
            bestInds[3] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "5"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 200, 5, 10);
        results[4][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[4][i] < bestIndsEval[4])
        {
            bestIndsEval[4] = results[4][i];
            bestInds[4] = esAlgorithm.getPopulation()[0]->toCSVString();
        }


        //temporização
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
        outputToFile("../results/" + experimentId +"/"+experimentPrefix+"-time.csv", to_string(duration.count()) + ",", true);
        cout << "Elapsed Time: " << duration.count() << "\n";
    }

    string csvOutput = "1+1,10+20-i,10+20-ni,5+10-i,5+10-ni\n";
    for (int j = 0; j < numRuns; j++)
    {
        csvOutput += to_string(results[0][j]) + "," + to_string(results[1][j]) + "," + to_string(results[2][j]) + "," + to_string(results[3][j]) + "," + to_string(results[4][j]) + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2] + "\n" + bestInds[3] + "\n" + bestInds[4];

    outputToFile("../results/" + experimentId +"/"+experimentPrefix+".csv", csvOutput, false);
    outputToFile("../results/" + experimentId +"/"+experimentPrefix+"-best_inds.csv", bestIndividuals, false);

    delete[] results[0];
    delete[] results[1];
    delete[] results[2];
    delete[] results[3];
    delete[] results[4];
    delete[] results;
    clearContext(&ctx);
}

void runGRN10ESComparisonExperiment()
{
    appContext ctx{};
    initializeGRN10Context(&ctx, TRAINING_MODE, 1);

    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(grn10EvaluationLSODA);
    esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);

    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        esAlgorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    int numRuns = 30;

    double **results = new double *[5];
    results[0] = new double[numRuns];
    results[1] = new double[numRuns];
    results[2] = new double[numRuns];
    results[3] = new double[numRuns];
    results[4] = new double[numRuns];

    vector<string> bestInds(5);
    vector<double> bestIndsEval(5, DBL_MAX);

    for (int i = 0; i < numRuns; i++)
    {
        cout << "Run " << to_string(i) << endl;

        cout << "0"
             << "\n";
        esAlgorithm.run1Plus1ES(i, 0.5, 0.817, 10, 2000);
        results[0][i] = esAlgorithm.getPopulation().back()->getEvaluation();
        if (results[0][i] < bestIndsEval[0])
        {
            bestIndsEval[0] = results[0][i];
            bestInds[0] = esAlgorithm.getPopulation().back()->toCSVString();
        }
        cout << "1"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 1000, 10, 20);
        results[1][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[1][i] < bestIndsEval[1])
        {
            bestIndsEval[1] = results[1][i];
            bestInds[1] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "2"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 1000, 10, 20);
        results[2][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[2][i] < bestIndsEval[2])
        {
            bestIndsEval[2] = results[2][i];
            bestInds[2] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "3"
             << "\n";
        esAlgorithm.runPopulationalIsotropicES(i, 0.5, 1000, 10, 20);
        results[3][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[3][i] < bestIndsEval[3])
        {
            bestIndsEval[3] = results[3][i];
            bestInds[3] = esAlgorithm.getPopulation()[0]->toCSVString();
        }

        cout << "4"
             << "\n";
        esAlgorithm.runPopulationalNonIsotropicES(i, 0.5, 1000, 10, 20);
        results[4][i] = esAlgorithm.getPopulation()[0]->getEvaluation();
        if (results[4][i] < bestIndsEval[4])
        {
            bestIndsEval[4] = results[4][i];
            bestInds[4] = esAlgorithm.getPopulation()[0]->toCSVString();
        }
    }

    string csvOutput = "1+1,10+20-i,10+20-ni,10+20-i,10+20-ni\n";
    for (int j = 0; j < numRuns; j++)
    {
        csvOutput += to_string(results[0][j]) + "," + to_string(results[1][j]) + "," + to_string(results[2][j]) + "," + to_string(results[3][j]) + "," + to_string(results[4][j]) + "\n";
    }

    string bestIndividuals =
            bestInds[0] + "\n" + bestInds[1] + "\n" + bestInds[2] + "\n" + bestInds[3] + "\n" + bestInds[4];

    outputToFile("../lsoda-comparison-GRN10-30runs-20000it.csv", csvOutput, false);
    outputToFile("../lsoda-best-individuals-GRN10.txt", bestIndividuals, false);

    delete[] results[0];
    delete[] results[1];
    delete[] results[2];
    delete[] results[3];
    delete[] results[4];
    delete[] results;

    clearContext(&ctx);
}


void test(){
    double ind_5[] = {1.25, 4, 1.02, 1.57, 3.43, 0.72, 0.5, 0.45, 0.51, 0.52, 13, 4, 3, 4, 16};
    double ind_5_2[] ={1.2163355099083872, 1.1264485098219865, 2.973714367061704, 2.952143123315177, 2.998260518457365, 0.5687249950503857, 0.4580723119903261, 0.46214892372246563, 0.6182568295500336, 0.5213082492659304, 0.7708877748759901, 0.1497642024548283, 4.254757908429968, 3.759370669969996, 4.784173526119725, 10.935884810737809, 24.595975874929724, 2.8109199678182635, 4.922623602327875};
    double ind_10[] = {1.73,2,0.81,0.11, 1.23, 1.78, 1.14, 1.04, 3.47, 3.21,
                       0.45, 0.56, 0.99, 0.77, 0.71, 0.66, 0.46, 0.48, 0.66, 0.99, 0.85, 0.61, 0.55, 0.46, 0.17,
                       20, 9, 24, 12, 2, 2, 6, 4, 7, 24, 2, 7, 21, 20, 3};
    appContext ctx{};

    initializeGRN5Context(&ctx, TRAINING_MODE,1);
    cout << to_string(grn5EvaluationRK4(ind_5_2, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN5Context(&ctx, TRAINING_MODE,1);
    cout << to_string(grn5EvaluationLSODA(ind_5_2, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN10Context(&ctx, TRAINING_MODE, 1);
    cout << to_string(grn10EvaluationRK4(ind_10, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN10Context(&ctx, TRAINING_MODE, 1);
    cout << to_string(grn10EvaluationLSODA(ind_10, &ctx)) << "\n";
    clearContext(&ctx);
}

void testGRN5LSODARK4(){
    appContext ctx{};

    //esse indivíduo deveria ter fitness ~26
    double ind0[19]={1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                     2.952143123315177, 2.998260518457365, 0.5687249950503857,
                     0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                     0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                     4.254757908429968, 3.759370669969996, 4.784173526119725,
                     10.935884810737809, 24.595975874929724, 2.8109199678182635,
                     4.922623602327875};

    double ind[19] = {4.052585, 1.136443, 0.852662, 3.120117, 0.104307, 0.754601,
                      0.475951, 0.712506, 0.852778, 0.845837, 0.173564, 0.666615,
                      6.125123, 10.331758, 23.741081, 7.461154, 7.816841, 19.812765,
                      7.771741};
    double ind2[19] ={3.829097,2.148080,5.000000,5.000000,0.100000,0.100000,
                      1.000000,0.100000,0.100000,1.000000,0.100000,0.100000,
                      15.128335,4.608234,25.000000,1.000000,21.855324,10.160282,
                      1.000000};

    double ind3[19] = { 0.1,5,0.851307,3.865580,0.1,
                        0.100000,1.000000,0.100000,0.100000,1.000000,
                        1.000000,0.100000,1.224451,6.558196,13.511094,
                        5.650953,22.882314,11.149282,1.000000};

    double ind4[19] = {2.697747,3.016997,0.100000,4.891896,0.100000,
                       0.100000,1.000000,0.100000,0.100000,
                       1.000000,0.100000,0.100000,1.000000,
                       13.671514,23.013432,8.502602,25.000000,
                       1.420259,1.262061};

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 10);
    cout << to_string(grn5EvaluationRK4(ind0, &ctx)) << "\n";
    clearContext(&ctx);

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    cout << to_string(grn5EvaluationLSODA(ind0, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 10);
    cout << to_string(grn5EvaluationRK4(ind3, &ctx)) << "\n";
    clearContext(&ctx);

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    cout << to_string(grn5EvaluationLSODA(ind3, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 10);
    cout << to_string(grn5EvaluationRK4(ind4, &ctx)) << "\n";
    clearContext(&ctx);

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    cout << to_string(grn5EvaluationLSODA(ind4, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 10);
    cout << to_string(grn5EvaluationRK4(ind, &ctx)) << "\n";
    clearContext(&ctx);

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    cout << to_string(grn5EvaluationLSODA(ind, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 10);
    cout << to_string(grn5EvaluationRK4(ind2, &ctx)) << "\n";
    clearContext(&ctx);

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    cout << to_string(grn5EvaluationLSODA(ind2, &ctx)) << "\n";
    clearContext(&ctx);
}

void testGRN10LSODARK4(){
    appContext ctx{};

    //deveria ter fitness ~56
    double ind[40] = {1.73,2,0.81,0.11, 1.23, 1.78,
                      1.14, 1.04, 3.47, 3.21, 0.45,
                      0.56, 0.99, 0.77, 0.71, 0.66,
                      0.46, 0.48, 0.66, 0.99, 0.85,
                      0.61, 0.55, 0.46, 0.17, 20,
                      9, 24, 12, 2, 2, 6, 4, 7,
                      24, 2, 7, 21, 20, 3 };


    initializeGRN10Context(&ctx, TRAINING_MODE, 20);
     cout << to_string(grn10EvaluationRK4(ind, &ctx)) << "\n";
    clearContext(&ctx);


    initializeGRN10Context(&ctx, TRAINING_MODE, 1);
    cout << to_string(grn10EvaluationLSODA(ind, &ctx)) << "\n";
    clearContext(&ctx);

}

void testCMAES(){
    appContext ctx{};
    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    GRNCoefProblem problem = GRNCoefProblem(&ctx);
    problem.setEvaluationFunction(grn5EvaluationLSODA);
    population pop = population(problem, 20, 0);
    cmaes alg = cmaes(10000, -1, -1, -1, -1, 0.5, 1e-6, 1e-6, true, true, 0);
    alg.set_verbosity(100);
    population newPop = alg.evolve(pop);
    cout << "Total evaluations: " <<newPop.get_problem().get_fevals()<< "\n";
    cout << "Best fitness: " << newPop.champion_f()[0]<< endl;
    cout << vectorToString(newPop.champion_x().data(),0, 18);
}

void testTestSet(){
    appContext context{};

    initializeGRN5Context(&context, TEST_MODE, 1);
    appContext* ctx = &context;
    //GRNEDOHelpers::setMode(&ctx, TEST_MODE);

    //appContext* ctx = (appContext*)(data);
    double ind[] = {  1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                      2.952143123315177, 2.998260518457365, 0.5687249950503857,
                      0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                      0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                      4.254757908429968, 3.759370669969996, 4.784173526119725,
                      10.935884810737809, 24.595975874929724, 2.8109199678182635,
                      4.922623602327875};
    //double* _ind = (double *)ind;
    /*ctx->individual = reinterpret_cast<double *>(&ind);
    lsodaWrapper(twoBody5VarLSODA, ctx, ctx->yout);
    double eval = difference(ctx->yout, ctx->expectedResult, ctx->nVariables, ctx->setStart, ctx->setEnd, ctx->nSteps);
*/

   double eval = grn5EvaluationLSODA(ind, ctx);
    //printGRNVector(ctx->expectedResult, 5, 50);
    //cout << "\n\n";
    cout << eval << endl;

}

void testCMAES2(){
    appContext ctx{};
    initializeGRN5Context(&ctx, TRAINING_MODE, 1);
    Algorithm esAlgorithm = Algorithm(ctx.IND_SIZE);
    esAlgorithm.setEvaluationFunction(grn5EvaluationLSODA);
    //esAlgorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);
    esAlgorithm.setContext(&ctx);
    esAlgorithm.runCMAES(0, 100000, 20);
    GRNEDOHelpers::setMode(&ctx, TEST_MODE);
    esAlgorithm.evaluate(esAlgorithm.getBestIndividual());

    cout << "Best fitness: " << esAlgorithm.getBestIndividual()->getEvaluation() << endl;
}

int main2(int argc, char** argv)
{
    string grnMode;
    string evalMode;
    string algName;

    if (argc == 6) {
        if(strcmp(argv[1], "grn5") == 0){
            grnMode = "grn5";
        }else if(strcmp(argv[1], "grn10") == 0) {
            grnMode = "grn10";
        }else{
            cout << "Invalid GRN mode" << endl;
        }

        if(strcmp(argv[2], "lsoda") == 0) {
            evalMode = "lsoda";
        }else if(strcmp(argv[2], "rk4") == 0) {
            evalMode = "rk4";
        }else{
            cout << "Invalid evaluation mode" << endl;
        }

        if(strcmp(argv[3], "cmaes") == 0 ||
        strcmp(argv[3], "es-i") == 0 ||
        strcmp(argv[3], "es-ni"  ) == 0 ||
        strcmp(argv[3], "de") == 0 ||
        strcmp(argv[3], "sade") == 0) {
            algName = argv[3];
        }else  {
            cout << "Invalid algorithm name" << endl;
        }

        int maxEvals = atoi(argv[4]);
        int seed = atoi(argv[5]);
        runExperimentRound(grnMode, evalMode, algName, maxEvals, seed);
    }

    //runExperimentRound("grn5", "lsoda", "cmaes", 1000, 0);
    return 0;


}

string findArg(int argc, char** argv, string argName){
    for(int i = 0; i < argc; i++){
        if(strcmp(argv[i], argName.c_str()) == 0){
            return string(argv[i+1]);
        }
    }
    return "";
}

vector<string> findArgList(int argc, char** argv, string argName, int listSize){
    vector<string> args;

    for(int i = 0; i < argc; i++){
        if(strcmp(argv[i], argName.c_str()) == 0){
            for(int j=0; j < listSize; j++){
                args.push_back(string(argv[i+j+1]));
            }
        }
    }
    return args;
}

bool argExists(int argc, char** argv, string argName){
    bool argFound = false;

    for(int i = 0; i < argc && !argFound; i++){
        if(strcmp(argv[i], argName.c_str()) == 0){
            argFound = true;
        }
    }
    return argFound;
}

std::map<string, string> parseArgs(int argc, char** argv){
    std::map<string, string> args;
    //if (argc == 13) {
        string inputFile = findArg(argc, argv, "-i");
        string grnModel = findArg(argc, argv, "-m");
        string evalMode = findArg(argc, argv, "-e");
        string algName = findArg(argc, argv, "-a");
        string maxEvals = findArg(argc, argv, "-n");
        string seed = findArg(argc, argv, "-s");

        if(argExists(argc, argv, "-ts")){
            args["testSet"] = findArg(argc, argv, "-ts");
        }
        if(argExists(argc, argv, "-sd")){
            vector<string> setIntervals = findArgList(argc, argv, "-sd", 4);
            args["trainingStart"] = setIntervals[0];
            args["trainingEnd"] = setIntervals[1];
            args["testStart"] = setIntervals[2];
            args["testEnd"] = setIntervals[3];
        }

        args["inputFile"] = inputFile;

        if(grnModel == "grn5"|| grnModel == "grn10"){
            args["grnModel"] = grnModel;
        }else{
            cout << "Invalid GRN model" << endl;
        }

        if(evalMode == "lsoda" || evalMode == "rk4") {
            args["evalMode"] = evalMode;
        }else{
            cout << "Invalid evaluation mode" << endl;
        }

        if(algName == "cmaes" || algName == "es-i" || algName == "es-ni"||
        algName == "de" ||algName == "sade" || algName == "1+1"){
            args["algName"] = algName;
        }
        else  {
            cout << "Invalid algorithm name" << endl;
        }

        args["maxEvals"] = maxEvals;
        args["seed"] = seed;
   // }else{
  //      cout << "Invalid number of arguments" << endl;
   // }

    return args;
}

int main(int argc, char** argv)
{
    std::map<string, string> args = parseArgs(argc, argv);
    string grnModel = args["grnModel"];
    string evalMode = args["evalMode"];
    string algName = args["algName"];
    string inputFile = args["inputFile"];
    int seed = stoi(args["seed"]);
    int maxEvals = stoi(args["maxEvals"]);

    //appContext ctx{};
    double (*func)(void*,void*);
    ProblemDescription ctx = GRNEDOHelpers::grn5ProblemDescription;

    //todo: usar novo initialize
    //todo: unir grn5EvaluationLSODA e grn10EvaluationLSODA em uma só função
    //todo: fazer o mesmo para RK4

    if(grnModel == "grn5"){
        ctx = {.IND_SIZE = 19, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 5,.N_SIZE = 7,.K_SIZE = 7};

        if(evalMode == "lsoda"){
            func = &grnEvaluationLSODATest;
        }
        else if(evalMode == "rk4"){
            func = &grnEvaluationRK4Test;
        }
    }
    else if(grnModel == "grn10"){
        ctx = {.IND_SIZE = 40, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 10,.N_SIZE = 15,.K_SIZE = 15};

        if(evalMode == "lsoda"){
            func = &grn10EvaluationLSODATest;
        }
        else if(evalMode == "rk4"){
            func = &grn10EvaluationRK4Test;
        }
    }


//    GRNSeries series = GRNSeries(inputFile);
//    GRNSeries trainingSeries = GRNSeries(series, 0, 34);
//    GRNSeries testSeries = GRNSeries(series, 35, 49);
    GRNSeries* inputSeries = new GRNSeries(inputFile);
    GRNSeries* trainingSeries;
    GRNSeries* testSeries;
   // GRNSeries trainingSeries;
    //GRNSeries testSeries;

    if(args.find("testSet") != args.end()){
        trainingSeries = new GRNSeries(inputFile);
        testSeries = new GRNSeries(args["testSet"]);

    }else if(args.find("trainingStart") != args.end()){
        trainingSeries = new GRNSeries(*inputSeries, stoi(args["trainingStart"]), stoi(args["trainingEnd"]));
        testSeries = new GRNSeries(*inputSeries, stoi(args["testStart"]), stoi(args["testEnd"]));
    }else {
        trainingSeries = new GRNSeries(inputFile);
        testSeries = new GRNSeries(inputFile);
    }

    Algorithm algorithm = Algorithm(*trainingSeries, *testSeries, func, ctx.IND_SIZE);
    algorithm.setContext(&ctx);
    algorithm.setSigmaBounds(ctx.MIN_STRATEGY, ctx.MAX_STRATEGY);

    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < ctx.TAU_SIZE; i++)
    {
        algorithm.setBounds(i, ctx.MIN_TAU, ctx.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE; i++)
    {
        algorithm.setBounds(i, ctx.MIN_K, ctx.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < ctx.TAU_SIZE + ctx.K_SIZE + ctx.N_SIZE; i++)
    {
        algorithm.setBounds(i, ctx.MIN_N, ctx.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    //todo: falta usar de fato os conjuntos de treino e teste
    runExperimentRoundTest(algorithm, algName, maxEvals, seed);

    delete inputSeries;
    delete trainingSeries;
    delete testSeries;
    return 0;


}

