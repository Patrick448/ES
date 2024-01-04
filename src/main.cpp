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
#include "GRN5Model.h"
#include "Algorithm.h"

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
    string grnModelName = args["grnModel"];
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

    if(grnModelName == "grn5"){
        ctx = {.IND_SIZE = 19, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 5,.N_SIZE = 7,.K_SIZE = 7, .modelFunction = GRNEDOHelpers::grn5Model};
    }
    else if(grnModelName == "grn10"){
        ctx = {.IND_SIZE = 40, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 10,.N_SIZE = 15,.K_SIZE = 15, .modelFunction = GRNEDOHelpers::grn10Model};
    }

    if(evalMode == "lsoda"){
        func = &grnEvaluationLSODA;
    }
    else if(evalMode == "rk4"){
        func = &grnEvaluationRK4;
    }

    GRNSeries* inputSeries = new GRNSeries(inputFile);
    GRNSeries* trainingSeries;
    GRNSeries* testSeries;

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

    Algorithm algorithm = Algorithm(*trainingSeries, *testSeries, func, ctx.IND_SIZE, &ctx);
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

