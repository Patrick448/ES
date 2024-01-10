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



void runExperimentRoundTest(Algorithm& algorithm, string algName, int maxEvals, int seed)
{
    //string resultCsv = "seed,eval,time,numEvals,ind\n";
    string resultCsv = "";
    auto beg = chrono::high_resolution_clock::now();
    Individual *bestInd = nullptr;
    double  bestEval = 0;

    if(algName=="cmaes"){
        algorithm.runCMAES(seed, maxEvals, 40);
    }else if(algName=="es-i"){
        algorithm.runPopulationalIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }else if(algName=="es-ni"){
        algorithm.runPopulationalNonIsotropicES(seed, 0.5, maxEvals, 15, 105);
    }
    else if(algName=="1+1"){
        algorithm.run1Plus1ES(seed, 0.5, 0.817, 10, maxEvals);
    }else if(algName=="de") {
        algorithm.runDE(seed, maxEvals, 40);
    }else if(algName=="sade") {
        algorithm.runSADE(seed, maxEvals, 40);
    }else{
        cout << "Invalid algorithm name: "<< algName << endl;
    }

    //todo: uma função que reavalia população segundo conjunto de teste
    //  esAlgorithm.evaluate(esAlgorithm.getBestIndividual(), test);

    //GRNEDOHelpers::setMode(&ctx, TEST_MODE);
    //esAlgorithm.evaluate(esAlgorithm.getBestIndividual());
    algorithm.reevaluateBestIndividualUsingTestSet();
    bestInd = algorithm.getBestIndividual();
    //GRNEDOHelpers::setMode(&ctx, TRAINING_MODE);


    //temporização
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(end - beg);
    resultCsv += to_string(seed) + ","
                 + to_string(bestInd->getEvaluation()) + ","
                 + to_string(duration.count()) + ","
                 + to_string(algorithm.getEvaluations()) + ",["
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

vector<double> csvLineToDoubleVector(string text)
{
    vector<double> result;
    string part;

    for (int i = 0; i < text.size(); i++)
    {
        if (text[i] == ',')
        {
            if(part != "") {

                result.push_back(stod(part));
                part = "";
            }

        }
        else
        {
            part += (text[i]);
        }
    }

    if(part != ""){
        result.push_back(stod(part));
    }
    return result;
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
        string individual = findArg(argc, argv, "-ind");
        string outputFile = findArg(argc, argv, "-o");

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
        args["grnModel"] = grnModel;
        args["evalMode"] = evalMode;
        args["algName"] = algName;
        args["maxEvals"] = maxEvals;
        args["seed"] = seed;
        args["individual"] = individual;
        args["outputFile"] = outputFile;
   // }else{
  //      cout << "Invalid number of arguments" << endl;
   // }

    return args;
}

void printGRNResultTest(){
    GRNSeries series = GRNSeries("DadosZeEduardo_DATA.txt");
    GRNSeries testSeries = GRNSeries(series, 114, 161);
    ProblemDescription desc =  {.IND_SIZE = 31, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
            .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 5,.N_SIZE = 13,.K_SIZE = 13, .modelFunction = GRNEDOHelpers::grn5NCYCModel};
    double ind0[31] = {0.509131859706313,4.705890160374611,0.100002011193567,
                       0.100000000000000,0.587427439207285,0.930307173541365,
                       0.271938661057075,0.995819658954900,0.440293368635028,
                       0.294531357448824,0.643266490779396,0.100000206828619,
                       0.999511289572582,0.443006784717132,0.999923201337715,
                       0.999971822072094,0.230921764968310,0.999981525668626,
                       1.027467752938993,20.590994565269891,1.000000000000000,
                       24.932168591371969,1.000000000000000,1.809261972447725,
                       1.603055004597870,25.000000000000000,24.960964294175983,
                       12.005951470838969,25.000000000000000,15.878469782118251,
                       24.523305696638271};

    appContext ctx{.series = &testSeries, .description = &desc};
    printODEIntSeries(ind0, &ctx, "",0);
}

int main(int argc, char** argv)
{
    //printGRNResultTest();
    //return 0;
    std::map<string, string> args = parseArgs(argc, argv);
    string grnModelName = args["grnModel"];
    string evalMode = args["evalMode"];
    string algName = args["algName"];
    string inputFile = args["inputFile"];
    string individual = args["individual"];
    string outputFile = args["outputFile"];
    int seed;
    int maxEvals;
    double (*func)(void*,void*);
    ProblemDescription description = GRNEDOHelpers::grn5ProblemDescription;

    if(args["seed"] != ""){
        seed = stoi(args["seed"]);
    }
    if(args["maxEvals"] != ""){
        maxEvals = stoi(args["maxEvals"]);
    }

    if(grnModelName == "grn5"){
        description = {.IND_SIZE = 19, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 5,.N_SIZE = 7,.K_SIZE = 7, .modelFunction = GRNEDOHelpers::grn5Model};
    }
    else if(grnModelName == "grn10"){
        description = {.IND_SIZE = 40, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 10,.N_SIZE = 15,.K_SIZE = 15, .modelFunction = GRNEDOHelpers::grn10Model};
    }
    else if(grnModelName == "grn5new"){
        description = {.IND_SIZE = 21, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 5,.N_SIZE = 8,.K_SIZE = 8, .modelFunction = GRNEDOHelpers::grn5NewModel};
    }
    else if(grnModelName == "grn10new"){
        description = {.IND_SIZE = 48, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 10,.N_SIZE = 19,.K_SIZE = 19, .modelFunction = GRNEDOHelpers::grn10NewModel};
    }else if(grnModelName == "grn10new2"){
        description = {.IND_SIZE = 60, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 10,.N_SIZE = 25,.K_SIZE = 25, .modelFunction = GRNEDOHelpers::grn10New2Model};
    }
    else if(grnModelName == "grn5ncyc"){
        //Modelo do DadosZeEduardo.txt
        description = {.IND_SIZE = 31, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 5,.N_SIZE = 13,.K_SIZE = 13, .modelFunction = GRNEDOHelpers::grn5NCYCModel};
    }else if(grnModelName == "grn4ncyc"){
        description = {.IND_SIZE = 26, .MIN_K = 0.1,.MAX_K = 1,.MIN_N = 1,.MAX_N = 25,.MIN_TAU = 0.1,.MAX_TAU = 5,
                .MIN_STRATEGY = 0.1,.MAX_STRATEGY = 10,.TAU_SIZE = 4,.N_SIZE = 11,.K_SIZE = 11, .modelFunction = GRNEDOHelpers::grn4NCYCModel};
    }else{
        cout << "Invalid GRN model: "<< grnModelName << endl;
    }

    if(evalMode == "lsoda"){
        func = &grnEvaluationLSODA;
    }
    else if(evalMode == "rk4"){
        //todo: corrigir - por enquanto o rk4 só vai funcionar se os intervalos forem uniformes
        func = &grnEvaluationRK4;
    }else{
        cout << "Invalid ODE solver name: "<<evalMode << endl;
    }

    GRNSeries* inputSeries = new GRNSeries(inputFile);
    GRNSeries* trainingSeries;
    GRNSeries* testSeries;

    if(args.find("testSet") != args.end()){
        //todo: aqui vai ter um problema com o maxValues, pois ele é inicializado com o conjunto que é carregado
        //  talvez seja melhor usar o maxValues como parte do indivíduo
        trainingSeries = new GRNSeries(inputFile);
        testSeries = new GRNSeries(args["testSet"]);
        //todo: pensar em uma forma melhor de fazer isso. Para manter os max values do conjunto de treino
        testSeries->loadMaxValuesFrom(*trainingSeries);

    }else if(args.find("trainingStart") != args.end()){
        trainingSeries = new GRNSeries(*inputSeries, stoi(args["trainingStart"]), stoi(args["trainingEnd"]));
        testSeries = new GRNSeries(*inputSeries, stoi(args["testStart"]), stoi(args["testEnd"]));
        testSeries->loadMaxValuesFrom(*trainingSeries);

    }else {
        trainingSeries = new GRNSeries(inputFile);
        testSeries = new GRNSeries(inputFile);
        testSeries->loadMaxValuesFrom(*trainingSeries);
    }

    if(argExists(argc, argv, "-ind")){
        appContext ctx{.series = testSeries, .description = &description};
        vector<double> indVec = csvLineToDoubleVector(individual);
        double* indArray = indVec.data();
        printODEIntSeries(indArray, &ctx, outputFile, 0);
        return 0;
    }


    Algorithm algorithm = Algorithm(*trainingSeries, *testSeries, func, description.IND_SIZE, &description);
    algorithm.setSigmaBounds(description.MIN_STRATEGY, description.MAX_STRATEGY);

   //todo: passar essas inicializações para dentro da classe Algorithm usando o objeto ctx
    // inicializa limites de tau, k e n
    int cont = 0;
    for (int i = 0; i < description.TAU_SIZE; i++)
    {
        algorithm.setBounds(i, description.MIN_TAU, description.MAX_TAU, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < description.TAU_SIZE + description.K_SIZE; i++)
    {
        algorithm.setBounds(i, description.MIN_K, description.MAX_K, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    for (int i = cont + 1; i < description.TAU_SIZE + description.K_SIZE + description.N_SIZE; i++)
    {
        algorithm.setBounds(i, description.MIN_N, description.MAX_N, Algorithm::LOWER_CLOSED, Algorithm::UPPER_CLOSED);
        cont = i;
    }

    runExperimentRoundTest(algorithm, algName, maxEvals, seed);

    delete inputSeries;
    delete trainingSeries;
    delete testSeries;
    return 0;

}

