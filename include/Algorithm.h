//
// Created by patri on 27/01/2023.
//

#ifndef ES_ALGORITHM_H
#define ES_ALGORITHM_H
#include <vector>
#include "dependencies.h"
#include "Individual.h"
#include "GRNCoefProblem.h"
#include "GRNSeries.h"

using namespace std;
class Algorithm {

private:
    vector<Individual*> population;
    vector<double> upperBounds;
    vector<double> lowerBounds;
    vector<int> upperBoundTypes;
    vector<int> lowerBoundTypes;
    int numDimensions;
    double (*evaluationFunction)(void*, void*);
    int (*grnModel)(double t, double *y, double *ydot, void *context);
    int algorithmType;
    int evaluationsCounter;
    double minSigma;
    double maxSigma;
    void *context;
    Individual* bestIndividual;
    GRNCoefProblem* problem;
    GRNSeries* trainingSeries;
    GRNSeries* testSeries;

public:
    static int UPPER_OPEN;
    static int UPPER_CLOSED;
    static int LOWER_OPEN;
    static int LOWER_CLOSED;
    static int UPPER;
    static int LOWER;
    static int ONE_PLUS_ONE;
    static int ISOTROPIC;
    static int NON_ISOTROPIC;

    Algorithm(int numDimensions);
    Algorithm(GRNSeries &trainingSeries, GRNSeries &testSeries,
              double (*evaluationFunction)(void *, void *), int numDimensions);
    ~Algorithm();
    vector<Individual*> getPopulation();
    void addIndividual(Individual* individual);
    double evaluate(Individual* ind);
    void setNumDimensions(int val);
    void setBounds(int index, double lower, double upper, int lowerBoundType, int upperBoundType);
    int getNumDimensions();
    void setEvaluationFunction(double (*evaluationFunction)(void*, void*));
    void setContext(void* ctx);
    void createPopulation(int seed, int numIndividuals);
    double getBound(int index, int which);
    void setAlgorithmType(int type);
    void run1Plus1ES(int seed, double initialSigma, double c, int n,  int maxEvals);
    double getMinSigma();
    double getMaxSigma();
    void validate(Individual* ind);
    string populationToCSVString();
    void runPopulationalIsotropicES(int seed, double sigmaVariation, int maxEvals, int numParents, int numOffspring);
    void runPopulationalNonIsotropicES(int seed, double sigmaVariation, int maxEvals, int numParents, int numOffspring);
    void clear();
    void setSigmaBounds(double min, double max);
    void reevaluateAllNoCounter();
    double getReevaluationByIndexNoCounter(int i);
    int getEvaluations();
    void runCMAES(int seed,int maxEvals, int populationSizes);
    void runDE(int seed,int maxEvals, int populationSizes);
    void runSADE(int seed,int maxEvals, int populationSize);
    Individual* getBestIndividual();
    double evaluationIncrementCounterWrapper(void *ind, void * context);
    void reevaluateBestIndividualUsingTestSet();

    void setGrnModel(int (*grnModel)(double, double *, double *, void *));

};


#endif //ES_ALGORITHM_H
