//
// Created by patri on 27/01/2023.
//

#ifndef ES_ESALGORITHM_H
#define ES_ESALGORITHM_H
#include <vector>
#include "Individual.h"
using namespace std;
class ESAlgorithm {

private:
    vector<Individual*> population;
    vector<double> upperBounds;
    vector<double> lowerBounds;
    vector<int> upperBoundTypes;
    vector<int> lowerBoundTypes;
    int numDimensions;
    double (*evaluationFunction)(void*, void*);
    int algorithmType;
    int evaluationsCounter;
    double minSigma;
    double maxSigma;
    void *context;

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

    ESAlgorithm(int numDimensions);
    ~ESAlgorithm();
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
    void run1Plus1ES(int seed, double initialSigma, double c, int n,  int maxIterations);
    double getMinSigma();
    double getMaxSigma();
    void validate(Individual* ind);
    string populationToCSVString();
    void runPopulationalIsotropicES(int seed, double sigmaVariation, int maxIterations, int numParents, int numOffspring);
    void runPopulationalNonIsotropicES(int seed, double sigmaVariation, int maxIterations, int numParents, int numOffspring);
    void clear();
    double setSigmaBounds(double min, double max);

};


#endif //ES_ESALGORITHM_H
