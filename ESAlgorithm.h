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
    double (*evaluationFunction)(Individual*);
    int algorithmType;
    int numParents;
    int numOffspring;
    int evaluationsCounter;

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

    ESAlgorithm(int numDimensions, int numParents, int numOffspring);
    ~ESAlgorithm();
    vector<Individual*> getPopulation();
    void addIndividual(Individual* individual);
    double evaluate(Individual* ind);
    void setNumDimensions(int val);
    void setBounds(int index, double lower, double upper, int lowerBoundType, int upperBoundType);
    int getNumDimensions();
    void setEvaluationFunction(double (*evaluationFunction)(Individual*));
    void createPopulation(int seed, int numIndividuals);
    double getBound(int index, int which);
    void setAlgorithmType(int type);
    
};


#endif //ES_ESALGORITHM_H
