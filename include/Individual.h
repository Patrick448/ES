//
// Created by patri on 27/01/2023.
//

#ifndef ES_INDIVIDUAL_H
#define ES_INDIVIDUAL_H
#include <vector>
#include <string>

using namespace std;

class Individual {

private:
    double* parameters;
    double* sigmas;
    int numParameters;
    double fitnessTrain;
    double fitnessTest;
    double relativeFitnessTrain;
    double relativeFitnessTest;
    double globalSigma;
    double *maxValues;
    double *fullParameters;

public:
    Individual(int numParameters);
    Individual(int numDimensions,double *parameters);
    ~Individual();
    void setFitnessTrain(double val);
    double getFitnessTrain();
    void setRelativeFitnessTrain(double val);
    double getRelativeFitnessTrain();
    void setRelativeFitnessTest(double val);
    double getRelativeFitnessTest();
    void setFitnessTest(double val);
    double getFitnessTest();
    void setParameter(int index, double val);
    void setSigma(int index, double val);
    void setGlobalSigma(double val);
    double getGlobalSigma();
    double getDimension(int index);
    //todo: rename to something like getCoefficients, getParameters, etc.
    double *getParameters();
    double getSigma(int index);
    string toCSVString();
    void setMaxValues(double *maxValues);
    double *getMaxValues() const;
    void setParameters(double *parameters);
    double *getFullParameters() const;


};


#endif //ES_INDIVIDUAL_H
