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
    int numDimensions;
    double evaluation;
    double globalSigma;
    double *maxValues;

public:
    Individual(int numDimensions);
    Individual(int numDimensions, double *dimensions);
    ~Individual();
    void setEvaluation(double val);
    double getEvaluation();
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

};


#endif //ES_INDIVIDUAL_H
