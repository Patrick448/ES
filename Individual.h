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
    vector<double> dimensions;
    vector<double> sigmas;
    int numDimensions;
    double evaluation;
    double globalSigma;

public:
    Individual(int numDimensions);
    ~Individual();
    void setEvaluation(double val);
    double getEvaluation();
    void setDimension(int index, double val);
    void setSigma(int index, double val);
    void setGlobalSigma(double val);
    double getGlobalSigma();
    double getDimension(int index);
    double getSigma(int index);
    string toString();

};


#endif //ES_INDIVIDUAL_H
