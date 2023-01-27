//
// Created by patri on 27/01/2023.
//

#ifndef ES_INDIVIDUAL_H
#define ES_INDIVIDUAL_H
#include <vector>

using namespace std;

class Individual {

private:
    vector<double> dimensions;
    int* sigmas;
    int numDimensions;
    double evaluation;

public:
    Individual(int numDimensions);
    ~Individual();
    void setEvaluation(double val);
    double getEvaluation();
    void setDimension(int index, double val);

};


#endif //ES_INDIVIDUAL_H
