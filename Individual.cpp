//
// Created by patri on 27/01/2023.
//

#include "Individual.h"

Individual::Individual(int numDimensions) {
    this->dimensions.resize(numDimensions);
    this->numDimensions = numDimensions;
}

Individual::~Individual() {

}

void Individual::setEvaluation(double val) {
    this->evaluation = val;
}

double Individual::getEvaluation() {
    return this->evaluation;
}

void Individual::setDimension(int index, double val){
    this->dimensions[index] = val;
}
