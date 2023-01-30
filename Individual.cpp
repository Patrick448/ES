//
// Created by patri on 27/01/2023.
//

#include "Individual.h"
#include <math.h>

Individual::Individual(int numDimensions) {
    this->dimensions.resize(numDimensions);
    this->sigmas.resize(numDimensions);
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

void Individual::setSigma(int index, double val){
    this->sigmas[index] = val;

}


double Individual::getDimension(int index){
    return this->dimensions[index];
}

double Individual::getSigma(int index){
    return this->sigmas[index];
}

void Individual::setGlobalSigma(double val){
    this->globalSigma = val;
}
double Individual::getGlobalSigma(){
    return this->globalSigma;
}

string doubleToString(double val, int precision){
    int next =(int)val;
    string str = to_string(next) + ".";
    double remaining = fabs(val - next);

    for(int i=0; i< precision; i++){
        int next = (int)(remaining*10);
        remaining = remaining*10 - next;
        str += to_string(next);
    }

    return str;
}

string Individual::toCSVString(){
    string indString = "";
    for(int i=0; i< this->numDimensions; i++){
        indString+= doubleToString(this->dimensions[i], 15) + ",";
    }

    indString+= doubleToString(this->evaluation, 15);
    return indString;
}


