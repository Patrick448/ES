//
// Created by patri on 27/01/2023.
//

#include "Individual.h"
#include <math.h>

Individual::Individual(int numParameters) {
    this->numParameters = numParameters;
    this->parameters = new double[numParameters];
    this->sigmas = new double[numParameters];

}

Individual::Individual(int numDimensions, double *parameters) {
    this->numParameters = numDimensions;
    this->parameters = new double[numDimensions];

    for(int i=0; i< numDimensions; i++){
        this->parameters[i] = parameters[i];
    }
}

Individual::~Individual() {
    delete [] this->parameters;
    delete [] this->sigmas;
}

void Individual::setFitnessTrain(double val) {
    this->fitnessTrain = val;
}

double Individual::getFitnessTrain() {
    return this->fitnessTrain;
}

void Individual::setFitnessTest(double val) {
    this->fitnessTest = val;
}

double Individual::getFitnessTest() {
    return this->fitnessTest;
}

void Individual::setParameter(int index, double val){
    this->parameters[index] = val;
}

void Individual::setSigma(int index, double val){
    this->sigmas[index] = val;

}


double Individual::getDimension(int index){
    return this->parameters[index];
}

double* Individual::getParameters() {
    return this->parameters;
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
    for(int i=0; i< this->numParameters; i++){
        indString+= doubleToString(this->parameters[i], 15) + ",";
    }

    indString+= doubleToString(this->fitnessTrain, 15);
    return indString;
}

void Individual::setMaxValues(double *maxValues) {
    this->maxValues = maxValues;
}

double *Individual::getMaxValues() const {
    return maxValues;
}

void Individual::setParameters(double *parameters) {

    for(int i=0; i< this->numParameters; i++){
        this->parameters[i] = parameters[i];
    }
}

double *Individual::getFullParameters() const {
    //initialize fullParameters with maxValues then


    return nullptr;
}

void Individual::setRelativeFitnessTrain(double val) {
    this->relativeFitnessTrain = val;
}

double Individual::getRelativeFitnessTrain() {
    return this->relativeFitnessTrain;
}

void Individual::setRelativeFitnessTest(double val) {
    this->relativeFitnessTest = val;
}

double Individual::getRelativeFitnessTest() {
    return this->relativeFitnessTest;
}


