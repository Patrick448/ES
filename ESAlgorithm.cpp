//
// Created by patri on 27/01/2023.
//

#include <random>
#include "ESAlgorithm.h"
#include <float.h>

int ESAlgorithm::UPPER_OPEN=0;
int ESAlgorithm::UPPER_CLOSED=1;
int ESAlgorithm::LOWER_OPEN=2;
int ESAlgorithm::LOWER_CLOSED=3;
int ESAlgorithm::UPPER=4;
int ESAlgorithm::LOWER=5;
int ESAlgorithm::ONE_PLUS_ONE=6;
int ESAlgorithm::ISOTROPIC=7;
int ESAlgorithm::NON_ISOTROPIC=8;

ESAlgorithm::ESAlgorithm(int numDimensions, int numParents, int numOffspring){
    this->numDimensions = numDimensions;
    this->upperBounds.resize(numDimensions);
    this->lowerBounds.resize(numDimensions);
    this->lowerBoundTypes.resize(numDimensions);
    this->upperBoundTypes.resize(numDimensions);
    this->numParents = numParents;
    this->numOffspring = numOffspring;

}
ESAlgorithm::~ESAlgorithm(){
    for(Individual* i: this->population){
        delete i;
    }
}
vector<Individual*> ESAlgorithm::getPopulation(){
    return this->population;
}
void ESAlgorithm::addIndividual(Individual* individual){
    this->population.push_back(individual);
}
double ESAlgorithm::evaluate(Individual* ind){
    double eval =this->evaluationFunction(ind);
    ind->setEvaluation(eval);
    return eval;
}
void ESAlgorithm::setNumDimensions(int val){
    this->numDimensions = val;
}

int ESAlgorithm::getNumDimensions(){
    return this->numDimensions;
}

void ESAlgorithm::setBounds(int index, double lower, double upper, int lowerBoundType, int upperBoundType){
    this->lowerBounds[index] = lower;
    this->upperBounds[index] = upper;

    this->lowerBoundTypes[index] = lowerBoundType;
    this->upperBoundTypes[index] = upperBoundType;
}

void ESAlgorithm::setEvaluationFunction(double (*evaluationFunction)(Individual*)){
    this->evaluationFunction = evaluationFunction;
}

void ESAlgorithm::createPopulation(int seed, int numIndividuals) {

    default_random_engine re(seed);

    // Getting a random double value


    for(int i=0; i<numIndividuals; i++){
        Individual* ind = new Individual(this->numDimensions);

        for(int j=0; j<this->numDimensions; j++){
            uniform_real_distribution<double> unif(this->getBound(j, ESAlgorithm::LOWER),
                                                   this->getBound(j, ESAlgorithm::UPPER));
            double newDim = unif(re);
            ind->setDimension(j, newDim);
        }

        this->population.push_back(ind);

    }
}

double ESAlgorithm::getBound(int index, int which){
    if(which==ESAlgorithm::UPPER){
        if(this->upperBoundTypes[index]==ESAlgorithm::UPPER_CLOSED){
            return upperBounds[index];
        }else{
            return upperBounds[index]-DBL_MIN;
        }
    }else{
        if(this->lowerBoundTypes[index]==ESAlgorithm::LOWER_CLOSED){
            return lowerBounds[index];
        }else{
            return lowerBounds[index]+DBL_MIN;
        }
    }

}

void ESAlgorithm::setAlgorithmType(int type) {
    this->algorithmType=type;
}