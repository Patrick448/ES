//
// Created by patri on 27/01/2023.
//

#include <random>
#include "ESAlgorithm.h"
#include <float.h>
#include <iostream>

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
        this->validate(ind);
        this->evaluate(ind);
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

double ESAlgorithm::getMinSigma() {
    return this->minSigma;
}

double ESAlgorithm::getMaxSigma() {
    return this->maxSigma;
}

void ESAlgorithm::validate(Individual* ind){

    for(int i=0; i< this->numDimensions; i++){
        double uBound = this->getBound(i, ESAlgorithm::UPPER);
        double lBound = this->getBound(i, ESAlgorithm::LOWER);
        double value = ind->getDimension(i);
        bool uType = this->upperBoundTypes[i];
        bool lType = this->lowerBoundTypes[i];

        if(uType == UPPER_CLOSED && value > uBound){
            ind->setDimension(i, uBound);
        }else if(uType==UPPER_OPEN && value >= uBound){
            ind->setDimension(i, uBound - DBL_MIN);
        }else if(lType==LOWER_CLOSED && value < lBound){
            ind->setDimension(i, lBound);
        }else if(lType==LOWER_OPEN && value <= lBound){
            ind->setDimension(i, lBound + DBL_MIN);
        }
    }
}

void ESAlgorithm::run1Plus1ES(int seed, double initialSigma, double c, int n,  int maxIterations){
    double sigma = initialSigma;
    this->evaluationsCounter = 0;
    vector<int> successHistory;
    successHistory.reserve(maxIterations);
    default_random_engine re(seed);

    Individual* ind = new Individual(this->numDimensions);

    for(int j=0; j<this->numDimensions; j++){
        uniform_real_distribution<double> unif(this->getBound(j, ESAlgorithm::LOWER),
                                               this->getBound(j, ESAlgorithm::UPPER));
        double newDim = unif(re);
        ind->setDimension(j, newDim);
    }
    this->validate(ind);
    this->evaluate(ind);
    //cout << ind->toString() + "\n";

    this->population.push_back(ind);

    normal_distribution<double> normal(0, 1);

    //algorithm iterations
    for(int i=0; i < maxIterations; i++){
        Individual* newInd = new Individual(this->numDimensions);

        for(int j=0; j<this->numDimensions; j++){
            double newDim = ind->getDimension(j) + normal(re)*sigma;
            newInd->setDimension(j, newDim);
        }
        this->validate(newInd);
        this->evaluate(newInd);
        //cout << newInd->toString() + "\n";

        if(newInd->getEvaluation() < ind->getEvaluation()){
            ind = newInd;
            successHistory.push_back(1);
            this->population.push_back(newInd);

        }else{
            successHistory.push_back(0);
        }

        if((i+1)%n == 0){
            int windowSize = 10*n;
            if(successHistory.size()<10*n){
                windowSize = successHistory.size();
            }
            int successes = std::accumulate(successHistory.begin()+successHistory.size()-windowSize,
                                            successHistory.end(), decltype(successHistory)::value_type(0));
            int failures = windowSize - successes;
            float ps = (float)successes / (float)(successes + failures);

            if(ps < 1.0/5.0){
                sigma = sigma*c;
            }else if(ps > 1.0/5.0){
                sigma = sigma/c;
            }
        }

       // cout << "Iter: " + to_string(i) +" | pop: " + to_string(this->population.size())+ "\n";
    }

}

string ESAlgorithm::populationToString(){
    string popString = "";
    for(int i=0; i < this->population.size(); i++){
        popString += this->population[i]->toString() + "\n";
    }
    return popString;
}
