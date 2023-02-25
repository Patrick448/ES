//
// Created by patri on 27/01/2023.
//

//#include <random>
//#include "ESAlgorithm.h"
//#include <float.h>
//#include <iostream>
//#include <vector>
//#include <algorithm>
#include "dependencies.h"

int ESAlgorithm::UPPER_OPEN=0;
int ESAlgorithm::UPPER_CLOSED=1;
int ESAlgorithm::LOWER_OPEN=2;
int ESAlgorithm::LOWER_CLOSED=3;
int ESAlgorithm::UPPER=4;
int ESAlgorithm::LOWER=5;
int ESAlgorithm::ONE_PLUS_ONE=6;
int ESAlgorithm::ISOTROPIC=7;
int ESAlgorithm::NON_ISOTROPIC=8;

ESAlgorithm::ESAlgorithm(int numDimensions){
    this->numDimensions = numDimensions;
    this->upperBounds.resize(numDimensions);
    this->lowerBounds.resize(numDimensions);
    this->lowerBoundTypes.resize(numDimensions);
    this->upperBoundTypes.resize(numDimensions);
    //this->numParents = numParents;
    //this->numOffspring = numOffspring;
    this->maxSigma = DBL_MAX;
    this->minSigma = -DBL_MAX;

}
ESAlgorithm::~ESAlgorithm(){
    this->clear();
}

void ESAlgorithm::clear() {
    this->evaluationsCounter=0;
    for(Individual* i: this->population){
        delete i;
    }
    this->population.clear();
}

vector<Individual*> ESAlgorithm::getPopulation(){
    return this->population;
}
void ESAlgorithm::addIndividual(Individual* individual){
    this->population.push_back(individual);
}
double ESAlgorithm::evaluate(Individual* ind){
    double eval =this->evaluationFunction(ind->getDimensions());
    ind->setEvaluation(eval);
    this->evaluationsCounter++;
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

void ESAlgorithm::setEvaluationFunction(double (*evaluationFunction)(double*)){
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

double ESAlgorithm::setSigmaBounds(double min, double max) {
    this->minSigma = min;
    this->maxSigma = max;
}



void ESAlgorithm::validate(Individual* ind){

    for(int i=0; i< this->numDimensions; i++){
        double uBound = this->getBound(i, ESAlgorithm::UPPER);
        double lBound = this->getBound(i, ESAlgorithm::LOWER);
        double value = ind->getDimension(i);
        int uType = this->upperBoundTypes[i];
        int lType = this->lowerBoundTypes[i];

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

double validatedSigma(double actual, double min, double max){
    if(actual < min){
        return min;
    }else if(actual> max){
        return max;
    }

    return actual;
}

void ESAlgorithm::run1Plus1ES(int seed, double initialSigma, double c, int n,  int maxIterations){
    this->clear();
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


    //algorithm iterations
    for(int i=0; i < maxIterations; i++){
        normal_distribution<double> normal(0, pow(sigma, 2.0));

        Individual* newInd = new Individual(this->numDimensions);
        for(int j=0; j<this->numDimensions; j++){
            double variation = normal(re);
            double newDim = ind->getDimension(j) + variation;
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
                sigma = validatedSigma(sigma, minSigma, maxSigma);
            }else if(ps > 1.0/5.0){
                sigma = sigma/c;
                sigma = validatedSigma(sigma, minSigma, maxSigma);
            }
        }

       // cout << "Iter: " + to_string(i) +" | pop: " + to_string(this->population.size())+ "\n";
    }

}

bool compareIndividuals(Individual* a, Individual* b){
    return a->getEvaluation() < b->getEvaluation();
}

/// deletes individuals from start to end (all inclusive)
void deleteIndividuals(vector<Individual*> &vec, int start, int end){
    for(int i=start; i<=end; i++){
        delete vec[i];
    }
    vec.erase(vec.begin()+start, vec.begin()+end+1);
}

void ESAlgorithm::runPopulationalIsotropicES(int seed, double sigmaVariation, int maxIterations, int numParents, int numOffspring){
    this->clear();
    vector<int> successHistory;
    successHistory.reserve(maxIterations);
    default_random_engine re(seed);
    uniform_real_distribution<double> unifSigmaDistribution(this->getMinSigma(),this->getMaxSigma());

    for(int i=0; i<numParents; i++){
        Individual* ind = new Individual(this->numDimensions);
        double newSigma = unifSigmaDistribution(re);
        ind->setGlobalSigma(newSigma);

        for(int j=0; j<this->numDimensions; j++){
            uniform_real_distribution<double> unifDimDistribution(this->getBound(j, ESAlgorithm::LOWER),
                                                   this->getBound(j, ESAlgorithm::UPPER));

            double newDim = unifDimDistribution(re);
            ind->setDimension(j, newDim);

        }
        this->validate(ind);
        this->evaluate(ind);
        //cout << ind->toString() + "\n";

        this->population.push_back(ind);
    }


    normal_distribution<double> normal1(0, pow(sigmaVariation, 2.0));

    //algorithm iterations
    this->population.reserve(this->population.size() + numOffspring);
    for(int i=0; i < maxIterations; i++){

        //mutate parents and generate offspring
        for(int j=0; j<numParents; j++){
            for(int k=0; k< ceil((double)numOffspring/(double)numParents); k++){
                Individual* newInd = new Individual(this->numDimensions);
                double newSigma = population[j]->getGlobalSigma() * exp(normal1(re));
                newSigma = validatedSigma(newSigma, minSigma, maxSigma);
                newInd->setGlobalSigma(newSigma);
                normal_distribution<double> dimMutationDistribution(0, newSigma);

                for(int d=0; d<this->numDimensions; d++){
                    double newDim = population[j]->getDimension(d) + dimMutationDistribution(re);
                    newInd->setDimension(d, newDim);
                }

                this->validate(newInd);
                this->evaluate(newInd);
                this->population.push_back(newInd);

            }
        }

        sort(this->population.begin(), this->population.end(), compareIndividuals);
        deleteIndividuals(this->population, numParents, this->population.size()-1);
    }

}


void ESAlgorithm::runPopulationalNonIsotropicES(int seed, double sigmaVariation, int maxIterations, int numParents, int numOffspring){
    this->clear();
    vector<int> successHistory;
    successHistory.reserve(maxIterations);
    default_random_engine re(seed);
    uniform_real_distribution<double> unifSigmaDistribution(this->getMinSigma(),this->getMaxSigma());

    for(int i=0; i<numParents; i++){
        Individual* ind = new Individual(this->numDimensions);

        for(int j=0; j<this->numDimensions; j++){
            uniform_real_distribution<double> unifDimDistribution(this->getBound(j, ESAlgorithm::LOWER),
                                                                  this->getBound(j, ESAlgorithm::UPPER));

            double newDim = unifDimDistribution(re);
            double newSigma = unifSigmaDistribution(re);

            ind->setSigma(j, newSigma);
            ind->setDimension(j, newDim);

        }
        this->validate(ind);
        this->evaluate(ind);

        this->population.push_back(ind);
    }


    normal_distribution<double> normal1(0, pow(sigmaVariation, 2.0));

    //algorithm iterations
    this->population.reserve(this->population.size() + numOffspring);
    for(int i=0; i < maxIterations; i++){

        //mutate parents and generate offspring
        for(int j=0; j<numParents; j++){
            for(int k=0; k< ceil((double)numOffspring/(double)numParents); k++){
                Individual* newInd = new Individual(this->numDimensions);

                for(int d=0; d<this->numDimensions; d++){
                    double newSigma = population[j]->getSigma(d) * exp(normal1(re)) * exp(normal1(re));
                    newSigma = validatedSigma(newSigma, minSigma, maxSigma);
                    normal_distribution<double> dimMutationDistribution(0, newSigma);

                    double newDim = population[j]->getDimension(d) + dimMutationDistribution(re);

                    newInd->setSigma(d, newSigma);
                    newInd->setDimension(d, newDim);
                }

                this->validate(newInd);
                this->evaluate(newInd);
                this->population.push_back(newInd);

            }
        }

        sort(this->population.begin(), this->population.end(), compareIndividuals);
        deleteIndividuals(this->population, numParents, this->population.size()-1);
    }

}


string ESAlgorithm::populationToCSVString(){
    string popString = "";
    for(int i=0; i < this->population.size(); i++){
        popString += this->population[i]->toCSVString() + "\n";
    }
    return popString;
}
