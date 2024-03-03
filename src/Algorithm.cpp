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
#include "GRNCoefProblem.h"
#include "Algorithm.h"

#include <pagmo/algorithms/cmaes.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/utils/constrained.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/sade.hpp>

int Algorithm::UPPER_OPEN=0;
int Algorithm::UPPER_CLOSED=1;
int Algorithm::LOWER_OPEN=2;
int Algorithm::LOWER_CLOSED=3;
int Algorithm::UPPER=4;
int Algorithm::LOWER=5;
int Algorithm::ONE_PLUS_ONE=6;
int Algorithm::ISOTROPIC=7;
int Algorithm::NON_ISOTROPIC=8;

using namespace pagmo;

Algorithm::Algorithm(int numDimensions){
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

Algorithm::Algorithm(GRNSeries &trainingSeries, GRNSeries &testSeries,
                     double (*evaluationFunction)(void *, void *), int numDimensions, ProblemDescription *description) {
    this->trainingSeries = &trainingSeries;
    this->testSeries = &testSeries;
    this->evaluationFunction = evaluationFunction;
    this->numDimensions = numDimensions;
    this->upperBounds.resize(numDimensions);
    this->lowerBounds.resize(numDimensions);
    this->lowerBoundTypes.resize(numDimensions);
    this->upperBoundTypes.resize(numDimensions);
    this->maxSigma = DBL_MAX;
    this->minSigma = -DBL_MAX;
    this->modelDescription = description;

}

Algorithm::~Algorithm(){
    this->clear();
}


void Algorithm::clear() {
    this->evaluationsCounter=0;

    for(Individual* i: this->population){
        delete i;
    }
    this->population.clear();

}

vector<Individual*> Algorithm::getPopulation(){
    return this->population;
}
void Algorithm::addIndividual(Individual* individual){
    this->population.push_back(individual);
}

double Algorithm::evaluate(Individual* ind){
    appContext trainCtx{.series = this->trainingSeries, .description = this->modelDescription};
    //appContext testCtx{.series = this->testSeries, .description = this->modelDescription};

    double fitnessTrain =this->evaluationFunction(ind->getParameters(), &trainCtx);
   // double testEval =this->evaluationFunction(ind->getParameters(), &testCtx);

    ind->setFitnessTrain(fitnessTrain);
   // ind->setFitnessTest(testEval);
    ind->setRelativeFitnessTrain(fitnessTrain / this->trainingSeries->getNumTimeSteps());
   // ind->setRelativeFitnessTest(testEval);
    this->evaluationsCounter++;
    return fitnessTrain;
}

bool compareIndividuals(Individual* a, Individual* b){
    return a->getFitnessTrain() < b->getFitnessTrain();
}

void Algorithm::reevaluateAllNoCounter(){
    for(Individual* ind: this->population){
        double eval =this->evaluationFunction(ind->getParameters(), this->context);
        ind->setFitnessTrain(eval);
    }
    sort(this->population.begin(), this->population.end(), compareIndividuals);

}


double Algorithm::getReevaluationByIndexNoCounter(int i){
    Individual* ind = this->population[i];
    double eval =this->evaluationFunction(ind->getParameters(), this->context);

    return eval;
}

void Algorithm::setNumDimensions(int val){
    this->numDimensions = val;
}

int Algorithm::getNumDimensions(){
    return this->numDimensions;
}

void Algorithm::setBounds(int index, double lower, double upper, int lowerBoundType, int upperBoundType){
    this->lowerBounds[index] = lower;
    this->upperBounds[index] = upper;

    this->lowerBoundTypes[index] = lowerBoundType;
    this->upperBoundTypes[index] = upperBoundType;
}

void Algorithm::setEvaluationFunction(double (*evaluationFunction)(void*, void*)){
    this->evaluationFunction = evaluationFunction;
}

void Algorithm::createPopulation(int seed, int numIndividuals) {

    default_random_engine re(seed);

    // Getting a random double value


    for(int i=0; i<numIndividuals; i++){
        Individual* ind = new Individual(this->numDimensions);

        for(int j=0; j<this->numDimensions; j++){
            uniform_real_distribution<double> unif(this->getBound(j, Algorithm::LOWER),
                                                   this->getBound(j, Algorithm::UPPER));
            double newDim = unif(re);
            ind->setParameter(j, newDim);
        }
        this->validate(ind);
        this->evaluate(ind);
        this->population.push_back(ind);

    }
}

double Algorithm::getBound(int index, int which){
    if(which == Algorithm::UPPER){
        if(this->upperBoundTypes[index] == Algorithm::UPPER_CLOSED){
            return upperBounds[index];
        }else{
            return upperBounds[index]-DBL_MIN;
        }
    }else{
        if(this->lowerBoundTypes[index] == Algorithm::LOWER_CLOSED){
            return lowerBounds[index];
        }else{
            return lowerBounds[index]+DBL_MIN;
        }
    }

}

void Algorithm::setContext(void *ctx) {
    this->context = ctx;
}

void Algorithm::setAlgorithmType(int type) {
    this->algorithmType=type;
}

double Algorithm::getMinSigma() {
    return this->minSigma;
}

double Algorithm::getMaxSigma() {
    return this->maxSigma;
}

void Algorithm::setSigmaBounds(double min, double max) {
    this->minSigma = min;
    this->maxSigma = max;
}

int Algorithm::getEvaluations() {
    return this->evaluationsCounter;
}

void Algorithm::validate(Individual* ind){

    for(int i=0; i< this->numDimensions; i++){
        double uBound = this->getBound(i, Algorithm::UPPER);
        double lBound = this->getBound(i, Algorithm::LOWER);
        double value = ind->getDimension(i);
        int uType = this->upperBoundTypes[i];
        int lType = this->lowerBoundTypes[i];

        if(uType == UPPER_CLOSED && value > uBound){
            ind->setParameter(i, uBound);
        }else if(uType==UPPER_OPEN && value >= uBound){
            ind->setParameter(i, uBound - DBL_MIN);
        }else if(lType==LOWER_CLOSED && value < lBound){
            ind->setParameter(i, lBound);
        }else if(lType==LOWER_OPEN && value <= lBound){
            ind->setParameter(i, lBound + DBL_MIN);
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

void Algorithm::run1Plus1ES(int seed, double initialSigma, double c, int n, int maxEvals){
    this->clear();
    double sigma = initialSigma;
    this->evaluationsCounter = 0;
    vector<int> successHistory;
    successHistory.reserve(maxEvals);
    default_random_engine re(seed);

    Individual* ind = new Individual(this->numDimensions);

    for(int j=0; j<this->numDimensions; j++){
        uniform_real_distribution<double> unif(this->getBound(j, Algorithm::LOWER),
                                               this->getBound(j, Algorithm::UPPER));
        double newDim = unif(re);
        ind->setParameter(j, newDim);
    }
    this->validate(ind);
    this->evaluate(ind);
    //cout << ind->toString() + "\n";

    this->population.push_back(ind);


    //algorithm iterations
    for(int i=0; i < maxEvals-1; i++){
        normal_distribution<double> normal(0, pow(sigma, 2.0));

        Individual* newInd = new Individual(this->numDimensions);
        for(int j=0; j<this->numDimensions; j++){
            double variation = normal(re);
            double newDim = ind->getDimension(j) + variation;
            newInd->setParameter(j, newDim);
        }
        this->validate(newInd);
        this->evaluate(newInd);
        //cout << newInd->toString() + "\n";

        if(newInd->getFitnessTrain() < ind->getFitnessTrain()){
            ind = newInd;
            successHistory.push_back(1);
            this->population.push_back(newInd);
            this->bestIndividual = newInd;
           // this->fitnessTrainHistory.push_back(newInd->getFitnessTrain());
            this->record(newInd);
        }else{
            successHistory.push_back(0);
            //this->fitnessTrainHistory.push_back(ind->getFitnessTrain());
            this->record(ind);
            delete newInd;
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


/// deletes individuals from start to end (all inclusive)
void deleteIndividuals(vector<Individual*> &vec, int start, int end){

    if(start >= end){
        return;
    }

    for(int i=start; i<=end; i++){
        delete vec[i];
    }
    vec.erase(vec.begin()+start, vec.begin()+end+1);
}


void Algorithm::runPopulationalIsotropicES(int seed, double sigmaVariation, int maxEvals, int numParents, int numOffspring){
    this->clear();
    //vector<int> successHistory;
    //successHistory.reserve(maxIterations);
    default_random_engine re(seed);
    uniform_real_distribution<double> unifSigmaDistribution(this->getMinSigma(),this->getMaxSigma());
    bool stop = false;

    for(int i=0; i<numParents; i++){
        Individual* ind = new Individual(this->numDimensions);
        double newSigma = unifSigmaDistribution(re);
        ind->setGlobalSigma(newSigma);

        for(int j=0; j<this->numDimensions; j++){
            uniform_real_distribution<double> unifDimDistribution(this->getBound(j, Algorithm::LOWER),
                                                   this->getBound(j, Algorithm::UPPER));

            double newDim = unifDimDistribution(re);
            ind->setParameter(j, newDim);

        }
        this->validate(ind);
        this->evaluate(ind);
        //cout << ind->toString() + "\n";

        this->population.push_back(ind);
    }


    normal_distribution<double> normal1(0, pow(sigmaVariation, 2.0));

    //algorithm iterations
    this->population.reserve(this->population.size() + numOffspring);
    while(!stop){
        //mutate parents and generate offspring
        for(int j=0; j<numParents && !stop; j++){
            for(int k=0; k< ceil((double)numOffspring/(double)numParents); k++){
                Individual* newInd = new Individual(this->numDimensions);
                double newSigma = population[j]->getGlobalSigma() * exp(normal1(re));
                newSigma = validatedSigma(newSigma, minSigma, maxSigma);
                newInd->setGlobalSigma(newSigma);
                normal_distribution<double> dimMutationDistribution(0, newSigma);

                for(int d=0; d<this->numDimensions; d++){
                    double newDim = population[j]->getDimension(d) + dimMutationDistribution(re);
                    newInd->setParameter(d, newDim);
                }

                this->validate(newInd);
                this->evaluate(newInd);
                this->population.push_back(newInd);

                if(this->getEvaluations() >= maxEvals){
                    stop = true;
                    break;
                    //return;
                }

            }
        }

        sort(this->population.begin(), this->population.end(), compareIndividuals);
        deleteIndividuals(this->population, numParents, this->population.size()-1);

        //this->fitnessTrainHistory.push_back(this->population[0]->getFitnessTrain());
        this->record(this->population[0]);
    }

    this->bestIndividual = this->population[0];

}

void Algorithm::runPopulationalNonIsotropicES(int seed, double sigmaVariation, int maxEvals, int numParents, int numOffspring){
    this->clear();
    //vector<int> successHistory;
    //successHistory.reserve(maxIterations);
    default_random_engine re(seed);
    uniform_real_distribution<double> unifSigmaDistribution(this->getMinSigma(),this->getMaxSigma());
    bool stop = false;

    for(int i=0; i<numParents; i++){
        Individual* ind = new Individual(this->numDimensions);

        for(int j=0; j<this->numDimensions; j++){
            uniform_real_distribution<double> unifDimDistribution(this->getBound(j, Algorithm::LOWER),
                                                                  this->getBound(j, Algorithm::UPPER));

            double newDim = unifDimDistribution(re);
            double newSigma = unifSigmaDistribution(re);

            ind->setSigma(j, newSigma);
            ind->setParameter(j, newDim);

        }
        this->validate(ind);
        this->evaluate(ind);

        this->population.push_back(ind);
    }


    normal_distribution<double> normal1(0, pow(sigmaVariation, 2.0));

    //algorithm iterations
    this->population.reserve(this->population.size() + numOffspring);
    while(!stop){
        //mutate parents and generate offspring
        for(int j=0; j<numParents && !stop; j++){
            for(int k=0; k< ceil((double)numOffspring/(double)numParents); k++){
                Individual* newInd = new Individual(this->numDimensions);

                for(int d=0; d<this->numDimensions; d++){
                    double newSigma = population[j]->getSigma(d) * exp(normal1(re)) * exp(normal1(re));
                    newSigma = validatedSigma(newSigma, minSigma, maxSigma);
                    normal_distribution<double> dimMutationDistribution(0, newSigma);

                    double newDim = population[j]->getDimension(d) + dimMutationDistribution(re);

                    newInd->setSigma(d, newSigma);
                    newInd->setParameter(d, newDim);
                }

                this->validate(newInd);
                this->evaluate(newInd);
                this->population.push_back(newInd);

                if(this->getEvaluations() >= maxEvals){
                    stop = true;
                    break;
                    //return;
                }

            }
        }

        sort(this->population.begin(), this->population.end(), compareIndividuals);
        deleteIndividuals(this->population, numParents, this->population.size()-1);

        //this->fitnessTrainHistory.push_back(this->population[0]->getFitnessTrain());
        this->record(this->population[0]);
    }

    this->bestIndividual = this->population[0];

}


void Algorithm::runCMAES(int seed, int maxEvals, int populationSize){
    this->clear();
    this->bestIndividual = new Individual(this->numDimensions);
    this->bestIndividual->setFitnessTrain(DBL_MAX);
    int maxGenerations, newSeed;
    srand(seed);
    appContext evaluationContext = appContext{.series = this->trainingSeries, .description = this->modelDescription};
    //vector<double> evalProgression;
    while(this->evaluationsCounter < maxEvals){
        newSeed = rand();
        maxGenerations = (maxEvals - this->evaluationsCounter - populationSize)/populationSize;
        GRNCoefProblem problem = GRNCoefProblem(&evaluationContext);
        problem.setEvaluationFunction(this->evaluationFunction);
        pagmo::population pop = pagmo::population(problem, populationSize, newSeed);
        cmaes alg = cmaes(maxGenerations, -1, -1, -1, -1, 0.5, 1e-6, 1e-6, false, true, newSeed);
        alg.set_verbosity(1);
        pagmo::population newPop = alg.evolve(pop);


        auto logs = alg.get_log();
        for(int i=0; i<logs.size(); i++){
            this->fitnessTrainHistory.push_back(get<2>(logs[i]));
        }

        this->evaluationsCounter += pop.get_problem().get_fevals() + newPop.get_problem().get_fevals();
        //todo: ver o que fazer com esse vazamento de memória
        Individual *newIndividual = new Individual(this->numDimensions, newPop.champion_x().data());
        this->evaluate(newIndividual);

        if(newIndividual->getFitnessTrain() < this->bestIndividual->getFitnessTrain()){
            delete this->bestIndividual;
            this->bestIndividual = newIndividual;
        }

    }
}


void Algorithm::runDE(int seed, int maxEvals, int populationSize){
    this->clear();
    this->bestIndividual = new Individual(this->numDimensions);
    this->bestIndividual->setFitnessTrain(DBL_MAX);
    int maxGenerations, newSeed;
    srand(seed);
    appContext evaluationContext = appContext{.series = this->trainingSeries, .description = this->modelDescription};


    while(this->evaluationsCounter < maxEvals){
        newSeed = rand();
        maxGenerations = (maxEvals - this->evaluationsCounter - populationSize)/populationSize;
        GRNCoefProblem problem = GRNCoefProblem(&evaluationContext);
        problem.setEvaluationFunction(this->evaluationFunction);
        pagmo::population pop = pagmo::population(problem, populationSize, newSeed);
        pagmo::de alg = pagmo::de(maxGenerations, 0.8, 0.9, 2u, 1e-6, 1e-6, newSeed);
        alg.set_verbosity(1);
        pagmo::population newPop = alg.evolve(pop);

        auto logs = alg.get_log();
        for(int i=0; i<logs.size(); i++){
            this->fitnessTrainHistory.push_back(get<2>(logs[i]));
        }


        this->evaluationsCounter += pop.get_problem().get_fevals() + newPop.get_problem().get_fevals();
        //todo: ver o que fazer com esse vazamento de memória
        Individual *newIndividual = new Individual(this->numDimensions, newPop.champion_x().data());
        this->evaluate(newIndividual);

        if(newIndividual->getFitnessTrain() < this->bestIndividual->getFitnessTrain()){
            delete this->bestIndividual;
            this->bestIndividual = newIndividual;
        }
    }
}
void Algorithm::runSADE(int seed, int maxEvals, int populationSize){
    this->clear();
    this->bestIndividual = new Individual(this->numDimensions);
    this->bestIndividual->setFitnessTrain(DBL_MAX);
    int maxGenerations, newSeed;
    srand(seed);
    appContext ctx = appContext{.series = this->trainingSeries, .description = this->modelDescription};


    while(this->evaluationsCounter < maxEvals){
        newSeed = rand();
        maxGenerations = (maxEvals - this->evaluationsCounter - populationSize)/populationSize;
        GRNCoefProblem problem = GRNCoefProblem(&ctx);
        problem.setEvaluationFunction(this->evaluationFunction);
        pagmo::population pop = pagmo::population(problem, populationSize, newSeed);
        // Standard parameters from documentation:
        // sade(unsigned gen = 1u, unsigned variant = 2u, unsigned variant_adptv = 1u, double ftol = 1e-6, double xtol = 1e-6, bool memory = false, unsigned seed = pagmo::random_device::next())
        pagmo::sade alg = pagmo::sade(maxGenerations, 2u, 1u, 1e-6, 1e-6, false, newSeed);
        alg.set_verbosity(1);
        pagmo::population newPop = alg.evolve(pop);

        auto logs = alg.get_log();
        for(int i=0; i<logs.size(); i++){
            this->fitnessTrainHistory.push_back(get<2>(logs[i]));
        }


        this->evaluationsCounter += pop.get_problem().get_fevals() + newPop.get_problem().get_fevals();
        //cout<< "Generations left: " << (maxEvals - this->evaluationsCounter)/populationSize << endl;
        // cout << "Best fitness: " << newPop.champion_f()[0]<< endl;
        // cout << "Evaluations: " << this->evaluationsCounter << endl;
        // cout << "Seed: " << newSeed << endl;

        //todo: ver o que fazer com esse vazamento de memória
        Individual *newIndividual = new Individual(this->numDimensions, newPop.champion_x().data());
        this->evaluate(newIndividual);

        if(newIndividual->getFitnessTrain() < this->bestIndividual->getFitnessTrain()){
            delete this->bestIndividual;
            this->bestIndividual = newIndividual;
        }
    }
}

Individual* Algorithm::getBestIndividual(){
    return this->bestIndividual;
}

string Algorithm::populationToCSVString(){
    string popString = "";
    for(int i=0; i < this->population.size(); i++){
        popString += this->population[i]->toCSVString() + "\n";
    }
    return popString;
}



void Algorithm::reevaluateBestIndividualUsingTestSet() {
    appContext ctx{.series = this->testSeries, .description = this->modelDescription};
    double eval =this->evaluationFunction(this->bestIndividual->getParameters(), &ctx);
    this->bestIndividual->setFitnessTrain(eval);
    //return eval;
}

vector<double> Algorithm::getEvalsProgress() {
    return this->fitnessTrainHistory;
}

void Algorithm::record(Individual *ind) {
    this->fitnessTrainHistory.push_back(ind->getFitnessTrain());
    //this->fitnessTestHistory.push_back(ind->getFitnessTest());
}

