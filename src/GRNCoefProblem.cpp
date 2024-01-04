//
// Created by patrick on 23/07/23.
//

#include "GRNCoefProblem.h"


GRNCoefProblem::GRNCoefProblem(appContext *context) {
    this->evaluationContext = context;
    this->problemDescription = context->description;
}

GRNCoefProblem::GRNCoefProblem(appContext *evaluationContext, ProblemDescription *problemDescription) {
    this->evaluationContext = evaluationContext;
    this->problemDescription = problemDescription;
}

GRNCoefProblem::GRNCoefProblem() {

}


vector_double GRNCoefProblem::fitness(const vector_double &v) const {
    double val = evaluationFunction((void*)v.data(), (void*)evaluationContext);
    return {val};
}

std::pair<vector_double, vector_double> GRNCoefProblem::get_bounds() const {
    int cont = 0;
    auto bounds = std::make_pair(vector_double(problemDescription->IND_SIZE), vector_double(problemDescription->IND_SIZE));

    for (int i = 0; i < problemDescription->TAU_SIZE; i++)
    {
        bounds.first[i] = this->problemDescription->MIN_TAU;
        bounds.second[i] = this->problemDescription->MAX_TAU;
        cont = i;
    }

    for (int i = cont + 1; i < problemDescription->TAU_SIZE + problemDescription->K_SIZE; i++)
    {
        bounds.first[i] = this->problemDescription->MIN_K;
        bounds.second[i] = this->problemDescription->MAX_K;
        cont = i;
    }

    for (int i = cont + 1; i < problemDescription->TAU_SIZE + problemDescription->K_SIZE + problemDescription->N_SIZE; i++)
    {
        bounds.first[i] = this->problemDescription->MIN_N;
        bounds.second[i] = this->problemDescription->MAX_N;
        cont = i;
    }

    return bounds;
}


void GRNCoefProblem::setEvaluationFunction(double (*evaluationFunction)(void*, void*)){
    this->evaluationFunction = evaluationFunction;
}

GRNCoefProblem::GRNCoefProblem(int dim) {
    this->dimension = dim;
    this->upperBounds = vector_double(dim);
    this->lowerBounds = vector_double(dim);
}

void GRNCoefProblem::setBounds(int index, double lower, double upper) {
    this->upperBounds[index] = upper;
    this->lowerBounds[index] = lower;
}

