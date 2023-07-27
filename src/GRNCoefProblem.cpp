//
// Created by patrick on 23/07/23.
//

#include "GRNCoefProblem.h"

/*
GRNCoefProblem::GRNCoefProblem(appContext *context) {
    this->context = context;
}

GRNCoefProblem::GRNCoefProblem() {

}


vector_double GRNCoefProblem::fitness(const vector_double &v) const {
    double val = evaluationFunction((void*)v.data(), (void*)context);
    return {val};
}

std::pair<vector_double, vector_double> GRNCoefProblem::get_bounds() const {
    int cont = 0;
    auto bounds = std::make_pair(vector_double(context->IND_SIZE), vector_double(context->IND_SIZE));

    for (int i = 0; i < context->TAU_SIZE; i++)
    {
        bounds.first[i] = this->context->MIN_TAU;
        bounds.second[i] = this->context->MAX_TAU;
        cont = i;
    }

    for (int i = cont + 1; i < context->TAU_SIZE + context->K_SIZE; i++)
    {
        bounds.first[i] = this->context->MIN_K;
        bounds.second[i] = this->context->MAX_K;
        cont = i;
    }

    for (int i = cont + 1; i < context->TAU_SIZE + context->K_SIZE + context->N_SIZE; i++)
    {
        bounds.first[i] = this->context->MIN_N;
        bounds.second[i] = this->context->MAX_N;
        cont = i;
    }

    return bounds;
}


void GRNCoefProblem::setEvaluationFunction(double (*evaluationFunction)(void*, void*)){
    this->evaluationFunction = evaluationFunction;
}

*/