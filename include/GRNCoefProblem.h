//
// Created by patrick on 23/07/23.
//

#ifndef ES_GRNCOEFPROBLEM_H
#define ES_GRNCOEFPROBLEM_H

#include "dependencies.h"

using namespace pagmo;

/// Class that represents the Problem to be solved by Pagmo (as required by the library)
class GRNCoefProblem {
private:
    appContext *evaluationContext;
    ProblemDescription* problemDescription;
    double (*evaluationFunction)(void *, void *);
    int dimension;
    vector_double upperBounds;
    vector_double lowerBounds;

public:
    explicit GRNCoefProblem(appContext *context);
    explicit GRNCoefProblem(appContext *evaluationContext, ProblemDescription *problemDescription);
    explicit GRNCoefProblem(int dim);
    explicit GRNCoefProblem();
    ~GRNCoefProblem(){};
    vector_double fitness(const vector_double &v) const;
    std::pair<vector_double, vector_double> get_bounds() const;
    void setEvaluationFunction(double (*evaluationFunction)(void*, void*));
    void setBounds(int index, double lower, double upper);
};

#endif //ES_GRNCOEFPROBLEM_H
