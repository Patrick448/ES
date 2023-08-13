//
// Created by patrick on 23/07/23.
//

#ifndef ES_GRNCOEFPROBLEM_H
#define ES_GRNCOEFPROBLEM_H

#include <utility>
#include "pagmo/types.hpp"
#include "appCtx.h"

using namespace pagmo;

/// Class that represents the Problem to be solved by Pagmo (as required by the library)
class GRNCoefProblem {
private:
    appContext *context;
    double (*evaluationFunction)(void *, void *);
public:
    explicit GRNCoefProblem(appContext *context);
    explicit GRNCoefProblem();
    ~GRNCoefProblem(){};
    vector_double fitness(const vector_double &v) const;
    std::pair<vector_double, vector_double> get_bounds() const;
    void setEvaluationFunction(double (*evaluationFunction)(void*, void*));


};

#endif //ES_GRNCOEFPROBLEM_H
