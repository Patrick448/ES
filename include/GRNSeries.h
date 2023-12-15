//
// Created by patrick on 07/12/23.
//

#ifndef ES_GRNSERIES_H
#define ES_GRNSERIES_H
#include "dependencies.h"

using namespace std;

class GRNSeries {
private:
    double **vectors;
    int numTimeSteps;
    int numColumns;
    bool matrixInitialized;
    double *maxValues;
public:
    GRNSeries();
    GRNSeries(string filepath);
    GRNSeries(GRNSeries &grnSeries, int start, int end);
    ~GRNSeries();
    void loadFromFile(string filename);
    void initializeMatrix(string filepath);

    double **getVectors() const;

    void setVectors(double **vectors);

    int getNumTimeSteps() const;

    int getNumColumns() const;



    bool isMatrixInitialized() const;

    void setMatrixInitialized(bool matrixInitialized);

    double *getMaxValues() const;
    void loadMaxValues();

};


#endif //ES_GRNSERIES_H
