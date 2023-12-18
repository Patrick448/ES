//
// Created by patrick on 07/12/23.
//

#ifndef ES_GRNSERIES_H
#define ES_GRNSERIES_H
#include <string>

using namespace std;

class GRNSeries {
private:
    double **vectors;
    int numTimeSteps;
    int numColumns;
    bool matrixInitialized;
    double *maxValues;
    double *initialValues;
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

    int getNumVariables() const;

    double getStartTime() const;

    double getEndTime() const;

    double* getInitialValues() const;

    bool isMatrixInitialized() const;

    void setMatrixInitialized(bool matrixInitialized);

    void initializeInitialValues();

    double *getMaxValues() const;
    void loadMaxValues();

};


#endif //ES_GRNSERIES_H
