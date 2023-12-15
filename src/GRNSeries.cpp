//
// Created by patrick on 07/12/23.
//
#include "GRNSeries.h"

GRNSeries::GRNSeries() {

}


GRNSeries::GRNSeries(string filepath) {
    loadFromFile(filepath);
    loadMaxValues();
}

GRNSeries::GRNSeries(GRNSeries &grnSeries, int start, int end) {
    int timeSteps = end - start + 1;
    this->numTimeSteps = timeSteps;
    this->numColumns = grnSeries.getNumColumns();

    this->vectors = new double*[grnSeries.getNumColumns()];
    for(int i=0; i<grnSeries.getNumColumns(); i++){
        this->vectors[i] = new double[timeSteps];
    }
    matrixInitialized = true;

    for(int i=0; i<grnSeries.getNumColumns(); i++){
        for(int j=0; j<timeSteps; j++){
            this->vectors[i][j] = grnSeries.getVectors()[i][j+start];
        }
    }

    loadMaxValues();
}

vector<double> extractDoubleRow(string text)
{
    vector<double> result;
    string part;

    for (int i = 0; i < text.size(); i++)
    {
        if (text[i] == ' ' || text[i] == '\n' || text[i] == '\r')
        {
            result.push_back(stod(part));
            part = "";
        }
        else
        {
            part += (text[i]);
        }
    }

    result.push_back(stod(part));
    return result;
}

void GRNSeries::initializeMatrix(string filepath){
    int rowCount = 0;
    int columnCount = 0;
    ifstream input(filepath);
    string textAux;

    if(!input.eof()) {
        getline(input, textAux);
        vector<double> lineElements = extractDoubleRow(textAux);
        columnCount = (int) lineElements.size();
        rowCount++;
    }

    for(int i=1; !input.eof(); i++){
        getline(input, textAux);
        rowCount++;
    }

    input.close();
    this->numTimeSteps = rowCount;
    this->numColumns = columnCount;

    this->vectors = new double*[this->numColumns];
    for(int i=0; i<this->numColumns; i++){
        this->vectors[i] = new double[this->numTimeSteps];
    }
    matrixInitialized = true;

}

void GRNSeries::loadFromFile(string filepath)//(string path, int numVectors, double *vectors[])
{
    initializeMatrix(filepath);
    ifstream input(filepath);
    string textAux;
    vector<double> lineElements;

    //todo: ver se preciso tratar caso de linha vazia
    for(int i=0; !input.eof(); i++){
        getline(input, textAux);
        lineElements = extractDoubleRow(textAux);

        for(int j=0; j<this->numColumns; j++){
            this->vectors[j][i] = lineElements[j];
        }
    }

    input.close();

}

double getMaxValue(double *values, int start, int end)
{
    double maxValue = 0;
    for (int i = start; i <= end; i++)
    {
        if (values[i] > maxValue)
        {
            maxValue = values[i];
        }
    }

    return maxValue;
}

void GRNSeries::loadMaxValues()
{
    int start = 0;
    int end = this->numTimeSteps - 1;
    this->maxValues = new double[this->numColumns];

    for (int i = 1; i < this->numColumns; i++)
    {
        this->maxValues[i - 1] = getMaxValue(this->vectors[i], start, end);
    }
}

double **GRNSeries::getVectors() const {
    return vectors;
}

void GRNSeries::setVectors(double **vectors) {
    GRNSeries::vectors = vectors;
}

int GRNSeries::getNumTimeSteps() const {
    return numTimeSteps;
}


int GRNSeries::getNumColumns() const {
    return numColumns;
}


bool GRNSeries::isMatrixInitialized() const {
    return matrixInitialized;
}

void GRNSeries::setMatrixInitialized(bool matrixInitialized) {
    GRNSeries::matrixInitialized = matrixInitialized;
}



GRNSeries::~GRNSeries() {
    if(!this->matrixInitialized){
        for(int i=0; i<this->numColumns; i++){
            delete[] this->vectors[i];
        }
        delete[] this->vectors;
        delete[] this->maxValues;
    }

}


double *GRNSeries::getMaxValues() const {
    return maxValues;
}



