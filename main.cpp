//
// Created by patri on 27/01/2023.
//
#include <iostream>
#include <fstream>
#include <math.h>
#include "ESAlgorithm.h"
using namespace std;

void outputToFile(string path, string text, bool append){
    ofstream outputf;

    if(append){
        outputf.open(path, std::ios_base::app);
    }else{
        outputf.open(path);
    }

    outputf << text;
    outputf.close();
}

double testFunc(Individual* ind){
    double eval =0;
    double x = ind->getDimension(0);
    double y = ind->getDimension(1);
    eval = pow(x + 2*y -7, 2) + pow(2*x + y - 5, 2);

    return eval;
}

int main(){
    cout << "\nHello\n";
    int numDim = 2;
    ESAlgorithm esAlgorithm = ESAlgorithm(numDim, 5, 10);

    for(int i=0; i<numDim; i++){
       esAlgorithm.setBounds(i, -10, 10, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
    }
   // esAlgorithm.createPopulation(0, 10);

   esAlgorithm.setEvaluationFunction(testFunc);
   esAlgorithm.run1Plus1ES(1, 1.0, 0.817, 10, 250);
   cout << esAlgorithm.populationToCSVString() + "\n";

    esAlgorithm.setSigmaBounds(0.1, 10);
    esAlgorithm.runPopulationalIsotropicES(1, 0.5, 100, 5, 10);
    cout << esAlgorithm.populationToCSVString() + "\n";

    esAlgorithm.setSigmaBounds(0.1, 10);
    esAlgorithm.runPopulationalNonIsotropicES(1, 0.5, 100, 5, 10);
    cout << esAlgorithm.populationToCSVString() + "\n";

    cout << "\nFinished\n";
}
