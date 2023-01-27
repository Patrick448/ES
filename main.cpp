//
// Created by patri on 27/01/2023.
//
#include <iostream>
#include "ESAlgorithm.h"
using namespace std;

int main(){
    cout << "Hello";
    ESAlgorithm esAlgorithm = ESAlgorithm(10, 5, 10);

    for(int i=0; i<10; i++){
        esAlgorithm.setBounds(i, 0, 10, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
    }
    esAlgorithm.createPopulation(0, 10);
}
