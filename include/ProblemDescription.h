//
// Created by patrick on 05/06/23.
//

#ifndef ES_PROBLEM_DESC_H
#define ES_PROBLEM_DESC_H


 struct ProblemDescription {
     int IND_SIZE;        // Tamanho do indivíduo (quantidade de coeficientes)
     double MIN_K;        // Menor valor que K pode assumir
     double MAX_K;        // Maior valor que K pode assumir
     double MIN_N;        // Menor valor que N pode assumir
     double MAX_N;        // Maior valor que N pode assumir
     double MIN_TAU;      // Menor valor que TAU pode assumir
     double MAX_TAU;      // Maior valor que TAU pode assumir
     double MIN_STRATEGY; // Menor valor que a estratégia pode assumir
     double MAX_STRATEGY; // Maior valor que a estratégia pode assumir
     int TAU_SIZE;
     int N_SIZE;
     int K_SIZE;
     int (*modelFunction) (double t, double *y, double *ydot, void *_data);

 } ;
#endif //ES_PROBLEM_DESC_H
