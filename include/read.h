#ifndef READ_HPP
#define READ_HPP

#include "dependencies.h"
void readFile(std::string path,int numVectors);
void readGRNFileToVectors(string path, int numVectors, double *vectors[]);

#endif // READ_HPP