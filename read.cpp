#include "dependencies.h"

vector<double> doubleVectorFromText(string text)
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

    return result;
}

vector<int> intVectorFromText(string text)
{
    vector<int> result;
    string part;

    for (int i = 0; i < text.size(); i++)
    {

        if (text[i] == ' ' || text[i] == '\n' || text[i] == '\r')
        {
            result.push_back(stof(part));
            part = "";
        }
        else
        {
            part += (text[i]);
        }
    }

    return result;
}

// Teste de leitura da instância
void readFile(string path, int numVectors)
{
    ifstream input(path);
    string textAux;
    vector<double> lineElements;
    vector<vector<double>> allElements(numVectors + 1);


    while(!input.eof()){
        getline(input, textAux);
        lineElements = doubleVectorFromText(textAux);

        //allElements.push_back(lineElements);
        for(int i=0; i<lineElements.size(); i++){
            allElements[i].push_back(lineElements[i]);
        }
        //getline(input, textAux);
    }

    input.close();


}

// Teste de leitura da instância
void readFileToVectors(string path, int numVectors, double *vectors[])
{
    ifstream input(path);
    string textAux;
    vector<double> lineElements;
    //vector<vector<double>> allElements(numVectors);


    for(int i=0; !input.eof(); i++){
        getline(input, textAux);
        lineElements = doubleVectorFromText(textAux);

        //allElements.push_back(lineElements);
        for(int j=0; j<numVectors; j++){
            vectors[j][i] = lineElements[j];
        }
       // getline(input, textAux);
    }

    input.close();

}
