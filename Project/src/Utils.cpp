#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

bool ImportPolygons (const string& fileName, unsigned int const & numF )
{
    ifstream file(fileName);

    if (file.fail())
    {
        std::cerr << "Errore nell'apertura del file " << std::endl;
        return false;
    }

    string line;
    getline(file, line); //leggiamo #Number of Fractures e la saltiamo
    getline(file, line); //leggiamo il numero di fratture numF



    file.close();
    return true;
}


