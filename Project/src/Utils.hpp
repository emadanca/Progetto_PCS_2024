#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <fstream>
#include <sstream>


// Funzione per leggere i dati dal file DFN per importare:
// -il numero di fratture
// -Id frattura
// -numero vertici
// -elenco dei vertici (3D)
bool ImportPolygons(const std::string& fileName);


#endif
