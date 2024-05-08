#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

// Funzione per leggere i dati da un file
bool leggiDati(const std::string& fileName, double& S, int& n, double*& w, double*& r);

// Funzione per calcolare il tasso di rendimento e il valore finale del portafoglio
void calcolaRendimento(double S, int n, double* w, double* r, double& rateOfReturn, double& V);

// Funzione per esportare i risultati su un file
bool esportaRisultati(const std::string& fileName, double S, int n, double* w, double* r, double rateOfReturn, double V);

#endif
