#include "Utils.hpp"
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

int main()
{
    // Lettura dei dati DFN da un file
    vector<Fracture> fractures;
    readDFN("DFN/FR3_data.txt", fractures);

    // Calcolo delle tracce
    vector<Trace> traces = findTrace(fractures);

    // Classificazione delle tracce come passanti o non-passanti
    classifyTraces(traces, fractures);

    // Ordinamento delle tracce per lunghezza
    sortTracesByLength(traces);

    // Stampa delle tracce su un file
    printTracesToFile(traces, "traces.txt");

    return 0;
}
