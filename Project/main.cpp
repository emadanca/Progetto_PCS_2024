#include "Utils.hpp"
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace DFNLibrary;

int main()
{
    /*vector<Fracture> fractures;
    vector<Trace> traces;
    vector<Trace> passanteTraces;
    vector<Trace> nonPassanteTraces;*/

    // Leggi i dati DFN dal file
    //readDFN("DFN/FR82_data.txt", fractures);
    // Creazione delle fratture di esempio
    Fracture f1;
    f1.idFracture= 0;
    f1.numVertices = 4;
    f1.vertices[1] = Point3D(0.0, 0.0, 0.0);
    f1.vertices[2] = Point3D(1.0, 0.0, 0.0);
    f1.vertices[3] = Point3D(1.0, 1.0, 0.0);
    f1.vertices[4] = Point3D(0.0, 1.0, 0.0);

    Fracture f2;
    f2.idFracture= 0;
    f2.numVertices = 4;
    f2.vertices[1] = Point3D(0.80, 0.0, -0.10);
    f2.vertices[2] = Point3D(0.80, 0.0, 0.299999);
    f2.vertices[3] = Point3D(0.8, 1.0, 0.29999);
    f2.vertices[4] = Point3D(0.8, 1, -0.1);

    //Fracture f1 = {0, 4, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}}};
    //Fracture f2 = {1, 4, {{0.80, 0.0, -0.10}, {0.8, 0.0, 0.299999}, {0.8, 1.0, 0.29999}, {0.8, 1, -0.1}}};

    // Dichiarazione degli intersection points e del flag onThePlane
    vector<Point3D> intersectionPoints(4);
    array<bool, 2> onThePlane;

    // Definizione della tolleranza
    double tol = 1e-6;

    // Trova gli intersection points
    bool intersectionFound = findIntersection2(f1, f2, intersectionPoints, tol, onThePlane);

    // Stampa gli intersection points se trovati
    if (intersectionFound)
    {
        cout << "Intersection Points:" << endl;
        for (const auto& point : intersectionPoints)
        {
            cout << point << endl;
        }
    } else
    {
        cout << "No intersection points found." << endl;
    }

    return 0;

    // Trova le tracce
    // traces = findTrace(fractures, 2.2e-16);

    // Classifica le tracce
    //classifyTraces(fractures, traces, passanteTraces, nonPassanteTraces);

    // Ordina le tracce per lunghezza
    //sortTracesByLength(traces);

    // Stampa le tracce su file
    //printTracesToFile(traces, "traces.txt");

}

