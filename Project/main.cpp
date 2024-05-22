#include "Utils.hpp"
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace DFNLibrary;

int main()
{

    // Creazione delle fratture di esempio
    Fracture f1;
    f1.idFracture= 0;
    f1.numVertices = 4;
    Point3D A;
    A.x = 0.0;
    A.y = 0.0;
    A.z = 0.0;

    Point3D B;
    B.x = 1;
    B.y = 0;
    B.z = 0;

    Point3D vettorialeAB = crossProduct(A,B);
    cout << "vettoriale ab : " << vettorialeAB << endl;
    Point3D C;
    C.x = 1;
    C.y = 1;
    C.z = 0;

    Point3D D;
    D.x = 0;
    D.y = 1;
    D.z = 0;

    f1.vertices.push_back(A);
    f1.vertices.push_back(B);
    f1.vertices.push_back(C);
    f1.vertices.push_back(D);

    Fracture f2;
    f2.idFracture= 0;
    f2.numVertices = 4;

    Point3D E;
    E.x = 0.8;
    E.y = 0;
    E.z = -0.1;

    Point3D F;
    F.x = 0.8;
    F.y = 0;
    F.z = 0.29999999;

    Point3D G;
    G.x = 0.8;
    G.y = 1;
    G.z = 0.29999999;

    Point3D H;
    H.x = 0.8;
    H.y = 1;
    H.z = -0.1;

    f2.vertices.push_back(E);
    f2.vertices.push_back(F);
    f2.vertices.push_back(G);
    f2.vertices.push_back(H);

    Point3D normal1 = crossProduct(f1.vertices[0], f1.vertices[1], f1.vertices[2]);
    cout << "normale 1" << normal1 << endl;

    Point3D normal2 = crossProduct(f2.vertices[0], f2.vertices[1], f2.vertices[2]);
    cout << "normale 2" << normal2 << endl;
    double tol = 1.e-16;

    vector<Fracture> fractures;
    vector<Trace> traces;
    array<bool,2> onThePlane;

    // Leggi i dati DFN dal file
    readDFN("DFN/FR10_data.txt", fractures);
    cout << "size fractures " << fractures.size() << endl;

    for (size_t i = 0; i < fractures.size(); i++) {
        Fracture f = fractures[i];
        cout << "Fracture " << i << " vertices:" << endl;
        cout << f.vertices.size() << endl;
        for (const auto& vertex : f.vertices) {
            cout << vertex << endl;
        }

    for (size_t i = 0; i < fractures.size(); ++i) {
        for (size_t j = i +1; j < fractures.size(); ++j) {
            const Fracture& f1 = fractures[i];
            // cout << "VERTICI VERI : " << f1.vertices[0] << ";" << f1.vertices[1] << ";" << f1.vertices[2] << endl;
            const Fracture& f2 = fractures[j];

            // Vettore per memorizzare gli intersectionPoints per questa coppia di fratture
            vector<Point3D> intersectionPoints(4);
            // cout << "ciao" << endl;

            // Trova gli intersectionPoints per questa coppia di fratture
            bool intersectionFound = findIntersection2(f1, f2, intersectionPoints, tol, onThePlane);
            // cout << "ciao" << endl;
            /*
            cout << " intersection point delle fratture " << i << "e" <<j << endl;
            cout << intersectionPoints[0] << endl;
            cout << intersectionPoints[1] << endl;
            cout << intersectionPoints[2] << endl;
            cout << intersectionPoints[3] << endl;*/
        }
    }

    //list<Trace> tracesList; // Creo il vettore che temporaneo di tracce
    vector<Trace> traces; // Creo il vettore definitivo

    // Trova le tracce
    traces = findTrace(fractures, 2.2e-16);

    // Stampa le tracce su file
    printTracesToFile(traces, "traces.txt");
    sortAndDivideTracesByFracture(traces, "tracesSorted.txt");

    return 0;

}
}

