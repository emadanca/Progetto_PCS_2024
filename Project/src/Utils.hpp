#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <fstream>
#include <string>

using namespace std; // Utilizzo dello spazio dei nomi standard

// Definizione della struttura per rappresentare un punto tridimensionale
struct Point3D
{
    double x;
    double y;
    double z;
};

// Definizione della struttura per rappresentare una frattura
struct Fracture
{
    int id;
    vector<Point3D> vertices; // Vettore dei vertici
};

// Definizione della struttura per rappresentare una traccia
struct Trace
{
    int id;
    int fracture_id1;
    int fracture_id2;
    Point3D start; // Punto di inizio
    Point3D end;   // Punto di fine
};

// Dichiarazione delle funzioni

// Funzione per leggere i dati dal file DFN e creare un vettore fatto di fratture che sono caratterizzare da id e coordinate

void readDFN(const string& filename, vector<Fracture>& fractures); // Funzione per leggere i dati DFN da un file

bool segmentIntersection (const Point3D& p1, const Point3D& p2, const Point3D& q1, const Point3D& q2 );

vector<Trace> findTrace(const vector<Fracture>& fractures);

double calculateTraceLength(const Trace& trace);

void sortTracesByLength(vector<Trace>& traces); // Funzione per ordinare le tracce per lunghezza

void classifyTraces(vector<Trace>& traces, const vector<Fracture>& fractures); // Funzione per classificare le tracce come passanti o non-passanti


void printTracesToFile(const vector<Trace>& traces, const string& filename); // Funzione per stampare le tracce su un file

void printClassifiedTracesToFile(const vector<vector<Trace>>& passanteTraces, const vector<vector<Trace>>& nonPassanteTraces, const string& filename); // Funzione per stampare le tracce classificate su un file




#endif // UTILS_HPP

