#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <fstream>
#include <string>
#include <array>

using namespace std; // Utilizzo dello spazio dei nomi standarD

namespace DFNLibrary {

// Definizione della struttura per rappresentare un punto tridimensionale
struct Point3D
{
    double x;
    double y;
    double z;

    // Default constructor
    Point3D() : x(0.0), y(0.0), z(0.0) {}

    Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

};

ostream& operator << (ostream& os, const Point3D& c);
Point3D operator-(const Point3D& p1, const Point3D& p2);
Point3D operator+(const Point3D& p1, const Point3D& p2);
bool operator==(const Point3D& p1, const Point3D& p2);
Point3D crossProduct(const Point3D& A, const Point3D& B, const Point3D& C);
double calculateD(const Point3D& normal, const Point3D& pointOnPlane);
Point3D crossProduct(const Point3D& v1, const Point3D& v2);
double distance(const Point3D& A, const Point3D& B);
double dotProduct(const Point3D& u, const Point3D& v);

// Definizione della struttura per rappresentare una frattura
struct Fracture
{
    int idFracture;
    int numVertices;
    Point3D normal; // normale del piano contenente la frattura
    Point3D center; // centrometro della frattura
    vector<Point3D> vertices; // Vettore dei vertici
};

// Definizione della struttura per rappresentare una traccia
struct Trace
{
    int idTrace;
    int fracture_id1;
    int fracture_id2;
    array<bool,2> Tips; // Indica se la traccia è passante (false) o non-passante (true). La dimensione è due perchè si riferisce a due poligoni
    double length; // Lunghezza della traccia
    Point3D start; // Punto di inizio
    Point3D end;   // Punto di fine
    vector<Point3D> extremes; // Vettore di punti di intersezione
};
// Dichiarazione delle funzioni
void readDFN(const string& filename, vector<Fracture>& fractures);

//#################### BOUNDINGBOX ######################
bool areFracturesFarApart(const Fracture& fracture1, const Fracture& fracture2);

//############### INTERSEZIONI FRATTURE ##############
Point3D intersectionPlaneLine(const Point3D& coefficienti, const double d, const Point3D& A, const Point3D& B);
bool findIntersection2(const Fracture& f1, const Fracture& f2, vector<Point3D>& intersectionPoints, double tol, array<bool,2> onThePlane);
bool findInternalPoints(vector<Point3D>& intersectionPoints, double tol, vector<Point3D>& internalPoints, array<bool,2>& traceTips);
vector<Trace> findTrace(const vector<Fracture>& fractures, double tol);
void printTracesToFile(const vector<Trace>& traces, const string& filename);
void sortAndDivideTracesByFracture(const vector<Trace>& traces, const string& filename);
}
#endif // UTILS_HPP
