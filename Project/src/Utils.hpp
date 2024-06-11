
#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

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

inline ostream& operator << (ostream& os, const Point3D& c)
{
    os <<"{"<< c.x << ";"<< c.y << ";" << c.z << "}" << endl;
    return os;
}
inline Point3D operator*(double scalar, Point3D c)
{
    return Point3D(c.x * scalar, c.y * scalar, c.z * scalar);
}
inline Point3D operator-(const Point3D& p1, const Point3D& p2)
{
    return Point3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

inline Point3D operator+(const Point3D& p1, const Point3D& p2)
{
    return Point3D(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

inline bool operator==(const Point3D& p1, const Point3D& p2)
{
    // Definisci una tolleranza per il confronto
    double tolerance = 1e-6; // Ad esempio, 1e-6 per confronti con doppia precisione (double)

    // Effettua il confronto delle coordinate considerando la tolleranza
    return (abs(p1.x - p2.x) < tolerance) &&
           (abs(p1.y - p2.y) < tolerance) &&
           (abs(p1.z - p2.z) < tolerance);
}
inline Point3D crossProduct(const Point3D& A, const Point3D& B, const Point3D& C)
{
    Point3D result;
    result.x = (B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y);
    result.y = (B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z);
    result.z = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
    // Calcola il termine noto dell'equazione del piano
    return result;
}
inline double dotProduct(const Point3D& u, const Point3D& v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}
// Funzione per calcolare il termine noto dell'equazione del piano
inline double calculateD(const Point3D& normal, const Point3D& pointOnPlane)
{
    return -dotProduct(normal, pointOnPlane);
}
inline Point3D crossProduct(const Point3D& v1, const Point3D& v2)
{
    Point3D result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}
Point3D crossProduct(const Point3D& A, const Point3D& B, const Point3D& C);
double calculateD(const Point3D& normal, const Point3D& pointOnPlane);
Point3D crossProduct(const Point3D& v1, const Point3D& v2);
double dotProduct(const Point3D& u, const Point3D& v);


// Funzione per calcolare la distanza euclidea tra due punti
inline double distance(const Point3D& A, const Point3D& B)
{
    double dx = B.x - A.x;
    double dy = B.y - A.y;
    double dz = B.z - A.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}
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

// Definizione della struttura per rappresentare una frattura
struct Fracture
{
    int idFracture;
    int numVertices;
    Point3D normal; // normale del piano contenente la frattura
    Point3D center; // centrometro della frattura
    vector<Point3D> vertices; // Vettore dei vertici
    vector<unsigned int> IndiciVertices; //vector con gli id dei vertici
    vector<Trace> PassingTraces;
    vector<Trace> nonPassingTraces;
};

void readDFN(const string& filename, vector<Fracture>& fractures);

//#################### BOUNDINGBOX ######################
bool areFracturesFarApart(const Fracture& fracture1, const Fracture& fracture2);

//############### INTERSEZIONI FRATTURE ##############
Point3D intersectionPlaneLine(const Point3D& coefficienti, const double d, const Point3D& A, const Point3D& B);
bool findIntersection2(const Fracture& f1, const Fracture& f2, vector<Point3D>& intersectionPoints, double tol, array<bool,2> onThePlane);
bool findInternalPoints(vector<Point3D>& intersectionPoints, double tol, vector<Point3D>& internalPoints, array<bool,2>& traceTips);
vector<Trace> findTrace(const vector<Fracture>& fractures, double tol);
void printTracesToFile(const vector<Trace>& traces, const string& filename);
// Funzione per confrontare due tracce in base all'ordine di intersezione
inline bool compareByLengthDesc(const Trace& trace1, const Trace& trace2)
{
    // Ordina per il primo elemento della coppia
    return trace1.length < trace2.length;
}

void sortAndDivideTracesByFracture(const vector<Trace>& traces, vector<Fracture>& fracture, const string& filename);
}
#endif // UTILS_HPP

