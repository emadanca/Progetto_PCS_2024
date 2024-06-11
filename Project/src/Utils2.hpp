#pragma once

#include "Utils2.hpp"
#include "PolygonalMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath> // Per utilizzare NaN
#include <algorithm>

using namespace std;
using namespace PolygonalLibrary;

namespace PolygonalLibrary{

template<typename T>
inline void VerificaRaddoppio( vector<T>& vec)
{

    //vec è un vettore di lunghezza l con n elementi
    unsigned int l = vec.capacity();
    unsigned int n = vec.size();

    if( n == l)
    {
        vector<T> new_vec;
        new_vec.reserve(2*l);
        for (unsigned int i = 0; i < n ; i++)
        {
            new_vec.push_back(vec[i]);
        }
        vec = new_vec;
    }
}
Point3D intersectionLines(const Point3D& p1, const Point3D& p2, const Point3D& dir1, const Point3D& dir2);

inline bool pointAlreadyExists(Point3D& pto, vector<Point3D>& punti, unsigned int& id)
{    //l'id mi serve quando cerco un punto tra i vertici
    double tol = 1e-9;
    bool exist = false;

    if(punti.size() == 0)
    {
        return exist;
    }

    for (unsigned int i = 0; i < punti.size(); i++){
        Point3D element = punti[i];
        bool uguagl_x = (abs(pto.x- element.x) < tol);
        bool uguagl_y = (abs(pto.y - element.y) < tol);
        bool uguagl_z = (abs(pto.z - element.z) < tol);

        if (uguagl_x && uguagl_y && uguagl_z){exist = true; id = i; return exist;}
    }
    return exist;
}

bool isPointOnBoundary(const Point3D& point, const Fracture& fracture);

void Sottopoligonazione(Fracture& fracture, vector<Trace>& PassingTraces, vector<Trace>& nonPassingTraces, PolygonalMesh& mesh, const double& tol, list<unsigned int>& IDVertices);

PolygonalMesh SubPolygonMesh (const Fracture& fracture, vector<Trace>& PassingTraces, vector<Trace>& nonPassingTraces);

}
