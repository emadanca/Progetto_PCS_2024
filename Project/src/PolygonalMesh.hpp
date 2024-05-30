#ifndef __POLYGONALMESH_H
#define __POLYGONALMESH_H


#include <iostream>
//#include "Eigen/Eigen"
#include <vector>
#include <map>
#include <list>
#include "Utils.hpp"


using namespace std;
//using namespace Eigen;
using namespace DFNLibrary;

namespace PolygonalLibrary {

/*definire un oggetto di tipo PolygonalMesh per
ogni frattura che memorizzi:
                            - i vertici CELL0D
                            - i lati CELL1D
                            - i poligoni CELL2D
generati dal taglio delle frattura con le sue tracce*/

struct PolygonalMesh
{
    unsigned int NumberCell0D = 0; ///< numero delle celle 0D
    vector<unsigned int> Cell0DId = {}; ///< Contiene gli ID dei vertici
    vector<Point3D> Cell0DCoordinates = {}; ///< Cell0D coordinates in Point3D (X,Y,Z)

    unsigned int NumberCell1D = 0; ///< numero delle celle 1D
    vector<unsigned int> Cell1DId = {}; ///< Cell1D id, size 1 x NumberCell1D contiene gli ID dei LATI
    vector<pair<unsigned int, unsigned int>> Cell1DVertices = {}; ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId) SAREBBERO TESTA E CODA

    unsigned int NumberCell2D = 0; ///< numero delle celle 2D
    vector<unsigned int> Cell2DId = {}; ///< Cell2D id, size 1 x NumberCell2D
    vector<vector<unsigned int>> Cell2DVertices = {}; ///< Cell2D Vertices indices
                                                      ///< CONTIENE GLI ID DEI VERTICI DEL POLIGONO
    vector<vector<unsigned int>> Cell2DEdges = {}; ///< Cell2D Cell1D indices
                                                   ///< CONTIENE GLI ID DEI SEGMENTI DEL POLIGONO
};

}

#endif
