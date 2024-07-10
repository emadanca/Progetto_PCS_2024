#ifndef TEST_HPP
#define TEST_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "Utils.hpp"
#include "Utils2.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace PolygonalLibrary;
using namespace DFNLibrary;

TEST(Point3DTest, OperatorOutput)
{
    Point3D p(1.0, 2.0, 3.0);
    ostringstream oss;
    oss << p;
    EXPECT_EQ(oss.str(), "{1;2;3}\n");
}

TEST(Point3DTest, OperatorMultiply)
{
    Point3D p1(1.0, 2.0, 3.0);
    double scalar = 2.0;
    Point3D result = scalar * p1;
    Point3D expected(2.0, 4.0, 6.0);
    EXPECT_EQ(result, expected);
}

TEST(Point3DTest, OperatorSubtract)
{
    Point3D p1(1.0, 2.0, 3.0);
    Point3D p2(0.5, 1.0, 1.5);
    Point3D result = p1 - p2;
    Point3D expected(0.5, 1.0, 1.5);
    EXPECT_EQ(result, expected);
}

TEST(Point3DTest, OperatorAdd)
{
    Point3D p1(1.0, 2.0, 3.0);
    Point3D p2(0.5, 1.0, 1.5);
    Point3D result = p1 + p2;
    Point3D expected(1.5, 3.0, 4.5);
    EXPECT_EQ(result, expected);
}

TEST(Point3DTest, DotProduct)
{
    Point3D u(1.0, 2.0, 3.0);
    Point3D v(2.0, 3.0, 4.0);
    double result = dotProduct(u, v);
    double expected = 20.0;
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST(Point3DTest, CrossProductTwoVectors)
{
    Point3D v1(1.0, 2.0, 3.0);
    Point3D v2(2.0, 3.0, 4.0);
    Point3D result = crossProduct(v1, v2);
    Point3D expected(-1.0, 2.0, -1.0); // Calculate the expected result manually or from a trusted source
    EXPECT_EQ(result, expected);
}

TEST(Point3DTest, CalculateD)
{
    Point3D normal(1.0, 1.0, 1.0);
    Point3D pointOnPlane(0.0, 0.0, 0.0);
    double result = calculateD(normal, pointOnPlane);
    double expected = 0.0;
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST(Point3DTest, DistanceCalculation)
{
    Point3D A(1.0, 2.0, 3.0);
    Point3D B(4.0, 5.0, 6.0);
    double result = distance(A, B);
    double expected = sqrt(3.0 * 3.0 + 3.0 * 3.0 + 3.0 * 3.0); // Distanza tra A e B calcolata manualmente
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST(Point3DTest, DistanceSamePoint)
{
    Point3D A(1.0, 2.0, 3.0);
    double result = distance(A, A); // Distanza tra lo stesso punto
    EXPECT_DOUBLE_EQ(result, 0.0);
}

TEST(Point3DTest, DistanceNegativeCoordinates)
{
    Point3D A(-1.0, -2.0, -3.0);
    Point3D B(-4.0, -5.0, -6.0);
    double result = distance(A, B);
    double expected = sqrt(3.0 * 3.0 + 3.0 * 3.0 + 3.0 * 3.0); // Distanza tra A e B con coordinate negative
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST(TraceTest, CompareByLengthDesc)
{
    Trace trace1;
    Trace trace2;
    Trace trace3;
    trace1.length = 3.0;
    trace2.length = 5.0;
    trace3.length = 1.0;

    // Verifica che trace1 sia maggiore di trace2 (in ordine decrescente per length)
    EXPECT_TRUE(compareByLengthDesc(trace1, trace2));

    // Verifica che trace2 sia maggiore di trace3
    EXPECT_FALSE(compareByLengthDesc(trace2, trace3));

    // Verifica che trace1 non sia maggiore di trace3
    EXPECT_FALSE(compareByLengthDesc(trace1, trace3));

    // Verifica che trace3 non sia maggiore di trace2
    EXPECT_TRUE(compareByLengthDesc(trace3, trace2));
}


// Test per la funzione intersectionPlaneLine
TEST(IntersectionTest, PlaneLineIntersection) {
    // Coefficienti del piano: ad esempio, piano x-y con normale (0, 0, 1) e d = 0
    Point3D planeCoefficients = {0.0, 0.0, 1.0};
    double d = 0.0;

    // Punto A sulla linea
    Point3D A = {0.0, 0.0, 0.0};

    // Punto B sulla linea
    Point3D B = {1.0, 1.0, 1.0};

    // Calcola l'intersezione della linea con il piano
    Point3D intersection = intersectionPlaneLine(planeCoefficients, d, A, B);

    // Punto di intersezione atteso
    Point3D expectedIntersection = {0.0, 0.0, 0.0};
    double tolerance = 1e-6;

    // Verifica che il punto di intersezione calcolato sia vicino al punto atteso
    // Verifica che le coordinate del punto di intersezione calcolato siano vicine a quelle del punto atteso
    ASSERT_NEAR(intersection.x, expectedIntersection.x, tolerance);
    ASSERT_NEAR(intersection.y, expectedIntersection.y, tolerance);
    ASSERT_NEAR(intersection.z, expectedIntersection.z, tolerance);
    //ASSERT_TRUE(comparePoints(intersection, expectedIntersection));
}


// Test per la funzione areFracturesFarApart
TEST(FractureTest, FarApartTest)
{
    // Creiamo due fratture con centri e vertici noti
    Point3D center1 = {0.0, 0.0, 0.0};
    std::vector<Point3D> vertices1 = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    Fracture fracture1;
    fracture1.vertices = vertices1;
    fracture1.center = center1;

    Point3D center2 = {10.0, 10.0, 10.0};
    std::vector<Point3D> vertices2 = {{11.0, 10.0, 10.0}, {10.0, 11.0, 10.0}};
    Fracture fracture2;
    fracture2.vertices= vertices2;
    fracture2.center = center2;

    // Chiamiamo la funzione per verificare se le fratture sono lontane
    bool result = areFracturesFarApart(fracture1, fracture2);

    // Verifichiamo che il risultato sia true (le fratture sono lontane)
    ASSERT_TRUE(result);
}


TEST(CheckVerticesTest, VerticesIntersectionTest)
{
    // Creiamo una frattura con vertici noti
    Fracture fracture;
    fracture.numVertices = 4;
    fracture.vertices = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};

    // Definiamo i coefficienti del piano (ad esempio, un piano z = 0)
    Point3D planeCoeff = {0.0, 0.0, 1.0};
    double planeD = 0.0;

    // Variabili per memorizzare i punti di intersezione
    vector<Point3D> intersectionPoints(2);
    bool onPlane = false;

    // Chiamiamo la funzione da testare
    bool result = DFNLibrary::checkVertices(fracture, planeCoeff, planeD, intersectionPoints, onPlane);

    // Verifichiamo che la funzione restituisca true (vertici attraversano il piano)
    ASSERT_TRUE(result);

    // Verifichiamo che onPlane sia false (i vertici non giacciono tutti sul piano)
    ASSERT_FALSE(onPlane);

    // Verifichiamo che i punti di intersezione siano corretti (considerando tolleranza)
    ASSERT_NEAR(intersectionPoints[0].x, 0.0, 1e-6);
    ASSERT_NEAR(intersectionPoints[0].y, 0.0, 1e-6);
    ASSERT_NEAR(intersectionPoints[0].z, 0.0, 1e-6);

    ASSERT_NEAR(intersectionPoints[1].x, 0.0, 1e-6);
    ASSERT_NEAR(intersectionPoints[1].y, 0.0, 1e-6);
    ASSERT_NEAR(intersectionPoints[1].z, 0.0, 1e-6);
}


// Test per verificare l'intersezione tra due fratture
TEST(FindIntersectionTest, IntersectionTest)
{
    // Creiamo due fratture con vertici noti
    Fracture f1, f2;
    f1.numVertices = 4;
    f1.vertices = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};

    f2.numVertices = 4;
    f2.vertices = {{0.5, 0.5, -0.5}, {1.5, 0.5, -0.5}, {1.5, 1.5, -0.5}, {0.5, 1.5, -0.5}};

    // Definiamo la tolleranza
    double tol = 1e-6;

    // Variabili per memorizzare i punti di intersezione
    vector<Point3D> intersectionPoints(4);
    array<bool, 2> onThePlane = {false, false};

    // Chiamiamo la funzione da testare
    bool result = findIntersection2(f1, f2, intersectionPoints, tol, onThePlane);

    // Verifichiamo che la funzione restituisca true (intersezione trovata)
    ASSERT_TRUE(result);

    // Verifichiamo che almeno una delle fratture abbia i vertici sul piano
    ASSERT_TRUE(onThePlane[0] || onThePlane[1]);

    // Verifichiamo che i punti di intersezione siano corretti
    ASSERT_NEAR(intersectionPoints[0].x, 0.0, tol);
    ASSERT_NEAR(intersectionPoints[0].y, 0.0, tol);
    ASSERT_NEAR(intersectionPoints[0].z, 0.0, tol);

    ASSERT_NEAR(intersectionPoints[1].x, 1.0, tol);
    ASSERT_NEAR(intersectionPoints[1].y, 0.0, tol);
    ASSERT_NEAR(intersectionPoints[1].z, 0.0, tol);

    ASSERT_NEAR(intersectionPoints[2].x, 0.5, tol);
    ASSERT_NEAR(intersectionPoints[2].y, 0.5, tol);
    ASSERT_NEAR(intersectionPoints[2].z, -0.5, tol);

    ASSERT_NEAR(intersectionPoints[3].x, 0.5, tol);
    ASSERT_NEAR(intersectionPoints[3].y, 1.0, tol);
    ASSERT_NEAR(intersectionPoints[3].z, -0.5, tol);
}

// Test case for findInternalPoints function
TEST(FindInternalPointsTest, IntersectionAndInternalPointsTest)
{
    // Definizione dei punti di intersezione (esempi di punti)
    vector<Point3D> intersectionPoints = {
        Point3D(0.0, 0.0, 0.0),   // Punto 0
        Point3D(1.0, 0.0, 0.0),   // Punto 1
        Point3D(0.5, 1.0, 0.0),   // Punto 2
        Point3D(1.0, 1.0, 0.0)    // Punto 3
    };

    // Definizione della tolleranza
    double tolerance = 1e-6;

    // Definizione di strutture per i risultati
    vector<Point3D> internalPoints;
    array<bool, 2> traceTips;

    // Chiamata alla funzione
    bool intersection = findInternalPoints(intersectionPoints, tolerance, internalPoints, traceTips);

    // Verifica dei risultati attesi
    ASSERT_TRUE(intersection);  // Ci si aspetta che ci sia un'intersezione valida
    ASSERT_EQ(internalPoints.size(), 2);  // Ci si aspetta esattamente 2 punti interni
    ASSERT_EQ(traceTips.size(), 2);  // Ci si aspetta esattamente 2 flag per i tipi di traccia

    // Verifica specifica dei risultati
    // Esempi di possibili asserzioni in base alla logica della tua funzione
    if (traceTips[0] == false && traceTips[1] == false) {
        // Ci si aspetta che entrambi i poligoni siano passanti
        ASSERT_NEAR(internalPoints[0].x, 0.0, tolerance);
        ASSERT_NEAR(internalPoints[0].y, 0.0, tolerance);
        ASSERT_NEAR(internalPoints[1].x, 1.0, tolerance);
        ASSERT_NEAR(internalPoints[1].y, 0.0, tolerance);
    } else if (traceTips[0] == false && traceTips[1] == true) {
        // Ci si aspetta che il primo poligono sia passante e il secondo non passante
        // Verifica specifica sui punti interni
        ASSERT_NEAR(internalPoints[0].x, 0.0, tolerance);
        ASSERT_NEAR(internalPoints[0].y, 0.0, tolerance);
        ASSERT_NEAR(internalPoints[1].x, 0.5, tolerance);
        ASSERT_NEAR(internalPoints[1].y, 1.0, tolerance);
    }
}

TEST(FRACTURETEST, TestFindTraces) {
    double tol = 10 * std::numeric_limits<double>::epsilon();

    // Definiamo manualmente le fratture
    Fracture f0;
    f0.idFracture = 0;
    f0.vertices = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0}
    };

    Fracture f1;
    f1.idFracture = 1;
    f1.vertices = {
        {0.8, 0.0, -0.1},
        {0.8, 0.0, 0.3},
        {0.8, 1.0, 0.3},
        {0.8, 1.0, -0.1}
    };

    Fracture f2;
    f2.idFracture = 2;
    f2.vertices = {
        {-0.237778, 0.5, -0.34444},
        {0.3161837, 0.5, -0.34444},
        {0.3161837, 0.5, 0.4528389},
        {-0.237778, 0.5, 0.4528389}
    };

    // Inseriamo le fratture in un vettore
    std::vector<Fracture> fractures = {f0, f1, f2};

    // Trova le tracce
    std::vector<Trace> vecTraces = findTrace(fractures, tol);

    std::array<bool, 2> traces_ok;

    // Verifica il numero di tracce trovate
    EXPECT_EQ(vecTraces.size(), 4);

    // Verifica i Tips per ciascuna traccia
    traces_ok = {false, false};
    EXPECT_EQ(vecTraces[0].Tips, traces_ok);

    traces_ok = {true, false};
    EXPECT_EQ(vecTraces[1].Tips, traces_ok);
    EXPECT_EQ(vecTraces[2].Tips, traces_ok);
    EXPECT_EQ(vecTraces[3].Tips, traces_ok);

    // Altri controlli possono essere aggiunti qui per verificare le estremità e altri dettagli se necessario.
}


#endif



