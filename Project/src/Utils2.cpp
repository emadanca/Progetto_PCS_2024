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

// Funzione per calcolare l'intersezione tra due rette in 3D
Point3D intersectionLines(const Point3D& p1, const Point3D& p2, const Point3D& dir1, const Point3D& dir2)
{
    // Calcola il prodotto vettoriale dei due vettori direzione
    Point3D cross = crossProduct(dir1, dir2);

    // Calcola il prodotto scalare del prodotto vettoriale con la differenza tra i punti delle rette
    // double dot = dotProduct(cross, {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z});

    // Calcola il denominatore per il coefficiente t
    double denominator = dotProduct(cross, cross);

    double tol = 1e-9;

    // Controlla se le rette sono parallele (il denominatore è zero)
    if (denominator < tol)
    {
        // Le rette sono parallele, non c'è intersezione
        cerr << "Le rette sono parallele, non c'è intersezione." << endl;
        return {NAN, NAN, NAN};
    }

    Point3D p21 = p2 - p1;
    Point3D crossValue = crossProduct(p21,dir2);
    double t = dotProduct(cross,crossValue) / denominator;

    cout << "Parametro t: " << t << endl;

    // Calcola il coefficiente t per la prima retta
    // double t = dot / denominator;
    // Calcola il punto di intersezione
    // Point3D intersectionPoint = {p1.x + t * dir1.x, p1.y + t * dir1.y, p1.z + t * dir1.z};

    Point3D intersectionPoint = p1 + t*dir1;
    return intersectionPoint;
}

// Funzione per determinare se un punto giace sul bordo di una frattura o è interno al lato
bool isPointOnBoundary(const Point3D& point, const Fracture& fracture)
{
    for (size_t i = 0; i < fracture.vertices.size(); ++i)
    {
        // Prendiamo i vertici consecutivi
        Point3D vertex1 = fracture.vertices[i];
        Point3D vertex2 = fracture.vertices[(i + 1) % fracture.vertices.size()]; // Vertice successivo, considerando il wrap-around

        // Calcoliamo il vettore da vertex1 a vertex2 e da vertex1 a point
        Point3D vector1 = {vertex2.x - vertex1.x, vertex2.y - vertex1.y, vertex2.z - vertex1.z};
        Point3D vector2 = {point.x - vertex1.x, point.y - vertex1.y, point.z - vertex1.z};

        // Calcoliamo il prodotto vettoriale tra i due vettori
        Point3D crossProductvalue = crossProduct(vector1, vector2);
        double epsilon = 2.2e-16;
        // Se la norma due del prodotto vettoriale è inferiore alla soglia epsilon, il punto giace sul lato
        if (dotProduct(crossProductvalue, crossProductvalue) < (epsilon * epsilon))
        {
            // Ora controlliamo se il punto è compreso tra i due vertici
            double dotProductValue = dotProduct(vector1, vector2);
            double lengthSquared = dotProduct(vector1, vector1);
            if (dotProductValue >= 0 && dotProductValue <= lengthSquared) {
                return true; // Il punto è interno al lato
            }
        }
    }
    return false; // Il punto non è su alcun bordo della frattura
}

void Sottopoligonazione(Fracture& fracture, vector<Trace>& PassingTraces, vector<Trace>& nonPassingTraces, PolygonalMesh& mesh, const double& tol, list<unsigned int>& IDVertices)
{
    /*Definiamo ricorsivamente il metodo:
            - prima divido la frattura in due in base alla prima traccia che dobbiamo considerare
            - vado a dividere le tracce restanti in base a quale poligono appartengono
            - se le tracce del poligono sono finite salvo il poligono nella mesh
            - else rifaccio*/

    Trace cut_trace; // Inizializziamo la traccia che taglierà il poligono

    if (!PassingTraces.empty()) // Verifica se ci sono tracce passanti
    {
        cout << "Ci sono tracce passanti" << endl;
        cut_trace = PassingTraces.front(); // Assegna la prima traccia
        PassingTraces.erase(PassingTraces.begin());// Rimuove la traccia e riordina il vettore e la memoria
    }
    else if (!nonPassingTraces.empty())
    {
        cout << "Ci sono tracce non passanti" << endl;
        cut_trace = nonPassingTraces.front();
        nonPassingTraces.erase(nonPassingTraces.begin());
    }
    else
    {
        // Se entrambi i vettori sono vuoti, la funzione termina
        return;
    }

    // Estrazione delle coordinate della traccia POINT3D e calcolo della retta
    Point3D p1 = cut_trace.start;
    Point3D p2 = cut_trace.end;

    // Calcola le differenze tra le coordinate dei punti e ricava la direttrice della traccia
    Point3D dir = p2-p1;

    // Questi vettori saranno utilizzati per memorizzare i nuovi punti di intersezione
    vector<Point3D> newPoints; // --> coordinate
    vector<unsigned int> id_newPoints; // --> indici dei punti
    newPoints.reserve(10);
    id_newPoints.reserve(10);

    for (int i = 0; i < fracture.numVertices; i++)
    {
        // Definisco i segmenti del poligono
        unsigned int id_start = fracture.IndiciVertices[i];  // Ottieniamo l'ID del vertice iniziale del segmento
        unsigned int id_end;  // Inizializzo l'ID del vertice di fine del segmento

        // Determiniamo l'ID del vertice finale
        // Caso particolare : l'ultimo vertice si collega al primo
        if(i==fracture.numVertices-1)  // Questo significa che stiamo esaminando l'ultimo vertice della frattura durante l'iterazione
        {
            id_end = fracture.IndiciVertices[0]; // Viene assegnato l'ID del primo vertice della frattura.
                // Questo perché, quando si arriva all'ultimo vertice della frattura,
                // il vertice successivo a esso è il primo vertice della frattura in ordine ciclico.
        }
        else
        {
            id_end = fracture.IndiciVertices[i+1]; // L'ID del vertice di fine del segmento viene assegnato all'ID del vertice successivo a quello corrente
        }

        Point3D origin = mesh.Cell0DCoordinates[id_start];
        Point3D end = mesh.Cell0DCoordinates[id_end];
        Point3D dir_f = end -origin; // Direzione frattura dati i punti origin e end
        Point3D crossNormals = crossProduct(dir_f, dir);
        double crossNormSquared = dotProduct(crossNormals, crossNormals);

        bool notParallel = !(crossNormSquared < tol*tol);
        if(notParallel)
        {
            Point3D punto = intersectionLines(origin, p1, dir_f, dir);
            cout << "Punto calcolato: " << punto << endl;

            // Verifica se il punto di intersezione si trova sul bordo di una delle fratture
            bool onBoundary = isPointOnBoundary(punto, fracture);

            if (onBoundary)
            {
                unsigned int a; // Contatore degli id
                cout << "Verifico se il punto esiste in newPoints." << endl;
                if(pointAlreadyExists(punto, newPoints, a)) // Vedo se questo punto appena trovato coincide con i punti di intersezione già trovati
                {
                    cout << "Il punto esiste in newPoints." << endl;

                    //vector<Point3D> coordFracture; // contiene le coordinate dei vertici della frattura corrente
                    //coordFracture.reserve(fracture.numVertices);
                    /*
                    for(unsigned int i = 0; i<fracture.numVertices; i++)
                    {
                        cout << "3 Aggiungo vertice alla frattura: " << mesh.Cell0DCoordinates[fracture.IndiciVertices[i]] << endl;
                        //coordFracture.push_back(mesh.Cell0DCoordinates[fracture.IndiciVertices[i]]);
                        fracture.vertices.push_back(mesh.Cell0DCoordinates[fracture.IndiciVertices[i]]);
                    }*/

                    unsigned int id_new_punto = 0;
                    if (!pointAlreadyExists(punto, fracture.vertices, id_new_punto))
                    {
                        cout << "Il punto " << punto << "non esiste in fracture.vertices, aggiungo." << endl;

                        // Caso in cui il prolungamento della traccia cade su un vertice della frattura
                        newPoints.push_back(punto);  // Appendo il punto di intersezione
                        id_newPoints.push_back(id_new_punto); // Appendo anche l'ID nuovo che coincide con l'ID del vertice
                    }
                    else if(!pointAlreadyExists(punto, mesh.Cell0DCoordinates, id_new_punto))
                    {
                        cout << "Il punto" << punto << " non esiste in mesh.Cell0DCoordinates, aggiungo." << endl;
                        // Caso in cui il prolungamento della traccia tagliante tocca una cella 0D già creata
                        newPoints.push_back(punto);  // Appendo il punto di intersezione
                        id_newPoints.push_back(id_new_punto);  // Push dell'id dei nuovi punti

                        // Inserisco l'id del punto nella lista
                        auto pos = find(IDVertices.begin(), IDVertices.end(), id_end); // Cerchiamo la posizione di id_end
                        IDVertices.insert(pos, id_new_punto); // Inseriamo il nuovo punto nella posizione corretta
                    }
                    else
                    {
                        cout << "Punto nuovo, lo aggiungo." << endl;
                        newPoints.push_back(punto); // Metto il punto di intersezione
                        id_new_punto = mesh.NumberCell0D; // Id diventa equivalente al numero di vertici del poligono
                        id_newPoints.push_back(id_new_punto); // Metto anche l'id

                        VerificaRaddoppio(mesh.Cell0DId); // Se c'è bisogno raddoppio la dimensione del vettore
                        mesh.Cell0DId.push_back(id_new_punto);

                        VerificaRaddoppio(mesh.Cell0DCoordinates);
                        mesh.Cell0DCoordinates.push_back(punto); // Aggiungo le coordinate

                        mesh.NumberCell0D += 1; // Aggiungo 1 al numero di vertici

                        auto pos = find(IDVertices.begin(), IDVertices.end(), id_end);
                        IDVertices.insert(pos, id_new_punto);
                    }
                }
            }
        }
    }

    // Se è stato fatto il taglio vado avanti e mi creo i due sottopoligoni
    Fracture fracture1;
    fracture1.idFracture = fracture.idFracture;
    fracture1.IndiciVertices.reserve(10);  // Riservo spazio

    Fracture fracture2;
    fracture2.idFracture = fracture.idFracture;
    fracture2.IndiciVertices.reserve(10);

    // Devo anche sapere quali sono i vertici delle due sottofratture
    // Mi baso sul fatto che so i due punti appena creati
    bool flag = true; // true --> fracture1, false --> fracture2
    for(unsigned int indice : IDVertices )
    {
        if (indice == id_newPoints[0] || indice == id_newPoints[1])// Questa condizione verifica se l'elemento attuale nella lista dei vertici è uguale a uno dei nuovi punti creati
        {
            // Ovviamente i nuovi punti fanno parte di entrambe le fratture
            fracture1.IndiciVertices.push_back(indice);
            // VerificaRaddoppio(fracture1.IndiciVertices);
            fracture1.vertices.push_back(newPoints[0]);
            fracture1.vertices.push_back(newPoints[1]);
            // VerificaRaddoppio(fracture1.vertices);

            fracture2.IndiciVertices.push_back(indice);
            VerificaRaddoppio(fracture2.IndiciVertices);
            fracture2.vertices.push_back(newPoints[0]);
            fracture2.vertices.push_back(newPoints[1]);
            VerificaRaddoppio(fracture.vertices);

            flag = !flag;  // Così che da ora in poi i vertici vengano assegnati alla seconda sotto-frattura
        }
        else
        {
            if (flag)
            {
                fracture1.IndiciVertices.push_back(indice);
                VerificaRaddoppio(fracture1.IndiciVertices);
                fracture1.IndiciVertices.push_back(indice);
            }
            else
            {
                fracture2.IndiciVertices.push_back(indice);
                VerificaRaddoppio(fracture2.IndiciVertices);
                fracture2.IndiciVertices.push_back(indice);
            }
        }
    }

    fracture1.IndiciVertices.shrink_to_fit(); // Tolgo quello che è vuoto
    fracture2.IndiciVertices.shrink_to_fit();
    fracture1.numVertices = fracture1.IndiciVertices.size();
    fracture2.numVertices= fracture2.IndiciVertices.size();


    // Determino le tracce passanti e quelle non per ogni sotto-frattura
    // creiamo i vettori contenenti le tracce passanti e non pasanti                                                                                   // per il futuro: provare a farlo con le list
    vector<Trace> traces1_Passing;
    vector<Trace> traces1_nonPassing;

    vector<Trace> traces2_Passing;
    vector<Trace> traces2_nonPassig;

    // TRACCE PASSANTI


    cout << "PASSING TRACE " << PassingTraces.empty() << endl;
    for (Trace& trace : PassingTraces)  // ciclo sulle tracce
    {

        // OBBIETTIVO: capire se le tracce appartengono a fracture1 o fracture2
        // caso particolare: appartengono a entrambe

        // passo 1: ricavo le coordinate degli estremi della traccia
        Point3D normal = fracture.normal;
        Point3D origin = trace.start;
        Point3D end = trace.end;

        Point3D separetorPlane = crossProduct(normal, p2-p1); // si fa attenzione a prendere il piano della traccia CUT
        // testiamo le tracce: dobbiamo capire se gli estremi della traccia stanno da una parte o dall'altra
        // prendiamo le distanze dal piano
        double d1 = dotProduct(end-p1,normal);
        double d2 = dotProduct(origin-p1,normal);

        int position;
        if (d1 > 0 && d2 > 0)
            position = 1; // Above the plane
        else if (d1 < 0 && d2 < 0)
            position = -1; // Below the plane
        else
            position = 0; // Crossing the plane

        if (dotProduct(separetorPlane, fracture.vertices[0] - cut_trace.start) < 0 ) position *= -1; // per capire se contiamo dal basso o dall'alto
        if (position == 1)
            fracture1.PassingTraces.push_back(trace);
        else if (position == -1)
            fracture2.PassingTraces.push_back(trace);
        else
        {
            fracture1.PassingTraces.push_back(trace);
            fracture2.PassingTraces.push_back(trace);
        }
    }

    // TRACCE NON PASSANTI

    for (Trace& trace : nonPassingTraces)  // ciclo sulle tracce
    {
        Point3D normal = fracture.normal;
        Point3D origin = trace.start;
        Point3D end = trace.end;

        Point3D separetorPlane = crossProduct(normal, p2-p1); // fai attenzione a prendere il piano della traccia CUT
        // testiamo le tracce: dobbiamo capire se gli estremi della traccia stanno da una parte o dall'altra
        // prendiamo le distanze dal piano
        double d1 = dotProduct(end-p1,normal);
        double d2 = dotProduct(origin-p1,normal);

        int position;
        if (d1 > 0 && d2 > 0) {
            position = 1; // Sopra il piano
        } else if (d1 < 0 && d2 < 0) {
            position = -1; //Sotto il piano
        } else {
            position = 0; // attraversa il piano
        }

        if (dotProduct(separetorPlane, fracture.vertices[0] - cut_trace.start) < 0 ) position *= -1; // per capire se contiamo dal basso o dall'alto
        if (position == 1)
            fracture1.nonPassingTraces.push_back(trace);
        else if (position == -1)
            fracture2.nonPassingTraces.push_back(trace);
        else
        {
            fracture1.nonPassingTraces.push_back(trace);
            fracture2.nonPassingTraces.push_back(trace);
        }
    }

    // GESTIONE DELLA MESH

    list<unsigned int> IDVertices1; // Lo popolo con gli indici dei vertici della fracture1
    for (int i = 0; i< fracture1.numVertices; i++)
    {
        IDVertices1.push_back(fracture1.IndiciVertices[i]);
    }

    list<unsigned int> IDVertices2; // Lo popolo con i vertici della fracture2
    for (int i = 0; i< fracture2.numVertices; i++)
    {
        IDVertices2.push_back(fracture2.IndiciVertices[i]);
    }

    // Se non ci sono tracce benissimo, non è necessario fare la sottopoligonazione, salvo e vado avanti facendo la mesh
    if (traces1_Passing.size() + traces1_nonPassing.size() == 0)
    {
        unsigned int num_cell1d = mesh.NumberCell1D;
        mesh.NumberCell1D += fracture1.numVertices;
        unsigned int num_celle2d = mesh.NumberCell2D;
        mesh.NumberCell2D +=1;

        // Creazione di due vettori ausiliari:
        // vengono utilizzati per memorizzare temporaneamente informazioni durante la creazione delle celle della mesh
        vector<unsigned int> vertices1; // memorizza gli identificatori dei vertici della frattura fracture1.
        vertices1.reserve(fracture1.numVertices);
        vector<unsigned int> edges1; // memorizza gli identificatori delle celle 1D corrispondenti ai segmenti della fracture1.
        edges1.reserve(fracture1.numVertices);

        for (int v = 0; v < fracture1.numVertices; v++)
        {
            unsigned int id_origin = fracture1.IndiciVertices[v];
            unsigned int id_end;
            if(v == fracture1.numVertices-1)
            {
                id_end = fracture1.IndiciVertices[0];
            }
            else
            {
                id_end = fracture1.IndiciVertices[v+1];
            }

            // Celle 1D
            unsigned int id_segm = num_cell1d+v;
            VerificaRaddoppio(mesh.Cell1DId);
            mesh.Cell1DId.push_back(id_segm);
            VerificaRaddoppio(mesh.Cell1DVertices);
            mesh.Cell1DVertices.push_back({id_origin, id_end});

            // Celle 2D
            vertices1.push_back(id_origin);
            edges1.push_back(id_segm);
        }

        VerificaRaddoppio(mesh.Cell2DId);
        mesh.Cell2DId.push_back(num_celle2d);
        VerificaRaddoppio(mesh.Cell2DEdges);
        mesh.Cell2DEdges.push_back(edges1);
        VerificaRaddoppio(mesh.Cell2DVertices);
        mesh.Cell2DVertices.push_back(vertices1);
    }
    else
    {
        Sottopoligonazione(fracture1, traces1_Passing, traces1_nonPassing, mesh, tol, IDVertices1);
    }

    // Si fa  la medesima cosa per la frattura2
    if (traces2_Passing.size() + traces2_nonPassig.size() == 0)
    {
        // salvo le celle 1D e 2D
        unsigned int num_cell1d = mesh.NumberCell1D;
        mesh.NumberCell1D += fracture2.numVertices;
        unsigned int num_celle2d = mesh.NumberCell2D;
        mesh.NumberCell2D +=1;

        // mi creo due vettori ausiliari
        vector<unsigned int> vertices2;
        vertices2.reserve(fracture2.numVertices);
        vector<unsigned int> edges2;
        edges2.reserve(fracture2.numVertices);

        for (int v = 0; v < fracture2.numVertices; v++)
        {
            unsigned int id_origin = fracture2.IndiciVertices[v];
            unsigned int id_end;
            if(v == fracture2.numVertices-1)
            {
                id_end = fracture2.IndiciVertices[0];
            }
            else{
                id_end = fracture2.IndiciVertices[v+1];
            }

            // Celle 1D
            unsigned int id_segm = num_cell1d+v;
            VerificaRaddoppio(mesh.Cell1DId);
            mesh.Cell1DId.push_back(id_segm);
            VerificaRaddoppio(mesh.Cell1DVertices);
            mesh.Cell1DVertices.push_back({id_origin, id_end});

            // Celle 2D
            vertices2.push_back(id_origin);
            edges2.push_back(id_segm);
        }

        VerificaRaddoppio(mesh.Cell2DId);
        mesh.Cell2DId.push_back(num_celle2d);
        VerificaRaddoppio(mesh.Cell2DEdges);
        mesh.Cell2DEdges.push_back(edges2);
        VerificaRaddoppio(mesh.Cell2DVertices);
        mesh.Cell2DVertices.push_back(vertices2);
    }
    else
    {
        Sottopoligonazione(fracture2, traces2_Passing, traces2_nonPassig, mesh, tol, IDVertices2);
    }
}

// La funzione SubPoligonMesh serve a creare una mesh poligonale partendo da una
//  frattura e da due serie di tracce (pass e nonPass).
// In pratica, questa funzione divide una frattura  in sottopoligoni più piccoli
// e gestisce l'organizzazione delle celle 0D, 1D e 2D all'interno della struttura di dati PolygonalMesh definita in hpp.

PolygonalMesh SubPolygonMesh (const Fracture& fracture, vector<Trace>& PassingTraces, vector<Trace>& nonPassingTraces)
{
    PolygonalMesh mesh;
    // ridefinisco la frattura in modo da iniziare l'indicizzazione delle celle 0D da 0
    Fracture f;  // definisco un nuovo oggetto frattura vuoto da riempire
    f.vertices = fracture.vertices;
    // vector<Point3D> vertices_f; // conterrà solo i vertici di questa frattura

    f.numVertices = fracture.numVertices; // questo nuovo oggetto avrà lo stesso numero di vertici della frattura che gli passo in input
    f.IndiciVertices.reserve(f.numVertices);  // inizializzo la memoria per ID
    f.IndiciVertices = fracture.IndiciVertices;
    //vertices_f.reserve(fracture.numVertices);  // memoria per le coordinate
    f.idFracture = fracture.idFracture;

    mesh.NumberCell0D = mesh.NumberCell0D + f.numVertices;  // il numero di vertici è lo stesso della frattura f
    mesh.Cell0DId.reserve(20);      // ATTENZIONE: non so quanti ID mi generano i tagli
    mesh.Cell0DCoordinates.reserve(20);   // non so nemmeno quante coordinate nuove ho

    // inizio a salvare le coordinate 0D
    for (int i = 0; i< fracture.numVertices; i++)  // ciclo sul numero di vertici della frattura che passo in input
    {
        Point3D pto = fracture.vertices[i];
        mesh.Cell0DCoordinates.push_back(pto); // aggiungo le coordinate del punto alla mesh
        mesh.Cell0DId.push_back(i);// aggiungo l'ID del punto alla mesh
    }

    if (PassingTraces.size() + nonPassingTraces.size() == 0)
    {
        // salvo le celle 1D
        mesh.NumberCell1D = mesh.NumberCell1D + fracture.numVertices;// numVertices == numEDGES
        mesh.Cell1DId.reserve(fracture.numVertices);  // memoria per i lati

        for (int v = 0; v < fracture.numVertices; v++)// ciclo sul numero di LATI
        {
            unsigned int id_origin = fracture.IndiciVertices[v]; // salvo gli ID dei vertici di origine
            unsigned int id_end; // e di fine
            if(v == fracture.numVertices)
            {
                id_end = fracture.IndiciVertices[0];
            } // l'ultimo vertice è collegato tramite un segmento al primo
            else
            {
                id_end = fracture.IndiciVertices[v+1];
            }
            mesh.Cell1DId.push_back(v);// aggiungo v agli ID dei segmenti
            mesh.Cell1DVertices.push_back({id_origin, id_end}); // aggiungo origine e fine
        }

        //salvo le celle 2D --> in questo caso ne ho solo una, la frattura stessa
        mesh.NumberCell2D = 1;   // perché ne ho una
        mesh.Cell2DId.reserve(1);  // quindi ho una cella di memoria
        mesh.Cell2DId.push_back(0);
        mesh.Cell2DVertices.reserve(1);
        mesh.Cell2DVertices.push_back(fracture.IndiciVertices);  // memorizzo gli ID dei vertici
        mesh.Cell2DEdges.reserve(1);
        mesh.Cell2DEdges.push_back(mesh.Cell1DId);  // metto gli ID dei segmenti
    }
    else
    {
        // Mi creo una lista di vertici per poterli ordinare
        list<unsigned int> vertici_f;
        for (int i = 0; i< f.numVertices; i++)
        {
            vertici_f.push_back(f.IndiciVertices[i]);
        }
        double tol = 1e-9;
        Sottopoligonazione(f, PassingTraces, nonPassingTraces, mesh, tol, vertici_f);
    }
    return mesh;
}
}
