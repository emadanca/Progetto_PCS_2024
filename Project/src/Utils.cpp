#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>
#include <map>
#include "Utils.hpp"

using namespace std;
using namespace DFNLibrary;

namespace DFNLibrary {

// Funzione per leggere i dati DFN da un file
void readDFN(const string& filename, vector<Fracture>& fractures)
{

    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    fractures.clear();

    string line;
    getline(file, line); // salta l'header # Number of Fractures

    int numFractures;
    getline(file, line);
    istringstream iss(line);
    iss >> numFractures;
    fractures.reserve(numFractures);
    cout << "il numero è : " << numFractures << endl;

    while (getline(file, line))
    {
        Fracture fracture;
        getline(file, line);
        istringstream iss(line);
        char delimiter;

        iss >> fracture.idFracture >> delimiter >> fracture.numVertices;
        // fracture.vertices.resize(fracture.numVertices);
        //cout << "la size dei vertici è : " << fracture.vertices.size() << endl;
        //cout << "ID: " << fracture.idFracture << ", delimiter:" << delimiter << ", NumVertices: " << fracture.numVertices << endl;

        getline(file, line); // salta # Vertices

        vector<double> x_values, y_values, z_values;
        getline(file, line); //leggiamo la riga contentente le x
        istringstream issx(line);

        for(int j = 0; j < fracture.numVertices; ++j)
        {

            double x;
            char delimiter;
            issx >> x >> delimiter;
            x_values.push_back(x);


        }

        getline(file, line); //leggiamo la riga contentente le y
        istringstream issy(line);

        for(int j = 0; j < fracture.numVertices; ++j)
        {

            double y;
            issy >> y>> delimiter;
            y_values.push_back(y);
        }

        getline(file, line); //leggiamo la riga contentente le z
        istringstream issz(line);

        for(int j = 0; j < fracture.numVertices; ++j)
        {

            double z;
            char delimiter;
            issz >> z >> delimiter;
            z_values.push_back(z);
        }
        unsigned int count_vertices = 0;
        // Creazione e idicizzazione dei punti
        vector<Point3D> vertices;
        for (size_t ii = 0; ii < x_values.size(); ++ii)
        {
            Point3D p;
            p.x = x_values[ii];
            p.y = y_values[ii];
            p.z = z_values[ii];
            // cout << "il punto è: " << p << endl;
            fracture.vertices.push_back(p);
            //unsigned int count_vertices = 0;
            fracture.IndiciVertices.push_back(count_vertices);
            cout << "indici dei vertici: " << fracture.IndiciVertices[ii] << endl;
            count_vertices++;
        }

        // Calcolo della normale per ogni frattura
        Point3D normal = crossProduct(fracture.vertices[0], fracture.vertices[1], fracture.vertices[2]);
        fracture.normal = normal;

        //calcolo del centrometro della frattura
        // Inizializziamo il centro come (0, 0, 0)
        Point3D center = {0.0, 0.0, 0.0};


        // Sommiamo le coordinate di tutti i vertici
        for (const auto& vertex : fracture.vertices)
        {
            center.x += vertex.x;
            center.y += vertex.y;
            center.z += vertex.z;
        }

        // Calcoliamo la media delle coordinate per ottenere il centro
        center.x /= fracture.numVertices;
        center.y /= fracture.numVertices;
        center.z /= fracture.numVertices;

        fracture.center = center;

        fractures.push_back(fracture);
    }

    file.close();
}

Point3D intersectionPlaneLine(const Point3D& coefficienti, const double d, const Point3D& A, const Point3D& B)
{

    // Calcoliamo il prodotto scalare tra i coefficienti del piano e il vettore direzione della retta (B - A)
    double dotProductvaue = dotProduct(coefficienti, B - A);

    // Definiamo una tolleranza per il confronto con zero
    const double epsilon = 1e-9;

    // Se il prodotto scalare è zero (o molto vicino a zero), la retta è parallela al piano e non ci sono intersezioni
    if (abs(dotProductvaue) < epsilon)
    {
        // Ritorniamo un punto che rappresenta l'assenza di intersezione
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }

    // Calcoliamo il parametro t della retta
    double t = -(A.x * coefficienti.x + A.y * coefficienti.y + A.z * coefficienti.z + d) / dotProductvaue;

    // Calcoliamo le coordinate del punto di intersezione
    Point3D intersection;
    intersection.x = A.x + t * (B.x - A.x);
    intersection.y = A.y + t * (B.y - A.y);
    intersection.z = A.z + t * (B.z - A.z);

    return intersection;
}

// Funzione per verificare se due fratture sono molto distanti
bool areFracturesFarApart(const Fracture& fracture1, const Fracture& fracture2)
{
    // Troviamo la distanza tra i centri delle due fratture
    double distanceCenters = distance(fracture1.center, fracture2.center);

    // Inizializziamo la somma dei raggi come zero
    double Radius1 = 0.0;

    // Per ciascuna frattura, troviamo il raggio massimo tra il centro e i vertici
    for (const auto& vertex : fracture1.vertices)
    {
        // Calcoliamo la distanza tra il vertice e il centro della frattura
        double distToCenter = distance(vertex, fracture1.center);

        // Se questa distanza supera il raggio attuale, la impostiamo come nuovo raggio massimo
        if (distToCenter > Radius1)
            Radius1 = distToCenter;
    }

    // Ripetiamo il processo per la seconda frattura

    double Radius2 = 0.0;
    for (const auto& vertex : fracture2.vertices)
    {
        double distToCenter = distance(vertex, fracture2.center);
        if (distToCenter > Radius2)
            Radius2 = distToCenter;
    }

    double sumOfRadii= Radius1 + Radius2; // somma dei raggi

    // Se la distanza tra i centri è maggiore della somma dei raggi, le fratture sono lontane
    if (distanceCenters > sumOfRadii)
        return true;
    else
        return false;
}

bool checkVertices(const Fracture& f, const Point3D& planeCoeff, double planeD, vector<Point3D>& ip, bool& onPlane)
{
    double tol = 1e-6; // Valore di tolleranza
    bool potentialIntersection = true;
    vector<Point3D> ipTemp(2); // inizializzo il vector
    ip = ipTemp;
    onPlane = false;
    unsigned int count = 0; // conta quanti vertici attraversano il piano.
    int firstOtherSide = -1; // memorizza l'indice del primo vertice che si trova dall'altra parte del piano rispetto al primo vertice
    // Inizialiazziamo le variabili con il PRIMO VERTICE che ha indice 0
    bool firstPositive = (dotProduct(f.vertices[0], planeCoeff) + planeD) > tol; // verifica se il primo vertice è sopra (true) o sotto (false) il piano, considerando la tolleranza tol
    //cout << "vertice 0 prod scalare normale piano " << dotProduct(f.vertices[0], planeCoeff) + planeD << endl;
    bool previousPositive = firstPositive; // inizializza previousPositive con il valore di firstPositive.

    for (unsigned int i = 1; i < f.numVertices; ++i)
    {
        double dotProductvalue = dotProduct(f.vertices[i],planeCoeff) + planeD; // determina se il vertice è sopra (true) o sotto (false) il piano
        bool positive = (dotProductvalue > tol); // se è maggiore di tol è true altrimenti è false
        //cout << "positive " << dotProductvalue << endl;

        if (abs(dotProductvalue) < tol) // se il vertice i è sul piano
        {
            // Se il vertice corrente è sul piano (abs(dotProduct) < tol),
            // verifica se anche i vertici adiacenti sono sul piano.
            // Se sì, imposta onPlane a true e memorizza i vertici nel vector ip.
            if (abs(dotProduct(f.vertices[(i + 1) % f.numVertices], planeCoeff) + planeD) < tol)
            {
                onPlane = true;
                ip[0] = f.vertices[i];
                ip[1] = f.vertices[(i + 1) % f.numVertices];
                return true;
            }
            else if (abs(dotProduct(f.vertices[(i - 1) % f.numVertices], (planeCoeff)) + planeD) < tol)
            {
                onPlane = true;
                ip[0] = f.vertices[(i - 1) % f.numVertices];
                ip[1] = f.vertices[i];
                return true;
            }
        }
        // Verifichiamo ora se i vertici passano ATTRAVERSO il piano
        if (positive != previousPositive) // significa che un vertice ha attraversato il piano
        {
            if (count == 0) // siamo al primo cambio
            {
                firstOtherSide = i; // memorizza l'indice
                count++;
            }
        }
        if ((positive != previousPositive) && count!=0)
            potentialIntersection =  true ;

        if ((previousPositive == positive) && count!=0)
            count++;

        previousPositive = positive;
    }

    if (count > 0)
    {
        int lastIndex = (firstOtherSide + count - 1) % f.numVertices; // calcola l'indice dell'ultimo vertice che si trova dall'altra parte del piano
        cout << "LAST INDEX del firstotherside: " << firstOtherSide << " del count : " << count << " e' " << lastIndex << endl;
        Point3D intPoint0 = intersectionPlaneLine(planeCoeff, planeD, f.vertices[firstOtherSide - 1], f.vertices[firstOtherSide]); // viene calcolato il punto di intersezione tra il piano e la linea che collega il vertice precedente all'indice firstOtherSide e il vertice all'indice firstOtherSide
        // cout <<"################"<< intPoint0 << endl;
        ip[0] = (intPoint0);
        Point3D intPoint1 = intersectionPlaneLine(planeCoeff, planeD, f.vertices[lastIndex], f.vertices[(lastIndex + 1) % f.numVertices]); // viene calcolato il punto di intersezione tra il piano e la linea che collega il vertice precedente all'indice firstOtherSide e il vertice all'indice firstOtherSide
        ip[1] = (intPoint1);
        return true;
    }
    else{
        return false;
    }

}

bool findIntersection2(const Fracture& f1, const Fracture& f2, vector<Point3D>& intersectionPoints, double tol, array<bool,2> onThePlane)
{
    //cout << "i vertici sono :" << f1.vertices[0] << " " << f1.vertices[1] << " " << f1.vertices[2] << endl;
    // Point3D normal1 = crossProduct(f1.vertices[0], f1.vertices[1], f1.vertices[2]);
    //cout << "normale 1" << normal1 << endl;
    double d1 = calculateD(f1.normal, f1.vertices[0]);
    //cout << "d1" << d1 << endl;

    // Point3D normal2 = crossProduct(f2.vertices[0], f2.vertices[1], f2.vertices[2]);
    //cout << "normale 2" << normal2 << endl;
    double d2 = calculateD(f2.normal, f2.vertices[0]);
    //cout << "d2" << d2 << endl;

    bool potentialIntersection = false; // Flag per indicare se ci sarà intersezione
    onThePlane[0] = onThePlane[1] = false;

    Point3D crossNormals = crossProduct(f2.normal, f1.normal);
    //cout << crossNormals << endl;
    double crossNormSquared = dotProduct(crossNormals, crossNormals);
    //cout << "crossNormSquared: " << crossNormSquared << endl;

    if (crossNormSquared < tol * tol)
    {
        // i piani sono paralleli quindi non si intersecano
        cout << "Planes are parallel, no intersection." << endl;
        potentialIntersection = false;
        return false;
    }

    vector<Point3D> ip1, ip2, intersectionPointsTemp(4);
    intersectionPoints = intersectionPointsTemp; // Copia i valori dalla variabile temporanea

    if (checkVertices(f2, f1.normal, d1, ip2, onThePlane[1]))
    {
        potentialIntersection = true;
    }

    if (potentialIntersection && checkVertices(f1, f2.normal, d2, ip1, onThePlane[0]))
    {
        cout << "INTERSECTION POINTS " << endl;
        intersectionPoints[0] = ip1[0];
        cout << intersectionPoints[0] << endl;
        intersectionPoints[1] = ip1[1];
        cout << intersectionPoints[1] << endl;
        intersectionPoints[2] = ip2[0];
        cout << intersectionPoints[2] << endl;
        intersectionPoints[3] = ip2[1];
        cout << intersectionPoints[3] << endl;
    }
    return potentialIntersection;
}

bool findInternalPoints(vector<Point3D>& intersectionPoints, double tol, vector<Point3D>& internalPoints, array<bool,2>& traceTips)
{
    // OBBIETTIVO: trovare gli estremi veri della traccia scoorendo su 4 possibili punti
    // Per farlo devo controllare la posizione reciproca dei punti per stabilire i due più interni (i punti stanno sulla stessa riga quindi posta confrontare i prodotti scalari
    // Inoltre a seconda della loro posizione e secondo gli estremi del poligono si può capire che tipo di traccia è (passante o non passante)
    // ATTENZIONE: la traccia è passante o non passante in relazione al poligo. quindi può essere passante per il primo e non passante per l'altro, o viceversa
    // - se i due interni appartengono allo stesso poligono, la traccia è passante per quel poligono
    // - se sono uno di un poligono e uno di un altro è non passante per entrambi
    // - se inoltre i punti interni e esterni coincidono a due a due è passante per entrambi
    // ricordando che in intersectionPoints ci sono prima due punti del primo poligono (0,1), poi due punti del secondo poligono (2,3)
    // escludo casi di punti coincidenti:

    vector<Point3D> internalPointsTemp(2);
    internalPoints = internalPointsTemp;

    bool intersection = true; //vede se c'è effettivamente intersezione

    double tol2 = tol*tol;

    // caso facile: i punti coincidono --> prodotto scalare uguale a zero
    // la traccia è passante per entrambi
    if((dotProduct(intersectionPoints[0]-intersectionPoints[2], intersectionPoints[0]-intersectionPoints[2])<tol2 &&
         dotProduct(intersectionPoints[1]-intersectionPoints[3], intersectionPoints[1]-intersectionPoints[3])<tol2) ||

        (dotProduct(intersectionPoints[1]-intersectionPoints[2], intersectionPoints[1]-intersectionPoints[2])<tol2 &&
         dotProduct(intersectionPoints[0]-intersectionPoints[3], intersectionPoints[0]-intersectionPoints[3])<tol2))
    {
        traceTips={false,false}; // RICORDA : PASSANTE = FALSE NONPASSANTE = TRUE
        internalPoints[0] = intersectionPoints[0];
        internalPoints[1] = intersectionPoints[1];
        //di uno dei due poligoni, va bene sia il primo (0,1) che il secondo (2,3)
    }

    //vedo i casi in cui solo uno coincide
    //caso 2 coincide con 0, i punti 1 e 3 possono essere:
    //-entrami a destra 031 013
    //uno a destra e uno a sinistra 103 301
    else if ((dotProduct(intersectionPoints[0]-intersectionPoints[2], intersectionPoints[0]-intersectionPoints[2]))<tol2)
    {
        // se 0 coincide con due facciamo questi controlli sul prodotti scalari

        if(dotProduct(intersectionPoints[0]-intersectionPoints[1], intersectionPoints[0]-intersectionPoints[3])>tol)
        {
            //vuol dire che lo zero sta a lato : i casi sono 013 031
            if(dotProduct(intersectionPoints[0]-intersectionPoints[1], intersectionPoints[3]-intersectionPoints[1])>tol)
            {
                //031

                traceTips = {true,false};
                internalPoints[0] = intersectionPoints[0];
                internalPoints[1] = intersectionPoints[3];

            }
            else
            {
                //013
                traceTips = {false,true};
                internalPoints[0] = intersectionPoints[0];
                internalPoints[1] = intersectionPoints[1];
            }
        }
        else{// se il prodotto scalare è negativo lo zero sta in mezzo
            //103
            internalPoints[0] = intersectionPoints[0];
            internalPoints[1] = intersectionPoints[2];
            //si vedrà dopo che ha lunghezza nulla
        }
    }

    // caso in cui 1 coincide con 2
    else if (dotProduct(intersectionPoints[1]-intersectionPoints[2], intersectionPoints[1]-intersectionPoints[2])<tol2)
    {
        // se 1 e 2 coincidono verifichiamo i vari prodotti scalri
        if(dotProduct(intersectionPoints[1]-intersectionPoints[0],intersectionPoints[1]-intersectionPoints[3])>tol)
        {
            // se 1 sta a lato i casi sono: 130 103 altrimenti 1 sta in mezzo
            if(dotProduct(intersectionPoints[1]-intersectionPoints[0],intersectionPoints[3]-intersectionPoints[0])>tol)
            {
                //130
                traceTips={true,false};
                internalPoints[0] = intersectionPoints[0];
                internalPoints[1] = intersectionPoints[3];
            }
            else
            {
                //103
                traceTips={false,true};
                internalPoints[0] = intersectionPoints[1];
                internalPoints[1] = intersectionPoints[0];

            }
        }
        else
        {
            //013
            internalPoints[0] = intersectionPoints[1];
            internalPoints[1] = intersectionPoints[2];
            //lunghezza nulla
        }
    }

    // caso in cui 0 coincide con 3
    else if (dotProduct(intersectionPoints[0]-intersectionPoints[3], intersectionPoints[0]-intersectionPoints[3])<tol2)
    {
        // se il  prodotto scalare è nullo 0 coincide con 3, facciamo le seguenti verifiche per capire dove si trovano i punti rispetto a 0
        if(dotProduct(intersectionPoints[0]-intersectionPoints[1], intersectionPoints[0]-intersectionPoints[2])>tol)
        {
            // se è positivo lo zero sta a lato, ci manca capire la posizione di 2 e 1
            if(dotProduct(intersectionPoints[0]-intersectionPoints[1],intersectionPoints[2]-intersectionPoints[1])>tol)
            {
                //021
                traceTips={true,false};
                internalPoints[0] = intersectionPoints[0];
                internalPoints[1] = intersectionPoints[2];
            }
            else{
                //012
                traceTips={false,true};
                internalPoints[0] = intersectionPoints[0];
                internalPoints[1] = intersectionPoints[1];

            }
        }
        else{
            //102
            internalPoints[0] = intersectionPoints[0];
            internalPoints[1] = intersectionPoints[3];
            //lunghezza nulla
        }
    }

    // caso in cui 1 coincide con 3
    else if (dotProduct(intersectionPoints[1]-intersectionPoints[3], intersectionPoints[1]-intersectionPoints[3])<tol2)
    {
        if(dotProduct(intersectionPoints[1]-intersectionPoints[0], intersectionPoints[1]-intersectionPoints[2])>tol)
        {
            if(dotProduct(intersectionPoints[1]-intersectionPoints[0], intersectionPoints[2]-intersectionPoints[0])>tol)
            {
                //120
                traceTips={true,false};
                internalPoints[0] = intersectionPoints[1];
                internalPoints[1] = intersectionPoints[2];
            }
            else
            {
                //102
                traceTips={false,true};
                internalPoints[0] = intersectionPoints[1];
                internalPoints[1] = intersectionPoints[0];

            }
        }
        else{
            //012
            internalPoints[0] = intersectionPoints[1];
            internalPoints[1] = intersectionPoints[3];
            //lunghezza nulla
        }
    }
    // se uno nessuno coincide dobbiamo controllare la posizione reciproca di tutti e 4
    else
    {
        //se non è passante per entrambi, vedo la posizione reciproca con i prodotti scalari:
        //confronto posizione di 2 rispetto a 0
        if(dotProduct(intersectionPoints[1]-intersectionPoints[0], intersectionPoints[2]-intersectionPoints[0])>tol)
        {
            //confronto la posizione 2 rispetto a 1 se il prodotto scalare è positivo si trova a destra di 2
            if(dotProduct(intersectionPoints[0]-intersectionPoints[1], intersectionPoints[2]-intersectionPoints[1])>tol)
            {
                //confronto la posizione 3 rispetto a 2
                if(dotProduct(intersectionPoints[0]-intersectionPoints[2],intersectionPoints[3]-intersectionPoints[2])>tol)
                {
                    //confronto posizione 3 rispetto a 0
                    if(dotProduct(intersectionPoints[1]-intersectionPoints[0],intersectionPoints[3]-intersectionPoints[0])>tol)
                    {
                        //0321
                        internalPoints[0] = intersectionPoints[3];
                        internalPoints[1] = intersectionPoints[2];
                        traceTips = {true, false};
                        // per il secondo è passante
                        // per il primo è non passante
                    }
                    else{
                        //3021
                        internalPoints[0] = intersectionPoints[0];
                        internalPoints[1] = intersectionPoints[2];
                        traceTips = {true, true};
                        // non passante per entrambi

                    }
                }
                else
                {
                    // controllo la posizione 3 rispetto 1
                    if(dotProduct(intersectionPoints[0]-intersectionPoints[1], intersectionPoints[3]-intersectionPoints[1])>tol)
                    {
                        //0231
                        internalPoints[0] = intersectionPoints[3];
                        internalPoints[1] = intersectionPoints[2];
                        traceTips = {true, false};

                    }
                    else{
                        //0213
                        internalPoints[0] = intersectionPoints[1];
                        internalPoints[1] = intersectionPoints[2];
                        // non passante per entrambi
                        traceTips = {true, true};

                    }
                }
            }
            else
            {
                // posizione 3 rispetto a 1
                if(dotProduct(intersectionPoints[0]-intersectionPoints[1],intersectionPoints[3]-intersectionPoints[1])>tol)
                {
                    // posizione 3 rispetto posizione 0
                    if(dotProduct(intersectionPoints[1]-intersectionPoints[0], intersectionPoints[3]-intersectionPoints[0])>tol)
                    {
                        //0312
                        internalPoints[0] = intersectionPoints[3];
                        internalPoints[1] = intersectionPoints[1];
                        traceTips = {true, true};
                    }
                    else
                    {
                        //3012
                        internalPoints[0] = intersectionPoints[0];
                        internalPoints[1] = intersectionPoints[1];
                        traceTips = {false, true};
                    }
                }
                else
                {
                    // controllo la posizione 3 rispetto a 2
                    if(dotProduct(intersectionPoints[0]-intersectionPoints[2],intersectionPoints[3]-intersectionPoints[2])>tol)
                    {
                        //0132
                        internalPoints[0] = intersectionPoints[3];
                        internalPoints[1] = intersectionPoints[1];

                        intersection=false;
                    }
                    else
                    {
                        //0123
                        internalPoints[0] = intersectionPoints[1];
                        internalPoints[1] = intersectionPoints[2];

                        intersection=false;
                    }
                }
            }
        }
        else{
            // controllo la posizione 3 rispetto a 0
            if(dotProduct(intersectionPoints[1]-intersectionPoints[0],intersectionPoints[3]-intersectionPoints[0])>tol)
            {
                // controllo la posizione 3 rispetto a 1
                if(dotProduct(intersectionPoints[0]-intersectionPoints[1], intersectionPoints[3]-intersectionPoints[1])>tol)
                {
                    //2031
                    internalPoints[0] = intersectionPoints[3];
                    internalPoints[1] = intersectionPoints[0];
                    traceTips = {true, true};
                }
                else
                {
                    //2013
                    internalPoints[0] = intersectionPoints[0];
                    internalPoints[1] = intersectionPoints[1];

                    traceTips = {false, true};

                }
            }
            else
            {
                // controllo la posizione 3 rispetto a 2
                if(dotProduct(intersectionPoints[0]-intersectionPoints[2], intersectionPoints[3]-intersectionPoints[2])>tol)
                {
                    //2301
                    internalPoints[0] = intersectionPoints[3];
                    internalPoints[1] = intersectionPoints[0];

                    intersection=false;
                }
                else
                {
                    //3201
                    internalPoints[0] = intersectionPoints[0];
                    internalPoints[1] = intersectionPoints[2];

                    intersection=false;
                }
            }
        }

    }
    return intersection;

}
vector<Trace> findTrace(const vector<Fracture>& fractures, double tol)
{
    //list<Trace> tracesList; // Creo il vettore che temporaneo di tracce
    vector<Trace> traces; // Creo il vettore definitivo

    int traceID = 0; // Inizializziamo l'ID della traccia a 0

    for (size_t i = 0; i < fractures.size(); ++i) //Scorre su tutte le fratture
    {
        //Fracture fracture1 = fractures[i];

        // Un altro ciclo for annidato che scorre le fratture rimanenti,
        //iniziando da quella successiva a quella corrente (i),
        //per evitare di ripetere le stesse coppie di fratture
        for (size_t j = i + 1; j < fractures.size(); ++j)
        {
            Fracture fracture1 = fractures[i];
            Fracture fracture2 = fractures[j % fractures.size()];
            // se le fratture sono molto lontante non è neceessario controllarle
            // passo alla coppia successiva
            bool distant = areFracturesFarApart(fractures[i],fractures[j]);
            if (!distant)
            {
                vector<Point3D> intersectionPoints(4);
                array<bool,2> onThePlane;
                bool intersectionPossibility = findIntersection2(fractures[i],fractures[j],intersectionPoints,tol,onThePlane);
                if(intersectionPossibility)
                {
                    // se c'è intersezione
                    // cerchiamo tra i 4 potenziali chi sono i punti di intersezione
                    vector<Point3D> extremities; //qui salverò i due punti estremi della traccia
                    array<bool,2> tips = {true,true}; //se resta così è non passante per entrambi
                    bool intersection = findInternalPoints(intersectionPoints,tol,extremities,tips);

                    //calcolo la lunghezza ed escludo il caso in cui sia un unico punto a toccare il poligono
                    double len = distance(extremities[0],extremities[1]);

                    if (len > tol && intersection)
                    {
                        Trace trace;
                        trace.idTrace= traceID++; // Assegna l'ID e incrementa traceID per la prossima traccia
                        trace.fracture_id1 = fractures[i].idFracture; // Assegna l'ID della frattura
                        trace.fracture_id2 = fractures[j].idFracture; // Assegna l'ID della frattura
                        trace.length = len;
                        trace.Tips = tips;
                        trace.start = extremities[0];
                        trace.end = extremities[1];

                        traces.push_back(trace);
                    }
                }
            }

        }
    }

    return traces;
}

// Funzioni per scrivere le tracce su un file
void printTracesToFile(const vector<Trace>& traces, const string& filename)
{
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "# Number of Traces\n" << traces.size() << "\n";
    file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<< endl;
    for (const Trace& trace : traces)
    {
        file << trace.idTrace << "; " << trace.fracture_id1 << "; " << trace.fracture_id2 << "; "
             << trace.start.x << "; " << trace.start.y << "; " << trace.start.z << "; "
             << trace.end.x << "; " << trace.end.y << "; " << trace.end.z << "\n";
    }

    file.close();
}
void sortAndDivideTracesByFracture(const vector<Trace>& traces, vector<Fracture>& fractures, const string& filename)
{

    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    map<int, vector<Trace>> fractureTraces;

    // Raggruppa le tracce per frattura
    for (const auto& trace : traces)
    {
        fractureTraces[trace.fracture_id1].push_back(trace);
        fractureTraces[trace.fracture_id2].push_back(trace);
    }

    // Elabora ogni frattura separatamente
    for (const auto& [fractureId, fractureTraceList] : fractureTraces)
    {
        vector<Trace> passingTraces;
        vector<Trace> nonPassingTraces;

        for (const auto& trace : fractureTraceList)
        {
            if ((trace.fracture_id1 == fractureId && trace.Tips[0] == false) ||
                (trace.fracture_id2 == fractureId && trace.Tips[1] == false))
            {
                passingTraces.push_back(trace);
            }
            else
            {
                nonPassingTraces.push_back(trace);
            }
        }

        // Ordina le tracce passanti e non passanti per lunghezza decrescente
        sort(passingTraces.begin(), passingTraces.end(), compareByLengthDesc);
        sort(nonPassingTraces.begin(), nonPassingTraces.end(), compareByLengthDesc);

        // Trova la frattura corrispondente nel vettore fractures e aggiorna i suoi campi
        for (auto& fracture : fractures)
        {
            if (fracture.idFracture == fractureId)
            {
                fracture.PassingTraces = passingTraces;
                fracture.nonPassingTraces = nonPassingTraces;
                break; // Esci dal ciclo una volta trovata la frattura corrispondente
            }
        }

        // Stampa i risultati per questa frattura
        file << "# FractureId; NumTraces\n" << fractureId << "; " << passingTraces.size() + nonPassingTraces.size() << "\n";

        file << "# TraceId; Tips; Length\n";
        for (const auto& trace : passingTraces) {
            file << trace.idTrace << "; " << false << "; " << trace.length << "\n";
        }
        for (const auto& trace : nonPassingTraces) {
            file << trace.idTrace << "; " << true << "; " << trace.length << "\n";
        }
    }

    file.close();
}

}
