#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

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
    int id; // Identificativo della traccia
    int fracture_id1; // Identificativo della prima frattura
    int fracture_id2; // Identificativo della seconda frattura
    bool Tips; // Indica se la traccia è passante (false) o non-passante (true)
    double length; // Lunghezza della traccia
    Point3D start; // Punto di inizio
    Point3D end;   // Punto di fine
};

// Funzione per leggere i dati DFN da un file
void readDFN(const string& filename, vector<Fracture>& fractures)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    fractures.clear();

    string line;
    int numFractures;
    getline(file, line);
    stringstream(line) >> numFractures;

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        Fracture fracture;
        int id, numVertices;
        stringstream ss(line);
        char delimiter;
        ss >> id >> delimiter >> numVertices;

        for (int i = 0; i < numVertices; ++i) {
            double x, y, z;
            getline(file, line);
            stringstream vertex_ss(line);
            vertex_ss >> x >> y >> z;
            fracture.vertices.push_back({x, y, z});
        }

        fracture.id = id;
        fractures.push_back(fracture);
    }

    file.close();
}
//
bool segmentIntersection (const Point3D& p1, const Point3D& p2, const Point3D& q1, const Point3D& q2 )
{
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dqx = q2.x - q1.x;
    double dqy = q2.y - q1.y;
    double d = dx * dqy - dy * dqx;
    //double d = (p2.x - p1.x) * (q2.y - q1.y) - (p2.y - p1.y) * (q2.x - q1.x); // calcolo del determinate di due punti

    if (d == 0) // Segmenti paralleli
        return false;

    //calcolo dei coefficienti che descrivono la posizione relativa di un punto di intersezione lungo i segmenti rispetto ai loro estremi
    //t e u ci dicono a che distanza si trovano i punti di intersezione lungo i rispettivi segmenti.
    double t = ((q1.x - p1.x) * dqy - (q1.y - p1.y) * dqx) / d;
    double u = ((q1.x - p1.x) * dy - (q1.y - p1.y) * dx) / d;

    return (t >= 0 && t <= 1 && u >= 0 && u <= 1);

}

vector<Trace> findTrace(const vector<Fracture>& fractures)
{
    vector<Trace> traces;

    int traceID = 0; // Inizializziamo l'ID della traccia a 0

    for (size_t i = 0; i < fractures.size(); ++i)
    {
        const vector<Point3D>& fracture1 = fractures[i].vertices;
        for (size_t j = i + 1; j < fractures.size(); ++j)
        {
            const vector<Point3D>& fracture2 = fractures[j].vertices;
            for (size_t k = 0; k < fracture1.size() - 1; ++k)
            {
                for (size_t l = 0; l < fracture2.size() - 1; ++l)
                {
                    if (segmentIntersection(fracture1[k], fracture1[k + 1], fracture2[l], fracture2[l + 1]))
                    {
                        Trace trace;
                        trace.id = traceID++; // Assegna l'ID e incrementa traceID per la prossima traccia
                        trace.fracture_id1 = i + 1; // Assegna l'ID della frattura
                        trace.fracture_id2 = j + 1; // Assegna l'ID della frattura
                        trace.start = fracture1[k];
                        trace.end = fracture1[k + 1];
                        traces.push_back(trace);
                    }
                }
            }
        }
    }

    return traces;
}

double calculateTraceLength(const Trace& trace)
{
    double dx = trace.end.x - trace.start.x;
    double dy = trace.end.y - trace.start.y;
    double dz = trace.end.z - trace.start.z;
    return sqrt(dx * dx + dy * dy + dz * dz); // Calcola la lunghezza usando il teorema di Pitagora
}

// Funzione per ordinare le tracce per lunghezza
void sortTracesByLength(vector<Trace>& traces)
{
    sort(traces.begin(), traces.end(), [](const Trace& a, const Trace& b) {
        return a.length > b.length; // Ordina in modo decrescente per lunghezza
    });
}

// Funzione per calcolare e classificare le tracce come passanti o non-passanti
void classifyTraces(vector<Trace>& traces, const vector<Fracture>& fractures)
{
    for (Trace& trace : traces)
    {
        const Fracture& fracture1 = fractures[trace.fracture_id1 - 1];
        const Fracture& fracture2 = fractures[trace.fracture_id2 - 1];

        // Verifica se la traccia è passante o non-passante
        bool passante = false;
        for (size_t i = 0; i < fracture1.vertices.size() - 1; ++i)
        {
            for (size_t j = 0; j < fracture2.vertices.size() - 1; ++j)
            {
                if (segmentIntersection(fracture1.vertices[i], fracture1.vertices[i + 1], fracture2.vertices[j], fracture2.vertices[j + 1]))
                {
                    passante = true;
                    break;
                }
            }
            if (passante)
                break;
        }
        trace.Tips = !passante; // Se non è passante, allora è non-passante

        // Calcola la lunghezza della traccia
        trace.length = calculateTraceLength(trace);
    }
}

// Funzione per scrivere le tracce su un file, includendo la classificazione e l'ordinamento per lunghezza

// Funzione per scrivere le tracce su un file
void printTracesToFile(const vector<Trace>& traces, const string& filename)
{
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "# Number of Traces\n" << traces.size() << "\n";
    for (const Trace& trace : traces) {
        file << trace.id << "; " << trace.fracture_id1 << "; " << trace.fracture_id2 << "; "
             << trace.start.x << "; " << trace.start.y << "; " << trace.start.z << "; "
             << trace.end.x << "; " << trace.end.y << "; " << trace.end.z << "\n";
    }

    file.close();
}
