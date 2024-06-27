#ifndef __FRACTURES_H
#define __FRACTURES_H

#include "Eigen/Eigen"
#include <iostream>
#include <utility>
#include <map>
#include <vector>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
struct Fracture{
    unsigned int numVertici;
    Matrix3Xd Vertici;
    vector<unsigned int> traccePassanti;
    vector<unsigned int> tracceNonPassanti;

    Fracture() = default; // costruttore di default

    Fracture(unsigned int& numVertici,
              Matrix3Xd& Vertici
              ):
        numVertici(numVertici),
        Vertici(Vertici)
        {} // prende in input le coordinate dei vertici e la lista dei vertici e inizializza i membri corrispondenti VerticesCoordinates e ListVertices con i valori passati come argomenti

    Vector3d Baricentro(Matrix3Xd &Poligono); //metodo in quanto proprietà della frattura
    double Raggio(Vector3d &centro, Matrix3Xd &Poligono);
    Vector4d TrovaPiano(Matrix3Xd &Poligono);
};

inline double DistanzaEuclidea(Vector3d &centro1, Vector3d &centro2) {
    double distanza = 0;
    for (unsigned i = 0; i < 3; i++) {
        distanza += pow(centro1[i] - centro2[i], 2);
    }
    return distanza;
}

bool IntersezioneSfere(Fracture& polygon1, Fracture& polygon2);


struct Trace{
    unsigned int id;
    pair<Vector3d, Vector3d> Vertices;
    double length;
    unsigned int id1;
    unsigned int id2;

    Trace() = default;
};

/*struct Compare {
    bool operator()(const double& a, const double& b) const {
        return a > b; // Ordine decrescente
    }
};*/ // struct utilizzata dalla multimappa per ordinare gli elementi in ordine decrescente



Matrix<double,2,3> IntersezionePiani(Fracture &polygon1, Fracture& polygon2);


inline Vector2d ParametriRetta (const Vector3d& P0, const Vector3d& P1, const Vector3d& Q, const Vector3d& dir_retta){
    Vector2d solution;
    MatrixXd A(3,2);
    A.col(0) = (P1 - P0); // Colonna per il parametro t
    A.col(1) = (-dir_retta); // Colonna per il parametro s
    Vector3d b = (Q - P0);
    if (((P1-P0).cross(dir_retta)).norm() != 0){ //controllo parallelismo
        solution = A.colPivHouseholderQr().solve(b);
    }
    // Risoluzione del sistema lineare
    else {
        solution(0) = -1;
        solution(1) = -1; //assegnamo soluzioni che scarta
    }

    return solution;
}

vector<Vector3d> Intersection_Point(Matrix<double,2,3> &retta, Matrix3Xd &vertici, const unsigned int& numVert);

inline bool isLess(Vector3d p1, Vector3d p0, Matrix<double,2,3> retta_inters) {
    return (p1[0] - p0[0]) * retta_inters(1,0) + (p1[1] - p0[1]) * retta_inters(1,1) + (p1[2] - p0[2]) * retta_inters(1,2) < 0;
}

pair<Vector3d, Vector3d> Traccia(vector<Vector3d> &intersezioni1, vector<Vector3d> &intersezioni2, Matrix<double,2,3> retta_inters);

bool Tips (vector<Vector3d>& intersezioni, pair<Vector3d,Vector3d>& verticiTraccia);

bool Find_Trace(Trace& trace, unsigned int& idT, Fracture& polygon1, Fracture& polygon2);

void ImportFracturesList(const string& filepath, Fracture& fracture, unordered_map<unsigned int, Fracture>& CollectionFractures);

inline bool compare(const pair<unsigned int, double>& coppia1, const pair<unsigned int, double>& coppia2){
    return coppia1.second > coppia2.second;
}
void OutputSort(vector<unsigned int>& IdTrace, unordered_map<unsigned int, Trace>& elencoTracce, ofstream& FileFracture, bool& tips);
}

#endif
