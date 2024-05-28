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
    MatrixXd Vertici;

    Fracture() = default; // costruttore di default

    Fracture(unsigned int& numVertici,
              MatrixXd& Vertici
              ):
        numVertici(numVertici),
        Vertici(Vertici)
        {} // prende in input le coordinate dei vertici e la lista dei vertici e inizializza i membri corrispondenti VerticesCoordinates e ListVertices con i valori passati come argomenti

    Vector3d Baricentro(MatrixXd &Poligono); //metodo in quanto proprietà della frattura
    double Raggio(Vector3d &centro, MatrixXd &Poligono);
    Vector4d TrovaPiano(MatrixXd &poligono);
};

inline double DistanzaEuclidea(Vector3d &centro1, Vector3d &centro2) {
    double distanza = 0;
    for (unsigned i = 0; i < 3; i++) {
        distanza += pow(centro1[i] - centro2[i], 2);
    }
    return distanza;
}
bool IntersezioneSfere(Fracture& polygons, MatrixXd& poly_1, MatrixXd& poly_2);


struct Trace{
    unsigned int id;
    pair<Vector3d, Vector3d> Vertices;
    unsigned int id_fract1;
    unsigned int id_fract2;
};


Matrix<double,2,3> IntersezionePiani(Fracture &polygons, MatrixXd &poly_1, MatrixXd &poly_2);


inline Vector2d ParametriRette (const Vector3d& P0, const Vector3d& P1, const Vector3d& Q, const Vector3d& dir_retta){

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

vector<Vector3d> Intersection_Point(Matrix<double,2,3> &retta, MatrixXd &vertici, const unsigned int& numVert);

inline bool isLessOrEqual(Vector3d p1, Vector3d p2, MatrixXd retta_int) {
    return (p1[0] - p2[0]) * retta_int(0,0) + (p1[1] - p2[1]) * retta_int(0,1) + (p1[2] - p2[2]) * retta_int(0,2) <= 0;
}

inline pair<Vector3d, Vector3d> Traccia(vector<Vector3d> &intersezioni1, vector<Vector3d> &intersezioni2){
    //trova l'inizio e la fine dell'intersezione
    Vector3d intersection_start;
    Vector3d intersection_end ;
    bool overlap = true;

    for (int i = 0; i < 3; ++i) {
        intersection_start[i] = max(intersezioni1[0][i], intersezioni2[0][i]);
        intersection_end[i] = min(intersezioni1[1][i], intersezioni2[1][i]);

        if (intersection_start[i] > intersection_end[i]) {
            overlap = false;
        }
        if (!overlap) {
            // Se non c'è sovrapposizione, possiamo restituire un valore indicativo
            // di nessuna intersezione, come due punti uguali o una coppia di zero.
            return {Vector3d::Zero(), Vector3d::Zero()};
        }

        return {intersection_start, intersection_end};
    }
}

void Find_Trace(Fracture& polygon, MatrixXd& vert_1, MatrixXd& vert_2, const unsigned int& numVert_1, const unsigned int& numVert_2);

void ImportFracturesList(const string& filepath, Fracture& fracture, unordered_map<unsigned int, Fracture>& CollectionFractures);

}

#endif
