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
struct Fractures{
    unsigned int numFractures;
    map<unsigned int, MatrixXd> FracturesMap;

    Fractures() = default; // costruttore di default

    Fractures(unsigned int& numFractures,
              map<unsigned int, MatrixXd>& FracturesMap
              ):
        numFractures(numFractures),
        FracturesMap(FracturesMap)
    {} // prende in input le coordinate dei vertici e la lista dei vertici e inizializza i membri corrispondenti VerticesCoordinates e ListVertices con i valori passati come argomenti

    Vector3d Baricentro(MatrixXd &Poligono);
    double Raggio(Vector3d &centro, MatrixXd &Poligono);
    Vector4d TrovaPiano(MatrixXd &poligono);
};

void ImportFracturesList(const string& filepath, Fractures& fractures);
double DistanzaEuclidea(Vector3d &centro1, Vector3d &centro2);
bool IntersezioneSfere(Fractures &polygons, MatrixXd &poly_1, MatrixXd &poly_2);

void ImportFracturesList(const string& filepath, Fractures& fractures);


struct Trace{
    unsigned int id;
    Matrix2d Vertices;
    unsigned int id_fract1;
    unsigned int id_fract2;
};

vector<Vector3d> Intersection_Point(Matrix<double, 2, 3> &retta, MatrixXd &vertici);

MatrixXd IntersezionePiani(Fractures &polygons, vector<Vector3d> intersezioni1, vector<Vector3d> intersezioni2,
                           MatrixXd &poly_1, MatrixXd &poly_2);

inline bool isLessOrEqual(Vector3d p1, Vector3d p2, Vector3d t) {
    return (p1[0] - p2[0]) * t[0] + (p1[1] - p2[1]) * t[1] + (p1[2] - p2[2]) * t[2] <= 0;
}


inline pair<Vector3d, Vector3d> Traccia(vector<Vector3d> &intersezioni1,
                                        vector<Vector3d> &intersezioni2){

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


void Find_Trace(vector<Vector3d> &intersezioni1, vector<Vector3d> &intersezioni2);

}


#endif
