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
    MatrixXd VerticesCoordinates;
    map<unsigned int, MatrixXd> FracturesMap;

    Fractures() = default; // costruttore di default

    Fractures(const MatrixXd& VerticesCoordinates,
              unsigned int& numFractures,
              map<unsigned int, MatrixXd>& FracturesMap
              ):
        VerticesCoordinates(VerticesCoordinates),
        numFractures(numFractures),
        FracturesMap(FracturesMap)
    {} // prende in input le coordinate dei vertici e la lista dei vertici e inizializza i membri corrispondenti VerticesCoordinates e ListVertices con i valori passati come argomenti

    vector<vector<unsigned int>> SpheresIntersection(); // Portarla fuori dalla struct
    pair<map<unsigned int, vector<unsigned int>>, map<unsigned int, vector<unsigned int>>> Passante_NonPassante(ofstream& outputFile);
    vector<double> Baricentro(MatrixXd Poligono);
    double Raggio(vector<double> centro, MatrixXd Poligono);
    double DistanzaEuclidea(vector<double> centro1, vector<double> centro2);
    Vector4d TrovaPiano(MatrixXd Poligono);
};
void ImportFracturesList(const string& filepath,
                         Fractures& fractures);
bool IntersezioneSfere(Fractures &polygons, unsigned int id1, unsigned int id2);
#endif
}
