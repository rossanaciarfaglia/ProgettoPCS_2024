#ifndef __FRACTURES_H
#define __FRACTURES_H

#include "Eigen/Eigen"
#include <iostream>
#include <utility>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
struct Fractures{
    map<unsigned int, MatrixXd> FracturesMap;
    MatrixXd VerticesCoordinates;
    vector<vector<unsigned int>> listVertices;
    unsigned int numFractures;

    Fractures() = default; //costruttore di default

    Fractures(const MatrixXd& VerticesCoordinates,
              const vector<vector<unsigned int>>& listVertices,
              unsigned int& numFractures,
              map<unsigned int, MatrixXd>& FracturesMap
              ):
        VerticesCoordinates(VerticesCoordinates),
        listVertices(listVertices),
        numFractures(numFractures)
    {} //prende in input le coordinate dei vertici e la lista dei vertici e inizializza i memebri corrispondenti VerticesCoordinates e ListVertices con i valori passati come argomenti

    vector<vector<unsigned int>> SpheresIntersection();
    //creiamo una funzione che ci restituisca una mappa con le tracce passanti e una con quelli non passanti
    pair<map<unsigned int, vector<unsigned int>>, map<unsigned int, vector<unsigned int>>> Passante_NonPassante();

//LEI QUI INSERISCE IL tRIANGULATEpOLYGONS
};
void ImportFracturesList(const string& filepath,
                     Fractures& fractures);
}
#endif
