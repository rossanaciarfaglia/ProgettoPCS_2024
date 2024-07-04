#pragma once

#include <ostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary {

struct PolygonalMesh {
    vector<unsigned int> IdCell0D;     // tutti gli id dei vertici
    vector<Vector3d> CoordinatesCell0D;    // le coordinate per ciascun vertice

    vector<unsigned int> IdCell1D;     // tutti gli id dei lati
    vector<Vector2i> VerticesCell1D;    // gli id di start ed end per ciascun lato

    vector<unsigned int> NumberElements2D;      // il num dei vertici (= num dei lati) che lo compongono
    vector<VectorXi> VerticesCell2D;        // gli id dei vertici in senso antiorario
    vector<VectorXi> EdgesCell2D;       // gli id dei lati in senso antiorario
};

inline void Add_Vert_to_Mesh(PolygonalMesh& mesh, pair<unsigned int, Vector3d> vertex){
    mesh.IdCell0D.push_back(vertex.first);
    mesh.CoordinatesCell0D.push_back(vertex.second);
}

inline void Add_Edge_to_Mesh(PolygonalMesh& mesh){
    mesh.IdCell1D.push_back(0);
    mesh.VerticesCell1D.push_back({0,0});
}

inline void Ordina_Punti(vector<unsigned int>& padre, vector<Vector3d> punti, vector<unsigned int>& ordinati){

}

}
