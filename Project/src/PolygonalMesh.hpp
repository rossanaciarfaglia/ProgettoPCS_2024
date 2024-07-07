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
    vector<vector<unsigned int>> VerticesCell2D;        // gli id dei vertici in senso antiorario
    vector<vector<unsigned int>> EdgesCell2D;       // gli id dei lati in senso antiorario
};

inline void Add_Vert_to_Mesh(PolygonalMesh& mesh, pair<unsigned int, Vector3d> vertex){
    mesh.IdCell0D.push_back(vertex.first);
    mesh.CoordinatesCell0D.push_back(vertex.second);
}


inline void Ordina_Punti(vector<unsigned int>& padre, vector<Vector3d> punti, vector<unsigned int>& ordinamento){
    vector<double> coordinate(padre.size(),0);
    Matrix3Xd A = punti[padre[1]] - punti[padre[0]];
    Vector3d b;
    VectorXd alpha;
    for(unsigned int pt=0; pt<padre.size(); pt++){
        b = punti[padre[pt]] - punti[padre[0]];
        alpha = A.colPivHouseholderQr().solve(b);
        coordinate[pt] = alpha[0];
    }

    for (unsigned int i=0; i<ordinamento.size(); i++){
        ordinamento[i] = i;
    }
    sort(ordinamento.begin(), ordinamento.end(), [&](const unsigned int& a, const unsigned int& b){return(coordinate[a] < coordinate[b]);});

}

}
