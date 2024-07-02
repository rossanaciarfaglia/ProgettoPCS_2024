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
}
