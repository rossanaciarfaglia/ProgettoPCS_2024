#pragma once
#include "Fractures.hpp"
#include "PolygonalMesh.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <deque>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
struct SottoPoligoni{
    unsigned int id;
    vector<pair<unsigned int, Vector3d>> Vertici;
    unsigned int numVertici;
    vector<unsigned int> Passanti;
    vector<unsigned int> NonPassanti;
    map<unsigned int, pair<pair<unsigned int,Vector3d>,pair<unsigned int,Vector3d>>> estremi; //gli id di tutte le tracce associati ai valori

    SottoPoligoni() = default;
};
//prendiamo l elenco delle tracce passanti
void DividiPoligono(unsigned int& id_tr, SottoPoligoni& frattura, unsigned int& id_sott, map<unsigned int, SottoPoligoni>& Sotto_poligoni,
                    map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni,const string& flag, unsigned int& idSP, unsigned int& idV,
                    PolygonalLibrary::PolygonalMesh& mesh);
void AnalizzaTraccia(Vector3d& start_taglio, Vector3d& end_taglio, SottoPoligoni& taglio, unsigned int& id_traccia, SottoPoligoni& uscente,
                     SottoPoligoni& entrante, Vector3d& VettoreEntrante, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni, unsigned int& idV);


inline unsigned int Regola_Mano_Destra(Vector3d u, Vector3d v, Vector3d controllo_uscente){
    /* ora effettuo il prodotto scalare tra il vettore che ho stabilito entrante e il vettore che viene dal prodotto vettoriale
       se è negativo allora la direzione del vettore considerato è opposta a quello dato, per cui il vettore risulta entrante */
    double prodotto_misto = (u[1]*v[2]-v[1]*u[2])*controllo_uscente[0] + (v[0]*u[2]-u[0]*v[2])*controllo_uscente[1] + (u[0]*v[1]-v[0]*u[1])*controllo_uscente[2];
    if(prodotto_misto > 0)      // return 0 corrisponde a entrante, return 1 uscente, return 2 a = 0
        return 1; //uscente
    else if(prodotto_misto < 0)
        return 0; //entrante
    else
        return 2; //uguale a 0
}
}

// namespace PolygonalLibrary{
// inline void Add_Vert_to_Mesh(PolygonalMesh& mesh, unsigned int& idV, Vector3d& Coordinates){
//     if(mesh.IdCell0D[idV] = idV){
//         return;     //c'è già
//     }
//     mesh.IdCell0D[idV] = idV;
//     mesh.CoordinatesCell0D[idV] = Coordinates;
// }
// }
