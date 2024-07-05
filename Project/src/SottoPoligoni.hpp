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
    vector<pair<unsigned int, pair<unsigned int, unsigned int>>> Lati;

    SottoPoligoni() = default;
};
//prendiamo l elenco delle tracce passanti
void DividiPoligono(unsigned int& id_tr, SottoPoligoni& frattura, unsigned int& id_sott, map<unsigned int, SottoPoligoni>& Sotto_poligoni,
                    map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni,const string& flag, unsigned int& idSP, unsigned int& idV,
                    PolygonalLibrary::PolygonalMesh& mesh, map<unsigned int, vector<unsigned int>>& mappaLati);
void AnalizzaTraccia(Vector3d& start_taglio, Vector3d& end_taglio, SottoPoligoni& taglio, unsigned int& id_traccia, SottoPoligoni& uscente,
                     SottoPoligoni& entrante, Vector3d& VettoreEntrante, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni, unsigned int& idV);


inline unsigned int Regola_Mano_Destra(Vector3d u, Vector3d v, Vector3d controllo_uscente){
    /* ora effettuo il prodotto scalare tra il vettore che ho stabilito entrante e il vettore che viene dal prodotto vettoriale
       se è negativo allora la direzione del vettore considerato è opposta a quello dato, per cui il vettore risulta entrante */
    Vector3d prod_vett = ProdottoVettoriale(u,v);
    double prodotto_misto = prod_vett[0]*controllo_uscente[0] + prod_vett[1]*controllo_uscente[1] + prod_vett[2]*controllo_uscente[2];
    if(prodotto_misto > 0)      // return 0 corrisponde a entrante, return 1 uscente, return 2 a = 0
        return 1; //uscente
    else if(prodotto_misto < 0)
        return 0; //entrante
    else
        return 2; //uguale a 0
}

inline bool Punto_su_Lato(Vector3d p1, Vector3d p2, Vector3d q){
    Matrix3Xd A = p2 - p1;
    Vector3d b = q - p1;
    VectorXd alpha = A.colPivHouseholderQr().solve(b);
    double alfa = alpha[0];

    if ((alfa < 0) || (alfa > 1))
        return false;
    if (abs(q[0] - p1[0] - alfa*(p2[0] - p1[0])) > 1e-8 || abs(q[1] - p1[1] - alfa*(p2[1] - p1[1])) > 1e-8 || abs(q[2] - p1[2] - alfa*(p2[2] - p1[2])) > 1e-8)
        return false;

    return true;
}

}

