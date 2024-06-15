#pragma once
#include "Fractures.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
struct SottoPoligoni{
    unsigned int id;
    Matrix3Xd Vertici;
    vector<unsigned int> Passanti;
    vector<unsigned int> NonPassanti;
    map<unsigned int, pair<Vector3d, Vector3d>> estremi; //tutti i punti che fanno la traccia associate anche ai valori
};
//prendiamo l elenco delle tracce passanti
void DividiPoligono(Trace& traccia, Fracture& frattura,list<unsigned int>& Sotto_poligoni, map<unsigned int, Trace>& elenco_traccia);
void AnalizzaTraccia(Trace& TracciaTaglio, Trace& Traccia, list<unsigned int>& id_sottopoligoni, Vector3d& VettoreEntrante, Fracture& frattura);


inline unsigned int Regola_Mano_Destra(Vector3d u, Vector3d v, Vector3d controllo_entrante){
    //ora effettuo il prodotto scalare tra il vettore che ho stabilito entrante e il vettore che viene dal prodotto vettoriale
    //se è negativo allora la direzione del vettore considerato è opposta a quello dato, per cui il vettore risulta uscente
    double prodotto_scalare = (u[1]*v[2]-v[1]*u[2])*controllo_entrante[0] +(v[0]*u[2]-u[0]*v[2])*controllo_entrante[1] + (u[0]*v[1]-v[0]*u[1])*controllo_entrante[2];
    if(prodotto_scalare < 0)      // return 0 corrisponde a entrante, return 1 uscente, return 2 a = 0
        return 1; //uscente
    else if(prodotto_scalare > 0)
        return 0; //entrante
    else
        return 2; //uguale a 0
}
}
