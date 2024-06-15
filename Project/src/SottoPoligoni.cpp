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
    map<unsigned int, Vector2d> estremi;
};


bool Pto_Retta (Vector3d p, Vector3d v1, Vector3d v2){ //punto, vertice1, vertice2
    double tol = 1e-9;
    if (fabs((p[2]-v1[2])/(v2[2]-v1[2])-(p[1]-v1[1])/(v2[1]-v1[1]))<tol &&
        fabs((p[2]-v1[2])/(v2[2]-v1[2])-(p[0]-v1[0])/(v2[0]-v1[0]))<tol &&
        fabs((p[1]-v1[1])/(v2[1]-v1[1])-(p[0]-v1[0])/(v2[0]-v1[0]))<tol){
        return true;
    }
    return false;

}

void DividiPoligono(Trace traccia, Fracture& frattura, list<unsigned int> Sotto_poligoni, vector<Trace> elenco_tracce){
    //creiamo il sottopoligono 1 (USCENTE) e sottopoligono2 (ENTRANTE)
    SottoPoligoni uscente;
    SottoPoligoni entrante;
    array<unsigned int,2> estremi_entrante;
    //assegnamo 0 al primo di uscente
    //assegnamo 0 e 1 al primo e al secondo di entrante
    Vector3d start = traccia[frattura.traccePassanti[1]].Vertices.first;        //è per questo che serve "estremi"
    Vector3d end = traccia[frattura.traccePassanti[1]].Vertices.second;
    entrante.Vertici.col(0) = start;
    uscente.Vertici.col(0) = start;     uscente.Vertici.col(1) = end;
    for (unsigned int i=0; i<frattura.numVertici; i++){
        Vector2d parametri = ParametriRetta(frattura.Vertici.col(i),frattura.Vertici.col((i+1)%frattura.numVertici), start, end-start);
        if (parametri[0] >= 0 && parametri[0] <= 1 && parametri[1] >= 0 && parametri[1] <= 1) { //controlliamo che c'è l'intersezione in quel lato
            //se è lo start
            if(Pto_Retta(start, frattura.Vertici.col(i), frattura.Vertici.col((i+1)%frattura.numVertici))){
                estremi_entrante[0] = (i+1)%frattura.numVertici; //prendo il più grande
            }
            else{estremi_entrante[1] = i; //prendo il più piccolo
            }
        }
    }
    unsigned int numVEntr = ((estremi_entrante[1]-estremi_entrante[0]+frattura.numVertici)%frattura.numVertici)+1;  //sono i vertici all'infuori dei due nuovi
    for (unsigned int e=0; e<numVEntr; e++){
        entrante.Vertici.col(e+1) = frattura.Vertici.col((estremi_entrante[0]+e)%frattura.numVertici);
    }
    entrante.Vertici.col(numVEntr+1) = end;
    for (unsigned int u=0; u<frattura.numVertici-numVEntr; u++){
        uscente.Vertici.col(u+2) = frattura.Vertici.col((estremi_entrante[1]+1+u)%frattura.numVertici);
    }
}




// void Taglio(Fracture& frattura, vector<Trace>& tracce){
//     tracce[frattura.traccePassanti[1]].Vertices;
//     frattura.Vertici;
//     Fracture Fratt_dx;
//     Fracture Fratt_sx;
//     Vector2d estremi_sx;
//     for (unsigned int i=0; i<frattura.numVertici; i++){
//         Vector2d parametri = ParametriRetta(frattura.Vertici.col(i),frattura.Vertici.col((i+1)%frattura.numVertici),tracce[frattura.traccePassanti[1]].Vertices.first, (tracce[frattura.traccePassanti[1]].Vertices.second-tracce[frattura.traccePassanti[1]].Vertices.first));
//         if (parametri[0] >= 0 && parametri[0] <= 1 && parametri[1] >= 0 && parametri[1] <= 1) { //controlliamo che c'è l'intersezione in quel lato
//             //se è lo start
//             if(Pto_Retta(tracce[frattura.traccePassanti[1]].Vertices.first,frattura.Vertici.col(i),frattura.Vertici.col((i+1)%frattura.numVertici))){
//                 estremi_sx[0] = frattura.Vertici.col(i); //prendo il piu piccolo
//             }
//             else{estremi_sx[1] = frattura.Vertici.col((i+1)%frattura.numVertici); //prendo il piu grande}
//         }
//     }




    // unsigned int np= frattura.traccePassanti.size();
    // for (unsigned int i=1; i<np; i++){
    //     frattura.traccePassanti[i];
    // }
}

