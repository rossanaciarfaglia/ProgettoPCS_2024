#include "Fractures.hpp"
#include "SottoPoligoni.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
bool Pto_Retta (Vector3d p, Vector3d v1, Vector3d v2){ //punto, vertice1, vertice2
    double tol = 1e-9;
    if (fabs((p[2]-v1[2])/(v2[2]-v1[2])-(p[1]-v1[1])/(v2[1]-v1[1]))<tol &&
        fabs((p[2]-v1[2])/(v2[2]-v1[2])-(p[0]-v1[0])/(v2[0]-v1[0]))<tol &&
        fabs((p[1]-v1[1])/(v2[1]-v1[1])-(p[0]-v1[0])/(v2[0]-v1[0]))<tol){
        return true;
    }
    return false;

}


void AnalizzaTraccia(Trace& TracciaTaglio, Trace& traccia, SottoPoligoni uscente, SottoPoligoni entrante, Vector3d VettoreEntrante, Fracture& frattura){
    //prendo gli estremi della traccia
    Vector3d start = traccia.Vertices.first;
    Vector3d end = traccia.Vertices.second;

    unsigned int count = 0; //per valutare se la traccia appartiene a un solo sottopoligono o a entrambi

    if(Regola_Mano_Destra(TracciaTaglio.Vertices.second - TracciaTaglio.Vertices.first, start - TracciaTaglio.Vertices.first, VettoreEntrante) == 1) {       //uscente
        //uscente.estremi.first e lo posiamo uguale a start
        uscente.estremi.first = start;
        count ++;}

    else if(Regola_Mano_Destra(TracciaTaglio.Vertices.second - TracciaTaglio.Vertices.first, start - TracciaTaglio.Vertices.first, VettoreEntrante) == 0)  //entrante
        entrante.estremi.first = start;
    else
        cout<<"qui è uguale a zero ma non so che farci"<<endl;

    //lo analizzo con l end
    if(Regola_Mano_Destra(TracciaTaglio.Vertices.second - TracciaTaglio.Vertices.first, end - TracciaTaglio.Vertices.first, VettoreEntrante) == 1) {       //uscente
        //uscente.estremi.first e lo posiamo uguale a start
        uscente.estremi.second = end;
        count ++;}

    else if(Regola_Mano_Destra(TracciaTaglio.Vertices.second - TracciaTaglio.Vertices.first, end - TracciaTaglio.Vertices.first, VettoreEntrante) == 0)  //entrante
        entrante.estremi.second = end;
    else
        cout<<"qui è uguale a zero ma non so che farci"<<endl;


    switch (count) {
    case 0:
        bool trovato; //analizza se troviamo l id in tracce passanti o non passanti della frattura grande
        for (int i = 0; i < frattura.traccePassanti.size(); i++) {
            if (frattura.traccePassanti.size()[i] == traccia.id)
                trovato = true;
            else
                trovato = false;
            break;
        }


        break;
    case 2:
        cout << "Hai scelto l'opzione 2" << std::endl;
        break;
    case 1:
        cout << "Hai scelto l'opzione 3" << std::endl;
        break;
    }
}


//non va bene il metodo della lista
void DividiPoligono(Trace traccia, Fracture& frattura, list<unsigned int>& Sotto_poligoni, vector<Trace> elenco_tracce){
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

    //iteriamo sui vertici dei poligoni
    for (unsigned int i=0; i < frattura.numVertici; i++){
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
    for (unsigned int u=0; u < frattura.numVertici-numVEntr; u++){
        uscente.Vertici.col(u+2) = frattura.Vertici.col((estremi_entrante[1]+1+u)%frattura.numVertici);
    }

    //stabilinamo una direzione rispetto alla quale definiremo la regola della mano desctra restituisce un vettore uscente o entrante(positivo sarà uscente, negativo entrante)
    Vector3d Direzione_Entrante = uscente.Vertici.col(1) - start;

    for(unsigned int k = 1; k < frattura.traccePassanti.size() ; k++ ){  //itero su tutte le passanti
        AnalizzaTraccia(elenco_tracce[frattura.traccePassanti[1]], elenco_tracce[frattura.traccePassanti[k]], uscente, entrante, Direzione_Entrante, frattura); //struttura trace delle tracce che corrispondono a questo id
    }

}

}




