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


void AnalizzaTraccia(Vector3d& start_taglio, Vector3d& end_taglio, SottoPoligoni& taglio, unsigned int& id_traccia, SottoPoligoni& uscente, SottoPoligoni& entrante, Vector3d& VettoreEntrante, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni){
    //prendo gli estremi della traccia
    Vector3d start = taglio.estremi[id_traccia].first;
    Vector3d end = taglio.estremi[id_traccia].second;

    unsigned int count = 0; //per valutare se la traccia appartiene a un solo sottopoligono o a entrambi

    if(Regola_Mano_Destra(end_taglio - start_taglio, start - start_taglio, VettoreEntrante) == 1) {       //uscente
        //uscente.estremi.first e lo posiamo uguale a start
        uscente.estremi[id_traccia].first = start;
        count ++;}

    else if(Regola_Mano_Destra(end_taglio - start_taglio, start - start_taglio, VettoreEntrante) == 0)  //entrante
        entrante.estremi[id_traccia].first = start;
    else
        cout<<"qui è uguale a zero ma non so che farci"<<endl;

    //lo analizzo con l'end
    if(Regola_Mano_Destra(end_taglio - start_taglio, end - start_taglio, VettoreEntrante) == 1) {       //uscente
        //uscente.estremi.first e lo posiamo uguale a start
        uscente.estremi[id_traccia].second = end;
        count ++;}
    else if(Regola_Mano_Destra(end_taglio - start_taglio, end - start_taglio, VettoreEntrante) == 0)  //entrante
        entrante.estremi[id_traccia].second = end;
    else {
        cout<<"le due tracce si sovrappongono"<<endl;
        return;}

    bool trovato = false; //analizza se troviamo l'id in tracce passanti o non passanti della frattura grande
    for (int i = 0; i < taglio.Passanti.size(); i++) {
        if (taglio.Passanti[i] == id_traccia)
            trovato = true;
        break;
    }
    switch (count) {
    case 0:
        if (trovato){
            entrante.Passanti.push_back(id_traccia);
        }
        else {
            entrante.NonPassanti.push_back(id_traccia);
        }
        break;

    case 2:
        if (trovato){
            uscente.Passanti.push_back(id_traccia);
        }
        else {
            uscente.NonPassanti.push_back(id_traccia);
        }
        break;


    case 1:
        //le tracce si intersecano
        //quindi dobbiamo trovare l'intersezione
        Vector2d solution = ParametriRetta(start, end, start_taglio, end_taglio-start_taglio);
        Vector3d intersezione = start + solution[0]*(end - start); //abbiamo trovato l'intersezione
        if (trovato){
            //devo inserire il punto di intersezione
            //non so se viene inserirlo come start o end --> lo passo alla funzione isLessorEqual e lo swap
            if ((uscente.estremi[id_traccia].first - start).norm() <= 1e-14){
                //significa che lo start sta in uscente
                uscente.estremi[id_traccia].second = intersezione;
                entrante.estremi[id_traccia].first = intersezione;
            }
            else {
                uscente.estremi[id_traccia].first = intersezione;
                entrante.estremi[id_traccia].second = intersezione;
            }
            //aggiungo alla traccia l'id dei sottopoligoni ad essa associati
            Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
            Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
            //adesso aggiungiamo la traccia ai sottopoligoni passanti
            uscente.Passanti.push_back(id_traccia);
            entrante.Passanti.push_back(id_traccia);
        }
        else { //siamo in non passanti
            bool flag = false; //controllo che sia passante per uno dei due
            //se il prodotto vettoriale tra i lati del sottopoligono e il vettore dato dalla distanza tra uno dei punti tra start ed end e uno dei vertici del lato
            //è 0 --> start o end (a seconda di quello scelto) appartiene al poligono
            for(unsigned int i = 0; i < uscente.numVertici; i++){ //qui vedo se è passante per il sottopoligono uscente
                if (Regola_Mano_Destra(uscente.Vertici.col((i+1)%uscente.numVertici)-uscente.Vertici.col(i),start-uscente.Vertici.col(i),uscente.Vertici.col((i+1)%uscente.numVertici)-uscente.Vertici.col(i)) == 2 ||
                    Regola_Mano_Destra(uscente.Vertici.col((i+1)%uscente.numVertici)-uscente.Vertici.col(i),end-uscente.Vertici.col(i),uscente.Vertici.col((i+1)%uscente.numVertici)-uscente.Vertici.col(i)) == 2){
                    if ((uscente.estremi[id_traccia].first - start).norm() <= 1e-14){
                        //significa che lo start sta in uscente
                        uscente.estremi[id_traccia].second = intersezione;
                        entrante.estremi[id_traccia].first = intersezione;
                    }
                    else {
                        uscente.estremi[id_traccia].first = intersezione;
                        entrante.estremi[id_traccia].second = intersezione;
                    }
                    Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
                    uscente.Passanti.push_back(id_traccia);
                    Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
                    entrante.NonPassanti.push_back(id_traccia);
                    flag = true;
                }
            }
            //ora controllo se è passante per il sottopoligono entrante
            if (flag == false){
                for(unsigned int i = 0; i < entrante.numVertici; i++){ //qui vedo se è passante per il sottopoligono uscente
                    if (Regola_Mano_Destra(entrante.Vertici.col((i+1)%entrante.numVertici)-entrante.Vertici.col(i),start-entrante.Vertici.col(i),entrante.Vertici.col((i+1)%entrante.numVertici)-entrante.Vertici.col(i)) == 2 ||
                        Regola_Mano_Destra(entrante.Vertici.col((i+1)%entrante.numVertici)-entrante.Vertici.col(i),end-entrante.Vertici.col(i),entrante.Vertici.col((i+1)%entrante.numVertici)-entrante.Vertici.col(i)) == 2){
                        if ((entrante.estremi[id_traccia].first - start).norm() <= 1e-14){
                            //significa che lo start sta in uscente
                            entrante.estremi[id_traccia].second = intersezione;
                            uscente.estremi[id_traccia].first = intersezione;
                        }
                        else {
                            entrante.estremi[id_traccia].first = intersezione;
                            uscente.estremi[id_traccia].second = intersezione;
                        }
                        Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
                        entrante.Passanti.push_back(id_traccia);
                        Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
                        uscente.NonPassanti.push_back(id_traccia);
                        flag = true;
                    }
                }
            }

            if(flag == false){
                if ((entrante.estremi[id_traccia].first - start).norm() <= 1e-14){
                    //significa che lo start sta in uscente
                    entrante.estremi[id_traccia].second = intersezione;
                    uscente.estremi[id_traccia].first = intersezione;
                }
                else {
                    entrante.estremi[id_traccia].first = intersezione;
                    uscente.estremi[id_traccia].second = intersezione;
                }
                Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
                entrante.NonPassanti.push_back(id_traccia);
                Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
                uscente.NonPassanti.push_back(id_traccia);
                flag = true;
            }
        }
        break;
    }
}


// void PrimoDividiPoligono(Trace traccia, Fracture& frattura, list<unsigned int>& Sotto_poligoni, vector<Trace> elenco_tracce){
//     //creiamo il sottopoligono 1 (USCENTE) e sottopoligono2 (ENTRANTE)
//     SottoPoligoni uscente;
//     SottoPoligoni entrante;
//     array<unsigned int,2> estremi_entrante;
//     //assegnamo 0 al primo di uscente
//     //assegnamo 0 e 1 al primo e al secondo di entrante
//     Vector3d start = traccia.Vertices.first;        //è per questo che serve "estremi"
//     Vector3d end = traccia[frattura.traccePassanti[1]].Vertices.second;
//     entrante.Vertici.col(0) = start;
//     uscente.Vertici.col(0) = start;     uscente.Vertici.col(1) = end;

//     //iteriamo sui vertici dei poligoni
//     for (unsigned int i=0; i < frattura.numVertici; i++){
//         Vector2d parametri = ParametriRetta(frattura.Vertici.col(i),frattura.Vertici.col((i+1)%frattura.numVertici), start, end-start);
//         if (parametri[0] >= 0 && parametri[0] <= 1 && parametri[1] >= 0 && parametri[1] <= 1) { //controlliamo che c'è l'intersezione in quel lato
//             //se è lo start
//             if(Pto_Retta(start, frattura.Vertici.col(i), frattura.Vertici.col((i+1)%frattura.numVertici))){
//                 estremi_entrante[0] = (i+1)%frattura.numVertici; //prendo il più grande
//             }
//             else{estremi_entrante[1] = i; //prendo il più piccolo
//             }
//         }
//     }

//     unsigned int numVEntr = ((estremi_entrante[1]-estremi_entrante[0]+frattura.numVertici)%frattura.numVertici)+1;  //sono i vertici all'infuori dei due nuovi
//     for (unsigned int e=0; e<numVEntr; e++){
//         entrante.Vertici.col(e+1) = frattura.Vertici.col((estremi_entrante[0]+e)%frattura.numVertici);
//     }
//     entrante.Vertici.col(numVEntr+1) = end;
//     entrante.numVertici = numVEntr + 2;
//     for (unsigned int u=0; u < frattura.numVertici-numVEntr; u++){
//         uscente.Vertici.col(u+2) = frattura.Vertici.col((estremi_entrante[1]+1+u)%frattura.numVertici);
//     }
//     uscente.numVertici = frattura.numVertici-numVEntr + 2;

//     //stabilinamo una direzione rispetto alla quale definiremo la regola della mano desctra restituisce un vettore uscente o entrante(positivo sarà uscente, negativo entrante)
//     Vector3d Direzione_Entrante = uscente.Vertici.col(1) - start;

//     for(unsigned int k = 1; k < frattura.traccePassanti.size() ; k++ ){  //itero su tutte le passanti
//         AnalizzaTraccia(elenco_tracce[frattura.traccePassanti[1]], elenco_tracce[frattura.traccePassanti[k]], uscente, entrante, Direzione_Entrante, frattura); //struttura trace delle tracce che corrispondono a questo id
//     }

// }



void DividiPoligono(unsigned int& id_tr, SottoPoligoni& frattura, list<SottoPoligoni>& Sotto_poligoni, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni){
    //creiamo il sottopoligono 1 (USCENTE) e sottopoligono2 (ENTRANTE)
    SottoPoligoni uscente;
    SottoPoligoni entrante;
    Vector3d start;
    Vector3d end;
    if (flag= passante){
        start = frattura.estremi[id_tr].first;
        end = frattura.estremi[id_tr].second;
    }
    else{
        //allunga traccia
        // start ed end sono i punti di intersezione
    }
    array<unsigned int,2> estremi_entrante;
    //assegnamo 0 al primo di uscente
    //assegnamo 0 e 1 al primo e al secondo di entrante
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
    entrante.numVertici = numVEntr + 2;
    for (unsigned int u=0; u < frattura.numVertici-numVEntr; u++){
        uscente.Vertici.col(u+2) = frattura.Vertici.col((estremi_entrante[1]+1+u)%frattura.numVertici);
    }
    uscente.numVertici = frattura.numVertici-numVEntr + 2;

    //stabilinamo una direzione rispetto alla quale definiremo la regola della mano destra restituisce un vettore uscente o entrante(positivo sarà uscente, negativo entrante)
    Vector3d Direzione_Entrante = uscente.Vertici.col(1) - start;

    for(unsigned int k = 0; k < frattura.estremi.size() ; k++ ){  //itero su tutte le tracce
        if(k == id_tr){
            break; }
        AnalizzaTraccia(start, end, frattura, k, uscente, entrante, Direzione_Entrante, Tracce_SottoPoligoni); //struttura trace delle tracce che corrispondono a questo id
    }
    Sotto_poligoni.push_back(uscente);
    Sotto_poligoni.push_back(entrante);
    Sotto_poligoni.pop_front();
}

}




