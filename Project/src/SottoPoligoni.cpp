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
    bool sovrapp_start = false;
    bool sovrapp_end = false;

    if(Regola_Mano_Destra(end_taglio - start_taglio, start - start_taglio, VettoreEntrante) == 1) {       //uscente
        uscente.estremi[id_traccia].first = start;
        count ++;}
    else if(Regola_Mano_Destra(end_taglio - start_taglio, start - start_taglio, VettoreEntrante) == 0)  //entrante
        entrante.estremi[id_traccia].first = start;
    else {
        //start è sulla traccia di taglio
        sovrapp_start = true;
    }

    //lo analizzo con l'end
    if(Regola_Mano_Destra(end_taglio - start_taglio, end - start_taglio, VettoreEntrante) == 1) {       //uscente
        uscente.estremi[id_traccia].second = end;
        count ++;}
    else if(Regola_Mano_Destra(end_taglio - start_taglio, end - start_taglio, VettoreEntrante) == 0)  //entrante
        entrante.estremi[id_traccia].second = end;
    else {
        //end è sulla traccia di taglio
        sovrapp_end = true;
    }

    bool trovato = false; //analizza se troviamo l'id in tracce passanti della frattura grande (altrimenti è in non passanti)
    for (int i = 0; i < taglio.Passanti.size(); i++) {
        if (taglio.Passanti[i] == id_traccia)
            trovato = true;
        break;
    }

    if (sovrapp_start = true){
        switch (count) {
        case 0:
            entrante.estremi[id_traccia].first = start;
            if (trovato){
                entrante.Passanti.push_back(id_traccia);
            }
            else {
                entrante.NonPassanti.push_back(id_traccia);
            }
            break;
        case 1:
            uscente.estremi[id_traccia].first = start;
            if (trovato){
                uscente.Passanti.push_back(id_traccia);
            }
            else {
                uscente.NonPassanti.push_back(id_traccia);
            }
            break;
        }
    }
    else if (sovrapp_end = true){   // Per assunzione della correttezza dei file forniti, qui possiamo scrivere "else if" invece di "if" perchè il caso di sovrapposizione completa di due tracce non esiste
        switch (count) {
        case 0:
            entrante.estremi[id_traccia].second = end;
            if (trovato){
                entrante.Passanti.push_back(id_traccia);
            }
            else {
                entrante.NonPassanti.push_back(id_traccia);
            }
            break;
        case 1:
            uscente.estremi[id_traccia].second = end;
            if (trovato){
                uscente.Passanti.push_back(id_traccia);
            }
            else {
                uscente.NonPassanti.push_back(id_traccia);
            }
            break;
        }
    }
    else{
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
            /*Le tracce si intersecano, quindi dobbiamo trovare l'intersezione
            In questo caso allo stesso id_traccia saranno associate due tracce "diverse" a seconda che sono in uscente o entrante (NON nello stesso) */
            Vector2d coeff = CoefficientiRette(start, end, start_taglio, end_taglio-start_taglio);
            Vector3d intersezione = start + coeff[0]*(end - start); //abbiamo trovato l'intersezione
            if (trovato){
                //devo inserire il punto di intersezione
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
                    if (Regola_Mano_Destra(uscente.Vertici[(i+1)%uscente.numVertici].second-uscente.Vertici[i].second, start-uscente.Vertici[i].second, VettoreEntrante) == 2 ||
                        Regola_Mano_Destra(uscente.Vertici[(i+1)%uscente.numVertici].second-uscente.Vertici[i].second, end-uscente.Vertici[i].second, VettoreEntrante) == 2){
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
                        if (Regola_Mano_Destra(entrante.Vertici[(i+1)%entrante.numVertici].second-entrante.Vertici[i].second, start-entrante.Vertici[i].second, VettoreEntrante) == 2 ||
                            Regola_Mano_Destra(entrante.Vertici[(i+1)%entrante.numVertici].second-entrante.Vertici[i].second, end-entrante.Vertici[i].second, VettoreEntrante) == 2){
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
}



void DividiPoligono(unsigned int& id_tr, SottoPoligoni& frattura, map<unsigned int, SottoPoligoni>& Sotto_poligoni, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni, const string& flag, unsigned int& idSP, unsigned int& idV){
    //creiamo il sottopoligono 1 (USCENTE) e sottopoligono2 (ENTRANTE)
    SottoPoligoni uscente;
    SottoPoligoni entrante;
    uscente.Vertici.reserve(5);
    entrante.Vertici.reserve(5);
    Vector3d start;
    Vector3d end;
    if (flag == "passanti"){
        start = frattura.estremi[id_tr].first;
        end = frattura.estremi[id_tr].second;
    }
    else if(flag == "nonpassanti"){
        //allunga traccia
        // start ed end sono i punti di intersezione
        Matrix<double,2,3> traccia;
        traccia.row(0) = frattura.estremi[id_tr].first;
        traccia.row(1) = frattura.estremi[id_tr].second - frattura.estremi[id_tr].first;
        Matrix3Xd vertici;
        vertici.resize(3, frattura.numVertici);
        for (unsigned int m=0; m<frattura.numVertici; m++){
            vertici.col(m) = frattura.Vertici[m].second;
        }
        vector<Vector3d> intersezioni = Intersection_Point(traccia, vertici, frattura.numVertici);

        if (isLess(intersezioni[1], intersezioni[0], traccia)) swap(intersezioni[1], intersezioni[0]);
        start = intersezioni[0];
        end = intersezioni[1];
    }
    array<unsigned int,2> estremi_entrante;
    //assegnamo 0 al primo di uscente
    //assegnamo 0 e 1 al primo e al secondo di entrante
    entrante.Vertici.push_back({idV, start});
    uscente.Vertici.push_back({idV, start});
    idV++;
    uscente.Vertici.push_back({idV, end});

    //iteriamo sui vertici dei poligoni
    for (unsigned int i=0; i < frattura.numVertici; i++){
        Vector2d parametri = CoefficientiRette(frattura.Vertici[i].second,frattura.Vertici[(i+1)%frattura.numVertici].second, start, end-start);
        if (parametri[0] >= 0 && parametri[0] <= 1 && parametri[1] >= 0 && parametri[1] <= 1) { //controlliamo che c'è l'intersezione in quel lato
            //se è lo start
            if(Pto_Retta(start, frattura.Vertici[i].second, frattura.Vertici[(i+1)%frattura.numVertici].second)){
                estremi_entrante[0] = (i+1)%frattura.numVertici; //prendo il più grande
            }
            else{estremi_entrante[1] = i; //prendo il più piccolo
            }
        }
    }

    unsigned int numVEntr = ((estremi_entrante[1]-estremi_entrante[0]+frattura.numVertici)%frattura.numVertici)+1;  //sono i vertici all'infuori dei due nuovi
    for (unsigned int e=0; e<numVEntr; e++){
        entrante.Vertici.push_back({idV-frattura.numVertici+estremi_entrante[0]+e, frattura.Vertici[(estremi_entrante[0]+e)%frattura.numVertici].second});
    }
    entrante.Vertici.push_back({idV, end});
    entrante.numVertici = numVEntr + 2;
    for (unsigned int u=0; u < frattura.numVertici-numVEntr; u++){
        uscente.Vertici.push_back({idV-frattura.numVertici+estremi_entrante[1]+1+u, frattura.Vertici[(estremi_entrante[1]+1+u)%frattura.numVertici].second});
    }
    uscente.numVertici = frattura.numVertici-numVEntr + 2;
    cout << "e: " << entrante.numVertici << endl << "u: " << uscente.numVertici << endl;

    //stabilinamo una direzione rispetto alla quale definiremo la regola della mano destra restituisce un vettore uscente o entrante(positivo sarà entrante, negativo uscente)
    Vector3d Direzione_Entrante = entrante.Vertici[1].second - start;

    for(unsigned int k = 0; k < frattura.estremi.size() ; k++){  //itero su tutte le tracce
        if(k == id_tr){
            continue; }
        AnalizzaTraccia(start, end, frattura, k, uscente, entrante, Direzione_Entrante, Tracce_SottoPoligoni); //struttura trace delle tracce che corrispondono a questo id
    }
    idV++;
    Sotto_poligoni.erase(idSP);
    idSP ++;
    Sotto_poligoni.insert({idSP, uscente});
    idSP ++;
    Sotto_poligoni.insert({idSP, entrante});
}

}




