#include "Fractures.hpp"
#include "SottoPoligoni.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <utility>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {

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
            Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
            break;
        case 1:
            uscente.estremi[id_traccia].first = start;
            if (trovato){
                uscente.Passanti.push_back(id_traccia);
            }
            else {
                uscente.NonPassanti.push_back(id_traccia);
            }
            Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
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
            Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
            break;
        case 1:
            uscente.estremi[id_traccia].second = end;
            if (trovato){
                uscente.Passanti.push_back(id_traccia);
            }
            else {
                uscente.NonPassanti.push_back(id_traccia);
            }
            Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
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
            Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
            break;

        case 2:
            if (trovato){
                uscente.Passanti.push_back(id_traccia);
            }
            else {
                uscente.NonPassanti.push_back(id_traccia);
            }
            Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
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
                    if (Regola_Mano_Destra(uscente.Vertici[(i+1)%uscente.numVertici].second - uscente.Vertici[i].second, start - uscente.Vertici[i].second, VettoreEntrante) == 2 ||
                        Regola_Mano_Destra(uscente.Vertici[(i+1)%uscente.numVertici].second - uscente.Vertici[i].second, end - uscente.Vertici[i].second, VettoreEntrante) == 2){
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
                    for(unsigned int i = 0; i < entrante.numVertici; i++){ //qui vedo se è passante per il sottopoligono entrante
                        if (Regola_Mano_Destra(entrante.Vertici[(i+1)%entrante.numVertici].second - entrante.Vertici[i].second, start - entrante.Vertici[i].second, VettoreEntrante) == 2 ||
                            Regola_Mano_Destra(entrante.Vertici[(i+1)%entrante.numVertici].second - entrante.Vertici[i].second, end - entrante.Vertici[i].second, VettoreEntrante) == 2){
                            if ((entrante.estremi[id_traccia].first - start).norm() <= 1e-14){
                                //significa che lo start sta in entrante
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



void DividiPoligono(unsigned int& id_tr, SottoPoligoni& frattura, unsigned int& id_sott, map<unsigned int, SottoPoligoni>& Sotto_poligoni, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni, const string& flag, unsigned int& idSP, unsigned int& idV, PolygonalLibrary::PolygonalMesh& mesh){
    //stabiliamo una direzione rispetto alla quale definiremo la regola della mano destra
    Vector3d Direzione_Uscente = ProdottoVettoriale(frattura.Vertici[1].second - frattura.Vertici[0].second,
                                                    frattura.Vertici[2].second - frattura.Vertici[1].second);   // Per assunzione della correttezza dei file, i poligoni hanno almeno 3 vertici
    //creiamo il sottopoligono 1 (USCENTE) e sottopoligono2 (ENTRANTE)
    SottoPoligoni uscente;
    SottoPoligoni entrante;
    uscente.id = idSP+1;
    entrante.id = idSP+2;
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
    if (){
        mesh.IdCell0D.push_back(idV);
        mesh.CoordinatesCell0D.push_back(start);
        mesh.IdCell0D.push_back(idV+1);
        mesh.CoordinatesCell0D.push_back(end);
    }


    bool segn; // tiene segno dell'ultimo prodotto misto per sapere se sono passato dal sottopoligono uscente a quello entrante (o viceversa)

    //iteriamo sui vertici dei poligoni
    cout << "vertici sottopoligono: ";
    if (Regola_Mano_Destra(end-start, frattura.Vertici[0].second-start, Direzione_Uscente) == 1){
        uscente.Vertici.push_back({frattura.Vertici[0].first, frattura.Vertici[0].second});
        segn = true;
    }
    else {
        entrante.Vertici.push_back({frattura.Vertici[0].first, frattura.Vertici[0].second});
        segn = false;
    }
    cout << frattura.Vertici[0].first << " ";
    for (unsigned int i=1; i < frattura.numVertici; i++){
        if (Regola_Mano_Destra(end-start, frattura.Vertici[i].second-start, Direzione_Uscente) == 1){
            if (segn == false){
                uscente.Vertici.push_back({idV + 1,end});
                entrante.Vertici.push_back({idV + 1, end});
                // PolygonalLibrary::Add_Vert_to_Mesh(mesh, idV+1, end);
                cout << idV+1 << " ";
            }
            uscente.Vertici.push_back({frattura.Vertici[i].first, frattura.Vertici[i].second});
            segn = true;
        }
        else {
            if (segn == true){
                uscente.Vertici.push_back({idV,start});
                entrante.Vertici.push_back({idV,start});
                // PolygonalLibrary::Add_Vert_to_Mesh(mesh, idV, start);
                cout << idV << " ";
            }
            entrante.Vertici.push_back({frattura.Vertici[i].first, frattura.Vertici[i].second});
            segn = false;
        }
        cout << frattura.Vertici[i].first << " ";
    }
    if (Regola_Mano_Destra(end-start, frattura.Vertici[0].second-start, Direzione_Uscente) == 1){
        if (segn == false){
            uscente.Vertici.push_back({idV + 1,end});
            entrante.Vertici.push_back({idV + 1, end});
            // PolygonalLibrary::Add_Vert_to_Mesh(mesh, idV+1, end);
            cout << idV+1 << " ";
        }
    }
    else {
        if (segn == true){
            uscente.Vertici.push_back({idV,start});
            entrante.Vertici.push_back({idV,start});
            // PolygonalLibrary::Add_Vert_to_Mesh(mesh, idV, start);
            cout << idV << " ";
        }
    }
    cout << endl;
    uscente.numVertici = uscente.Vertici.size();
    entrante.numVertici = entrante.Vertici.size();
    cout << "e: " << entrante.numVertici << "       u: " << uscente.numVertici << endl;
    cout << "uscente: ";
    for(unsigned int l=0; l<uscente.numVertici; l++){
        cout << uscente.Vertici[l].first << " ";
    }
    cout << endl;
    cout << "entrante: ";
    for(unsigned int l=0; l<entrante.numVertici; l++){
        cout << entrante.Vertici[l].first << " ";
    }
    cout << endl;


    for(auto k = frattura.estremi.begin(); k != frattura.estremi.end(); k++){  //itero su tutte le tracce
        unsigned int key = k->first;
        if(key == id_tr){
            continue; }
        AnalizzaTraccia(start, end, frattura, key, uscente, entrante, Direzione_Uscente, Tracce_SottoPoligoni); //struttura trace delle tracce che corrispondono a questo id
    }
    idV += 2;
    Sotto_poligoni.erase(id_sott);
    idSP ++;
    Sotto_poligoni.insert({idSP, uscente});
    idSP ++;
    Sotto_poligoni.insert({idSP, entrante});
    idSP ++;
}

}




