#include "Fractures.hpp"
#include "SottoPoligoni.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <utility>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {

void AnalizzaTraccia(pair<unsigned int,Vector3d>& start_taglio, pair<unsigned int,Vector3d>& end_taglio, SottoPoligoni& taglio, unsigned int& id_traccia, SottoPoligoni& uscente, SottoPoligoni& entrante, Vector3d& VettoreUscente, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni, unsigned int& idV){
    //prendo gli estremi della traccia
    pair<unsigned int,Vector3d> start = taglio.estremi[id_traccia].first;
    pair<unsigned int,Vector3d> end = taglio.estremi[id_traccia].second;

    Tracce_SottoPoligoni[id_traccia].remove(taglio.id);

    unsigned int count = 0; //per valutare se la traccia appartiene a un solo sottopoligono o a entrambi
    bool sovrapp_start = false;
    bool sovrapp_end = false;

    if(Regola_Mano_Destra(end_taglio.second - start_taglio.second, start.second - start_taglio.second, VettoreUscente) == 1) {       //uscente
        uscente.estremi[id_traccia].first = start;
        count ++;}
    else if(Regola_Mano_Destra(end_taglio.second - start_taglio.second, start.second - start_taglio.second, VettoreUscente) == 0)  //entrante
        entrante.estremi[id_traccia].first = start;
    else {
        //start è sulla traccia di taglio
        sovrapp_start = true;
    }

    //lo analizzo con l'end
    if(Regola_Mano_Destra(end_taglio.second - start_taglio.second, end.second - start_taglio.second, VettoreUscente) == 1) {       //uscente
        uscente.estremi[id_traccia].second = end;
        count ++;}
    else if(Regola_Mano_Destra(end_taglio.second - start_taglio.second, end.second - start_taglio.second, VettoreUscente) == 0)  //entrante
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
            /* Le tracce si intersecano, quindi dobbiamo trovare l'intersezione
            In questo caso allo stesso id_traccia saranno associate due tracce "diverse" a seconda che sono in uscente o entrante (NON nello stesso) */
            Vector2d coeff = CoefficientiRette(start.second, end.second, start_taglio.second, end_taglio.second-start_taglio.second);
            Vector3d intersezione = start.second + coeff[0]*(end.second - start.second); //abbiamo trovato l'intersezione
            if (trovato){
                //devo inserire il punto di intersezione
                if ((uscente.estremi[id_traccia].first.second - start.second).norm() <= 1e-14){
                    //significa che lo start sta in uscente
                    uscente.estremi[id_traccia].second = {idV,intersezione};
                    entrante.estremi[id_traccia].first = {idV,intersezione};
                }
                else {
                    uscente.estremi[id_traccia].first = {idV,intersezione};
                    entrante.estremi[id_traccia].second = {idV,intersezione};
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
                    if (Regola_Mano_Destra(uscente.Vertici[(i+1)%uscente.numVertici].second - uscente.Vertici[i].second, start.second - uscente.Vertici[i].second, VettoreUscente) == 2 ||
                        Regola_Mano_Destra(uscente.Vertici[(i+1)%uscente.numVertici].second - uscente.Vertici[i].second, end.second - uscente.Vertici[i].second, VettoreUscente) == 2){
                        if ((uscente.estremi[id_traccia].first.second - start.second).norm() <= 1e-14){
                            //significa che lo start sta in uscente
                            uscente.estremi[id_traccia].second = {idV,intersezione};
                            entrante.estremi[id_traccia].first = {idV,intersezione};
                        }
                        else {
                            uscente.estremi[id_traccia].first = {idV,intersezione};
                            entrante.estremi[id_traccia].second = {idV,intersezione};
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
                        if (Regola_Mano_Destra(entrante.Vertici[(i+1)%entrante.numVertici].second - entrante.Vertici[i].second, start.second - entrante.Vertici[i].second, VettoreUscente) == 2 ||
                            Regola_Mano_Destra(entrante.Vertici[(i+1)%entrante.numVertici].second - entrante.Vertici[i].second, end.second - entrante.Vertici[i].second, VettoreUscente) == 2){
                            if ((entrante.estremi[id_traccia].first.second - start.second).norm() <= 1e-14){
                                //significa che lo start sta in entrante
                                entrante.estremi[id_traccia].second = {idV,intersezione};
                                uscente.estremi[id_traccia].first = {idV,intersezione};
                            }
                            else {
                                entrante.estremi[id_traccia].first = {idV,intersezione};
                                uscente.estremi[id_traccia].second = {idV,intersezione};
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
                    if ((entrante.estremi[id_traccia].first.second - start.second).norm() <= 1e-14){
                        //significa che lo start sta in uscente
                        entrante.estremi[id_traccia].second = {idV,intersezione};
                        uscente.estremi[id_traccia].first = {idV,intersezione};
                    }
                    else {
                        entrante.estremi[id_traccia].first = {idV,intersezione};
                        uscente.estremi[id_traccia].second = {idV,intersezione};
                    }
                    Tracce_SottoPoligoni[id_traccia].push_back(entrante.id);
                    entrante.NonPassanti.push_back(id_traccia);
                    Tracce_SottoPoligoni[id_traccia].push_back(uscente.id);
                    uscente.NonPassanti.push_back(id_traccia);
                    flag = true;
                }
            }
            idV++;
            break;
        }
    }
}



void DividiPoligono(unsigned int& id_tr, SottoPoligoni& frattura, unsigned int& id_sott, map<unsigned int, SottoPoligoni>& Sotto_poligoni, map<unsigned int, list<unsigned int>>& Tracce_SottoPoligoni, const string& flag, unsigned int& idSP, unsigned int& idV, PolygonalLibrary::PolygonalMesh& mesh, map<unsigned int, vector<unsigned int>>& mappaLati){
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
    pair<unsigned int, Vector3d> start;
    pair<unsigned int, Vector3d> end;
    cout << "Nuovo poligono" << endl;
    if (flag == "passanti"){
        start = frattura.estremi[id_tr].first;
        end = frattura.estremi[id_tr].second;
    }
    else if(flag == "nonpassanti"){
        //allunga traccia
        // start ed end sono i punti di intersezione
        Matrix<double,2,3> traccia;
        traccia.row(0) << frattura.estremi[id_tr].first.second[0], frattura.estremi[id_tr].first.second[1], frattura.estremi[id_tr].first.second[2];
        traccia.row(1) = frattura.estremi[id_tr].second.second - frattura.estremi[id_tr].first.second;
        Matrix3Xd vertici;
        vertici.resize(3, frattura.numVertici);
        for (unsigned int m=0; m<frattura.numVertici; m++){
            vertici.col(m) = frattura.Vertici[m].second;
        }
        vector<Vector3d> intersezioni = Intersection_Point(traccia, vertici, frattura.numVertici);

        if (isLess(intersezioni[1], intersezioni[0], traccia)) swap(intersezioni[1], intersezioni[0]);

        if(intersezioni[0] == frattura.estremi[id_tr].first.second || intersezioni[0] == frattura.estremi[id_tr].second.second){
            start = frattura.estremi[id_tr].first;
            end = {idV, intersezioni[1]};
            PolygonalLibrary::Add_Vert_to_Mesh(mesh, end);
            mappaLati[id_tr].push_back(idV);
            idV++;
        }
        else if (intersezioni[1] == frattura.estremi[id_tr].first.second || intersezioni[1] == frattura.estremi[id_tr].second.second){
            start = {idV, intersezioni[0]};
            end = frattura.estremi[id_tr].second;
            PolygonalLibrary::Add_Vert_to_Mesh(mesh, start);
            mappaLati[id_tr].push_back(idV);
            idV++;
        }
        else {
            start = {idV, intersezioni[0]};
            end = {idV+1, intersezioni[1]};
            PolygonalLibrary::Add_Vert_to_Mesh(mesh, start);
            PolygonalLibrary::Add_Vert_to_Mesh(mesh, end);
            mappaLati[id_tr].push_back(idV);
            mappaLati[id_tr].push_back(idV+1);
            idV += 2;
        }
    }

    for (auto lato : frattura.Lati){
        if (Punto_su_Lato(mesh.CoordinatesCell0D[lato.second.first], mesh.CoordinatesCell0D[lato.second.second], start.second)){
            mappaLati[lato.first].push_back(start.first);

        }
        else if (Punto_su_Lato(mesh.CoordinatesCell0D[lato.second.first], mesh.CoordinatesCell0D[lato.second.second], end.second)){
            mappaLati[lato.first].push_back(end.first);

        }
    }

    bool segno; // tiene segno dell'ultimo prodotto misto per sapere se sono passato dal sottopoligono uscente a quello entrante (o viceversa)
    bool check_trace = false;
    //iteriamo sui lati dei poligoni
    if (Regola_Mano_Destra(end.second - start.second, mesh.CoordinatesCell0D[frattura.Lati[0].second.first] - start.second, Direzione_Uscente) == 1){
        segno = true;
    }
    else if (Regola_Mano_Destra(end.second - start.second, mesh.CoordinatesCell0D[frattura.Lati[0].second.first] - start.second, Direzione_Uscente) == 0){
        segno = false;
    }
    else {

    }
    for (unsigned int l=0; l<frattura.numVertici; l++){
        if (Regola_Mano_Destra(end.second - start.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.second] - start.second, Direzione_Uscente) == 1){
            if (segno == false){
                uscente.Vertici.push_back(end);
                entrante.Vertici.push_back(end);
                if (check_trace == false){
                    uscente.Lati.push_back({id_tr, {start.first, end.first}});
                    uscente.Lati.push_back({frattura.Lati[l].first, {end.first, frattura.Lati[l].second.second}});
                    entrante.Lati.push_back({frattura.Lati[l].first, {frattura.Lati[l].second.first, end.first}});
                    entrante.Lati.push_back({id_tr, {end.first, start.first}});
                    check_trace = true;
                }
                else{
                    uscente.Lati.push_back({frattura.Lati[l].first, {end.first, frattura.Lati[l].second.second}});
                    entrante.Lati.push_back({frattura.Lati[l].first, {frattura.Lati[l].second.first, end.first}});
                }
            }
            else {
                uscente.Lati.push_back(frattura.Lati[l]);
            }
            uscente.Vertici.push_back({frattura.Lati[l].second.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.second]});
            segno = true;
        }

        else if(Regola_Mano_Destra(end.second - start.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.second] - start.second, Direzione_Uscente) == 0) {
            if (segno == true){
                uscente.Vertici.push_back(start);
                entrante.Vertici.push_back(start);
                if (check_trace == false){
                    uscente.Lati.push_back({frattura.Lati[l].first, {frattura.Lati[l].second.first, start.first}});
                    uscente.Lati.push_back({id_tr, {start.first, end.first}});
                    entrante.Lati.push_back({id_tr, {end.first, start.first}});
                    entrante.Lati.push_back({frattura.Lati[l].first, {start.first, frattura.Lati[l].second.second}});
                    check_trace = true;
                }
                else{
                    uscente.Lati.push_back({frattura.Lati[l].first, {frattura.Lati[l].second.first, start.first}});
                    entrante.Lati.push_back({frattura.Lati[l].first, {start.first, frattura.Lati[l].second.second}});
                }
            }
            else {
                entrante.Lati.push_back(frattura.Lati[l]);
            }
            entrante.Vertici.push_back({frattura.Lati[l].second.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.second]});
            segno = false;
        }

        else {
            if(Regola_Mano_Destra(end.second - start.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.first] - start.second, Direzione_Uscente) == 2) {
                cout << "allineati" << endl;
                for (unsigned int i=0; i<frattura.numVertici; i++){
                    cout << frattura.Lati[i].second.first << ": " << mesh.CoordinatesCell0D[frattura.Lati[i].second.first][0] << " " << mesh.CoordinatesCell0D[frattura.Lati[i].second.first][1] << " " << mesh.CoordinatesCell0D[frattura.Lati[i].second.first][2] << "         ";
                }
                cout << endl;
                cout << "problematici: " << frattura.Lati[l].second.first << "   " << frattura.Lati[l].second.second << endl;
                cout << "traccia: " << start.second[0] << " " << start.second[1] << " " << start.second[2] << "       " << end.second[0] << " " << end.second[1] << " " << end.second[2] << " " << endl;
                // Tracce_SottoPoligoni[id_tr].remove(frattura.id);
                return;
            }
            cout << "vertice su traccia" << endl;
            uscente.Vertici.push_back({frattura.Lati[l].second.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.second]});
            entrante.Vertici.push_back({frattura.Lati[l].second.second, mesh.CoordinatesCell0D[frattura.Lati[l].second.second]});
            if(segno == true){
                uscente.Lati.push_back(frattura.Lati[l]);
                uscente.Lati.push_back({id_tr, {start.first, end.first}});
                entrante.Lati.push_back({id_tr, {end.first, start.first}});
            }
            else{
                entrante.Lati.push_back(frattura.Lati[l]);
                entrante.Lati.push_back({id_tr, {end.first, start.first}});
                uscente.Lati.push_back({id_tr, {start.first, end.first}});
            }
            segno = !segno;
        }
    }
    uscente.numVertici = uscente.Vertici.size();
    entrante.numVertici = entrante.Vertici.size();
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
        AnalizzaTraccia(start, end, frattura, key, uscente, entrante, Direzione_Uscente, Tracce_SottoPoligoni, idV); //struttura trace delle tracce che corrispondono a questo id
    }
    Sotto_poligoni.erase(id_sott);
    idSP ++;
    Sotto_poligoni.insert({idSP, uscente});
    idSP ++;
    Sotto_poligoni.insert({idSP, entrante});
    idSP ++;
}


void Convertitore_struct(SottoPoligoni& primo, unsigned int& idSP, unsigned int& idstart, unsigned int& idV, unsigned int& idL, unordered_map<unsigned int, Trace>& elenco_tracce, unordered_map<unsigned int, Fracture>& CollectionFractures, unsigned int& idP, map<unsigned int, SottoPoligoni>& Sotto_poligoni, map<unsigned int, vector<unsigned int>>& mappaLati, PolygonalLibrary::PolygonalMesh& mesh){
    primo.id = idSP;
    idstart = idV;
    for (unsigned int v=0; v<CollectionFractures[idP].numVertici; v++){
        primo.Vertici.push_back({idV, CollectionFractures[idP].Vertici.col(v)});
        primo.Lati.push_back({idL,{idV, idstart + (idV+1-idstart)%CollectionFractures[idP].numVertici}});
        PolygonalLibrary::Add_Vert_to_Mesh(mesh, {idV, CollectionFractures[idP].Vertici.col(v)});
        mappaLati.insert({idL,{idV, idstart + (idV+1-idstart)%CollectionFractures[idP].numVertici}});
        idL++;
        idV++;
    }
    primo.Passanti = CollectionFractures[idP].traccePassanti;
    primo.NonPassanti = CollectionFractures[idP].tracceNonPassanti;
    primo.numVertici = CollectionFractures[idP].numVertici;

    for(unsigned int i = 0; i < primo.Passanti.size(); i++){
        primo.estremi.insert({primo.Passanti[i], elenco_tracce[primo.Passanti[i]].Vertices});
    }
    for(unsigned int i = 0; i < primo.NonPassanti.size(); i++){
        primo.estremi.insert({primo.NonPassanti[i], elenco_tracce[primo.NonPassanti[i]].Vertices});
    }
    Sotto_poligoni.insert({idSP, primo});    //aggiungiamo alla lista dei sottopoligoni primo
}

}




