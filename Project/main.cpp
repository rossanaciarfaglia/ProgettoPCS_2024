#include "Fractures.hpp"
#include "SottoPoligoni.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include "PolygonalMesh.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace GeometryLibrary;
using namespace PolygonalLibrary;


int main() {
    bool tips;
    string filepath = "./DFN/FR3_data.txt";
    Fracture fracture;
    unordered_map<unsigned int, Fracture> CollectionFractures; // Il costo computazionale è O(1), non O(logn)


    ImportFracturesList(filepath, fracture, CollectionFractures);
    Trace trace;
    ofstream FileTrace("FileTrace.txt");
    if (!FileTrace.is_open()) {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    PolygonalMesh mesh;
    unsigned int idV = 0;
    unsigned int idL = 0;
    map<unsigned int, vector<unsigned int>> mappaLati;

    unordered_map<unsigned int, Trace> elenco_tracce;
    unsigned int n_key = CollectionFractures.size();
    for (unsigned int idF1 = 0; idF1 < n_key; idF1++){
        for(unsigned int idF2 = (idF1+1); idF2<n_key;idF2++){
            if(IntersezioneSfere(CollectionFractures[idF1], CollectionFractures[idF2])){
                cout << "Le sfere delle fratture " << idF1 << " e " << idF2 << " si intersecano" << endl;
                if(Find_Trace(trace, idL, CollectionFractures[idF1], CollectionFractures[idF2], idV, mesh)){
                    trace.id1 = idF1;
                    trace.id2 = idF2;
                    elenco_tracce[idL] = trace;
                    mappaLati.insert({idL,{trace.Vertices.first.first, trace.Vertices.second.first}});
                    idL++;
                }
            }
        }
    }


    FileTrace<<"#Number of Traces"<<endl;
    FileTrace<< idL <<endl;
    FileTrace<<"#TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;
    for(unsigned int i = 0; i < idL ; i++ ){
        FileTrace<< i << " "<<elenco_tracce[i].id1<< " "<<elenco_tracce[i].id2<< " "<<elenco_tracce[i].Vertices.first.second[0]<<" "<<elenco_tracce[i].Vertices.first.second[1]<<" "<<
            elenco_tracce[i].Vertices.first.second[2]<<" "<<elenco_tracce[i].Vertices.second.second[0]<<" "<<elenco_tracce[i].Vertices.second.second[1]<<" "<<elenco_tracce[i].Vertices.second.second[2]<<endl;
    }

    //SECONDO FILE OUTPUT
    ofstream FileFracture("FileFratture.txt");
    if (!FileFracture.is_open()) {
        cerr << "Error opening output file." << endl;
    }
    for(unsigned int numFratture = 0; numFratture < n_key; numFratture ++){
        FileFracture << "#FractureId; NumTraces" << endl;
        /* contiamo il numero delle tracce associate alla frattura
        accediamo alla struttura Fratture e alle mappe di passanti e non passanti */
        FileFracture << numFratture << " " << CollectionFractures[numFratture].traccePassanti.size() + CollectionFractures[numFratture].tracceNonPassanti.size() << endl;
        FileFracture << "#TraceId; Tips; Length" << endl;
        //vediamo se sta in passanti o non passanti
        tips = true;
        OutputSort(CollectionFractures[numFratture].traccePassanti, elenco_tracce, FileFracture, tips);
        tips = false;
        OutputSort(CollectionFractures[numFratture].tracceNonPassanti, elenco_tracce, FileFracture, tips);
    }
    cout << "seconda parte" << endl;


    /* per la seconda parte
    cicliamo sulle fratture */
    map<unsigned int, SottoPoligoni> Sotto_poligoni;
    unsigned int idSP = 0;
    map<unsigned int, list<unsigned int>> Tracce_SottoPoligoni; //lista che associa ad ogni traccia i sottopoligoni che la toccano; verrà aggiornata dopo AnalizzaTraccia

    unsigned int idstart;
    for (unsigned int idP=0; idP < CollectionFractures.size(); idP++){
        SottoPoligoni primo;
        primo.id = 0;
        idstart = idV;
        for (unsigned int v=0; v<CollectionFractures[idP].numVertici; v++){
            primo.Vertici.push_back({idV, CollectionFractures[idP].Vertici.col(v)});
            primo.Lati.push_back({idL,{idV, idstart + (idV+1-idstart)%CollectionFractures[idP].numVertici}});
            Add_Vert_to_Mesh(mesh, {idV, CollectionFractures[idP].Vertici.col(v)});
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

        string flag_p = "passanti";
        string flag_np = "nonpassanti";

        cout << "ex sottopol: " << idSP << endl;
        for (unsigned int h=0; h<primo.numVertici; h++){
            cout << primo.Vertici[h].first << " ";
        }
        cout << endl;

        if (primo.Passanti.size() != 0){
            DividiPoligono(primo.Passanti[0], primo, idSP, Sotto_poligoni, Tracce_SottoPoligoni, flag_p, idSP, idV, mesh, mappaLati);     //DividiPoligono per la prima frattura

            for (unsigned int j=1; j<primo.Passanti.size(); j++){
                for (auto& id_sott : Tracce_SottoPoligoni[primo.Passanti[j]]){
                    DividiPoligono(primo.Passanti[j], Sotto_poligoni[id_sott], id_sott, Sotto_poligoni, Tracce_SottoPoligoni, flag_p, idSP, idV, mesh, mappaLati);
                }
            }
            for (unsigned int k=0; k<primo.NonPassanti.size(); k++){
                for (auto& id_sott : Tracce_SottoPoligoni[primo.NonPassanti[k]]){
                    DividiPoligono(primo.NonPassanti[k], Sotto_poligoni[id_sott], id_sott, Sotto_poligoni, Tracce_SottoPoligoni, flag_np, idSP, idV, mesh, mappaLati);
                }
            }
        }
        else if (primo.NonPassanti.size() != 0) {
            DividiPoligono(primo.NonPassanti[0], primo, idSP, Sotto_poligoni, Tracce_SottoPoligoni, flag_np, idSP, idV, mesh, mappaLati);       //DividiPoligono per la prima frattura

            for (unsigned int k=1; k<primo.NonPassanti.size(); k++){
                for (auto& id_sott : Tracce_SottoPoligoni[primo.NonPassanti[k]]){
                    DividiPoligono(primo.NonPassanti[k], Sotto_poligoni[id_sott], id_sott, Sotto_poligoni, Tracce_SottoPoligoni, flag_np, idSP, idV, mesh, mappaLati);
                }
            }
        }
        else{
            cout << "No tracce per poligono " << idSP << endl;
            idSP++;
        }
    }
    cout << "!" << endl;

    map<unsigned int, vector<pair<unsigned int, pair<unsigned int, unsigned int>>>> listaFigli;
    unsigned int idPa;
    vector<unsigned int> Pa;
    for (auto padre : mappaLati){
        idPa = padre.first;
        Pa = padre.second;

    }



    cout << "IdCell0D:  ";
    for (unsigned int i=0; i<mesh.IdCell0D.size(); i++){
        cout << mesh.IdCell0D[i] << " ";
    }
    cout << endl;

    for (auto k : mappaLati){
        cout << k.first << ": ";
        for (auto kk : k.second){
            cout << kk << " ";
        }
        cout << endl;
    }



    // Close the output file stream
    FileTrace.close();
    FileFracture.close();
    return 0;
}
