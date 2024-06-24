#include "Fractures.hpp"
#include "SottoPoligoni.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace GeometryLibrary;


int main()
{
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

    unsigned int idT = 0;
    unordered_map<unsigned int, Trace> elenco_tracce;
    unsigned int n_key = CollectionFractures.size();
    for (unsigned int idF1 = 0; idF1 < n_key; idF1++){
        for(unsigned int idF2 = (idF1+1); idF2<n_key;idF2++){
            if(IntersezioneSfere(CollectionFractures[idF1], CollectionFractures[idF2])){
                cout << "Le sfere delle fratture " << idF1 << " e " << idF2 << " si intersecano" << endl;
                if(Find_Trace(trace, idT, CollectionFractures[idF1], CollectionFractures[idF2])){
                    trace.id1 = idF1;
                    trace.id2 = idF2;
                elenco_tracce[idT] = trace;
                idT += 1;}
            }
        }
    }

    FileTrace<<"#Number of Traces"<<endl;
    FileTrace<< idT <<endl;
    FileTrace<<"#TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;
    for(unsigned int i = 0; i < idT ; i++ ){
        FileTrace<< i << " "<<elenco_tracce[i].id1<< " "<<elenco_tracce[i].id2<< " "<<elenco_tracce[i].Vertices.first[0]<<" "<<elenco_tracce[i].Vertices.first[1]<<" "<<
            elenco_tracce[i].Vertices.first[2]<<" "<<elenco_tracce[i].Vertices.second[0]<<" "<<elenco_tracce[i].Vertices.second[1]<<" "<<elenco_tracce[i].Vertices.second[2]<<endl;
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

    // Matrix3Xd matrice = (MatrixXd(3,4) << 3,2,2,3,
    //                      2.6,2.6,2.6,2.6,
    //                      1,1,2,2).finished();
    // Vector3d vettore = {1,1,1};
    // VectorXd risult = matrice.col(1)-vettore;
    // for(auto& i:risult){
    //     cout << i << endl;
    // }





    // !MODO 1!
    /* per la seconda parte
    cicliamo sulle fratture */
    map<unsigned int, SottoPoligoni> Sotto_poligoni2;
    unsigned int idSP = 0;

    for (unsigned int i = 0; i < CollectionFractures.size() ; i++){
        map<unsigned int, list<unsigned int>> Tracce_SottoPoligoni; //lista che associa ad ogni traccia i sottopoligoni che la toccano; verrà aggiornata dopo AnalizzaTraccia
        //prima divisione
        SottoPoligoni primo; //adattiamo la funzione dividi poligono anche per il primo taglio
        primo.id = 0;
        primo.Vertici = CollectionFractures[i].Vertici;
        primo.Passanti = CollectionFractures[i].traccePassanti;
        primo.NonPassanti = CollectionFractures[i].tracceNonPassanti;
        primo.numVertici = CollectionFractures[i].numVertici;
        //manca la mappa di estremi e per poter accedere agli elementi usiamo Trace
        cout << "xd " << primo.Passanti.size() << endl;

        for(unsigned int i = 0; i < primo.Passanti.size(); i++){
            primo.estremi.insert({primo.Passanti[i], elenco_tracce[primo.Passanti.at(i)].Vertices});
        }
        for(unsigned int i = 0; i < primo.NonPassanti.size(); i++){
            primo.estremi.insert({primo.NonPassanti[i], elenco_tracce[primo.NonPassanti.at(i)].Vertices});
        }

        Sotto_poligoni2.insert({idSP, primo});    //aggiungiamo alla lista dei sottopoligoni primo

        string flag_p = "passanti";
        string flag_np = "nonpassanti";

        //FARE DIVIDIPOLIGONO PER PRIMO
        DividiPoligono(primo.Passanti[0], primo, Sotto_poligoni2, Tracce_SottoPoligoni, flag_p, idSP);

        for (unsigned int i=1; i<primo.Passanti.size(); i++){
            for (auto& id_sott : Tracce_SottoPoligoni[primo.Passanti[i]]){
                DividiPoligono(i, Sotto_poligoni2[id_sott], Sotto_poligoni2, Tracce_SottoPoligoni, flag_p, idSP);
            }
        }
        for (unsigned int i=0; i<primo.NonPassanti.size(); i++){
            for (auto& id_sott : Tracce_SottoPoligoni[primo.NonPassanti[i]]){
                DividiPoligono(i, Sotto_poligoni2[id_sott], Sotto_poligoni2, Tracce_SottoPoligoni, flag_np, idSP);
            }
        }

    }


    // // !MODO 2!
    // /* per la seconda parte
    // cicliamo sulle fratture */
    // list<SottoPoligoni> Sotto_poligoni;
    // map<unsigned int, list<unsigned int>> Tracce_SottoPoligoni; //lista che associa ad ogni traccia i sottopoligoni che la toccano; verrà aggiornata dopo AnalizzaTraccia

    // for (unsigned int i = 0; i < CollectionFractures.size() ; i++){
    //     //prima divisione
    //     SottoPoligoni primo; //adattiamo la funzione dividi poligono anche per il primo taglio
    //     primo.id = 0;
    //     primo.Vertici = CollectionFractures[i].Vertici;
    //     primo.Passanti = CollectionFractures[i].traccePassanti;
    //     primo.NonPassanti = CollectionFractures[i].tracceNonPassanti;
    //     primo.numVertici = CollectionFractures[i].numVertici;
    //     //manca la mappa di estremi e per poter accedere agli elementi usiamo Trace
    //     for(unsigned int i = 0; i < primo.Passanti.size(); i++){
    //         primo.estremi.insert({primo.Passanti[i], elenco_tracce[primo.Passanti.at(i)].Vertices});
    //     }
    //     for(unsigned int i = 0; i < primo.NonPassanti.size(); i++){
    //         primo.estremi.insert({primo.NonPassanti[i], elenco_tracce[primo.Passanti.at(i)].Vertices});
    //     }

    //     Sotto_poligoni.push_back(primo);    //aggiungiamo alla lista dei sottopoligoni primo
    // }

    // string flag_p = "passanti";
    // string flag_np = "nonpassanti";
    // for (SottoPoligoni sott = Sotto_poligoni.begin(); sott != Sotto_poligoni.end(); sott++){
    //     if(!sott.Passanti.empty()){
    //         DividiPoligono(sott.Passanti[0], sott, Sotto_poligoni, Tracce_SottoPoligoni, flag_p);
    //     }
    //     else if (!sott.NonPassanti.empty()){
    //         DividiPoligono(sott.NonPassanti[0], sott, Sotto_poligoni, Tracce_SottoPoligoni, flag_np);
    //     }
    //     else {
    //         // Non ha più tracce -> lo salvo nella mesh
    //     }
    // }





    // Close the output file stream
    FileTrace.close();
    FileFracture.close();
    return 0;
}
