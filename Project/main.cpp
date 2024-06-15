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
    vector<Trace> elenco_tracce;
    unsigned int n_key = CollectionFractures.size();
    for (unsigned int idF1 = 0; idF1 < n_key; idF1++){
        for(unsigned int idF2 = (idF1+1); idF2<n_key;idF2++){
            if(IntersezioneSfere(CollectionFractures[idF1], CollectionFractures[idF2])){
                cout << "Le sfere delle fratture " << idF1 << " e " << idF2 << " si intersecano" << endl;
                if(Find_Trace(trace, idT, CollectionFractures[idF1], CollectionFractures[idF2])){
                    trace.id1 = idF1;
                    trace.id2 = idF2;
                elenco_tracce.push_back(trace);
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
        //contiamo il numero delle tracce associate alla frattura
        //accediamo alla struttura Fratture e alle mappe di passanti e non passanti
        FileFracture << numFratture << " " << CollectionFractures[numFratture].traccePassanti.size() + CollectionFractures[numFratture].tracceNonPassanti.size() << endl;
        FileFracture << "#TraceId; Tips; Length" << endl;
        //vediamo se sta in passanti o non passanti
        tips = true;
        OutputSort(CollectionFractures[numFratture].traccePassanti, elenco_tracce, FileFracture, tips);
        tips = false;
        OutputSort(CollectionFractures[numFratture].tracceNonPassanti, elenco_tracce, FileFracture, tips);
    }



    /* per la seconda parte
    cicliamo sulle fratture */
    list<SottoPoligoni> Sotto_poligoni;
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
        for(auto iteratore = primo.Passanti.begin(); iteratore != primo.Passanti.end(); iteratore ++){
            primo.estremi[iteratore] = elenco_tracce[iteratore].Vertices;
        }
        for(auto iteratore = primo.NonPassanti.begin(); iteratore != primo.NonPassanti.end(); iteratore ++){
            primo.estremi[iteratore] = elenco_tracce[iteratore].Vertices;
        }
        //dobbiamo aggiungere alla lista dei sottopoligoni primo
        Sotto_poligoni.push_back(primo);
        //facciamo una funzione a parte per il primo taglio


        //FARE DIVIDIPOLIGONO PER PRIMO



        for (unsigned int i=0; i<primo.Passanti.size(); i++){
            for (auto s_pol = Tracce_SottoPoligoni[i].begin(); s_pol != Tracce_SottoPoligoni[i].end(); s_pol++){
                //flag       NB da passare anche in DividiPoligono
                DividiPoligono(i, primo, Sotto_poligoni, Tracce_SottoPoligoni);
            }
        }


    }


    // Close the output file stream
    FileTrace.close();
    FileFracture.close();
    return 0;
}
