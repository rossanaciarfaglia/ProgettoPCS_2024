#include "Fractures.hpp"
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
            if(IntersezioneSfere(fracture, CollectionFractures[idF1].Vertici, CollectionFractures[idF2].Vertici)){
                cout << "Le sfere delle fratture " << idF1 << " e " << idF2 << " si intersecano" << endl;
                if(Find_Trace(fracture, trace, idT, CollectionFractures[idF1], CollectionFractures[idF2])){
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
        FileFracture<<"#FractureId; NumTraces"<<endl;
        //contiamo il numero delle tracce associate alla frattura
        //accediamo alla struttura Fratture e alle mappe di passanti e non passanti
        FileFracture<<numFratture<<" "<<CollectionFractures[numFratture].traccePassanti.size() + CollectionFractures[numFratture].tracceNonPassanti.size()<<endl;
        FileFracture<<"#TraceId; Tips; Length"<<endl;
        //vediamo se sta in passanti o non passanti
        tips = false;
        OutputSort(CollectionFractures[numFratture].traccePassanti, elenco_tracce, FileFracture, tips);
        tips = true;
        OutputSort(CollectionFractures[numFratture].tracceNonPassanti, elenco_tracce, FileFracture, tips);
    }


    // Call the Passante_NonPassante function
    //auto result = fractures.Passante_NonPassante(outputFile);

    // Close the output file stream
    FileTrace.close();
    FileFracture.close();
    return 0;
}
