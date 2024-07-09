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
    string filepath = "./DFN/FR10_data.txt";
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


    /* per la seconda parte
    cicliamo sulle fratture */
    map<unsigned int, SottoPoligoni> Sotto_poligoni;
    unsigned int idSP = 0;
    map<unsigned int, list<unsigned int>> Tracce_SottoPoligoni; //lista che associa ad ogni traccia i sottopoligoni che la toccano; verrà aggiornata dopo AnalizzaTraccia

    unsigned int idstart;
    for (unsigned int idP=0; idP < CollectionFractures.size(); idP++){
        SottoPoligoni primo;
        Convertitore_struct(primo, idSP, idstart, idV, idL, elenco_tracce, CollectionFractures, idP, Sotto_poligoni, mappaLati, mesh);

        string flag_p = "passanti";
        string flag_np = "nonpassanti";       

        if (primo.Passanti.size() != 0){
            Tracce_SottoPoligoni[primo.Passanti[0]].push_back(primo.id);
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
            idSP++;
        }
        Tracce_SottoPoligoni = {};
    }


    map<unsigned int, vector<pair<unsigned int, pair<unsigned int, unsigned int>>>> listaFigli;
    unsigned int idF = 0;
    unsigned int idPa;
    vector<unsigned int> Pa;
    for(auto padre : mappaLati){
        idPa = padre.first;
        listaFigli.insert({idPa, {}});
        Pa = padre.second;
        vector<unsigned int> ordinamento(Pa.size(),0);
        Ordina_Punti(Pa, mesh.CoordinatesCell0D, ordinamento);

        for(unsigned int k = 0; k < ordinamento.size()-1; k++){
            listaFigli[idPa].push_back({idF, {Pa[ordinamento[k]], Pa[ordinamento[k+1]]}});
            idF++;
        }
    }

    for(auto p : listaFigli){
        for(auto f : p.second){
            mesh.IdCell1D.push_back(f.first);
            mesh.VerticesCell1D.push_back({f.second.first, f.second.second});
        }
    }

    vector<unsigned int> edges;
    vector<unsigned int> vertices;
    unsigned int lid, start, end;
    bool start_trovato;

    for(auto spol : Sotto_poligoni){
        SottoPoligoni polig = spol.second;

        for(auto lato : polig.Lati){
            lid = lato.first;
            start = lato.second.first;
            end = lato.second.second;
            start_trovato = false;

            vector<pair<unsigned int, pair<unsigned int, unsigned int>>> figli = listaFigli[lid];
            for(int f_ind = 0; f_ind < figli.size(); f_ind++){

                if(start_trovato == false && figli[f_ind].second.first == start){
                    start_trovato = true;
                }
                else if(start_trovato == false && figli[f_ind].second.first == end){
                    for(int f_ind2 = figli.size()-1; f_ind2 >= f_ind; f_ind2--){
                        if(start_trovato == false && figli[f_ind2].second.second == start){
                            start_trovato = true;
                        }
                        if(start_trovato == true){
                            vertices.push_back(figli[f_ind2].second.second);
                            edges.push_back(figli[f_ind2].first);
                        }
                    }
                    break;
                }

                if(start_trovato == true){
                    vertices.push_back(figli[f_ind].second.first);
                    edges.push_back(figli[f_ind].first);
                    if(figli[f_ind].second.second == end)
                        break;
                }
            }
        }
        mesh.EdgesCell2D.push_back(edges);
        mesh.VerticesCell2D.push_back(vertices);
        mesh.NumberElements2D.push_back(vertices.size());
        edges = {};
        vertices = {};
    }


    // Close the output file stream
    FileTrace.close();
    FileFracture.close();
    return 0;


<<<<<<< HEAD
    // Per esportare in Paraview
    string exportFolder = "./";
=======
     // Per esportare in Paraview
    string exportFolder = "C:/Users/User/OneDrive/Desktop/ProgettoPCS_2024/debug/";
>>>>>>> fb0f8bdbd0108d3af9a02158139dcc5eaa0767bc
    Gedim::UCDUtilities exporter;

    MatrixXd points(3, mesh.CoordinatesCell0D.size());
    for (unsigned int p = 0; p < mesh.CoordinatesCell0D.size(); p++) {
        points(0,p) = mesh.CoordinatesCell0D[p][0] ;
        points(1,p) = mesh.CoordinatesCell0D[p][1] ;
        points(2,p) = mesh.CoordinatesCell0D[p][2];
    };

    exporter.ExportPoints(exportFolder + "/Celle0Ds.inp",
                          points);

    MatrixXi edge(2,mesh.VerticesCell1D.size());
    for (unsigned int  i = 0; i < mesh.VerticesCell1D.size(); i++) {
        edge(0,i) = mesh.VerticesCell1D[i][0];
        edge(1,i) = mesh.VerticesCell1D[i][1];
    };

    exporter.ExportSegments(exportFolder + "/Geometry1Ds.inp",
                            points,
                            edge,
                            {},
                            {},
                            {});


}
