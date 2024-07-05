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
        // (primo,idstart, idV, CollectionFractures, idP, )
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

    // for(auto spol : Sotto_poligoni){
    //     SottoPoligoni polig = spol.second;
    //     for(auto lato : polig.Lati){
    //         unsigned int lid = lato.first;
    //         unsigned int start = lato.second.first;
    //         unsigned int end = lato.second.second;
    //         bool start_trovato = false;
    //         vector<pair<unsigned int, pair<unsigned int, unsigned int>>> figli = listaFigli[lid];
    //         for(unsigned int f_ind = 0; f_ind < figli.size(); f_ind++){
    //             if(figli[f_ind].second.first == start){
    //                 vertices.push_back(figli[f_ind].second.first);
    //                 edges.push_back(figli[f_ind].first);
    //                 for(unsigned int f_ind2 = f_ind+1; f_ind2 < figli.size(); f_ind++){
    //                     if(figli[f_ind2].second.first != end){
    //                         vertices.push_back(figli[f_ind2].second.first);
    //                         edges.push_back(figli[f_ind2].first);
    //                     }
    //                     else{
    //                         break;
    //                     }
    //                 }
    //                 break;
    //             }
    //             else if(figli[f_ind].second.first == end){
    //                 for(unsigned int f_ind2 = figli.size()-1; f_ind2 >= f_ind; f_ind2--){
    //                     if(start_trovato == false && figli[f_ind2].second.second == start){
    //                         start_trovato = true;
    //                     }
    //                     else if(start_trovato == true){
    //                         vertices.push_back(figli[f_ind2].second.second);
    //                         edges.push_back(figli[f_ind2].first);
    //                     }
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }
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




















    // for(auto spol : Sotto_poligoni){
    //     SottoPoligoni pol = spol.second;
    //     for(auto l : pol.Lati){
    //         // id lato = l.first
    //         //listaFigli[l.first] = vector
    //         for(unsigned int idfi=0; idfi<listaFigli[l.first].size(); idfi++){
    //             if(listaFigli[l.first][idfi].second.first == l.second.first){
    //                 edges.push_back(listaFigli[l.first][idfi].first);
    //                 vertices.push_back(listaFigli[l.first][idfi].second.first);
    //                 if(listaFigli[l.first][idfi].second.second == l.second.second){
    //                     break;
    //                 }
    //                 else{
    //                     for(unsigned int h = idfi+1; h<listaFigli[l.first].size(); h++){
    //                         edges.push_back(listaFigli[l.first][h].first);
    //                         vertices.push_back(listaFigli[l.first][h].second.first);
    //                         if(listaFigli[l.first][h].second.second == l.second.second){
    //                             break;
    //                         }
    //                     }
    //                 }
    //             }
    //             else if(listaFigli[l.first][idfi].second.second == l.second.second){
    //                 for(unsigned int opp = listaFigli[l.first].size()-1; opp>idfi; opp--){
    //                     if(listaFigli[l.first][opp].second.second == l.second.first){
    //                         edges.push_back(listaFigli[l.first][opp].first);
    //                         vertices.push_back(listaFigli[l.first][opp].second.second);
    //                         if(listaFigli[l.first][opp].second.first == l.second.second){
    //                             break;
    //                         }
    //                         else{
    //                             for(unsigned int h = opp-1; h<idfi; h++){
    //                                 edges.push_back(listaFigli[l.first][h].first);
    //                                 vertices.push_back(listaFigli[l.first][h].second.second);
    //                                 if(listaFigli[l.first][h].second.first == l.second.second){
    //                                     break;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     mesh.EdgesCell2D.push_back(edges);
    //     mesh.VerticesCell2D.push_back(vertices);
    //     mesh.NumberElements2D.push_back(vertices.size());
    //     edges = {};
    //     vertices = {};
    // }







    // VectorXi edges;
    // VectorXi vertices;
    // for(auto spol = Sotto_poligoni.begin(); spol != Sotto_poligoni.end(); spol++){
    //     SottoPoligoni pol = spol -> second;
    //     for(unsigned int posizIdL=0; posizIdL<pol.Lati.size(); posizIdL++){     // Scorro su Lati del poligono
    //         for(unsigned int j=0; j<listaFigli[pol.Lati[posizIdL].first].size(); j++){  // Scorro sul vettore di listaFigli associato all'idLato studiato
    //             if (listaFigli[pol.Lati[posizIdL].first][j].second.first == pol.Lati[posizIdL].second.first){
    //                 edges[posizIdL + j] = listaFigli[pol.Lati[posizIdL].first][j].first;
    //                 vertices[posizIdL + j] = listaFigli[pol.Lati[posizIdL].first][j].second.first;
    //                 if(listaFigli[pol.Lati[posizIdL].first][j].second.second == pol.Lati[posizIdL].second.second){
    //                     vertices[posizIdL + j] = listaFigli[pol.Lati[posizIdL].first][j].second.second;
    //                     break;
    //                 }
    //                 else{
    //                     for(unsigned int k=j+1; k<listaFigli[pol.Lati[posizIdL].first].size(); k++){
    //                         edges[posizIdL + j] = listaFigli[pol.Lati[posizIdL].first][k].first;
    //                         vertices[posizIdL + j] = listaFigli[pol.Lati[posizIdL].first][k].second.first;
    //                         if(listaFigli[pol.Lati[posizIdL].first][k].second.second == pol.Lati[posizIdL].second.second){
    //                             vertices[posizIdL + j] = listaFigli[pol.Lati[posizIdL].first][k].second.second;
    //                             break;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     edges = {};
    //     vertices = {};
    // }



    cout << "IdCell0D:  ";
    for (unsigned int i=0; i<mesh.IdCell0D.size(); i++){
        cout << mesh.IdCell0D[i] << " ";
    }
    cout << endl;

    cout << "IdCell1D:  ";
    for (unsigned int i=0; i<mesh.IdCell1D.size(); i++){
        cout << mesh.IdCell1D[i] << " ";
    }
    cout << endl;

    cout << "NumbElemCell2d:    ";
    for(unsigned int i=0; i<mesh.NumberElements2D.size(); i++){
        cout << mesh.NumberElements2D[i] << " ";
    }
    cout << endl;

    cout << "VerticesCell2d" << endl;
    for(unsigned int i=0; i<mesh.VerticesCell2D.size(); i++){
        for(unsigned int j=0; j<mesh.VerticesCell2D[i].size(); j++){
            cout << mesh.VerticesCell2D[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "EdgesCell2d" << endl;
    for(unsigned int i=0; i<mesh.EdgesCell2D.size(); i++){
        for(unsigned int j=0; j<mesh.EdgesCell2D[i].size(); j++){
            cout << mesh.EdgesCell2D[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    for (auto k : mappaLati){
        cout << k.first << ": ";
        for (auto kk : k.second){
            cout << kk << " ";
        }
        cout << endl;
    }



    cout << "ListaFigli" << endl;
    for(auto k : listaFigli){
        cout << k.first;
        for(auto v : k.second){
            cout << "     " << v.first << ":  " << v.second.first << ", " << v.second.second << endl;
        }
        cout << endl;
    }
    cout << endl;



    // Close the output file stream
    FileTrace.close();
    FileFracture.close();
    return 0;
}
