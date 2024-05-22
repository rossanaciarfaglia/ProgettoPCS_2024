#include "Fractures.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>

using namespace std;
using namespace Eigen;

namespace GeometryLibrary {

// I dati del baricentro sono in un vettore lungo 3
Vector3d Fractures::Baricentro(MatrixXd Poligono) {
    Vector3d baricentro;
    for(unsigned int riga = 0; riga < 3; riga++) {
        double sum = 0;
        for(unsigned int colonna = 0; colonna < Poligono.cols(); colonna++) {
            sum += Poligono(riga, colonna);
        }
        baricentro[riga] = sum / Poligono.cols();      
    }
    return baricentro;
}


// Raggio delle sfere !!AL QUADRATO!!
double Fractures::Raggio(Vector3d baricentro, MatrixXd Poligono) {
    double R = 0;
    double dist;
    for (unsigned int numV=0; numV<Poligono.cols(); numV++){
        dist = 0;
        for (unsigned int i=0; i<3; i++){
            dist += (baricentro[i]-Poligono(i,numV))*(baricentro[i]-Poligono(i,numV));
        }
        if (dist > R){
            R = dist;
        }
    }
    return R;
}


Vector4d Fractures::TrovaPiano(MatrixXd poligono){
    Vector3d u;
    Vector3d v;
    for (unsigned int ax=0; ax<3; ax++){           // (Sono x,y e z)
        u[ax] = poligono(ax,2) - poligono(ax,0);   // u = P2-P0
        v[ax] = poligono(ax,1) - poligono(ax,0);   // v = P1-P0
    }
    double u_norm = u.norm();
    cout << "norma di u: " << u_norm << endl;
    double v_norm = v.norm();

    Vector4d n_d;   // è il vettore normale n + la costante d
    n_d[0] = (u[1]*v[2]-v[1]*u[2])/(u_norm*v_norm);
    n_d[1] = (v[0]*u[2]-u[0]*v[2])/(u_norm*v_norm);
    n_d[2] = (u[0]*v[1]-v[0]*u[1])/(u_norm*v_norm);

    n_d[3] = n_d[0]*poligono(0,0) + n_d[1]*poligono(1,0) + n_d[2]*poligono(2,0);

    return n_d;
}


bool IntersezioneSfere(Fractures &polygons, unsigned int id1, unsigned int id2){
    double tol_quad = 100 * numeric_limits<double>::epsilon() * numeric_limits<double>::epsilon();

    MatrixXd Poly1 = polygons.FracturesMap[id1];
    Vector3d baricentro1 = polygons.Baricentro(Poly1);
    double R1 = polygons.Raggio(baricentro1, Poly1);
    MatrixXd Poly2 = polygons.FracturesMap[id2];
    Vector3d baricentro2 = polygons.Baricentro(Poly2);
    double R2 = polygons.Raggio(baricentro2, Poly2);
    if (DistanzaEuclidea(baricentro1, baricentro2) - (R1+R2+2*sqrt(R1*R2)) < tol_quad){
        return true;
    }
    return false;
}


double DistanzaEuclidea(Vector3d centro1, Vector3d centro2) {
    double distanza = 0;
    for (unsigned i = 0; i < 3; i++) {
        distanza += pow(centro1[i] - centro2[i], 2);
    }
    return distanza;
}

pair<map<unsigned int, vector<unsigned int>>, map<unsigned int, vector<unsigned int>>> Fractures::Passante_NonPassante(ofstream& outputFile) {
    double tol = 10 * numeric_limits<double>::epsilon();
    unsigned int id_Traccia = 0;   //mi serve per il numero dell'id Traccia
    vector<vector<unsigned int>> Lista_Intersezioni_Sfere = SpheresIntersection(); // contiene id delle fratture le cui sfere si intersecano
    map<unsigned int, vector<unsigned int>> Passanti = {};
    map<unsigned int, vector<unsigned int>> NonPassanti = {};

    for(unsigned int i = 0; i < Lista_Intersezioni_Sfere.size(); i++) {
        unsigned int id1 = Lista_Intersezioni_Sfere[i][0];
        unsigned int id2 = Lista_Intersezioni_Sfere[i][1];
        Vector4d Piano1 = TrovaPiano(FracturesMap[id1]);
        Vector4d Piano2 = TrovaPiano(FracturesMap[id2]);

        if (fabs(Piano1[0] - Piano2[0]) < tol && fabs(Piano1[1] - Piano2[1]) < tol && fabs(Piano1[2] - Piano2[2]) < tol) {
            if(fabs(Piano1[3] - Piano2[3]) < tol) {
                // coincidenti
                cout << "Piani coincidenti tra le fratture " << id1 << " e " << id2 << endl;
                outputFile << "I PIANI SONO COINCIDENTI QUINDI NON POSSONO ESSERCI INTERSEZIONI CHE GENERANO TRACCE" << endl;
                outputFile << id1 << " e " << id2 << endl;
            } else {
                // paralleli
                outputFile << "I PIANI SONO PARALLELI QUINDI NON POSSONO ESSERCI INTERSEZIONI" << endl;
                outputFile << id1 << " e " << id2 << endl;
            }
        } else {
            // incidenti
            cout << "Piani incidenti tra le fratture " << id1 << " e " << id2 << endl;
            vector<Vector3d> intersezioni1;
            vector<Vector3d> intersezioni2;
            pair<Vector3d, Vector3d> traccia;
            map<unsigned int, pair<Vector3d, Vector3d>> mappa_traccia;
            Find_Trace(Piano1, Piano2, id1, id2, FracturesMap[id1], FracturesMap[id2], id_Traccia, mappa_traccia);
            id_Traccia+=1;
        }
    }
    return {Passanti, NonPassanti};
}




void ImportFracturesList(const string& filepath, Fractures& fractures) {
    ifstream file(filepath);

    if (file.fail()) {
        cerr << "Errore nell'apertura del file" << endl;
        return;
    } else {
        string line;
        getline(file, line); // legge la prima riga e la salta
        getline(file, line);
        istringstream convertInt(line); // converte la linea in intero
        convertInt >> fractures.numFractures;
        unsigned int IdFracture;
        unsigned int numeroVertici;
        char delimiter;
        for (unsigned int i = 0; i < fractures.numFractures; i++) {
            getline(file, line); // salta la linea "FractureId; NumVertices"
            getline(file, line);
            istringstream convertLine(line);
            convertLine >> IdFracture >> delimiter >> numeroVertici;
            fractures.VerticesCoordinates = MatrixXd::Zero(3, numeroVertici);
            getline(file, line); // salta la linea "Vertices"
            for (unsigned int riga = 0; riga < 3; riga++) {
                getline(file, line);
                istringstream convertCoordinates(line);
                for (unsigned int colonna = 0; colonna < numeroVertici; colonna++) {
                    convertCoordinates >> fractures.VerticesCoordinates(riga, colonna); // per ogni colonna un vertice diverso : riga 1 = x1,x2,x3,x4
                    if (colonna != numeroVertici - 1)
                        convertCoordinates >> delimiter; // l'ultima colonna non ha il delimitatore
                }
            }
            fractures.FracturesMap[IdFracture] = fractures.VerticesCoordinates; // aggiorna la mappa delle fratture
        }
    }

    file.close();
}
}
