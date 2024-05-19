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
vector<double> Fractures::Baricentro(MatrixXd Poligono) {
    vector<double> baricentro(3);
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
double Fractures::Raggio(vector<double> centro, MatrixXd Poligono) {
    double R = 0;
    double dist;
    for (unsigned int numV=0; numV<Poligono.cols(); numV++){
        dist = 0;
        for (unsigned int i=0; i<3; i++){
            dist += (centro[i]-Poligono(i,numV))*(centro[i]-Poligono(i,numV));
        }
        if (dist > R){
            R = dist;
        }
    }
    return R;
}

double DistanzaEuclidea(vector<double> centro1, vector<double> centro2) {
    double distanza = 0;
    for (unsigned i = 0; i < 3; i++) {
        distanza += pow(centro1[i] - centro2[i], 2);        // NB è il quadrato della dist
    }
    return distanza;
}

Vector4d Fractures::TrovaPiano(MatrixXd Poligono) {
    vector<double> P1 = {};
    vector<double> P2 = {};
    for (unsigned int colonna = 0; colonna < 3; colonna++) {
        P1.push_back(Poligono(0, colonna) - Poligono(1, colonna));
        P2.push_back(Poligono(0, colonna) - Poligono(2, colonna));
    }
    double a = P1[1] * P2[2] - P1[2] * P2[1];
    double b = -(P1[0] * P2[2] - P1[2] * P2[0]);
    double c = P1[0] * P2[1] - P1[1] * P2[0];
    double d = a * Poligono(Poligono.rows() - 1, 0) + b * Poligono(Poligono.rows() - 1, 1) + c * Poligono(Poligono.rows() - 1, 2);
    return {a, b, c, d};
}

// vector<vector<unsigned int>> SpheresIntersection() {
//     vector<vector<unsigned int>> SphereIntersection = {};
//     ofstream outputFile("No_Intersection.txt");

//     for(unsigned int key = 0; key < numFractures; key++) {
//         vector<double> baricentro1 = Baricentro(FracturesMap[key]);
//         double R1 = Distanza(baricentro1, FracturesMap[key]);
//         for(unsigned int chiavi = key + 1; chiavi < numFractures; chiavi++) {
//             vector<double> baricentro2 = Baricentro(FracturesMap[chiavi]);
//             double R2 = Distanza(baricentro2, FracturesMap[chiavi]);
//             if(R1 + R2 <= DistanzaEuclidea(baricentro1, baricentro2)) {
//                 outputFile << "NON CI SONO INTERSEZIONI TRA LE SFERE DELLE SEGUENTI FRATTURE" << endl;
//                 outputFile << key << " e " << chiavi << endl;
//             } else {
//                 SphereIntersection.push_back({key, chiavi});
//             }
//         }
//     }
//     return SphereIntersection;
// }
bool IntersezioneSfere(Fractures &polygons, unsigned int id1, unsigned int id2){
    double tol = 10 * numeric_limits<double>::epsilon();

    MatrixXd Poly1 = polygons.FracturesMap[id1];
    vector<double> Baricentro1 = polygons.Baricentro(Poly1);
    double R1 = polygons.Raggio(Baricentro1, Poly1);
    MatrixXd Poly2 = polygons.FracturesMap[id2];
    vector<double> Baricentro2 = polygons.Baricentro(Poly2);
    double R2 = polygons.Raggio(Baricentro2, Poly2);
    cout << "b1: " << Baricentro1[0] << " " << Baricentro1[1] << " "<< Baricentro1[2] << " "<< endl;
    cout << "b2: " << Baricentro2[0] << " " << Baricentro2[1] << " "<< Baricentro2[2] << " "<< endl;
    cout << "r1: " << sqrt(R1) << endl;
    cout << "r2: " << sqrt(R2) << endl;
    cout << "distanza: " << sqrt(DistanzaEuclidea(Baricentro1, Baricentro2)) << endl;
    cout << sqrt(DistanzaEuclidea(Baricentro1, Baricentro2)) - (sqrt(R1)+sqrt(R2)) << endl;
    cout << "tol: " << tol << endl;
    if (sqrt(DistanzaEuclidea(Baricentro1, Baricentro2)) - (sqrt(R1)+sqrt(R2)) < tol){
        return true;
    }
    return false;
}


// Vector3d ProiezioneVertice(Vector3d &Retta, unsigned int &id){
//     vector<vector<double>> Proiez_vectors = {};
//     Vector3d Proiezione;
//     for (unsigned int vert_i=0; vert_i< Fractures::FracturesMap[id].rows(); vert_i++){
//         double dist = abs(Retta[0]*FracturesMap[id][vert_i][0] + Retta[1]*FracturesMap[id][vert_i][1] + Retta[2]*FracturesMap[id][vert_i][2])/
//                       sqrt(Retta[0]*Retta[0] + Retta[1]*Retta[1] + Retta[2]*Retta[2]);
//         Proiezione[0] = FracturesMap[id][vert_i][0] - dist*(1/Retta[0])/sqrt(1/(Retta[0]*Retta[0]) + 1/(Retta[1]*Retta[1]) + 1/(Retta[2]*Retta[2]));
//         Proiezione[1] = FracturesMap[id][vert_i][1] - dist*(1/Retta[1]);
//         Proiezione[2] = FracturesMap[id][vert_i][2] - dist*(1/Retta[2]);

//         Proiez_vectors.push_back(Proiezione);
//     }
// }

// pair<map<unsigned int, vector<unsigned int>>, map<unsigned int, vector<unsigned int>>> Fractures::Passante_NonPassante(ofstream& outputFile) {
//     cout << "ciao" << endl;
//     double tol = 10 * numeric_limits<double>::epsilon();
//     vector<vector<unsigned int>> Lista_Intersezioni_Sfere = SpheresIntersection();
//     map<unsigned int, vector<unsigned int>> Passanti;
//     map<unsigned int, vector<unsigned int>> NonPassanti;

//     for(unsigned int i = 0; i < Lista_Intersezioni_Sfere.size(); i++) {
//         unsigned int id1 = Lista_Intersezioni_Sfere[i][0];
//         unsigned int id2 = Lista_Intersezioni_Sfere[i][1];
//         Vector4d Piano1 = TrovaPiano(FracturesMap[id1]);
//         Vector4d Piano2 = TrovaPiano(FracturesMap[id2]);

//         if (fabs(Piano1[0] - Piano2[0]) < tol && fabs(Piano1[1] - Piano2[1]) < tol && fabs(Piano1[2] - Piano2[2]) < tol) {
//             if(fabs(Piano1[3] - Piano2[3]) < tol) {
//                 // coincidenti
//                 cout << "Piani coincidenti tra le fratture " << id1 << " e " << id2 << endl;
//             }
//             else {
//                 // paralleli
//                 outputFile << "I PIANI SONO PARALLELI QUINDI NON POSSONO ESSERCI INTERSEZIONI" << endl;
//                 outputFile << id1 << " e " << id2 << endl;
//             }
//         }
//         else {
//             // incidenti
//             cout << "Piani incidenti tra le fratture " << id1 << " e " << id2 << endl;
//             Vector3d Retta;
//             Retta [0] = (Piano1[1]-Piano2[2])*(Piano2[1]-Piano1[2]);
//             Retta [1] = (Piano1[0]-Piano2[2])*(Piano2[0]-Piano1[2]);
//             Retta [2] = (Piano1[0]-Piano2[1])*(Piano2[0]-Piano1[1]);


//         }
//     }
//     return {Passanti, NonPassanti};
// }

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
                    convertCoordinates >> fractures.VerticesCoordinates(riga, colonna);
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
