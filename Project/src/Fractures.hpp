#ifndef __FRACTURES_H
#define __FRACTURES_H

#include "Eigen/Eigen"
#include <iostream>
#include <utility>
#include <map>
#include <vector>
using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
struct Fractures{
    unsigned int numFractures;
    MatrixXd VerticesCoordinates;
    map<unsigned int, MatrixXd> FracturesMap;

    Fractures() = default; // costruttore di default

    Fractures(const MatrixXd& VerticesCoordinates,
              unsigned int& numFractures,
              map<unsigned int, MatrixXd>& FracturesMap
              ):
        VerticesCoordinates(VerticesCoordinates),
        numFractures(numFractures),
        FracturesMap(FracturesMap)
    {} // prende in input le coordinate dei vertici e la lista dei vertici e inizializza i membri corrispondenti VerticesCoordinates e ListVertices con i valori passati come argomenti

    Vector3d Baricentro(MatrixXd Poligono);
    double Raggio(Vector3d centro, MatrixXd Poligono);
    Vector4d TrovaPiano(MatrixXd poligono);
};

void ImportFracturesList(const string& filepath, Fractures& fractures);
double DistanzaEuclidea(Vector3d centro1, Vector3d centro2);
bool IntersezioneSfere(Fractures &polygons, unsigned int id1, unsigned int id2);


void ImportFracturesList(const string& filepath,
                         Fractures& fractures);

inline bool isLessOrEqual(Vector3d p1, Vector3d p2, Vector3d t) {
    return (p1[0] - p2[0]) * t[0] + (p1[1] - p2[1]) * t[1] + (p1[2] - p2[2]) * t[2] <= 0;
}


inline pair<Vector3d, Vector3d>  Traccia(vector<Vector3d> punti1, vector<Vector3d> punti2, Vector3d t){

    //trova l'inizio e la fine dell'intersezione
    Vector3d intersection_start;
    Vector3d intersection_end ;
    bool overlap = true;

    for (int i = 0; i < 3; ++i) {
        intersection_start[i] = max(punti1[0][i], punti2[0][i]);
        intersection_end[i] = min(punti1[1][i], punti2[1][i]);

        if (intersection_start[i] > intersection_end[i]) {
            overlap = false;
        }
        if (!overlap) {
            // Se non c'è sovrapposizione, possiamo restituire un valore indicativo
            // di nessuna intersezione, come due punti uguali o una coppia di zero.
            return {Vector3d::Zero(), Vector3d::Zero()};
        }

        return {intersection_start, intersection_end};
    }
}


inline vector<Vector3d> Intersection_Point(Matrix<double, 2, 3> retta, MatrixXd vertici){
    Vector3d punto_intersezione ;
    vector<Vector3d> intersezioni;
    for(unsigned int c = 0; c < vertici.cols(); c ++){
        double t_line ;
        double t_segment ;
        if(c == vertici.cols() - 1){
            t_line = ((vertici(0,0) - vertici(0,c)) * (retta(1,0) - vertici(1,c)) -
                      (vertici(1,0) - vertici(1,c)) * (retta(0,0) - vertici(0,c))) /
                     ((vertici(1,0) - vertici(1,c)) * retta(0,1) -
                      (vertici(0,0) - vertici(1,c)) * retta(1,1));
            t_segment = (retta(0,0) + retta(1,0) * t_line - vertici(0,c))/(vertici(0,0) - vertici(0,c));
        }
        else{
            t_line =((vertici(0,c+1) - vertici(0,c)) * (retta(1,0) - vertici(1,c)) -
                      (vertici(1,c+1) - vertici(1,c)) * (retta(0,0) - vertici(0,c))) /
                     ((vertici(1,c+1) - vertici(1,c)) * retta(0,1) -
                      (vertici(0,c+1) - vertici(1,c)) * retta(1,1));
            t_segment = (retta(0,0) + retta(1,0) * t_line - vertici(0,c))/(vertici(0,c+1) - vertici(0,c));

        }

        cout<<"t_line: "<<t_line;
        cout<<"t_segment: "<<t_segment;
        //La prima parte della condizione if verifica se il parametro t_line è compreso tra 0 e 1, assicurandosi che il punto di intersezione sia all'interno della porzione della retta considerata. La seconda parte della condizione if verifica se il parametro t_segment è compreso tra 0 e 1, assicurandosi che il punto di intersezione si trovi all'interno del segmento del poligono considerato. Solo se entrambe queste condizioni sono soddisfatte, il punto di intersezione è considerato valido.
        //verifica se il punto di intersezione si trova nel segmento
        if (t_line >= 0 && t_line <= 1 && t_segment >= 0 && t_segment <= 1) {
            // Calcola le coordinate del punto di intersezione
            punto_intersezione[0] = retta(0,0) + retta(1,0) * t_line;
            punto_intersezione[1] = retta(0,1) + retta(1,1) * t_line;
            punto_intersezione[2] = retta(0,2) + retta(1,2) * t_line;
            intersezioni.push_back(punto_intersezione);

        }
        //altrimenti non ci sono intersezioni lato, retta di intersezione piani
        //non ci serve memorizzarli
    }
    return intersezioni;
}


//typedef Matrix<double, 3, 1> Vector3d;
inline void Find_Trace(Vector4d& Piano1, Vector4d& Piano2, unsigned int& id1, unsigned int& id2, MatrixXd& Poligono1, MatrixXd Poligono2, unsigned int id_Traccia,
                       map<unsigned int, pair<Vector3d, Vector3d>>&  mappa_traccia) {
    // Poligono1 e Poligono2   riga1 = x1 x2 x3 x4, riga 2 = y1 y2 y3 y4
    Matrix<double, 2, 3> retta_intersezione = {};
    Vector3d n1 = {Piano1[0], Piano1[1], Piano1[2]};
    Vector3d n2 = {Piano2[0], Piano2[1], Piano2[2]};
    Vector3d t;
    // facciamo il prodotto vettoriale tra n1 e n2 e troviamo t (vettore direzione della retta)
    t[0] = n1[1] * n2[2] - n1[2] * n2[1];
    t[1] = -(n1[0] * n2[2] - n1[2] * n2[0]);
    t[2] = n1[0] * n2[1] - n1[1] * n2[0];

    Matrix3d A;
    A.row(0) = n1.transpose();
    A.row(1) = n2.transpose();
    A.row(2) = t.transpose();

    Vector3d P= {};

    Vector3d b ={};
    b<<Piano1[3], Piano2[3],0; //vettore delle d
        // Risoluzione del sistema con il metodo QR
    P = A.colPivHouseholderQr().solve(b);

    for(unsigned int colonna = 0; colonna < 3 ; colonna ++ ){
        retta_intersezione(0,colonna) = P[colonna];  //la x0 y0 z0
        retta_intersezione(1,colonna) = t[colonna];  //la a, b, c --> direzioni lungo i piani
        //tali che la retta r(t) = (x0, y0, z0) + t (a,b,c)
    }
    //return retta_intersezione;
    vector<Vector3d> intersezioni1;
    vector<Vector3d> intersezioni2;
    intersezioni2 = Intersection_Point(retta_intersezione, Poligono2);
    intersezioni1 = Intersection_Point(retta_intersezione, Poligono1);



    for(unsigned int i  = 0; i < intersezioni1.size(); i++){
        for(unsigned int j=0; j < 3; j++){
            cout<<intersezioni1[i][j]<<"  ";
        }
        cout<<endl;
    }
    //analizzo i punti di intersezione e trovo la traccia
    //sicuramente è passante o non passante, dobbiamo capire quale dei due è

    //cerchiamo l'intersezione tra i segmenti individuati da intersezioni1 e inetrsezioni2
    //verifichiamo se i segmenti sono allineati
    //se il primo punto non è piu piccolo del punto 2
    if (intersezioni1.size() < 2 || intersezioni2.size() < 2) {
        cerr << "Errore: non abbastanza punti di intersezione trovati." << endl;
        return;
    }

    if (!isLessOrEqual(intersezioni1[0], intersezioni1[1], t)) swap(intersezioni1[1], intersezioni1[0]);
    if (!isLessOrEqual(intersezioni2[0], intersezioni2[1], t)) swap(intersezioni2[1], intersezioni2[0]);


    pair<Vector3d, Vector3d>  traccia = Traccia(intersezioni1, intersezioni2, t);

    mappa_traccia[id_Traccia] = traccia;
}
}


#endif
