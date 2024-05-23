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
Vector3d Fractures::Baricentro(MatrixXd &Poligono) {
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
double Fractures::Raggio(Vector3d &baricentro, MatrixXd &Poligono) {
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


Vector4d Fractures::TrovaPiano(MatrixXd &poligono){
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


bool IntersezioneSfere(Fractures &polygons, MatrixXd &poly_1, MatrixXd &poly_2){
    double tol_quad = 100 * numeric_limits<double>::epsilon() * numeric_limits<double>::epsilon();

    Vector3d baricentro1 = polygons.Baricentro(poly_1);
    double R1 = polygons.Raggio(baricentro1, poly_1);
    Vector3d baricentro2 = polygons.Baricentro(poly_2);
    double R2 = polygons.Raggio(baricentro2, poly_2);
    if (DistanzaEuclidea(baricentro1, baricentro2) - (R1+R2+2*sqrt(R1*R2)) < tol_quad){
        return true;
    }
    return false;
}


double DistanzaEuclidea(Vector3d &centro1, Vector3d &centro2) {
    double distanza = 0;
    for (unsigned i = 0; i < 3; i++) {
        distanza += pow(centro1[i] - centro2[i], 2);
    }
    return distanza;
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
            MatrixXd VerticesCoordinates = MatrixXd::Zero(3, numeroVertici);
            getline(file, line); // salta la linea "Vertices"
            for (unsigned int riga = 0; riga < 3; riga++) {
                getline(file, line);
                istringstream convertCoordinates(line);
                for (unsigned int colonna = 0; colonna < numeroVertici; colonna++) {
                    convertCoordinates >> VerticesCoordinates(riga, colonna); // per ogni colonna un vertice diverso : riga 1 = x1,x2,x3,x4
                    if (colonna != numeroVertici - 1)
                        convertCoordinates >> delimiter; // l'ultima colonna non ha il delimitatore
                }
            }
            fractures.FracturesMap[IdFracture] = VerticesCoordinates; // aggiorna la mappa delle fratture
        }
    }

    file.close();
}


vector<Vector3d> Intersection_Point(Matrix<double, 2, 3> &retta, MatrixXd &vertici){
    Vector3d punto_intersezione ;
    vector<Vector3d> intersezioni;
    intersezioni.reserve(2);
    int num_c = vertici.cols();
    for(unsigned int c = 0; c < num_c; c ++){
        double t_line ;
        double t_segment ;
        t_line =((vertici(0,(c+1)%(num_c-1)) - vertici(0,c)) * (retta(1,0) - vertici(1,c)) -
                  (vertici(1,(c+1)%(num_c-1)) - vertici(1,c)) * (retta(0,0) - vertici(0,c))) /
                 ((vertici(1,(c+1)%(num_c-1)) - vertici(1,c)) * retta(0,1) -
                  (vertici(0,(c+1)%(num_c-1)) - vertici(1,c)) * retta(1,1));
        t_segment = (retta(0,0) + retta(1,0) * t_line - vertici(0,c))/(vertici(0,(c+1)%(num_c-1)) - vertici(0,c));

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


Matrix<double,2,3> IntersezionePiani(Fractures &polygons, vector<Vector3d> intersezioni1, vector<Vector3d> intersezioni2,
                           MatrixXd &poly_1, MatrixXd &poly_2) {

    double tol = 10 * numeric_limits<double>::epsilon();

    Vector4d n1 = polygons.TrovaPiano(poly_1);
    Vector4d n2 = polygons.TrovaPiano(poly_2);

    Vector3d t;
    t[0] = n1[1]*n2[2]-n1[2]*n2[1];
    t[1] = n1[2]*n2[0]-n1[0]*n2[2];
    t[2] = n1[0]*n2[1]-n1[1]*n2[0];

    if(t[0]-0 >= tol && t[1]-0 >= tol && t[2]-0 >= tol) {
        // Non sono paralleli o complanari

        Vector3d n1_reduced;
        n1_reduced << n1[0], n1[1], n1[2];
        Vector3d n2_reduced;
        n2_reduced << n2[0], n2[1], n2[2];
        Matrix3d A;
        A.row(0) = n1_reduced.transpose();
        A.row(1) = n2_reduced.transpose();
        A.row(2) = t.transpose();

        Vector3d P;

        Vector3d b;
        b << n1[3], n2[3], 0; //vettore delle d
            // Risoluzione del sistema con il metodo QR
        P = A.colPivHouseholderQr().solve(b);

        Matrix<double,2,3> retta_intersezione;
        for(unsigned int colonna = 0; colonna < 3 ; colonna ++ ){
            retta_intersezione(0,colonna) = P[colonna];  //la x0 y0 z0
            retta_intersezione(1,colonna) = t[colonna];  //la a, b, c --> direzioni lungo i piani
            //tali che la retta r(t) = (x0, y0, z0) + t (a,b,c)
        }

        intersezioni1 = Intersection_Point(retta_intersezione, poly_1);
        intersezioni2 = Intersection_Point(retta_intersezione, poly_2);

        if (intersezioni1.size() == 2 && intersezioni2.size() == 2) {
            return retta_intersezione;
        }
    }
}


//typedef Matrix<double, 3, 1> Vector3d;
void Find_Trace(vector<Vector3d> &intersezioni1, vector<Vector3d> &intersezioni2,MatrixXd retta) {

    if (!isLessOrEqual(intersezioni1[0], intersezioni1[1], retta)) swap(intersezioni1[1], intersezioni1[0]);
    if (!isLessOrEqual(intersezioni2[0], intersezioni2[1], retta)) swap(intersezioni2[1], intersezioni2[0]);


    pair<Vector3d, Vector3d>  traccia = Traccia(intersezioni1, intersezioni2);

    //mappa_traccia[id_Traccia] = traccia;
}
}
