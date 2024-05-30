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
Vector3d Fracture::Baricentro(MatrixXd &poligono) {
    Vector3d baricentro;
    for(unsigned int riga = 0; riga < 3; riga++) {
        double sum = 0;
        for(unsigned int colonna = 0; colonna < numVertici; colonna++) {
            sum += poligono(riga, colonna);
        }
        baricentro[riga] = sum / numVertici;
    }
    return baricentro;
}


// Raggio delle sfere !!AL QUADRATO!!
double Fracture::Raggio(Vector3d &baricentro, MatrixXd &poligono) {
    double R = 0;
    double dist;
    for (unsigned int colonna=0; colonna<numVertici; colonna++){
        dist = 0;   //passando alla colonna successiva devo riazzerare la distanza
        for (unsigned int i=0; i<3; i++){
            dist += (baricentro[i]-poligono(i,colonna))*(baricentro[i]-poligono(i,colonna));
        }
        if (dist > R){    //salvo solo se la ristanza è maggiore della distanza
            R = dist;
        }
    }
    return R;
}


Vector4d Fracture::TrovaPiano(MatrixXd &poligono){
    Vector3d u;
    Vector3d v;
    for (unsigned int coordinate=0; coordinate<3; coordinate++){           // (Sono x,y e z)
        u[coordinate] = poligono(coordinate,2) - poligono(coordinate,0);   // u = P2-P0
        v[coordinate] = poligono(coordinate,1) - poligono(coordinate,0);   // v = P1-P0
    }
    // double u_norm = u.norm();
    // double v_norm = v.norm();

    Vector4d piano;   // è il vettore normale n + la costante d
    piano[0] = (u[1]*v[2]-v[1]*u[2]);
    piano[1] = (v[0]*u[2]-u[0]*v[2]);
    piano[2] = (u[0]*v[1]-v[0]*u[1]);

    piano[3] = piano[0]*poligono(0,0) + piano[1]*poligono(1,0) + piano[2]*poligono(2,0);

    return piano;
}


bool IntersezioneSfere(Fracture& polygons, MatrixXd& poly_1, MatrixXd& poly_2){
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



Matrix<double,2,3> IntersezionePiani(Fracture &polygon, MatrixXd &poly_1, MatrixXd &poly_2) {

    Vector4d n1 = polygon.TrovaPiano(poly_1);
    Vector4d n2 = polygon.TrovaPiano(poly_2);
    Vector3d t;
    t[0] = n1[1]*n2[2]-n1[2]*n2[1];
    t[1] = n1[2]*n2[0]-n1[0]*n2[2];
    t[2] = n1[0]*n2[1]-n1[1]*n2[0];

    if(t.norm() != 0) {
    // Non sono paralleli o complanari

        Matrix3d A;
        A.row(0) = n1.head<3>();
        A.row(1) = n2.head<3>();
        A.row(2) = t;

        Vector3d Q;
        Vector3d b;
        b << n1[3], n2[3], 0;  //vettore delle d
        // Risoluzione del sistema con il metodo QR
        FullPivLU<Matrix3d> lu_decomp(A);
        Q = lu_decomp.solve(b);

        Matrix<double,2,3> retta_intersezione;

        retta_intersezione.row(0) = Q;  //la x0 y0 z0
        retta_intersezione.row(1) = t;  //la a, b, c --> direzioni lungo i piani
            //tali che la retta r(t) = (x0, y0, z0) + t (a,b,c)
        return retta_intersezione;
    }
}


vector<Vector3d> Intersection_Point(Matrix<double, 2, 3>& retta, MatrixXd& vertici, const unsigned int& numVert){
    Vector3d punto_intersezione;
    vector<Vector3d> intersezioni;
    intersezioni.reserve(2);
    Vector2d system_solution;
    for(unsigned int c = 0; c < numVert; c++){
        if (c == numVert - 1){ // l'ultimo vertice viene confrontato con il primo
            system_solution = ParametriRette(vertici.col(c), vertici.col(0), retta.row(0), retta.row(1));
            // controllo che alpha sia coerente anche con la coordinata z : (1-alpha)z0+alpha*z1 = qz+tdz
        }
        else {
            system_solution = ParametriRette(vertici.col(c), vertici.col(c+1), retta.row(0), retta.row(1));
        }

        if (system_solution[0] >= 0 && system_solution[0] <= 1) {   // Questo è il segmento
            // Calcola le coordinate del punto di intersezione
            punto_intersezione[0] = retta(0,0) + retta(1,0) * system_solution[1];       // Questa è la retta
            punto_intersezione[1] = retta(0,1) + retta(1,1) * system_solution[1];
            punto_intersezione[2] = retta(0,2) + retta(1,2) * system_solution[1];
            intersezioni.push_back(punto_intersezione);
        }
        //altrimenti non ci sono intersezioni lato, retta di intersezione piani
        //non ci serve memorizzarli
    }
    return intersezioni;
}


pair<Vector3d, Vector3d> Traccia(vector<Vector3d> &intersezioni1, vector<Vector3d> &intersezioni2, Matrix<double,2,3> retta_inters){
    //trova l'inizio e la fine dell'intersezione
    Vector3d intersection_start;
    Vector3d intersection_end ;

    double alfa_0 = 0;
    double alfa_1 = 0;
    double beta_0 = 0;
    double beta_1 = 0;
    for (unsigned int c=0; c<3; c++){
        alfa_0 += (intersezioni1[0][c]-retta_inters(0,c))*retta_inters(1,c);
        alfa_1 += (intersezioni1[1][c]-retta_inters(0,c))*retta_inters(1,c);
        beta_0 += (intersezioni2[0][c]-retta_inters(0,c))*retta_inters(1,c);
        beta_1 += (intersezioni2[1][c]-retta_inters(0,c))*retta_inters(1,c);
    }
    cout << "a0: " << alfa_0 << endl << "a1: " << alfa_1 << endl << "b0: " << beta_0 << endl << "b1: " << beta_1 << endl;

    double alfa_start = max(alfa_0, beta_0);
    double alfa_end = min(alfa_1, beta_1);
    cout << "a_s: " << alfa_start << endl << "a_e: " << alfa_end << endl;
    if (alfa_start > alfa_end){
        // Se non c'è sovrapposizione, possiamo restituire un valore indicativo
        // di nessuna intersezione, come due punti uguali o una coppia di zero.
        return {Vector3d::Zero(), Vector3d::Zero()};
    }

    double norma_r = retta_inters.row(1).norm();
    for (int i = 0; i < 3; i++) {
        intersection_start[i] = retta_inters(0,i) + (alfa_start/(norma_r*norma_r))*retta_inters(1,i);
        intersection_end[i] = retta_inters(0,i) + (alfa_end/(norma_r*norma_r))*retta_inters(1,i);
    }
    cout << "s:" << intersection_start << endl << "e: " << intersection_end << endl;


    return {intersection_start, intersection_end};
}

bool Tips (vector<Vector3d>& intersezioni, pair<Vector3d,Vector3d>& verticiTraccia){
    unsigned int passante = 0;
    // controlliamo il primo poligono
    if (verticiTraccia.first[0] == intersezioni[0][0] && verticiTraccia.first[1] == intersezioni[0][1] && verticiTraccia.first[2] == intersezioni[0][2])
    { passante += 1;}
    if(verticiTraccia.second[0] == intersezioni[1][0] && verticiTraccia.second[1] == intersezioni[1][1] && verticiTraccia.second[2] == intersezioni[1][2])
    { passante += 1;}

    if(passante == 2){
        return true;}
    else {return false;}

}

//typedef Matrix<double, 3, 1> Vector3d;
void Find_Trace(Fracture& polygon, Trace& trace, unsigned int& idT,Fracture& poligono1, Fracture& poligono2) {
    Matrix<double,2,3> retta_intersezione = IntersezionePiani(polygon, poligono1.Vertici, poligono2.Vertici);

    vector<Vector3d> intersezioni1 = Intersection_Point(retta_intersezione, poligono1.Vertici, poligono1.numVertici);
    vector<Vector3d> intersezioni2 = Intersection_Point(retta_intersezione, poligono2.Vertici, poligono2.numVertici);

    cout << "inter1[0]: " << intersezioni1[0] << endl;
    cout << "inter1[1]: " << intersezioni1[1] << endl;
    cout << "inter2[0]: " << intersezioni2[0] << endl;
    cout << "inter2[1]: " << intersezioni2[1] << endl;

    if (intersezioni1.size() < 2 || intersezioni2.size() < 2) {
        return;
    }

    if (isLess(intersezioni1[1], intersezioni1[0], retta_intersezione)) swap(intersezioni1[1], intersezioni1[0]);
    if (isLess(intersezioni2[1], intersezioni2[0], retta_intersezione)) swap(intersezioni2[1], intersezioni2[0]);

    cout << "i1[0]:" << intersezioni1[0] << endl;
    cout << "i1[1]: " << intersezioni1[1] << endl;
    cout << "i2[0]: " << intersezioni2[0] << endl;
    cout << "i2[1]: " << intersezioni2[1] << endl;

    trace.id = idT;
    trace.Vertices = Traccia(intersezioni1, intersezioni2, retta_intersezione);
    if(Tips(intersezioni1, trace.Vertices)){
        poligono1.traccePassanti.push_back(idT);}
    else {poligono1.tracceNonPassanti.push_back(idT);}
    if(Tips(intersezioni2, trace.Vertices)){
        poligono2.traccePassanti.push_back(idT);}
    else {poligono2.tracceNonPassanti.push_back(idT);}
}




void ImportFracturesList(const string& filepath, Fracture& fracture, unordered_map<unsigned int, Fracture>& CollectionFractures){//nome struttura
    ifstream file(filepath);
    unsigned int numFratture;
    unsigned int id;
    if (file.fail()) {
        cerr << "Errore nell'apertura del file" << endl;
        return;
    } else {
        string line;
        getline(file, line); // legge la prima riga e la salta
        getline(file, line);
        istringstream convertInt(line); // converte la linea in intero
        convertInt >> numFratture;
        char delimiter;
        for (unsigned int i = 0; i < numFratture; i++) {
            getline(file, line); // salta la linea "FractureId; NumVertices"
            getline(file, line);
            istringstream convertLine(line);
            convertLine >> id >> delimiter >> fracture.numVertici;
            fracture.Vertici.resize(3,fracture.numVertici);
            getline(file, line); // salta la linea "Vertices"
            for (unsigned int riga = 0; riga < 3; riga++) {
                getline(file, line);
                istringstream convertCoordinates(line);
                for (unsigned int colonna = 0; colonna < fracture.numVertici; colonna++) {
                    convertCoordinates >> fracture.Vertici(riga, colonna); // per ogni colonna un vertice diverso : riga 1 = x1,x2,x3,x4
                    if (colonna != fracture.numVertici - 1)
                        convertCoordinates >> delimiter; // l'ultima colonna non ha il delimitatore
                }
            }
            CollectionFractures[id] = fracture;
        }
    }

    file.close();
}



}

