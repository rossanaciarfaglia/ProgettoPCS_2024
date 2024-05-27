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
Vector3d Fracture::Baricentro(MatrixXd &Poligono) {
    Vector3d baricentro;
    for(unsigned int riga = 0; riga < 3; riga++) {
        double sum = 0;
        for(unsigned int colonna = 0; colonna < numVertici; colonna++) {
            sum += Poligono(riga, colonna);
        }
        baricentro[riga] = sum / numVertici;
    }
    cout << "baricentro ok" << endl;
    return baricentro;
}


// Raggio delle sfere !!AL QUADRATO!!
double Fracture::Raggio(Vector3d &baricentro, MatrixXd &Poligono) {
    double R = 0;
    double dist;
    for (unsigned int colonna=0; colonna<numVertici; colonna++){
        dist = 0;   //passando alla colonna successiva devo riazzerare la distanza
        for (unsigned int i=0; i<3; i++){
            dist += (baricentro[i]-Poligono(i,colonna))*(baricentro[i]-Poligono(i,colonna));
        }
        if (dist > R){    //salvo solo se la ristanza è maggiore della distanza
            R = dist;
        }
    }
    cout << "raggio ok" << endl;
    return R;
}


Vector4d Fracture::TrovaPiano(MatrixXd &poligono){
    Vector3d u;
    Vector3d v;
    for (unsigned int coordinate=0; coordinate<3; coordinate++){           // (Sono x,y e z)
        u[coordinate] = poligono(coordinate,2) - poligono(coordinate,0);   // u = P2-P0
        v[coordinate] = poligono(coordinate,1) - poligono(coordinate,0);   // v = P1-P0
    }
    double u_norm = u.norm();
    double v_norm = v.norm();

    Vector4d piano;   // è il vettore normale n + la costante d
    piano[0] = (u[1]*v[2]-v[1]*u[2])/(u_norm*v_norm);
    piano[1] = (v[0]*u[2]-u[0]*v[2])/(u_norm*v_norm);
    piano[2] = (u[0]*v[1]-v[0]*u[1])/(u_norm*v_norm);

    piano[3] = piano[0]*poligono(0,0) + piano[1]*poligono(1,0) + piano[2]*poligono(2,0);

    cout << "trova piano ok" << endl;
    return piano;
}


bool IntersezioneSfere(Fracture& polygons, MatrixXd& poly_1, MatrixXd& poly_2){
    double tol_quad = 100 * numeric_limits<double>::epsilon() * numeric_limits<double>::epsilon();

    Vector3d baricentro1 = polygons.Baricentro(poly_1);
    double R1 = polygons.Raggio(baricentro1, poly_1);
    Vector3d baricentro2 = polygons.Baricentro(poly_2);
    double R2 = polygons.Raggio(baricentro2, poly_2);
    cout << "inters sfere ok" << endl;
    if (DistanzaEuclidea(baricentro1, baricentro2) - (R1+R2+2*sqrt(R1*R2)) < tol_quad){
        return true;
    }
    return false;
}



Matrix<double,2,3> IntersezionePiani(Fracture &polygon, MatrixXd &poly_1, MatrixXd &poly_2) {

    double tol = 10 * numeric_limits<double>::epsilon();

    Vector4d n1 = polygon.TrovaPiano(poly_1);
    Vector4d n2 = polygon.TrovaPiano(poly_2);
    Vector3d t;
    t[0] = n1[1]*n2[2]-n1[2]*n2[1];
    t[1] = n1[2]*n2[0]-n1[0]*n2[2];
    t[2] = n1[0]*n2[1]-n1[1]*n2[0];

    if(t.norm() != 0) {
    // Non sono paralleli o complanari
        cout << "entra nell'if" << endl;

        // Vector3d n1_reduced;
        // n1_reduced << n1[0], n1[1], n1[2];
        // Vector3d n2_reduced;
        // n2_reduced << n2[0], n2[1], n2[2];
        Matrix3d A;
        A.row(0) = n1.head<3>();
        A.row(1) = n2.head<3>();
        A.row(2) = t;
        cout << "salva A" << endl;

        Vector3d P;

        Vector3d b;
        b << n1[3], n2[3], 0;  //vettore delle d
        // Risoluzione del sistema con il metodo QR
        FullPivLU<Matrix3d> lu_decomp(A);
        P = lu_decomp.solve(b);
        cout << "risolve la decomposizione" << endl;

        Matrix<double,2,3> retta_intersezione;

        for(unsigned int colonna = 0; colonna < 3 ; colonna ++ ){
            retta_intersezione(0,colonna) = P[colonna];  //la x0 y0 z0
            retta_intersezione(1,colonna) = t[colonna];  //la a, b, c --> direzioni lungo i piani
            //tali che la retta r(t) = (x0, y0, z0) + t (a,b,c)
        }
        cout << "inters piani ok" << endl;
        return retta_intersezione;
    }
}


vector<Vector3d> Intersection_Point(Matrix<double, 2, 3> &retta, MatrixXd &vertici){    //metrere qui la condizione di =2
    Vector3d punto_intersezione ;
    vector<Vector3d> intersezioni;
    intersezioni.reserve(2);
    int num_c = vertici.cols(); // usare numVert della structure
    Vector2d system_solution;
    for(unsigned int c = 0; c < num_c; c ++){
        system_solution = ParametriRette(vertici.col(c),vertici.col((c+1)%(num_c-1)),retta.row(0),retta.row(1));
        if (system_solution[0] >= 0 && system_solution[0] <= 1) {
            // Calcola le coordinate del punto di intersezione
            punto_intersezione[0] = retta(0,0) + retta(1,0) * system_solution[1];
            punto_intersezione[1] = retta(0,1) + retta(1,1) * system_solution[1];
            punto_intersezione[2] = retta(0,2) + retta(1,2) * system_solution[1];
            intersezioni.push_back(punto_intersezione);
        }
        //altrimenti non ci sono intersezioni lato, retta di intersezione piani
        //non ci serve memorizzarli
    }
    if (intersezioni.size() == 2){
        return intersezioni;
    }
    cout << "intersection point ok" << endl;
}


//typedef Matrix<double, 3, 1> Vector3d;
pair<Vector3d, Vector3d> Find_Trace(Fracture& polygon, MatrixXd& vert_1, MatrixXd& vert_2) {
    Matrix<double,2,3> retta_intersezione = IntersezionePiani(polygon, vert_1, vert_2);

    vector<Vector3d> intersezioni1 = Intersection_Point(retta_intersezione, vert_1);
    vector<Vector3d> intersezioni2 = Intersection_Point(retta_intersezione, vert_2);

    if (!isLessOrEqual(intersezioni1[0], intersezioni1[1], retta_intersezione)) swap(intersezioni1[1], intersezioni1[0]);
    if (!isLessOrEqual(intersezioni2[0], intersezioni2[1], retta_intersezione)) swap(intersezioni2[1], intersezioni2[0]);

    pair<Vector3d, Vector3d> traccia;
    traccia = Traccia(intersezioni1, intersezioni2);
    return traccia;
    cout << "find trace ok" << endl;
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

