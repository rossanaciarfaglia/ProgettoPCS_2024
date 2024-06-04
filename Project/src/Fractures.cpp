#include "Fractures.hpp"
#include "MergeSortAlgorithm.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <functional> //per usare greater per ordinare in ordine decrescente le chiavi della multimap
#include <algorithm>

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

    double alfa_start = max(alfa_0, beta_0);
    double alfa_end = min(alfa_1, beta_1);
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

    return {intersection_start, intersection_end};
}

bool Tips (vector<Vector3d>& intersezioni, pair<Vector3d,Vector3d>& verticiTraccia){
    cout<<"intersezioni"<<endl;
    for (unsigned int i = 0; i < 3; i++){
        for (unsigned int j = 0; j < 3; j++){
            cout<<intersezioni[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"vertici traccia"<<endl;
    for (unsigned int i = 0; i < 3; i++){
            cout<<verticiTraccia.first[i]<<" ";
        }
        cout<<endl;

        for (unsigned int i = 0; i < 3; i++){
            cout<<verticiTraccia.second[i]<<" ";
        }
        cout<<endl;

    unsigned int passante = 0;
        double tol = 1e-9;
    // controlliamo il primo poligono
        if (fabs(verticiTraccia.first[0] - intersezioni[0][0]) <= tol &&
            fabs(verticiTraccia.first[1] - intersezioni[0][1]) <= tol &&
            fabs(verticiTraccia.first[2] - intersezioni[0][2]) <= tol) {
        passante += 1;
    }

        if(fabs(verticiTraccia.second[0] - intersezioni[1][0]) <= tol &&
           fabs(verticiTraccia.second[1] - intersezioni[1][1]) <= tol &&
           fabs(verticiTraccia.second[2] - intersezioni[1][2]) <= tol) {
        passante += 1;
    }

    if(passante == 2){
        return true;}
    else {
            return false;}
}

//typedef Matrix<double, 3, 1> Vector3d;
bool Find_Trace(Fracture& polygon, Trace& trace, unsigned int& idT,Fracture& poligono1, Fracture& poligono2) {
    bool tips;
    Matrix<double,2,3> retta_intersezione = IntersezionePiani(polygon, poligono1.Vertici, poligono2.Vertici);

    vector<Vector3d> intersezioni1 = Intersection_Point(retta_intersezione, poligono1.Vertici, poligono1.numVertici);
    vector<Vector3d> intersezioni2 = Intersection_Point(retta_intersezione, poligono2.Vertici, poligono2.numVertici);
<<<<<<< HEAD

    cout << "inter1[0]: " << intersezioni1[0] << endl;
    cout << "inter1[1]: " << intersezioni1[1] << endl;
    cout << "inter2[0]: " << intersezioni2[0] << endl;
    cout << "inter2[1]: " << intersezioni2[1] << endl;
=======
>>>>>>> 2dee1827ee75f087d6261ed90e0bb71b793df925

    if (intersezioni1.size() < 2 || intersezioni2.size() < 2) {
        return false;
    }

    if (isLess(intersezioni1[1], intersezioni1[0], retta_intersezione)) swap(intersezioni1[1], intersezioni1[0]);
    if (isLess(intersezioni2[1], intersezioni2[0], retta_intersezione)) swap(intersezioni2[1], intersezioni2[0]);


    pair<Vector3d, Vector3d> a = Traccia(intersezioni1, intersezioni2, retta_intersezione);
    if(a.first == Vector3d::Zero() && a.second == Vector3d::Zero()){
        return false;
    }
    trace.Vertices = a;
    //if (isLess(trace.Vertices.first, trace.Vertices.second, retta_intersezione)) swap(trace.Vertices.second, trace.Vertices.first);
    trace.lenght = sqrt(DistanzaEuclidea(trace.Vertices.second, trace.Vertices.first));
    tips = Tips(intersezioni1, trace.Vertices);
    if(tips)
    {poligono1.traccePassanti.push_back(idT);}
    else {poligono1.tracceNonPassanti.push_back(idT);}

    tips =Tips(intersezioni2, trace.Vertices);
    if(tips){
        poligono2.traccePassanti.push_back(idT);}
    else {poligono2.tracceNonPassanti.push_back(idT);}
    trace.id = idT;
    return true;
}


void OutputSort (const vector<unsigned int>& IdTrace, const vector<Trace>& elencoTracce, ofstream& FileFracture, bool& tips){
    unordered_map<unsigned int, double> dizionario;
    for(unsigned int i = 0; i < IdTrace.size(); i++){
        dizionario.insert({IdTrace[i], elencoTracce[IdTrace[i]].lenght});
    }
    //trasferimento elementi del dizionario in un vettore di coppie per poter utilizzare sort
    vector<pair<unsigned int, double>> coppie_traccia(dizionario.begin(), dizionario.end());
    //ordinamento del vettore (in ordine decrescente)
    sort(coppie_traccia.begin(), coppie_traccia.end(), compare);
    for (const auto& elem : coppie_traccia) {
        // Fai qualcosa con elem
        FileFracture << elem.first << " " << tips << " " << elem.second << endl;
    }
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

