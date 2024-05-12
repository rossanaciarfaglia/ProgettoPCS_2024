#include "Fractures.hpp"
#include <Eigen/Dense>
#include<fstream>
#include<sstream>
#include<map>
#include<iomanip>  //per setprecision
#include<vector>
#include<cmath>
#include<utility>
using namespace std;
using namespace Eigen;
namespace GeometryLibrary{
//*************************************************************
void ImportFracturesList (const string& filepath,
                         Fractures &fractures)
{
//legge i dati relativi ai poligoni da un file di testo e popola l'oggetto fractures che abbiamo creato nella struct in Fractures.hpp con queste informazioni
    ifstream file(filepath);

if(file.fail()){
        cerr<<"errore nell'apertura del file"<<endl;}
else{
    string line;
    getline(file, line); //legge lap rima riga e la salta
    getline(file,line);
    istringstream convertInt(line); //converte la line in intero
    convertInt>> fractures.numFractures;                                        //ok
    unsigned int IdFracture;
    unsigned int numeroVertici;
    char delimiter = ';';
    for(unsigned int i=0; i< fractures.numFractures; i++){
        getline(file, line); //sto saltando la line = "FattureId; NumVertices"
        getline(file, line);
        istringstream convertLine(line);
        convertLine >> IdFracture >> delimiter >>numeroVertici;
        //cout<<IdFracture<<" "<<numeroVertici<<endl;
        fractures.VerticesCoordinates = MatrixXd::Zero(3, numeroVertici);
        getline(file, line); //sto saldando la linea "Vertices"
        //adesso analizziamo i vertici
        for (unsigned int colonna = 0; colonna < numeroVertici; colonna++){
        getline(file,line);
        istringstream convertCoordinates(line); //converte la line in intero
        for (unsigned int riga = 0; riga < 3; riga++){
            convertCoordinates>> fractures.VerticesCoordinates(riga, colonna); //stiamo riempiendo la matrice VerticesMatrix
            if (colonna != numeroVertici -1)
                convertCoordinates>> delimiter; //perche l'ultima colonna non ha il ;
        }
    }
    fractures.FracturesMap[IdFracture] = fractures.VerticesCoordinates;
}
}
}


vector<double> Baricentro(MatrixXd Poligono){
    vector<double> baricentro(3);
    for(unsigned int colonna = 0; colonna < Poligono.cols(); colonna++){
        double sum = 0;
        for(unsigned int riga = 0; riga < Poligono.rows(); riga++){
            sum += Poligono(riga,colonna);}
        baricentro[colonna] = sum/Poligono.rows();
    }
    return baricentro;
}

double Distanza(vector<double> centro, MatrixXd Poligono){
    double R = 0;
    //per velocità di calcolo usiamo il raggio al quadrato e lo chiamo R
    //abbiamo fissato la riga i=[0]
    for (unsigned int colonna = 0; colonna<Poligono.cols(); colonna++){
        R += pow(Poligono(0,colonna) - centro[colonna],2);
    }
    return R;
}


double DistanzaEuclidea(vector<double> centro1, vector<double> centro2){
    double distanza = 0;
    for (unsigned i = 0; i <3 ; i++){
        distanza += pow(centro1[i] - centro2[i],2);
    }
    return distanza;
}

Vector4d TrovaPiano(MatrixXd Poligono){
    //troviamo il primo vettore tra i primi due vertici
    vector<double> P1 = {};
    vector<double> P2 = {};
    for (unsigned int colonna = 0; colonna < 3; colonna ++){
        P1.push_back(Poligono(0,colonna) - Poligono(1,colonna));
    //troviamo il secondo vettore tra il primo e il terzo vertice
        P2.push_back(Poligono(0,colonna) - Poligono(2,colonna));
    }
    //ora troviamo il prodotto vettoriale tra i due, arrivando cosi alle componeneti della normale
    //che corrispondono ad a,b,c del piano ax+by+cz = d
    double a = P1[1]*P2[2] - P1[2]*P2[1];
    double b = -(P1[0]*P2[2] - P1[2]*P2[0]);
    double c = P1[0]*P2[1] - P1[1]*P2[0];
    //troviamo la d imponendo il passaggio per il vettore 4 che non abbiamo considerato nel piano
    double d = a*Poligono(Poligono.rows()-1, 0) + b*Poligono(Poligono.rows()-1, 1) + c*Poligono(Poligono.rows()-1, 2);
    return {a,b,c,d};
    }




//  QUESTA FUNZIONE SI OCCUPA DI DEFINIRE L'INTERSEZIONE TRA SFERE,
//  ELIMINANDO LE FRATTURE, LE CUE SFERE CIRCOSCRITTE, NON SI INTERSECANO
vector<vector<unsigned int>> Fractures::SpheresIntersection(){
    vector<vector<unsigned int>> SphereIntersection = {};
    //2--> apriamo un file di output in cui inseriamo le fratture che non hanno intersezioni reciproche con le sfere
    ofstream outputFile("No_Intersection.txt");
    //trovare il baricentro
    //1--> fissiamo una frattura e analizziamo l'eventuale intersezione con le altre in base alle sfere circoscritte
    for(unsigned int key = 0; key < numFractures; key++){
        vector<double> baricentro1 = {};
        baricentro1 = Baricentro(FracturesMap[key]);
        double R1 = Distanza(baricentro1, FracturesMap[key]);
        for(unsigned int chiavi = key +1; chiavi < numFractures; chiavi ++){ //prendo solo le chiavi con id maggiore
            vector<double> baricentro2 = {};
            baricentro2 = Baricentro(FracturesMap[chiavi]);
            double R2 = Distanza(baricentro2, FracturesMap[chiavi]);
            if(R1+R2 <= DistanzaEuclidea(baricentro1, baricentro2)){
                outputFile<<"NON CI SONO INTERSEZIONI TRA LE SFERE DELLE SEGUENTI FRATTURE"<<endl;
                outputFile<< key << " e " << chiavi<<endl;
            }
            else{
                SphereIntersection.push_back({key, chiavi});
            }
        }
      }
    return SphereIntersection;
    }

pair<map<unsigned int, vector<unsigned int>>, map<unsigned int, vector<unsigned int>>> Fractures::Passante_NonPassante(){
    vector<vector<unsigned int>> Lista_Intersezioni_Sfere = SpheresIntersection();
    for(unsigned int i = 0; i < Lista_Intersezioni_Sfere.size(); i++){
        //effettuo il calcolo del piano
        unsigned int id1 = Lista_Intersezioni_Sfere[i][0];
        unsigned int id2 = Lista_Intersezioni_Sfere[i][1];
        Vector4d Piano1 = TrovaPiano(FracturesMap[id1]);
        Vector4d Piano2 = TrovaPiano(FracturesMap[id2]);

        //controllo se i piani sono paralleli e coincidenti

    }
}

} //geometry library


