#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace GeometryLibrary;


int main()
{
    string filepath = "./DFN/FR3_data.txt";
    Fracture fracture;
    unordered_map<unsigned int, Fracture> CollectionFractures; // Il costo computazionale è O(1), non O(logn)
    ImportFracturesList(filepath, fracture, CollectionFractures);


//prendiamo ID vettore e ne prendiamo l'id
    unsigned int n_key = CollectionFractures.size();
    for (unsigned int idF1 = 0; idF1<n_key; idF1++){
        for(unsigned int idF2 = (idF1+1); idF2<n_key;idF2++){
            unsigned int idT = 0;

            if(IntersezioneSfere(fracture, CollectionFractures[idF1].Vertici, CollectionFractures[idF2].Vertici)){
                cout << "Le sfere delle fratture " << idF1 << " e " << idF2 << " si intersecano" << endl;
                //Studio del rapporto tra i due piani
                pair<Vector3d, Vector3d> traccia = Find_Trace(fracture, CollectionFractures[idF1].Vertici, CollectionFractures[idF2].Vertici);
                cout << "traccia.first[0]" << endl;
            }
        }
    }

    // Prepare the output file stream
    ofstream outputFile("No_Intersection.txt");
    if (!outputFile.is_open()) {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    // Call the Passante_NonPassante function
    //auto result = fractures.Passante_NonPassante(outputFile);

    // Close the output file stream
    outputFile.close();

    return 0;
}
