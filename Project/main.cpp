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
    Fractures fractures;
    ImportFracturesList(filepath, fractures);

    for (unsigned int idF=0; idF<fractures.numFractures-1; idF++){
        unsigned int idF2 = (idF+1)%fractures.numFractures;
        if(IntersezioneSfere(fractures,idF,idF2)){
            cout << "Le sfere delle fratture " << idF << " e " << idF2 << " si intersecano" << endl;
            // TrovaPiano per entrambe le fratture:
            Vector4d n1 = fractures.TrovaPiano(fractures.FracturesMap[idF]);
            Vector4d n2 = fractures.TrovaPiano(fractures.FracturesMap[idF2]);

            // Studio del rapporto tra i due piani
        }
    }

    // Prepare the output file stream
    ofstream outputFile("No_Intersection.txt");
    if (!outputFile.is_open()) {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    // Call the Passante_NonPassante function
    auto result = fractures.Passante_NonPassante(outputFile);

    // Close the output file stream
    outputFile.close();

    return 0;
}
