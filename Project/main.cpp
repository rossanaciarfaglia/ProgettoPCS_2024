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
        // cout << "hello" << endl;
        unsigned int idF2 = idF+1;
        if(IntersezioneSfere(fractures,idF,idF2)){
            cout << "Le sfere delle fratture " << idF << " e " << idF2 << " si intersecano" << endl;
            // TrovaPiano per entrambe le fratture
            // Studio del rapporto tra i due piani
        }
    }

    // // Prepare the output file stream
    // ofstream outputFile("output.txt");
    // if (!outputFile.is_open()) {
    //     cerr << "Error opening output file." << endl;
    //     return -1;
    // }

    // // Call the Passante_NonPassante function
    // auto result = fractures.Passante_NonPassante(outputFile);

    // // Close the output file stream
    // outputFile.close();

    return 0;
}
