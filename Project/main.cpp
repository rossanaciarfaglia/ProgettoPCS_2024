#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
using namespace std;
using namespace GeometryLibrary;
int main()
{
    string filepath = "./DFN/FR3_data.txt";
    Fractures fractures;
    ImportFracturesList(filepath, fractures);

    return 0;
}
