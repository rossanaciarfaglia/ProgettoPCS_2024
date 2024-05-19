#pragma once

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "UCDUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
TEST(PolygonsTests, TestSpheresIntersection){
    Fractures poligoniTest;

    MatrixXd Triangolo = (MatrixXd(3,3) << 3,3,2,
                                          2,1,1,
                                          1,1,1).finished();
    double a = 2.6;
    MatrixXd Quadrilatero = (MatrixXd(3,4) << 3,2,2,3,
                                              a,a,a,a,
                                              1,1,2,2).finished();
    unsigned int id1 = 0;
    unsigned int id2 = 1;
    poligoniTest.FracturesMap[id1] = Triangolo;
    poligoniTest.FracturesMap[id2] = Quadrilatero;

    if (!IntersezioneSfere(poligoniTest,id1,id2)){
        cout << "non si intersecano" << endl;
    }
    else {
        cout << "si intersecano" << endl;
    }

}
}
