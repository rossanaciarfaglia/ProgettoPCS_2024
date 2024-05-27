#pragma once

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "UCDUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace GeometryLibrary {
TEST(PolygonsTest, TestBaricentro_1){
    unsigned int n = 3;
    MatrixXd Triangolo = (MatrixXd(3,n) << 3,3,2,
                          2,1,1,
                          1,1,1).finished();

    Fracture poligoniTest = {n,Triangolo};
    Vector3d result = poligoniTest.Baricentro(Triangolo);
    Vector3d expected_res = {2.66667,1.33333,1};
    for (unsigned int i=0; i<3; i++){
        EXPECT_NEAR(result[i], expected_res[i], 0.00001);
    }
}
TEST(PolygonsTest, TestBaricentro_2){
    unsigned int n = 4;
    MatrixXd Quadrilatero = (MatrixXd(3,n) << 3,2,2,3,
                             2.6,2.6,2.6,2.6,
                             1,1,2,2).finished();
    Fracture poligoniTest = {n,Quadrilatero};
    Vector3d result = poligoniTest.Baricentro(Quadrilatero);
    Vector3d expected_res = {2.5,2.6,1.5};
    for (unsigned int i=0; i<3; i++){
        EXPECT_NEAR(result[i], expected_res[i], 0.00001);
    }
}


TEST(PolygonsTest, TestRaggio_1){
    unsigned int n = 3;
    MatrixXd Triangolo = (MatrixXd(3,n) << 3,3,2,
                          2,1,1,
                          1,1,1).finished();
    Fracture poligoniTest = {n,Triangolo};
    Vector3d baricentro = poligoniTest.Baricentro(Triangolo);

    double result = poligoniTest.Raggio(baricentro,Triangolo);

    EXPECT_NEAR(result, 0.55556, 0.00001);
}


TEST(PolygonsTest, TestRaggio_2){
    unsigned int n = 4;
    MatrixXd Quadrilatero = (MatrixXd(3,n) << 3,2,2,3,
                             2.6,2.6,2.6,2.6,
                             1,1,2,2).finished();
    Fracture poligoniTest = {n,Quadrilatero};
    Vector3d baricentro = poligoniTest.Baricentro(Quadrilatero);

    double result = poligoniTest.Raggio(baricentro,Quadrilatero);

    EXPECT_NEAR(result, 0.5, 0.00001);
}


TEST(PolygonsTest, TestTrovaPiano_1){
    unsigned int n = 3;
    MatrixXd Triangolo = (MatrixXd(3,n) << 3,3,2,
                          2,1,1,
                          1,1,1).finished();
    Fracture poligoniTest = {n,Triangolo};
    Vector4d result = poligoniTest.TrovaPiano(Triangolo);
    Vector4d expected_res = {0,0,1/1.41421,1/1.41421};

    for (unsigned int i=0; i<4; i++){
        cout << result[i] << " ";
        EXPECT_NEAR(result[i], expected_res[i], 0.00001);
    }
    cout << endl;
}


TEST(PolygonsTest, TestTrovaPiano_2){
    unsigned int n = 4;
    MatrixXd Quadrilatero = (MatrixXd(3,n) << 3,2,2,3,
                             2.6,2.6,2.6,2.6,
                             1,1,2,2).finished();
    Fracture poligoniTest = {n,Quadrilatero};
    Vector4d result = poligoniTest.TrovaPiano(Quadrilatero);
    Vector4d expected_res = {0,-1/1.41421,0,-1.3*1.41421};

    for (unsigned int i=0; i<4; i++){
        EXPECT_NEAR(result[i], expected_res[i], 0.00001);
    }
}

TEST(PolygonsTest, TestTrovaPiano_3){
    unsigned int n = 3;
    MatrixXd Triangolo = (MatrixXd(3,3) << 0,0.5,1.5,
                                       1,3,2,
                                       0,2,0).finished();
    Fracture poligoniTest = {n,Triangolo};
    Vector4d result = poligoniTest.TrovaPiano(Triangolo);
    Vector4d expected_res = {0.386244,-0.579365,0.482805,-0.579365};

    for (unsigned int i=0; i<4; i++){
        EXPECT_NEAR(result[i], expected_res[i], 0.00001);
    }
}


// TEST(IntersectionsTests, TestIntersezioneSfere){
//     map<int, Fracture> mappa;
//     MatrixXd Triangolo = (MatrixXd(3,3) << 3,3,2,
//                                           2,1,1,
//                                           1,1,1).finished();
//     MatrixXd Quadrilatero = (MatrixXd(3,4) << 3,2,2,3,
//                                               2.6,2.6,2.6,2.6,
//                                               1,1,2,2).finished();
//     Fracture poligoniTest1 = {};
//     mappa[0];

//     bool result = IntersezioneSfere(poligoniTest,Triangolo,Quadrilatero);

//     EXPECT_TRUE(result);
// }


// TEST(IntersectionsTests, TestDistanzaEuclidea){
//     Vector3d baricentro_1 = {2.66667,1.33333,1};
//     Vector3d baricentro_2 = {2.5,2.6,1.5};
//     double result = DistanzaEuclidea(baricentro_1, baricentro_2);

//     EXPECT_NEAR(result, 1.88223, 0.00001);
// }

}
