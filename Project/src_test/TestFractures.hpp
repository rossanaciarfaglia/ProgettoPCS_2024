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
    Vector4d expected_res = {0,0,1,1};

    for (unsigned int i=0; i<4; i++){
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
    Vector4d expected_res = {0,-1,0,-2.6};

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
    Vector4d expected_res = {2,-3,2.5,-3};

    for (unsigned int i=0; i<4; i++){
        EXPECT_NEAR(result[i], expected_res[i], 0.00001);
    }
}


TEST(IntersectionsTests, TestIntersezioneSfere){
    unsigned int t = 3;
    unsigned int q = 4;
    MatrixXd Triangolo = (MatrixXd(3,3) << 3,3,2,
                                          2,1,1,
                                          1,1,1).finished();
    MatrixXd Quadrilatero = (MatrixXd(3,4) << 3,2,2,3,
                                              2.6,2.6,2.6,2.6,
                                              1,1,2,2).finished();
    Fracture triangoloTest = {t,Triangolo};
    Fracture quadrilateroTest = {q,Quadrilatero};

    bool result = IntersezioneSfere(triangoloTest,quadrilateroTest);

    EXPECT_TRUE(result);
}


TEST(IntersectionsTests, TestDistanzaEuclidea){
    Vector3d baricentro_1 = {2.66667,1.33333,1};
    Vector3d baricentro_2 = {2.5,2.6,1.5};
    double result = DistanzaEuclidea(baricentro_1, baricentro_2);

    EXPECT_NEAR(result, 1.88223, 0.00001);
}


TEST(IntersectionsTests, TestIntersezionePiani){
    unsigned int t = 3;
    unsigned int q = 4;
    MatrixXd Triangolo = (MatrixXd(3,3) << 3,3,2,
                          2,1,1,
                          1,1,1).finished();
    MatrixXd Quadrilatero = (MatrixXd(3,4) << 3,2,2,3,
                             2.6,2.6,2.6,2.6,
                             1,1,2,2).finished();
    Fracture triangoloTest = {t,Triangolo};
    Fracture quadrilateroTest = {q,Quadrilatero};

    Matrix<double,2,3> result = IntersezionePiani(triangoloTest,quadrilateroTest);
    Matrix<double,2,3> expected_res;
    expected_res << 0,2.6,1,
                    1,0,0;

    for (unsigned int i=0; i<2; i++){
        for (unsigned int j=0; j<3; j++){
            EXPECT_NEAR(result(i,j), expected_res(i,j), 0.00001);
        }
    }
}


// TEST(IntersectionsTests, TestParametriRetta_1){     // Non ho capito
//     Vector3d P0 = {3,2,1};
//     Vector3d P1 = {3,1,1};
//     Vector3d P2 = {2,1,1};

//     Vector3d Q = {0,2.6,1};
//     Vector3d coeff_retta = {1,0,0};

//     Vector2d result_01 = ParametriRetta(P0,P1,Q,coeff_retta);
//     Vector2d expectedres_01 = {};
// }


TEST(IntersectionsTests, TestIntersectionPoint){
    unsigned int n = 3;
    MatrixXd Triangolo = (MatrixXd(3,3) << 3,3,2,
                          2,1,1,
                          1,1,1).finished();
    Matrix<double,2,3> retta;
    retta << 0,2.6,1,
            1,0,0;

    vector<Vector3d> result = Intersection_Point(retta,Triangolo,n);

    for (unsigned int i=0; i<2; i++){
        EXPECT_TRUE(isinf(result[0][i]));
    }
}



}
