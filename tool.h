#pragma once
#include <iostream>
#include <Eigen\Dense>
#include <fstream>
#include <string>
#include <queue>
#include <vector>
#include <stack>

const int INF = INF;

using namespace Eigen;
using namespace std;

typedef struct Edge {
    int start; int end;
    int f1 = -1; int f2 = -1;
    
}Edge;

typedef struct STvertex {
    int idx = 0;
    float Cal = 0;
}STvertex;

typedef struct f_patch {
    int idx;
    bool isDiscover=false;
    float weight = INF;
    friend bool operator< (f_patch a, f_patch b) {
        if (a.weight == INF) return true;
        else if (b.weight == INF) return false;
        else return a.weight > b.weight;
    }
}f_patch;


float GeodesicDistance(int f1, int f2, MatrixXf& vertixs, MatrixXi& faces);

float AngleDis(int f1, int f2, MatrixXf& vertixs, MatrixXi faces);

float CalWeight(MatrixXf& vertixs, MatrixXi& faces, MatrixXi& Dual, MatrixXf& weights, float delta = 0.8);

void CalDistMat(int rows, MatrixXi& Dual, MatrixXf& weights, MatrixXf& Dist, int& RepA, int& RepB);

void DistThread(int f1, int step, int rows, const MatrixXi& Dual, const MatrixXf& weights, MatrixXf& Dist);

void minDist(int f1, int rows, const MatrixXi& Dual, const MatrixXf& weights, MatrixXf& Dist);

void argMin(VectorXf& pA, VectorXf& pB, MatrixXf& Dist, int& RepA, int& RepB);

void STDiv(MatrixXf& Cal, vector<int>& C, vector<int>& A, vector<int>& B);