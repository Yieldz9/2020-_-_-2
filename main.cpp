#include <iostream>
#include <Eigen\Dense>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "tool.h"


using namespace Eigen;
using namespace std;


//未完成
//eta值的选择
//对偶图可以遍历边啊，你个大SB
void BuildGraph(ifstream& fin, MatrixXf& vertexs, MatrixXi& faces, MatrixXi& DualGraph, vector<Edge>& E, bool IsCreat) {
    
    char buffer[128];
    int v = 0;
    int f = 0;
    
    fin.open("C:/project/C++/dinosaur.2k.obj");
    fin.getline(buffer, 128);
    for (int i = 0; i < 128; i++) {
        if (buffer[i] == ' ') if (buffer[i + 1] >= '0' && buffer[i + 1] <= '9') if (v == 0) v = atoi(buffer + i + 1); else f = atoi(buffer + i + 1);
        if (f != 0) break;
    }
    faces.resize(f, 3);
    DualGraph.resize(f, 3);
    DualGraph.setConstant(-1);
    vertexs.resize(v, 3);
    //vertexs
    for (int i = 0; i < v;i++) {
        fin.getline(buffer, 128);
        int k = 0;
        for (int j = 0; j < 128; j++) {
            if (buffer[j] == ' ') {
                vertexs(i,k) = (float)atof(buffer + j + 1); k++;
            }
            if (k == 3) break;
        }
    }
    //faces
    for (int i = 0; i < f;i++) {
        int k = 0;
        fin.getline(buffer, 128);
        for (int j = 0; j < 128;j++) {
            if (buffer[j] == ' ') {
                faces(i, k) = atoi(buffer + j + 1) - 1 ; k++;
            }
            if (k == 3) break;
        }
    }
    //DualGraph
    Edge e[3];
    for (int i = 0; i < f; i++) {
            for (int j = 0; j < 3; j++) {
                e[j].start = faces(i, j); e[j].end = faces(i, (j + 1) % 3);
            }
            for (int j = 0; j < 3; j++) {
                int k;
                for (k = 0; k < E.size(); k++) {
                    if ((E[k].start == e[j].start && E[k].end == e[j].end) ||
                        (E[k].end == e[j].start && E[k].start == e[j].end)) {
                        if (E[k].f2 == -1) {
                            E[k].f2 = i;
                            for (int m = 0; m < 3; m++) {
                                if (DualGraph(E[k].f1, m) == -1) { DualGraph(E[k].f1, m) = E[k].f2; break; }

                            }
                            for (int m = 0; m < 3; m++) {
                                if (DualGraph(E[k].f2, m) == -1) {
                                    DualGraph(E[k].f2, m) = E[k].f1; break;
                                }
                            }
                            break;
                        }
                    }

                }
                if (k == E.size()) { e[j].f1 = i; E.push_back(e[j]); }
            }
        }
    
}



void output_final(ofstream& fout, MatrixXf& vertexs, MatrixXi& faces, vector<int>& A, vector<int>& B, int& RepA, int& RepB) {

    int flag = 0;
    int i, j=0, k=0;
    for (i = 0; i < vertexs.rows(); i++) {
        if (i == faces(RepA, 0) || i == faces(RepA, 1) || i == faces(RepA, 2)) { fout << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 255 255 255" << endl; }
        else if (i == faces(RepB, 0) || i == faces(RepB, 1) || i == faces(RepB, 2)) { fout << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 255 255 255" << endl; }
        else {
            flag = 0;
            for (j; j < A.size(); j++) {
                if (A[j] == i) { fout << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 255 0 0" << endl; flag = 1; }
                else if (A[j] > i) break;
            }
            if (!flag) for (k; k < B.size(); k++) {
                if (B[k] == i) { fout << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 0 255 0" << endl; flag = 1; }
                else if (B[k] > i) break;
            }
            if (!flag) fout << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 0 0 255" << endl;
        }
    }
    cout << "A:" << A.size() << " B:" << B.size() << endl;
    for (i = 0; i < faces.rows(); i++) {
        fout << "f " << faces(i, 0) + 1 << ' ' << faces(i, 1) + 1 << ' ' << faces(i, 2) + 1 << endl;
    }
}
int main() {
    
    float delta = 0.8;
    float Xi = 0.04;//0.02
    int epoch = 50;

    MatrixXi faces;
    MatrixXf vertexs;
    MatrixXi DualGraph;
    MatrixXf weights;
    MatrixXf Distance;
    MatrixXf STCap;
    //vector<vector<STvertex> > STCap;
    VectorXf pA, pB;
    vector<Edge> Edges;
    vector<int> FC, A, B, C;
    vector<int> S, T, STvertexs;
    int RepA, RepB;
    int temp_repa = 0, temp_repb = 0;
    float avgAngDis = 0;
    RepA = 0; RepB = 0;
    ifstream fin;

    BuildGraph(fin, vertexs,faces, DualGraph,Edges, false);
    pA.resize(faces.rows());
    pB.resize(faces.rows());
    
    avgAngDis = CalWeight(vertexs, faces, DualGraph, weights, delta);
    cout << "Weight over" << endl;

    CalDistMat(faces.rows(), DualGraph, weights, Distance, RepA, RepB);
    cout << "Dist over" << endl;
    cout <<0<< " RepA: " << RepA << "RepB: " << RepB << endl;
    //cout << Distance;
    //迭代聚类
    for (int j = 0; j < epoch;j++) {
        int A=0, B=0;
        for (int i = 0; i < faces.rows();i++) {
            pA(i) = Distance(i, RepA) / (Distance(i, RepA) + Distance(i, RepB));
            if (pA(i) > 0.5 + Xi) A++;
            pB(i) = Distance(i, RepB) / (Distance(i, RepA) + Distance(i, RepB));
            if (pB(i) > 0.5 + Xi) B++;
        }
        cout << "A: " << A << "B: " << B << endl;
        temp_repa = RepA; 
        temp_repb = RepB;
        argMin(pA, pB, Distance, RepA, RepB);
        //if (j != epoch-1 && Distance(RepA, RepB) < (faces.rows() / 100) * Distance.mean()) { RepA += rand()%50; RepB += rand()%50; }
        cout << j << " RepA: " << RepA << " RepB: " << RepB << endl;
        if ((temp_repb == RepB && temp_repa == RepA)||(temp_repa == RepB && temp_repb == RepA)) break;

    }
    //输出
    /*
    ofstream fout;
    fout.open("C:/project/C++/out_ini_dinasaur_0.05.obj");
    output_ini(Xi, fout, vertexs, faces, pA, pB);*/
    //区域点集A,B,C
    for (int i = 0; i < faces.rows(); i++) {
        int flag = 0;
        if (pA[i] > 0.5 + Xi) {
            //内部的点
            for (int j = 0; j < 3; j++) {
                if (pA[DualGraph(i, j)] <= 0.5 + Xi && pA[DualGraph(i, j)] >= 0.5 - Xi) { 
                    if (flag == 2) flag = 3;
                    else flag = 1;
                }
                else if (pA[DualGraph(i, j)] < 0.5 - Xi) {
                    if (flag == 1) flag = 3;
                    else flag = 2;
                }
            }
            if (flag == 2 || flag == 0) {
                for (int j = 0; j < 3; j++) {
                    A.push_back(faces(i, j));
                    
                }
            }

        }
        else if (pA[i] < 0.5 - Xi) {
            for (int j = 0; j < 3; j++) {
                if (pA[DualGraph(i, j)] >= 0.5 - Xi) { flag = 1; break; }
            }
            if (!flag) {
                for (int j = 0; j < 3; j++) {
                    B.push_back(faces(i, j));
                    
                }
            }
        }
        else { 
            for (int j = 0; j < 3; j++) { C.push_back(faces(i, j)); } 
            FC.push_back(i);
        }
    }
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    sort(C.begin(), C.end());
    A.erase(unique(A.begin(), A.end()), A.end());
    B.erase(unique(B.begin(), B.end()), B.end());
    C.erase(unique(C.begin(), C.end()), C.end());
    int i = 0;
    
    for (int j = 0; j < C.size(); j++) {
        i = 0;
        while (i < A.size()) {
            if (A[i] == C[j]) {
                A.erase(A.begin() + i);
                break;
            }
            else i++;
        }
    }
    
    for (int j = 0; j < C.size(); j++) {
        i = 0;
        while (i < B.size()) {
            if (B[i] == C[j]) {
                B.erase(B.begin() + i);
                break;
            }
            else i++;
        }
    }
    
    //S,T
    for (int i = 0; i < FC.size(); i++) {
        for (int j = 0; j < 3; j++) {
            if (pA[DualGraph(FC[i], j)] > 0.5 + Xi) {
                int f1 = DualGraph(FC[i], j), f2 = FC[i];
                for (int k = 0; k < 3; k++) {
                    if (faces(f1, k) != faces(f2, 0) && faces(f1, k) != faces(f2, 1) && faces(f1, k) != faces(f2, 2)) S.push_back(faces(f1, k));
                }
            }
            if (pA[DualGraph(FC[i], j)] < 0.5 - Xi) {
                int f1 = DualGraph(FC[i], j), f2 = FC[i];
                for (int k = 0; k < 3; k++) {
                    if (faces(f1, k) != faces(f2, 0) && faces(f1, k) != faces(f2, 1) && faces(f1, k) != faces(f2, 2)) T.push_back(faces(f1, k));
                }
            }
        }
    }
    FC.clear();
    sort(S.begin(), S.end());
    S.erase(unique(S.begin(), S.end()), S.end());
    sort(T.begin(), T.end());
    T.erase(unique(T.begin(), T.end()), T.end());
    int j = 0;
    for (int i = 0; i < C.size(); i++) {
        j = 0;
        while (j < S.size()) {
            if (C[i] == S[j]) {
                S.erase(S.begin() + j);
                break;
            }
            else j++;
        }
    }
    for (int i = 0; i < C.size(); i++) {
        j = 0;
        while (j < T.size()) {
            if (C[i] == T[j]) {
                T.erase(T.begin() + j);
                break;
            }
            else j++;
        }
    }

    //ST图邻接矩阵，还是用邻接表吧，虽然改值的时候会麻烦一点
    STCap.resize(C.size() + 2, C.size() + 2);
    STCap.setConstant(0);
    for (int i = 0; i < C.size() + 2; i++) {
        Edge temp;
        STvertex v;
        STCap(i, i) = 0;
        //S到各点距离
        if (i == 0) {
            
            for (int j = 0; j < S.size(); j++) {
                for (int k = 0; k < C.size(); k++) {
                    temp.end = S[j]; temp.start = C[k];
                    for (int m = 0; m < Edges.size(); m++) {
                        if ((temp.start == Edges[m].start && temp.end == Edges[m].end) ||
                            (temp.end == Edges[m].start && temp.start == Edges[m].end)) {
                            /*v.idx = k + 1;
                            v.Cal = 1 / (1 + AngleDis(Edges[m].f1, Edges[m].f2, vertexs, faces) / avgAngDis);
                            STCap[0].push_back(v);
                            v.idx = i;
                            STCap[k + 1].push_back(v);*/
                            STCap(0, k + 1) = -1;
                            STCap(k + 1, 0) = -1;
                            //STCap(0, k + 1) = 1 / (1 + AngleDis(Edges[m].f1, Edges[m].f2, vertexs, faces) / avgAngDis); 
                            //STCap(k + 1, 0) = STCap(0, k + 1);
                        }
                    }
                }
            }
        }
        //C点集之间的距离
        else if (i != C.size() + 1) {
            
            for (int j = i + 1; j < C.size() + 1; j++) {
                temp.start = C[i - 1]; temp.end = C[j - 1];
                for (int k = 0; k < Edges.size(); k++) {
                    if ((temp.start == Edges[k].start && temp.end == Edges[k].end) ||
                        (temp.end == Edges[k].start && temp.start == Edges[k].end)) {
                        /*v.idx = j + 1;
                        v.Cal = 1 / (1 + AngleDis(Edges[k].f1, Edges[k].f2, vertexs, faces) / avgAngDis);
                        STCap[i].push_back(v);
                        v.idx = i;
                        STCap[j + 1].push_back(v);*/
                        STCap(i, j) = 1 / (1 + AngleDis(Edges[k].f1, Edges[k].f2, vertexs, faces) / avgAngDis);
                        STCap(j, i) = STCap(i, j);
                    }
                }
                
            }
        }
        //T到各点距离
        else {
            for (int j = 0; j < T.size(); j++) {
                for (int k = 0; k < C.size(); k++) {
                    temp.end = T[j]; temp.start = C[k];
                    for (int m = 0; m < Edges.size(); m++) {
                        if ((temp.start == Edges[m].start && temp.end == Edges[m].end) ||
                            (temp.end == Edges[m].start && temp.start == Edges[m].end)) {
                            /*v.idx = k + 1;
                            v.Cal = 1 / (1 + AngleDis(Edges[m].f1, Edges[m].f2, vertexs, faces) / avgAngDis);
                            STCap[0].push_back(v);
                            v.idx = i;
                            STCap[k + 1].push_back(v);*/
                            STCap(i, k + 1) = -1;
                            STCap(k + 1, i) = -1;
                            //STCap(k + 1, i) = 1 / (1 + AngleDis(Edges[m].f1, Edges[m].f2, vertexs, faces) / avgAngDis);
                            //STCap(i, k + 1) = STCap(k + 1, i);
                        }
                    }
                }
            }
        }
    }
    vector<int> CA, CB,FA,FB;
    //cout << STCap;
    STDiv(STCap, C, CA, CB);
    FA.resize(A.size() + CA.size());
    FB.resize(B.size() + CB.size());
    merge(A.begin(), A.end(), CA.begin(), CA.end(), FA.begin());
    merge(B.begin(), B.end(), CB.begin(), CB.end(), FB.begin());
    
    
    ofstream temp("C:/project/C++/lines.obj");
    j = 0;int k = 0, p = 0;
    for (int i = 0; i < vertexs.rows(); i++) {
        if (p < C.size() && i == C[p]) {
            temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 0 0 1" << endl;
            p++;
        }
        else if (j < S.size() && i == S[j]) {
            temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 1 0 0" << endl;
            j++;
        }
        else if (k < T.size() && i == T[k]) {
            temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 0 1 0" << endl;
            k++;
        }
        else temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 1 1 1" << endl;
    }
    for (int i = 0; i < faces.rows(); i++) {
        temp << "f " << faces(i, 0) + 1 << ' ' << faces(i, 1) + 1 << ' ' << faces(i, 2) + 1 << endl;
    }
    temp.close();
    temp.open("C:/project/C++/line.obj");
    j = 0; k = 0;
    for (int i = 0; i < vertexs.rows(); i++) {
       
        if (j < CA.size() && i == CA[j]) {
            temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 1 0 0" << endl;
            j++;
        }
        else if (k < CB.size() && i == CB[k]) {
            temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 0 1 0" << endl;
            k++;
        }
        else temp << "v " << vertexs(i, 0) << ' ' << vertexs(i, 1) << ' ' << vertexs(i, 2) << " 1 1 1" << endl;
    }
    for (int i = 0; i < faces.rows(); i++) {
        temp << "f " << faces(i, 0) + 1 << ' ' << faces(i, 1) + 1 << ' ' << faces(i, 2) + 1 << endl;
    }
    temp.close();



    
    //补充确定区域的A，B点集,这个点集分割还是有一点问题
    
  
    
    ofstream fout;
    fout.open("C:/project/C++/out_dinasaur_0.05.obj");
    output_final(fout, vertexs, faces, A, B,RepA,RepB);
    fout.close();
    fout.open("C:/project/C++/finalout_dinasaur_0.05.obj");
    output_final(fout, vertexs, faces, FA, FB,RepA,RepB);
    fout.close();
    
    
    //最小割未实现
    //最后的输出也未实现
    return 0;
}