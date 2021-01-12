#include "tool.h"
#include <thread>

//#include "cuda_main.cuh"


float GeodesicDistance(int f1, int f2, MatrixXf& vertixs, MatrixXi& faces) {
    Vector3f M, N, e;
    int A = -1, B = -1, C = -1, D = -1;
    for (int i = 0; i < 3;i++) {
        for (int j = 0; j < 3;j++) {
            if (faces(f1, i) == faces(f2, j)) if (A == -1) A = faces(f1, i); else B = faces(f1, i);
            if (B != -1) break;
        }
    }
    for (int i = 0; i < 3; i++) {
        if (A != faces(f1, i) && B != faces(f1, i)) C = faces(f1, i);
        if (A != faces(f2, i) && B != faces(f2, i)) D = faces(f2, i);
    }
    if (A == -1 || B == -1) return -10;
    /*if (A == 2002 || B == 2002 || C == 2002 || D == 2002) {
        cout << endl;
    }*/
    else {
        for (int i = 0; i < 3; i++) {
            M(i) = (vertixs(A, i) + vertixs(C, i) - 2 * vertixs(B, i)) / 3;
            N(i) = (vertixs(A, i) + vertixs(D, i) - 2 * vertixs(B, i)) / 3;
            e(i) = vertixs(A, i) - vertixs(B, i);
        }
        e.normalize();
        float temp1 = M.cross(e).norm() + e.cross(N).norm();

        return sqrt(temp1 * temp1 + (M.dot(e) - N.dot(e)) * (M.dot(e) - N.dot(e)));
    }
}

float AngleDis(int f1,int f2, MatrixXf& vertixs, MatrixXi faces) {
    Vector3f t1, t2, t3, t4, n1, n2, flag;
    float eta;
    int i;
    int A = -1, B = -1, C = -1, D = -1;
    for (int i = 0; i < 3;i++) {
        for (int j = 0; j < 3;j++) {
            if (faces(f1, i) == faces(f2, j)) if (A == -1) A = faces(f1, i); else B = faces(f1, i);
            if (B != -1) break;
        }
    }
    for (int i = 0; i < 3; i++) {
        if (A != faces(f1, i) && B != faces(f1, i)) C = faces(f1, i);
        if (A != faces(f2, i) && B != faces(f2, i)) D = faces(f2, i);
    }
    if (A == -1 || B == -1) return -10;
    /*if (A == 2002 || B == 2002 || C == 2002 || D == 2002) {
        cout << endl;
    }*/
    for (i = 0; i < 3;i++) {
        t1(i) = vertixs(faces(f1,1), i) - vertixs(faces(f1,0), i);
        t2(i) = vertixs(faces(f1,2), i) - vertixs(faces(f1,1), i);
        t3(i) = vertixs(faces(f2,1), i) - vertixs(faces(f2,0), i);
        t4(i) = vertixs(faces(f2,2), i) - vertixs(faces(f2,1), i);
        flag(i) = vertixs(A, i) - (vertixs(C, i) + vertixs(D, i)) / 2;
    }
    n1 = t1.cross(t2);
    n2 = t3.cross(t4);
    if (flag.dot(n1) > 0) eta = 1;
    else eta = 0.8;
    return eta * (1 - n1.dot(n2) / (n1.norm() * n2.norm()));
}

//这里Dual里存的序号应该从0开始,权重矩阵,可以考虑用稀疏矩阵的方法处理？
float CalWeight(MatrixXf& vertixs, MatrixXi& faces, MatrixXi& Dual, MatrixXf& weights, float delta) {
    //计算平均距离
    int i, j, rows;
    float avgGeodis = 0;
    float avgAgdis = 0;
    rows = faces.rows();
    for (i = 0; i < rows; i++) {
        for (j = 0; j < 3; j++) {
            avgGeodis = avgGeodis + GeodesicDistance(i,Dual(i,j), vertixs,faces);
            avgAgdis = avgAgdis + AngleDis(i, Dual(i, j), vertixs, faces);
            /*if (GeodesicDistance(i, Dual(i, j), vertixs, faces) == -10) {
                cout << "Geo:" << i << ' ' << Dual(i, j) << endl;
                cout << "face i:" << faces(i, 0) << ' ' << faces(i, 1) << ' ' << faces(i, 2) << endl;
                cout << "face Di:" << faces(Dual(i, j), 0) << ' ' << faces(Dual(i, j), 1) << ' ' << faces(Dual(i, j), 2) << endl;
            }
            if (AngleDis(i, Dual(i, j), vertixs, faces) == -10) {
                cout << "Angel:" << i << ' ' << Dual(i, j) << endl;
            }*/
            
        }
    }
    avgGeodis /= 3 * rows;
    avgAgdis /= 3 * rows;
    weights.resize(rows, rows);
    weights.setConstant(-10);
    //这里可以优化,对称的特点
    for (i = 0; i < rows;i++) {
        for (j = 0; j < 3;j++) {
            /*if (i == 4000 || Dual(i, j) == 4000) {
                cout << "";
            }
            if (i == 0 && Dual(i, j) == 3) {
                cout << endl;
            }*/
            weights(i, Dual(i, j)) = delta * GeodesicDistance(i, Dual(i, j), vertixs, faces) / avgGeodis
                + (1 - delta) * AngleDis(i, Dual(i, j), vertixs, faces) / avgAgdis;
            
                    
        }
        weights(i, i) = 0;
    }
    return avgAgdis;
}


//计算f1,f2之间的最短距离

void minDist(int f1, int rows, const MatrixXi& Dual, const MatrixXf& weights, MatrixXf& Dist) {
    //mt.lock();
    vector<f_patch> Patches;
    int num = rows-f1;
    Patches.resize(rows);
    for (int i = 0; i < rows; i++) {
        Patches[i].idx = i;
        Patches[i].weight = weights(f1, i);
    }
    priority_queue<f_patch> Q;
    Q.push(Patches[f1]);
    f_patch x;
    while (!Q.empty()) {
        x = Q.top();
        Q.pop();
        if (Patches[x.idx].isDiscover) continue;
        if (x.idx >= f1) {
            Dist(f1, x.idx) = Patches[x.idx].weight;
            Dist(x.idx, f1) = Patches[x.idx].weight;
            num--;
            if (num == 0) return;
        }
        Patches[x.idx].isDiscover = true;
        for (int i = 0; i < 3;i++) {
            if (Patches[Dual(x.idx, i)].weight == -10) Patches[Dual(x.idx, i)].weight = Patches[x.idx].weight + weights(Dual(x.idx, i), x.idx);
            else {
                Patches[Dual(x.idx, i)].weight = Patches[Dual(x.idx, i)].weight <= Patches[x.idx].weight + weights(Dual(x.idx, i), x.idx)
                    ? Patches[Dual(x.idx, i)].weight : Patches[x.idx].weight + weights(Dual(x.idx, i), x.idx);
            }
            if (!Patches[Dual(x.idx, i)].isDiscover) Q.push(Patches[Dual(x.idx, i)]);
        } 
    }
    //mt.unlock();
}

void DistThread(int f1, int step, int rows, const MatrixXi& Dual, const MatrixXf& weights, MatrixXf& Dist) {
    for (int i = f1; i < rows; i += step) {
        minDist(i, rows, Dual, weights, Dist);
    }
}

//计算得到距离矩阵和相距最远的初始种子
//cuda加速
void CalDistMat( int rows, MatrixXi& Dual, MatrixXf& weights, MatrixXf& Dist, int& RepA, int& RepB) {
    float max = 0;
    float dist;
    Dist.resize(rows, rows);
    //GPUForDis(rows, rows, (float*)weights.data(), (int*)Dual.data(), (float*)Dist.data());
    
        int n = thread::hardware_concurrency();
        vector<thread> T(n - 1);
        for (int i = 0; i < n - 1; i++) {
            T[i] = thread(DistThread, i, n - 1, rows, Dual, weights, ref(Dist));
        }
        for (int i = 0; i < n - 1; i++) {
            T[i].join();
        }
        /*
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if (Dist(i, j) < 0) cout <<"false" << endl;
            }
        }*/
        //minDist(i, rows, Dual, weights, Dist);
        
        
        //cout << max << '-' << RepA <<'-'<< RepB << endl;
    
    Dist.maxCoeff(&RepA, &RepB);
}

//迭代更新聚类中心
void argMin(VectorXf& pA, VectorXf& pB, MatrixXf& Dist, int& RepA, int& RepB) {
    /*float sumA=0, sumB=0,maxA=0,maxB=0;
    for (int i = 0; i < pA.rows(); i++) {
        sumA = 0; sumB = 0;
        for (int j = 0; j < pA.rows(); j++) {
            sumA += pA(j) * Dist(i, j);
            sumB += pB(j) * Dist(i, j);
        }
        if (i == 0) { maxA = sumA; maxB = sumB; RepA = i, RepB = i; }
        else {
            if (maxA > sumA) { maxA = sumA; RepA = i; }
            if (maxB > sumB) { maxB = sumB; RepB = i; }
        }
    }*/
    VectorXf temp = pA.transpose() * Dist;
    temp.minCoeff(&RepA);
    temp = pB.transpose() * Dist;
    temp.minCoeff(&RepB);
}

//注意Cal的0行和最后一行是S和T
//未调试，复杂度有点高，使用邻接表
void STDiv(MatrixXf& Cal, vector<int>& C, vector<int>& A, vector<int>& B) {
    //DPS找到分割线，广搜找到A，B类的点
    queue<int> Q;
    VectorXi pre;
    VectorXf flow;
    int i, x, reg_x = 0;
    float maxFlow = -1;
    pre.resize(Cal.rows());
    flow.resize(Cal.rows());
    //不断BFS求最小割
    while (1) {
        while (!Q.empty()) {
            Q.pop();
        }
        pre.setConstant(-1);
        
        pre(0) = 0;//S点的前驱没有意义
        flow(0) = -1;
        Q.push(0);
        while (!Q.empty()) {
            int x = Q.front(); Q.pop();
            if (x == Cal.rows() - 1) break;
            for (i = 1; i < Cal.rows(); i++) {
                if (Cal(x, i) != 0 && pre(i) == -1) {
                    pre(i) = x;
                    if (flow(x) == -1) flow(i) = Cal(x, i);
                    else if (Cal(x, i) == -1) flow(i) = flow(x);
                    else flow(i) = min(flow(x), Cal(x, i));
                    Q.push(i);
                }
            }
        }
        
        if (pre(Cal.rows() - 1) == -1) break;
        else {
            int j = Cal.rows() - 1;
            while (j != 0) {
                if (flow(Cal.rows() - 1) == -1) {
                    Cal(pre(j), j) -= flow(Cal.rows() - 1);
                    Cal(j, pre(j)) -= flow(Cal.rows() - 1);
                }
                else if(Cal(pre(j), j)!=-1){
                    Cal(pre(j), j) -= flow(Cal.rows() - 1);
                    Cal(j, pre(j)) -= flow(Cal.rows() - 1);
                }
                
                
                j = pre(j);
            }
        }
        
    }
    while (!Q.empty()) {
        Q.pop();
    }
    pre.setConstant(-1);
    pre(0) = 0;//S点的前驱没有意义
    flow(0) = -1;
    Q.push(0);
    while (!Q.empty()) {
        int x = Q.front(); Q.pop();
        for (i = 1; i < Cal.rows(); i++) {
            if (Cal(x, i) != 0 && pre(i) == -1) {
                pre(i) = x;
                if (flow(x) == -1) flow(i) = Cal(x, i);
                else flow(i) = min(flow(x), Cal(x, i));
                Q.push(i);
            }
        }
    }
    for (i = 1; i < Cal.rows() - 1; i++) {
        if (pre(i) == -1) B.push_back(C[i - 1]);
        else A.push_back(C[i - 1]);
    }
    

    /*
    //反复深搜
    //0是未被访问，1是不能走的点，2是已经走过的点
    while (1) {
        state.setConstant(0);
        maxCal.clear();
        Q.push(0);
        while (!Q.empty()) {
            x = Q.top();
            if (x == Cal.rows() - 1) {
                reg_x = Cal.rows() - 1; Q.pop(); break;
            }
            
            for (int i = 0; i < Cal.cols(); i++) {
                if ((state(i)==0) && Cal(x, i) != 0) {
                    Q.push(i);
                    state(i) = 1;
                    maxCal.push_back(Cal(x, i));
                    break;
                }
                else if (i == Cal.cols() - 1) { 
                    Q.pop();
                    if (x != 0) maxCal.pop_back();
                    state(x) = 1;
                }
            }
        }
        if (maxCal.size() != 0) sort(maxCal.begin(), maxCal.end());
        if (reg_x != 0) while (!Q.empty()) {
            x = Q.top(); Q.pop();
            Cal(x, reg_x) -= maxCal[0];
            Cal(reg_x, x) = Cal(x, reg_x);
            reg_x = x;
        }
        else break;
    }
    //广搜
    //这里状态有访问,未被访问,已进入A，B
    
    P.push(0);
    state.setConstant(0);
    while (!P.empty()) {
        x = P.front(); P.pop();
        state(x) = 1;
        for (int i = 1; i < Cal.cols() - 1; i++) {
            if ((state(i) == 0) && Cal(x, i) != -10 && Cal(x, i) != 0) { A.push_back(C[i - 1]); state(i) = 2; P.push(i); }
        }
    }
    sort(A.begin(), A.end());
    int j = 0;
    for (int i = 0; i < C.size(); i++) {
        if (j<A.size() && C[i] == A[j]) j++;
        else B.push_back(C[i]);
    }*/

}



//测试
/*
int main() {
    MatrixXf vertixs;
    MatrixXi faces;
    vertixs.resize(4, 3);
    faces.resize(4, 3);
    vertixs << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    faces << 0, 1, 2, 0, 1, 3, 3, 2, 1, 0, 3, 2;
    //faces << 0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2;
    cout << GeodesicDistance(0, 2, vertixs, faces)<<endl;
    cout << AngleDis(0, 2, vertixs, faces)<<endl;
    return 0;
}*/