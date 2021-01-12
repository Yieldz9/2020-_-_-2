
#include "cuda_main.cuh"

//a*rows+x
void matrix_mul_cpu(int i, int j, float* M, int* N, float* P, int width)
{

    //if (j==0) printf("SP:%d %d\n", i, j);
    int k;
    int x;
    int front = -1, rear = -1, size = 0;
    f_patch* Patches = (f_patch*)malloc(sizeof(f_patch) * width);
    int* Q = (int*)malloc(sizeof(int) * width);
    for (k = 0; k < width; k++) {
        Patches[k].weight = M[i * width + k];
        Patches[k].IsDiscover = 0;
        Q[k] = INF;
    }
    Q[0] = i;
    front = 0;
    rear = 0;
    size = 1;
    while (size != 0) {
        x = Q[front]; front = (front + 1) % width; size--;
        if (Patches[x].IsDiscover) continue;
        if (x == j) { P[i * width + j] = Patches[j].weight; return; }
        Patches[x].IsDiscover = 1;
        for (int a = 0; a < 3; a++) {
            if (Patches[N[a * 4000 + x]].weight == -10) { Patches[N[a * 4000 + x]].weight = Patches[x].weight + M[x * width + N[a * 4000 + x]]; }
            else Patches[N[a * 4000 + x]].weight = Patches[N[a * 4000 + x]].weight <= Patches[x].weight + M[x * width + N[a * 4000 + x]]
                ? Patches[N[a * 4000 + x]].weight : Patches[x].weight + M[x * width + N[a * 4000 + x]];
        }
        int max = -1;
        for (int a = 0; a < width; a++) {
            if (!Patches[a].IsDiscover) {
                if (max == -1) max = a;
                else if (Patches[a].weight!=-10) max = Patches[max].weight <= Patches[a].weight ? max : a;
            }
        }
        rear = (rear + 1) % width; Q[rear] = max; size++;
    }
}

__global__ void matrix_mul_gpu(float* M, int* N, float* P, int width)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = threadIdx.y + blockDim.y * blockIdx.y;
    //if (j==0) printf("SP:%d %d\n", i, j);
    int k;
    int x;
    int front = -1, rear = -1, size = 0;
    f_patch* Patches = (f_patch*)malloc(sizeof(f_patch) * width);
    int* Q = (int*)malloc(sizeof(int) * width);
    for (k = 0; k < width; k++) {
        Patches[k].weight = M[i * width + k];
        Q[k] = INF;
    }
    Q[0] = i;
    front = 0;
    rear = 0;
    size = 1;
    while (size != 0) {
        x = Q[front]; front = (front + 1) % width; size--;
        if (Patches[x].IsDiscover) continue;
        if (x == j) { P[i * width + j] = Patches[j].weight; return; }
        Patches[x].IsDiscover = 1;
        for (int a = 0; a < 3; a++) {
            if (Patches[N[x * 3 + a]].weight == -10) { Patches[N[x * 3 + a]].weight = Patches[x].weight + M[x * width + N[x * 3 + a]]; }
            else Patches[N[x * 3 + a]].weight = Patches[N[x * 3 + a]].weight <= Patches[x].weight + M[x * width + N[x * 3 + a]]
                ? Patches[N[x * 3 + a]].weight : Patches[x].weight + M[x * width + N[x * 3 + a]];
        }
        int max = -1;
        for (int a = 0; a < width; a++) {
            if (!Patches[a].IsDiscover) {
                if (max == -1) max = a;
                else if (Patches[a].weight != -10) max = Patches[max].weight <= Patches[a].weight ? max : a;
            }
        }
        rear = (rear + 1) % width; Q[rear] = max; size++;
    }
}


void GPUForDis(int Row, int Col, float* weights , int* Dual, float* Dist)
{
    //struct timeval start, end;
    //gettimeofday(&start, NULL);
    /*
    float* weights = (float*)malloc(sizeof(float) * Row * Col);
    int* Dual = (int*)malloc(sizeof(int) * Row * 3);
    float* Dist = (float*)malloc(sizeof(float) * Row * Col);*/
    //malloc device memory
    float* d_dataWei, * d_dataDist;
    int* d_dataDual;
    /*cudaMalloc((void**)&d_dataWei, sizeof(float) * Row * Col);
    cudaMalloc((void**)&d_dataDual, sizeof(int) * Row * 3);
    cudaMalloc((void**)&d_dataDist, sizeof(float) * Row * Col);
    //set value
    

    cudaMemcpy(d_dataWei, weights, sizeof(float) * Row * Col, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dataDual, Dual, sizeof(int) * Row * 3, cudaMemcpyHostToDevice);
    //dim3 threadPerBlock(16, 16);
    dim3 threadPerBlock(16, 16);
    //dim3 blockNumber(2,2);
    dim3 blockNumber((Col + threadPerBlock.x - 1) / threadPerBlock.x, (Row + threadPerBlock.y - 1) / threadPerBlock.y);
    printf("Block(%d,%d)   Grid(%d,%d).\n", threadPerBlock.x, threadPerBlock.y, blockNumber.x, blockNumber.y);
    matrix_mul_gpu << <blockNumber, threadPerBlock >> > (d_dataWei, d_dataDual, d_dataDist, Col);
    //拷贝计算数据-一级数据指针
    cudaMemcpy(Dist, d_dataDist, sizeof(float) * Row * Col, cudaMemcpyDeviceToHost);*/
    for (int i = 0; i < 4000; i++) {
        for (int j = 0; j < 4000; j++) {
            matrix_mul_cpu(i, j, weights, Dual,Dist, 4000);
        }
    }
    for (int i = 0; i < 4; i++) {
        printf("%f %f %f %f\n", Dist[i * 4 + 0], Dist[i * 4 + 1], Dist[i * 4 + 2], Dist[i * 4 + 3]);
    }
    //释放内存
    
    cudaFree(d_dataWei);
    cudaFree(d_dataDual);
    cudaFree(d_dataDist);

    /*gettimeofday(&end, NULL);
    int timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("total time is %d ms\n", timeuse / 1000);*/

}