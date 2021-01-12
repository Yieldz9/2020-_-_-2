#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define INF -1;
typedef struct f_patch {
    float weight;
    int IsDiscover = 0;
}f_patch;

void matrix_mul_cpu(int i, int j, float* M, int* N, float* P, int width);
__global__ void matrix_mul_gpu(float* M, int* N, float* P, int width);
void GPUForDis(int Row, int Col, float* weights, int* Dual, float* Dist);