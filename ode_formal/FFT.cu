#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <cufft.h>

#define CUDA_CALL(err) \
    if ((err) != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
        exit(EXIT_FAILURE); \
    }

#define CUFFT_CALL(err) \
    if ((err) != CUFFT_SUCCESS) { \
        fprintf(stderr, "CUFFT error %d at %s:%d\n", (int)(err), __FILE__, __LINE__); \
        exit(EXIT_FAILURE); \
    }

typedef cufftComplex Complex;

const int NX = 128, NY = 128;
const int KX =  64, KY =  64;
const int PX = NX + KX - 1;
const int PY = NY + KY - 1;
const int SIZE = PX * PY;

__global__ void complexMul(Complex *H, const Complex *A, const Complex *B, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        float ar = A[idx].x, ai = A[idx].y;
        float br = B[idx].x, bi = B[idx].y;
        H[idx].x = ar * br - ai * bi;
        H[idx].y = ar * bi + ai * br;
    }
}

__global__ void complexMulAdd(Complex *H, const Complex *A, const Complex *B, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        float ar = A[idx].x, ai = A[idx].y;
        float br = B[idx].x, bi = B[idx].y;
        H[idx].x += ar * br - ai * bi;
        H[idx].y += ar * bi + ai * br;
    }
}

__global__ void normalize(Complex *data, int N, float scale) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        data[idx].x *= scale;
        data[idx].y *= scale;
    }
}

void padAndConvert(const float *in, Complex *out, int wx, int wy) {
    for (int i = 0; i < SIZE; ++i) out[i].x = out[i].y = 0.0f;
    for (int j = 0; j < wy; ++j)
        for (int i = 0; i < wx; ++i)
            out[i + j * PX].x = in[i + j * wx];
}

int main() {
    float *h_Mx = (float*)malloc(NX * NY * sizeof(float));
    float *h_My = (float*)malloc(NX * NY * sizeof(float));
    float *h_Mz = (float*)malloc(NX * NY * sizeof(float));
    float *h_D[9];
    for (int i = 0; i < 9; i++) h_D[i] = (float*)malloc(KX * KY * sizeof(float));

    for (int j = 0; j < NY; j++)
        for (int i = 0; i < NX; i++) {
            int idx = i + j * NX;
            h_Mx[idx] = 1.0f; h_My[idx] = 0.0f; h_Mz[idx] = 0.0f;
        }
    for (int j = 0; j < KY; j++)
        for (int i = 0; i < KX; i++) {
            int idx = i + j * KX;
            h_D[0][idx] = 1.0f; h_D[1][idx] = 0.0f; h_D[2][idx] = 0.0f;
            h_D[3][idx] = 0.0f; h_D[4][idx] = 1.0f; h_D[5][idx] = 0.0f;
            h_D[6][idx] = 0.0f; h_D[7][idx] = 0.0f; h_D[8][idx] = 1.0f;
        }

    Complex *h_Mx_p = (Complex*)malloc(SIZE * sizeof(Complex));
    Complex *h_My_p = (Complex*)malloc(SIZE * sizeof(Complex));
    Complex *h_Mz_p = (Complex*)malloc(SIZE * sizeof(Complex));
    Complex *h_D_p[9];
    for (int i = 0; i < 9; i++) h_D_p[i] = (Complex*)malloc(SIZE * sizeof(Complex));

    padAndConvert(h_Mx, h_Mx_p, NX, NY);
    padAndConvert(h_My, h_My_p, NX, NY);
    padAndConvert(h_Mz, h_Mz_p, NX, NY);
    for (int i = 0; i < 9; i++) padAndConvert(h_D[i], h_D_p[i], KX, KY);

    Complex *d_Mx, *d_My, *d_Mz, *d_Hx, *d_Hy, *d_Hz;
    Complex *d_D[9];
    CUDA_CALL(cudaMalloc(&d_Mx, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMalloc(&d_My, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMalloc(&d_Mz, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMalloc(&d_Hx, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMalloc(&d_Hy, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMalloc(&d_Hz, SIZE * sizeof(Complex)));
    for (int i = 0; i < 9; i++) CUDA_CALL(cudaMalloc(&d_D[i], SIZE * sizeof(Complex)));

    CUDA_CALL(cudaMemcpy(d_Mx, h_Mx_p, SIZE * sizeof(Complex), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_My, h_My_p, SIZE * sizeof(Complex), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Mz, h_Mz_p, SIZE * sizeof(Complex), cudaMemcpyHostToDevice));
    for (int i = 0; i < 9; i++)
        CUDA_CALL(cudaMemcpy(d_D[i], h_D_p[i], SIZE * sizeof(Complex), cudaMemcpyHostToDevice));

    cufftHandle plan;
    CUFFT_CALL(cufftPlan2d(&plan, PY, PX, CUFFT_C2C));

    for (int i = 0; i < 9; i++) CUFFT_CALL(cufftExecC2C(plan, d_D[i], d_D[i], CUFFT_FORWARD));
    CUFFT_CALL(cufftExecC2C(plan, d_Mx, d_Mx, CUFFT_FORWARD));
    CUFFT_CALL(cufftExecC2C(plan, d_My, d_My, CUFFT_FORWARD));
    CUFFT_CALL(cufftExecC2C(plan, d_Mz, d_Mz, CUFFT_FORWARD));

    CUDA_CALL(cudaMemset(d_Hx, 0, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMemset(d_Hy, 0, SIZE * sizeof(Complex)));
    CUDA_CALL(cudaMemset(d_Hz, 0, SIZE * sizeof(Complex)));

    int thr = 256, blk = (SIZE + thr - 1) / thr;
    complexMul <<<blk,thr>>>(d_Hx, d_Mx, d_D[0], SIZE);
    complexMulAdd<<<blk,thr>>>(d_Hx, d_My, d_D[1], SIZE);
    complexMulAdd<<<blk,thr>>>(d_Hx, d_Mz, d_D[2], SIZE);
    complexMul <<<blk,thr>>>(d_Hy, d_Mx, d_D[3], SIZE);
    complexMulAdd<<<blk,thr>>>(d_Hy, d_My, d_D[4], SIZE);
    complexMulAdd<<<blk,thr>>>(d_Hy, d_Mz, d_D[5], SIZE);
    complexMul <<<blk,thr>>>(d_Hz, d_Mx, d_D[6], SIZE);
    complexMulAdd<<<blk,thr>>>(d_Hz, d_My, d_D[7], SIZE);
    complexMulAdd<<<blk,thr>>>(d_Hz, d_Mz, d_D[8], SIZE);

    CUFFT_CALL(cufftExecC2C(plan, d_Hx, d_Hx, CUFFT_INVERSE));
    CUFFT_CALL(cufftExecC2C(plan, d_Hy, d_Hy, CUFFT_INVERSE));
    CUFFT_CALL(cufftExecC2C(plan, d_Hz, d_Hz, CUFFT_INVERSE));
    float scale = 1.0f / SIZE;
    normalize<<<blk,thr>>>(d_Hx, SIZE, scale);
    normalize<<<blk,thr>>>(d_Hy, SIZE, scale);
    normalize<<<blk,thr>>>(d_Hz, SIZE, scale);
    CUDA_CALL(cudaDeviceSynchronize());

    Complex *h_Hx = (Complex*)malloc(SIZE * sizeof(Complex));
    Complex *h_Hy = (Complex*)malloc(SIZE * sizeof(Complex));
    Complex *h_Hz = (Complex*)malloc(SIZE * sizeof(Complex));
    CUDA_CALL(cudaMemcpy(h_Hx, d_Hx, SIZE * sizeof(Complex), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_Hy, d_Hy, SIZE * sizeof(Complex), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_Hz, d_Hz, SIZE * sizeof(Complex), cudaMemcpyDeviceToHost));

    int cx = (PX/2) + (PY/2)*PX;
    printf("Hx(中心) = %f\n", h_Hx[cx].x);
    printf("Hy(中心) = %f\n", h_Hy[cx].x);
    printf("Hz(中心) = %f\n", h_Hz[cx].x);

    CUFFT_CALL(cufftDestroy(plan));
    cudaFree(d_Mx); cudaFree(d_My); cudaFree(d_Mz);
    cudaFree(d_Hx); cudaFree(d_Hy); cudaFree(d_Hz);
    for (int i = 0; i < 9; i++) cudaFree(d_D[i]);
    free(h_Mx); free(h_My); free(h_Mz);
    free(h_Hx); free(h_Hy); free(h_Hz);
    free(h_Mx_p); free(h_My_p); free(h_Mz_p);
    for (int i = 0; i < 9; i++) { free(h_D[i]); free(h_D_p[i]); }
    return 0;
}

