/* Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cufft.h>

typedef cufftComplex Complex; // cufftComplex is float2

// Pad data
int PadData(const Complex *signal,
            Complex      **padded_signal,
            int NX, int NY,
            const Complex *filter_kernel,
            Complex      **padded_filter_kernel,
            int KX, int KY, int &PX, int &PY)
{
    PX = NX + KX - 1;
    PY = NY + KY - 1;
    int Nfft = PX * PY;

    *padded_signal = (Complex*)malloc(sizeof(Complex) * Nfft);
    *padded_filter_kernel = (Complex*)malloc(sizeof(Complex) * Nfft);

    for (int y = 0; y < PY; y++) {
        for (int x = 0; x < PX; x++) {
            // pad signal
            if (x < NX && y < NY)
                (*padded_signal)[y * PX + x] = signal[y * NX + x];
            else
                (*padded_signal)[y * PX + x] = make_cuFloatComplex(0, 0);
            // pad kernel
            if (x < KX && y < KY)
                (*padded_filter_kernel)[y * PX + x] = filter_kernel[y * KX + x];
            else
                (*padded_filter_kernel)[y * PX + x] = make_cuFloatComplex(0, 0);
        }
    }
    return Nfft;
}

// CPU 端参考卷积
void Convolve(const Complex *signal, int NX, int NY,
              const Complex *filter, int KX, int KY,
              Complex *output)
{
    int ox = KX / 2, oy = KY / 2;
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            Complex acc = make_cuFloatComplex(0, 0);
            for (int ky = 0; ky < KY; ky++) {
                for (int kx = 0; kx < KX; kx++) {
                    int ix = x + kx - ox;
                    int iy = y + ky - oy;
                    if (ix >= 0 && ix < NX && iy >= 0 && iy < NY) {
                        Complex a = signal[iy * NX + ix];
                        Complex b = filter[(KY - 1 - ky) * KX + (KX - 1 - kx)];
                        acc.x += a.x * b.x - a.y * b.y;
                        acc.y += a.x * b.y + a.y * b.x;
                    }
                }
            }
            output[y * NX + x] = acc;
        }
    }
}

// CUDA 核函数：点乘 & 缩放
__global__ void ComplexPointwiseMulAndScale(Complex *a, const Complex *b, int size, float scale)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    for (int i = tid; i < size; i += stride) {
        Complex t;
        t.x = a[i].x * b[i].x - a[i].y * b[i].y;
        t.y = a[i].x * b[i].y + a[i].y * b[i].x;
        a[i].x = t.x * scale;
        a[i].y = t.y * scale;
    }
}

// 主函数
int main()
{
    int NX = 128;
    int NY = 128;
    int KX = 21;
    int KY = 21;
    // 1. 初始化输入信号和卷积核
    Complex *h_M = (Complex*)malloc(sizeof(Complex) * NX * NY);
    Complex *h_D = (Complex*)malloc(sizeof(Complex) * KX * KY);
    for (int i = 0; i < NX * NY; i++) {
        h_M[i].x = rand() / (float)RAND_MAX; h_M[i].y = 0;
    }
    for (int i = 0; i < KX * KY; i++) {
        h_D[i].x = rand() / (float)RAND_MAX; h_D[i].y = 0;
    }

    // 2. Pad to PX × PY
    Complex *h_pM, *h_pD;
    int PX, PY;
    int Nfft = PadData(h_M, &h_pM, NX, NY, h_D, &h_pD, KX, KY, PX, PY);
    size_t bytes = sizeof(Complex) * Nfft;

    // 3. Allocate and copy to device
    Complex *d_M, *d_D;
    cudaMalloc(&d_M, bytes);
    cudaMalloc(&d_D, bytes);
    cudaMemcpy(d_M, h_pM, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_D, h_pD, bytes, cudaMemcpyHostToDevice);

    // 4. 创建 FFT 计划
    cufftHandle plan;
    cufftPlan2d(&plan, PY, PX, CUFFT_C2C);

    // 5. Forward FFT
    cufftExecC2C(plan, d_M, d_M, CUFFT_FORWARD);
    cufftExecC2C(plan, d_D, d_D, CUFFT_FORWARD);

    // 6. 点乘 & 缩放
    float scale = 1.0f / Nfft;
    ComplexPointwiseMulAndScale<<<32, 256>>>(d_M, d_D, Nfft, scale);
    cudaDeviceSynchronize();

    // 7. IFFT 回空间域
    cufftExecC2C(plan, d_M, d_M, CUFFT_INVERSE);

    // 8. 拷回 host
    Complex *h_Hpad = (Complex*)malloc(bytes);
    cudaMemcpy(h_Hpad, d_M, bytes, cudaMemcpyDeviceToHost);

    // 9. 裁剪到 NX × NY 区域，并与 CPU 结果比较
    Complex *h_Href = (Complex*)malloc(sizeof(Complex) * NX * NY);
    Convolve(h_M, NX, NY, h_D, KX, KY, h_Href);

    int ox = KX / 2, oy = KY / 2;
    float error = 0;
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i_pad = (y + oy) * PX + (x + ox);
            int i_ref = y * NX + x;
            float dx = h_Hpad[i_pad].x - h_Href[i_ref].x;
            float dy = h_Hpad[i_pad].y - h_Href[i_ref].y;
            error += dx * dx + dy * dy;
        }
    }
    printf("L2 error: %.6f\n", sqrt(error));

    // 10. 清理资源
    cufftDestroy(plan);
    cudaFree(d_M); cudaFree(d_D);
    free(h_M); free(h_D); free(h_pM); free(h_pD);
    free(h_Hpad); free(h_Href);

    return 0;
}
