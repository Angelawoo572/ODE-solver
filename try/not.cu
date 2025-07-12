/**
 * Pure CUDA solver using fixed-step RK4 integration
 * for systems of the form:
 *   dm1/dt = m3*f2 - m2*f3 + g1 - m*g*m1
 *   dm2/dt = m1*f3 - m3*f1 + g2 - m*g*m2
 *   dm3/dt = m2*f1 - m1*f2 + g3 - m*g*m3
 *
 * This example runs 32 groups (3 vars each), total 96 equations.
 */

#include <stdio.h>
#include <cuda_runtime.h>

#define GROUPSIZE 3
#define NGROUPS 32
#define NEQ (GROUPSIZE * NGROUPS)
#define DT 1e-2f
#define T0 0.0f
#define TEND 0.4f
#define NSTEPS ((int)((TEND - T0) / DT))

__constant__ float che = 0.2f;
__constant__ float chk = 1.0f;
__constant__ float alpha = 0.02f;
__constant__ float msk[3] = {0.0f, 0.0f, 1.0f};

// Kernel: compute dy/dt
__global__ void f_kernel(const float* y, float* dydt, float* h, float* mh, int neq) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid > 2 && tid < neq - 3) {
        int iq = tid - 3;
        int ip = tid + 3;
        int ix = tid - tid % 3;
        int iy = ix + 1;
        int iz = iy + 1;
        int imsk = tid % 3;

        h[tid] = che * (y[iq] + y[ip]) + msk[imsk] * chk * y[iz];
    }
    __syncthreads();

    if (tid > 2 && tid < neq - 3) {
        int i = tid - tid % 3;
        int j = i + 1;
        int k = j + 1;
        mh[tid] = y[i]*h[i] + y[j]*h[j] + y[k]*h[k];

        j = tid + (tid + 1) % 3;
        k = tid + (tid + 2) % 3;

        dydt[tid] = y[k]*h[j] - y[j]*h[k] + alpha * (h[tid] - mh[tid] * y[tid]);
    } else if (tid < neq) {
        dydt[tid] = 0.0f;
    }
}

__global__ void copy_and_axpy(float* dst, const float* src, const float* dx, float a, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        dst[i] = src[i] + a * dx[i];
}

__global__ void rk4_update(float* y, const float* k1, const float* k2,
                           const float* k3, const float* k4, float dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < NEQ)
        y[i] += dt / 6.0f * (k1[i] + 2.0f*k2[i] + 2.0f*k3[i] + k4[i]);
}

int main() {
    float *y, *tmp, *k1, *k2, *k3, *k4, *h, *mh;
    cudaMalloc(&y,   NEQ * sizeof(float));
    cudaMalloc(&tmp, NEQ * sizeof(float));
    cudaMalloc(&k1,  NEQ * sizeof(float));
    cudaMalloc(&k2,  NEQ * sizeof(float));
    cudaMalloc(&k3,  NEQ * sizeof(float));
    cudaMalloc(&k4,  NEQ * sizeof(float));
    cudaMalloc(&h,   NEQ * sizeof(float));
    cudaMalloc(&mh,  NEQ * sizeof(float));

    float h_y[NEQ];
    for (int i = 0; i < NGROUPS; i++) {
        int ix = 3 * i;
        int iy = ix + 1;
        int iz = iy + 1;
        if (i == 0) {
            h_y[ix] = 0; h_y[iy] = 0; h_y[iz] = 1;
        } else if (i == NGROUPS - 1) {
            h_y[ix] = 0; h_y[iy] = 0; h_y[iz] = -1;
        } else if (i < NGROUPS / 2) {
            h_y[ix] = 0; h_y[iy] = 0.0175; h_y[iz] = 0.998;
        } else {
            h_y[ix] = 0; h_y[iy] = 0.0175; h_y[iz] = -0.998;
        }
    }
    cudaMemcpy(y, h_y, NEQ * sizeof(float), cudaMemcpyHostToDevice);

    dim3 blockDim(NEQ);
    dim3 gridDim((NEQ + blockDim.x - 1) / blockDim.x);

    for (int step = 0; step < NSTEPS; ++step) {
        f_kernel<<<gridDim, blockDim>>>(y, k1, h, mh, NEQ);
        cudaDeviceSynchronize();

        copy_and_axpy<<<gridDim, blockDim>>>(tmp, y, k1, 0.5f * DT, NEQ);
        f_kernel<<<gridDim, blockDim>>>(tmp, k2, h, mh, NEQ);
        cudaDeviceSynchronize();

        copy_and_axpy<<<gridDim, blockDim>>>(tmp, y, k2, 0.5f * DT, NEQ);
        f_kernel<<<gridDim, blockDim>>>(tmp, k3, h, mh, NEQ);
        cudaDeviceSynchronize();

        copy_and_axpy<<<gridDim, blockDim>>>(tmp, y, k3, DT, NEQ);
        f_kernel<<<gridDim, blockDim>>>(tmp, k4, h, mh, NEQ);
        cudaDeviceSynchronize();

        rk4_update<<<gridDim, blockDim>>>(y, k1, k2, k3, k4, DT);
        cudaDeviceSynchronize();

        if (step % 10 == 0) {
            cudaMemcpy(h_y, y, NEQ * sizeof(float), cudaMemcpyDeviceToHost);
            printf("t = %.3f | group 0 = [%g, %g, %g]\n", T0 + step * DT, h_y[0], h_y[1], h_y[2]);
        }
    }

    cudaFree(y); cudaFree(tmp); cudaFree(k1); cudaFree(k2);
    cudaFree(k3); cudaFree(k4); cudaFree(h); cudaFree(mh);
    return 0;
}

