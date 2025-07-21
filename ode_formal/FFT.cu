// vector_fft_convolution.cu
// 单 GPU + CUFFT：对 3 分量 M 做 3×3 核卷积，得到 H_x,H_y,H_z

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <cufft.h>

// error check
#define CUDA_CHECK(err) \
  if((err)!=cudaSuccess){ \
    fprintf(stderr,"CUDA Error %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(err)); \
    exit(EXIT_FAILURE); }

#define CUFFT_CHECK(err) \
  if((err)!=CUFFT_SUCCESS){ \
    fprintf(stderr,"CUFFT Error %s:%d: %d\n",__FILE__,__LINE__,(int)err); \
    exit(EXIT_FAILURE); }

// 复数类型
typedef cufftComplex Complex;

// 频域：H = A * B
__global__ void complexMul(Complex *H, const Complex *A, const Complex *B, int N){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if(idx<N){
    float ar=A[idx].x, ai=A[idx].y;
    float br=B[idx].x, bi=B[idx].y;
    H[idx].x = ar*br - ai*bi;
    H[idx].y = ar*bi + ai*br;
  }
}
// 频域累加：H += A * B
__global__ void complexMulAdd(Complex *H, const Complex *A, const Complex *B, int N){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if(idx<N){
    float ar=A[idx].x, ai=A[idx].y;
    float br=B[idx].x, bi=B[idx].y;
    H[idx].x += (ar*br - ai*bi);
    H[idx].y += (ar*bi + ai*br);
  }
}
// IFFT 后做 1/(PX*PY) 归一化
__global__ void normalize(Complex *H, int N, float scale){
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if(idx<N){
    H[idx].x *= scale;
    H[idx].y *= scale;
  }
}

int main(){
  // —— 用户可修改这几行 —— 
  const int NX = 128, NY = 128;   // M 的原始尺寸
  const int KX =  64, KY =  64;   // D 的原始尺寸
  // ————————————————
  const int PX = NX + KX - 1;
  const int PY = NY + KY - 1;
  const int SIZE = PX * PY;

  // 1) 在 Host 上分配 M 的 3 分量 和 9 个核
  float *h_Mx = (float*)malloc(sizeof(float)*NX*NY);
  float *h_My = (float*)malloc(sizeof(float)*NX*NY);
  float *h_Mz = (float*)malloc(sizeof(float)*NX*NY);
  float *h_Dxx = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dxy = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dxz = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dyx = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dyy = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dyz = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dzx = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dzy = (float*)malloc(sizeof(float)*KX*KY);
  float *h_Dzz = (float*)malloc(sizeof(float)*KX*KY);

  // — TODO: 填充 h_M* 和 h_D** 数组，例如：
  for(int j=0;j<NY;j++){
    for(int i=0;i<NX;i++){
      int idx = i + j*NX;
      h_Mx[idx] = 1.0f;   // 你的 Mx
      h_My[idx] = 0.0f;   // 你的 My
      h_Mz[idx] = 0.0f;   // 你的 Mz
    }
  }
  for(int j=0;j<KY;j++){
    for(int i=0;i<KX;i++){
      int idx = i + j*KX;
      // 举例：只用 Dxx，其它设 0
      h_Dxx[idx] = 1.0f; 
      h_Dxy[idx] = 0.0f;
      h_Dxz[idx] = 0.0f;
      h_Dyx[idx] = 0.0f;
      h_Dyy[idx] = 1.0f; // 也可以组合不同核
      h_Dyz[idx] = 0.0f;
      h_Dzx[idx] = 0.0f;
      h_Dzy[idx] = 0.0f;
      h_Dzz[idx] = 1.0f;
    }
  }

  // 2) 扩零填充到 PX×PY，并转成 Complex
  auto pad = [&](const float *in, Complex *out, int wx, int wy){
    // 清 0
    for(int i=0;i<SIZE;i++) out[i].x = out[i].y = 0.0f;
    // 填实部
    for(int j=0;j<wy;j++){
      for(int i=0;i<wx;i++){
        out[i + j*PX].x = in[i + j*wx];
      }
    }
  };
  Complex *h_Mx_p = (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_My_p = (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Mz_p = (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dxx_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dxy_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dxz_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dyx_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dyy_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dyz_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dzx_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dzy_p= (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Dzz_p= (Complex*)malloc(sizeof(Complex)*SIZE);

  pad(h_Mx, h_Mx_p, NX, NY);
  pad(h_My, h_My_p, NX, NY);
  pad(h_Mz, h_Mz_p, NX, NY);
  pad(h_Dxx,h_Dxx_p,KX,KY);
  pad(h_Dxy,h_Dxy_p,KX,KY);
  pad(h_Dxz,h_Dxz_p,KX,KY);
  pad(h_Dyx,h_Dyx_p,KX,KY);
  pad(h_Dyy,h_Dyy_p,KX,KY);
  pad(h_Dyz,h_Dyz_p,KX,KY);
  pad(h_Dzx,h_Dzx_p,KX,KY);
  pad(h_Dzy,h_Dzy_p,KX,KY);
  pad(h_Dzz,h_Dzz_p,KX,KY);

  // 3) 在 GPU 上开空间
  Complex *d_Mx, *d_My, *d_Mz;
  Complex *d_Dxx,*d_Dxy,*d_Dxz;
  Complex *d_Dyx,*d_Dyy,*d_Dyz;
  Complex *d_Dzx,*d_Dzy,*d_Dzz;
  Complex *d_HxH,*d_HyH,*d_HzH; // 频域 H^
  CUDA_CHECK(cudaMalloc(&d_Mx, sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_My, sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Mz, sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dxx,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dxy,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dxz,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dyx,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dyy,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dyz,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dzx,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dzy,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_Dzz,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_HxH,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_HyH,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMalloc(&d_HzH,sizeof(Complex)*SIZE));

  // 4) 拷 M 和 D 到 Device
  CUDA_CHECK(cudaMemcpy(d_Mx, h_Mx_p, sizeof(Complex)*SIZE, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_My, h_My_p, sizeof(Complex)*SIZE, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Mz, h_Mz_p, sizeof(Complex)*SIZE, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dxx,h_Dxx_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dxy,h_Dxy_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dxz,h_Dxz_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dyx,h_Dyx_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dyy,h_Dyy_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dyz,h_Dyz_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dzx,h_Dzx_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dzy,h_Dzy_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_Dzz,h_Dzz_p,sizeof(Complex)*SIZE,cudaMemcpyHostToDevice));

  // 5) 创建 CUFFT plan (尺寸 PY×PX)
  cufftHandle plan;
  CUFFT_CHECK(cufftPlan2d(&plan, PY, PX, CUFFT_C2C));

  // 6) 正向 FFT: Mx,My,Mz, D**
  Complex* allD[9] = {d_Dxx,d_Dxy,d_Dxz,d_Dyx,d_Dyy,d_Dyz,d_Dzx,d_Dzy,d_Dzz};
  for(int i=0;i<9;i++)
    CUFFT_CHECK(cufftExecC2C(plan, allD[i], allD[i], CUFFT_FORWARD));
  CUFFT_CHECK(cufftExecC2C(plan, d_Mx, d_Mx, CUFFT_FORWARD));
  CUFFT_CHECK(cufftExecC2C(plan, d_My, d_My, CUFFT_FORWARD));
  CUFFT_CHECK(cufftExecC2C(plan, d_Mz, d_Mz, CUFFT_FORWARD));

  // 7) 清零频域 H^
  CUDA_CHECK(cudaMemset(d_HxH,0,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMemset(d_HyH,0,sizeof(Complex)*SIZE));
  CUDA_CHECK(cudaMemset(d_HzH,0,sizeof(Complex)*SIZE));

  // 8) 在频域累加 3×3 乘积
  int thr=256, blk=(SIZE+thr-1)/thr;
  // Hx_hat = Dxx*Mx + Dxy*My + Dxz*Mz
  complexMul <<<blk,thr>>>(d_HxH, d_Mx, d_Dxx, SIZE);
  complexMulAdd<<<blk,thr>>>(d_HxH, d_My, d_Dxy, SIZE);
  complexMulAdd<<<blk,thr>>>(d_HxH, d_Mz, d_Dxz, SIZE);
  // Hy_hat = Dyx*Mx + Dyy*My + Dyz*Mz
  complexMul <<<blk,thr>>>(d_HyH, d_Mx, d_Dyx, SIZE);
  complexMulAdd<<<blk,thr>>>(d_HyH, d_My, d_Dyy, SIZE);
  complexMulAdd<<<blk,thr>>>(d_HyH, d_Mz, d_Dyz, SIZE);
  // Hz_hat = Dzx*Mx + Dzy*My + Dzz*Mz
  complexMul <<<blk,thr>>>(d_HzH, d_Mx, d_Dzx, SIZE);
  complexMulAdd<<<blk,thr>>>(d_HzH, d_My, d_Dzy, SIZE);
  complexMulAdd<<<blk,thr>>>(d_HzH, d_Mz, d_Dzz, SIZE);

  // 9) 逆 FFT + 归一化
  CUFFT_CHECK(cufftExecC2C(plan, d_HxH, d_HxH, CUFFT_INVERSE));
  CUFFT_CHECK(cufftExecC2C(plan, d_HyH, d_HyH, CUFFT_INVERSE));
  CUFFT_CHECK(cufftExecC2C(plan, d_HzH, d_HzH, CUFFT_INVERSE));
  float scale = 1.0f/(float)SIZE;
  normalize<<<blk,thr>>>(d_HxH, SIZE, scale);
  normalize<<<blk,thr>>>(d_HyH, SIZE, scale);
  normalize<<<blk,thr>>>(d_HzH, SIZE, scale);
  cudaDeviceSynchronize();

  // 10) 拷回 Host
  Complex *h_Hx = (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Hy = (Complex*)malloc(sizeof(Complex)*SIZE);
  Complex *h_Hz = (Complex*)malloc(sizeof(Complex)*SIZE);
  CUDA_CHECK(cudaMemcpy(h_Hx, d_HxH, sizeof(Complex)*SIZE, cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(h_Hy, d_HyH, sizeof(Complex)*SIZE, cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(h_Hz, d_HzH, sizeof(Complex)*SIZE, cudaMemcpyDeviceToHost));

  // 11) 输出检查（取中心点）
  int cx = (PX/2)+(PY/2)*PX;
  printf("Hx(中心) = %f\n", h_Hx[cx].x);
  printf("Hy(中心) = %f\n", h_Hy[cx].x);
  printf("Hz(中心) = %f\n", h_Hz[cx].x);

  // 12) 清理
  cufftDestroy(plan);
  cudaFree(d_Mx); cudaFree(d_My); cudaFree(d_Mz);
  cudaFree(d_Dxx);cudaFree(d_Dxy);cudaFree(d_Dxz);
  cudaFree(d_Dyx);cudaFree(d_Dyy);cudaFree(d_Dyz);
  cudaFree(d_Dzx);cudaFree(d_Dzy);cudaFree(d_Dzz);
  cudaFree(d_HxH);cudaFree(d_HyH);cudaFree(d_HzH);
  free(h_Mx); free(h_My); free(h_Mz);
  free(h_Dxx);free(h_Dxy);free(h_Dxz);
  free(h_Dyx);free(h_Dyy);free(h_Dyz);
  free(h_Dzx);free(h_Dzy);free(h_Dzz);
  free(h_Mx_p);free(h_My_p);free(h_Mz_p);
  free(h_Dxx_p);free(h_Dxy_p);free(h_Dxz_p);
  free(h_Dyx_p);free(h_Dyy_p);free(h_Dyz_p);
  free(h_Dzx_p);free(h_Dzy_p);free(h_Dzz_p);
  free(h_Hx); free(h_Hy); free(h_Hz);

  return 0;
}
