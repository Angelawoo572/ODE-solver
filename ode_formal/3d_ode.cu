/**
problem: three rate equations:
   dm1/dt = m3*f2 - m2*f3 + g1 - m*g*m1
   dm2/dt = m1*f3 - m3*f1 + g2 - m*g*m2
   dm3/dt = m2*f1 - m1*f2 + g3 - m*g*m3
on the interval from t = 0.0 to t = 4.e10, with
This program solves the problem with the BDF method
*/
#include <cvode/cvode.h> /* prototypes for CVODE fcts., consts.           */
#include <nvector/nvector_cuda.h> /* access to cuda N_Vector                       */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h> /* defs. of sunrealtype, int                        */
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <math.h>
/* Problem Constants */
#define GROUPSIZE 3               /* number of equations per group */
#define indexbound 2
#define ONE 1
#define TWO 2
#define RTOL      SUN_RCONST(1.0e-5) /* scalar relative tolerance            */
#define ATOL1     SUN_RCONST(1.0e-5) /* vector absolute tolerance components */
#define ATOL2     SUN_RCONST(1.0e-5)
#define ATOL3     SUN_RCONST(1.0e-5)
#define T0        SUN_RCONST(0.0)  /* initial time           */
#define T1        SUN_RCONST(0.1)  /* first output time      */
#define DT    ((T1 - T0) / NOUT)
// #define NOUT      120             /* number of output times */
#define ZERO SUN_RCONST(0.0)
#define BX 16
#define BY 16
// constant memory
__constant__ float msk[3]={0.0f,0.0f,1.0f};
__constant__ float nsk[3]={1.0f,0.0f,0.0f};
__constant__ float chk=1.0f;
__constant__ float che =0.0f;
__constant__ float alpha=0.02f; // 0.0f
__constant__ float chg = 1.0f;
__constant__ float cha = 1.5f; //0.2
__constant__ float chb = 0.0f;
/* user data structure for parallel*/
typedef struct
{
  // int ngroups; // number of groups
  int nx, ny;
  int neq; // number of equations
  // SOA
  sunrealtype *d_y0, *d_y1, *d_y2;   // m_x, m_y, m_z
  sunrealtype *d_yd0, *d_yd1, *d_yd2; // dm_x/dτ, dm_y/dτ, dm_z/dτ
} UserData;


/*
*-------------------------------
* Functions called by the solver
*-------------------------------
*/


/* Right hand side function evaluation kernel. */
__global__ static void f_kernel(
  const sunrealtype * __restrict__ y0,
  const sunrealtype * __restrict__ y1,
  const sunrealtype * __restrict__ y2,
  sunrealtype * __restrict__ yd0,
  sunrealtype * __restrict__ yd1,
  sunrealtype * __restrict__ yd2,
  int nx, int ny)
{
  // linear group index and base pointer for 3 components
  int ix = blockIdx.x * BX + threadIdx.x;
  int iy = blockIdx.y * BY + threadIdx.y;
  if (ix >= nx || iy >= ny) return;
  int gid = iy * nx + ix;
  // compute how many blocks in every row
  int blocks_x = (nx + BX - 1) / BX; // blockDim.s
  // blockIdx.x decides phase 0 = red, 1 = blue
  int phase = blockIdx.x / blocks_x;
  // checkerboard partition: red if (i+j) % 2 == 0
  bool is_red = (((ix+iy)&1) == 0);
  bool is_red_phase = (phase == 0);
  if (is_red != is_red_phase) return;

  __shared__ sunrealtype sx[BY+2][BX+2];
  __shared__ sunrealtype sy[BY+2][BX+2];
  __shared__ sunrealtype sz[BY+2][BX+2];

  int tx = threadIdx.x + 1;
  int ty = threadIdx.y + 1;
  sx[ty][tx] = y0[gid];
  sy[ty][tx] = y1[gid];
  sz[ty][tx] = y2[gid];

  if (threadIdx.x == 0 && ix > 0) {
    sx[ty][0] = y0[iy*nx + (ix-1)];
    sy[ty][0] = y1[iy*nx + (ix-1)];
    sz[ty][0] = y2[iy*nx + (ix-1)];
  }
  if (threadIdx.x == BX-1 && ix < nx-1) {
    sx[ty][BX+1] = y0[iy*nx + (ix+1)];
    sy[ty][BX+1] = y1[iy*nx + (ix+1)];
    sz[ty][BX+1] = y2[iy*nx + (ix+1)];
  }
  if (threadIdx.y == 0 && iy > 0) {
    sx[0][tx] = y0[(iy-1)*nx + ix];
    sy[0][tx] = y1[(iy-1)*nx + ix];
    sz[0][tx] = y2[(iy-1)*nx + ix];
  }
  if (threadIdx.y == BY-1 && iy < ny-1) {
    sx[BY+1][tx] = y0[(iy+1)*nx + ix];
    sy[BY+1][tx] = y1[(iy+1)*nx + ix];
    sz[BY+1][tx] = y2[(iy+1)*nx + ix];
  }
  __syncthreads();

  float mx = sx[ty][tx], my = sy[ty][tx], mz = sz[ty][tx];
  float mxl = sx[ty][tx-1], mxr = sx[ty][tx+1], mxu = sx[ty-1][tx], mxd = sx[ty+1][tx];
  float myl = sy[ty][tx-1], myr = sy[ty][tx+1], myu = sy[ty-1][tx], myd = sy[ty+1][tx];
  float mzl = sz[ty][tx-1], mzr = sz[ty][tx+1], mzu = sz[ty-1][tx], mzd = sz[ty+1][tx];

  // 1) 计算局部场 h
  float hx = che * (mxl + mxr + mxu + mxd)
             + msk[0] * (chk*mz + cha)
             + nsk[0] * (mxl + mxr) * chb;
  float hy = che * (myl + myr + myu + myd)
             + msk[1] * (chk*mz + cha)
             + nsk[1] * (myl + myr) * chb;
  float hz = che * (mzl + mzr + mzu + mzd)
             + msk[2] * (chk*mz + cha)
             + nsk[2] * (mzl + mzr) * chb;

  // 2) 计算 m*h
  float dot = mx*hx + my*hy + mz*hz;

  // 3) 计算 dm/dτ
  yd0[gid] = chg * (mz*hy - my*hz) + alpha * (hx - dot*mx);
  yd1[gid] = chg * (mx*hz - mz*hx) + alpha * (hy - dot*my);
  yd2[gid] = chg * (my*hx - mx*hy) + alpha * (hz - dot*mz);
}


/* Right hand side function. This just launches the CUDA kernel
  to do the actual computation. At the very least, doing this
  saves moving the vector data in y and ydot to/from the device
  every evaluation of f. */


static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData *udata = (UserData*)user_data;
  int nx = udata->nx, ny = udata->ny;
  size_t ng = (size_t)nx * ny;        // 每个分量的长度

  // 1) 拿到 y 和 ydot 在设备上的“基址”
  sunrealtype *y_base    = N_VGetDeviceArrayPointer_Cuda(y);
  sunrealtype *ydot_base = N_VGetDeviceArrayPointer_Cuda(ydot);

  // 2) 切分成 SoA
  sunrealtype *y0  = y_base;            // m_x
  sunrealtype *y1  = y_base +   ng;      // m_y
  sunrealtype *y2  = y_base + 2*ng;      // m_z

  sunrealtype *yd0 = ydot_base;         // dm_x/dτ
  sunrealtype *yd1 = ydot_base +   ng;   // dm_y/dτ
  sunrealtype *yd2 = ydot_base + 2*ng;   // dm_z/dτ

  // 3) 启动你的 kernel
  dim3 block(16,16);
  int blocks_x = (nx + block.x - 1) / block.x;
  int blocks_y = (ny + block.y - 1) / block.y;
  dim3 grid(2*blocks_x, blocks_y);

  f_kernel<<<grid, block>>>(
    /* input  */ y0,  y1,  y2,
    /* output */ yd0, yd1, yd2,
    /* dims   */ nx,  ny);
  cudaDeviceSynchronize();
  
  cudaError_t cuerr = cudaGetLastError();
  if (cuerr != cudaSuccess)
  {
       fprintf(stderr, ">>> ERROR in f: cudaGetLastError returned %s\n",
               cudaGetErrorName(cuerr));
       return (-1);
  }
  return (0);
}


/*
*-------------------------------
* Private helper functions
*-------------------------------
*/
static void PrintOutput(sunrealtype t, int nx, int ny,
                        const sunrealtype* y0,
                        const sunrealtype* y1,
                        const sunrealtype* y2)
{
    printf("\n==== At time t = %0.4e ====\n", t);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = j * nx + i;
#if defined(SUNDIALS_EXTENDED_PRECISION)
            printf("(%d,%d): m = %14.6Le  %14.6Le  %14.6Le\n",
                   i, j, y0[idx], y1[idx], y2[idx]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
            printf("(%d,%d): m = %14.6e  %14.6e  %14.6e\n",
                   i, j, y0[idx], y1[idx], y2[idx]);
#else
            printf("(%d,%d): m = %14.6e  %14.6e  %14.6e\n",
                   i, j, y0[idx], y1[idx], y2[idx]);
#endif
        }
    }
}


/*
* Get and print some final statistics
*/
static void PrintFinalStats(void* cvode_mem, SUNLinearSolver LS)
{
   long int nst, nfe, nsetups, nni, ncfn, netf, nge;


   CVodeGetNumSteps(cvode_mem, &nst);
   CVodeGetNumRhsEvals(cvode_mem, &nfe);
   CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
   CVodeGetNumErrTestFails(cvode_mem, &netf);
   CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
   CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
   CVodeGetNumGEvals(cvode_mem, &nge);


   printf("\nFinal Statistics:\n");
   printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld", nst, nfe,
         nsetups);
   printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n", nni, ncfn,
         netf, nge);
}


/*
*-------------------------------
* Main Program
*-------------------------------
*/
int main(int argc, char* argv[])
{
  SUNContext sunctx; // SUNDIALS context
  sunrealtype *abstol_data; // Host-side pointers to solution and tolerance data
  sunrealtype t;
  sunrealtype tout;
  N_Vector y, abstol; // SUNDIALS vector structures for solution and absolute tolerance
  SUNLinearSolver LS; // Linear solver object (cuSolverSp QR)
  SUNNonlinearSolver NLS;
  void* cvode_mem; // CVODE integrator memory
  int retval, iout; // return status and output counter
  int neq;// Problem size: number of equations, groups, and loop index
  UserData udata;
  

  /* Parse command-line to get number of groups */
  int nx = 128, ny = 128;
  neq     = nx * ny * 3;
  size_t ng = nx * ny;
  // 1) Allocate and initialize host‑side SoA arrays
  sunrealtype *h_y0 = (sunrealtype*)malloc(nx*ny*sizeof(sunrealtype));
  sunrealtype *h_y1 = (sunrealtype*)malloc(nx*ny*sizeof(sunrealtype));
  sunrealtype *h_y2 = (sunrealtype*)malloc(nx*ny*sizeof(sunrealtype));

  /* Fill user data */
  udata.nx  = nx;
  udata.ny  = ny;
  udata.neq = neq;
  
  /* Create SUNDIALS context */
  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  /* Allocate CUDA vectors for solution and tolerances */
  y     = N_VNew_Cuda(neq, sunctx);
  abstol= N_VNew_Cuda(neq, sunctx);
  abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

  sunrealtype *y_host = N_VGetHostArrayPointer_Cuda(y);
  for (size_t i = 0; i < ng; i++) {
    y_host[i        ] = h_y0[i];
    y_host[i +   ng ] = h_y1[i];
    y_host[i + 2*ng ] = h_y2[i];
  }
  // 这条才是真正把初始条件拷到 device
  N_VCopyToDevice_Cuda(y);


   /* Initialize y and abstol on host then copy to device */
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int idx = j*nx + i;
      if (j < ny/2 && i < nx / 2) {
        h_y0[idx] = 0.0;
        h_y1[idx] = 0.0175;
        h_y2[idx] = 0.998;
      } else {
        h_y0[idx] = 0.0;
        h_y1[idx] = 0.0175;
        h_y2[idx] = -0.998;
      }
      abstol_data[idx + 0] = ATOL1;
      abstol_data[idx + 1] = ATOL2;
      abstol_data[idx + 2] = ATOL3;
    }
  }
  N_VCopyToDevice_Cuda(y);
  N_VCopyToDevice_Cuda(abstol);

  // copy initial host data
  cudaMemcpy(udata.d_y0, h_y0, nx*ny*sizeof(sunrealtype), cudaMemcpyHostToDevice);
  cudaMemcpy(udata.d_y1, h_y1, nx*ny*sizeof(sunrealtype), cudaMemcpyHostToDevice);
  cudaMemcpy(udata.d_y2, h_y2, nx*ny*sizeof(sunrealtype), cudaMemcpyHostToDevice);
  /* Create and initialize CVODE solver memory */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  CVodeInit(cvode_mem, f, T0, y);
  CVodeSetUserData(cvode_mem, &udata);
  CVodeSVtolerances(cvode_mem, RTOL, abstol);

  /* Matrix-free GMRES linear solver (no Jacobian needed) */
  NLS = SUNNonlinSol_Newton(y, sunctx);
  CVodeSetNonlinearSolver(cvode_mem, NLS);
  LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
  CVodeSetLinearSolver(cvode_mem, LS, NULL);

  /* Print header */
  printf("\nGroup of independent 3-species kinetics problems\n\n");
  printf("number of groups = %d %d %d \n", nx, ny, nx * ny);


  /* Time-stepping loop */
  float ttotal=500.0f;
  iout = T0;
  tout = T1;
  int NOUT=ttotal/T1;
  while (iout < NOUT) {
    // &t cvode实际走到的地方
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    // copy solution back to host and print all groups
    if (retval == CV_SUCCESS) {
      iout++;
      tout += T1; // T0 + iout*T1
    }else {
      fprintf(stderr, "CVode error at output %d: retval = %d\n", iout, retval);
      break;
    }
    N_VCopyFromDevice_Cuda(y);
    sunrealtype *y_host2 = N_VGetHostArrayPointer_Cuda(y);
    // y_host2[0..ng-1] 是 m1, [ng..2ng-1] 是 m2, [2ng..3ng-1] 是 m3
    PrintOutput(t, nx, ny,
                y_host2,
                y_host2 +   ng,
                y_host2 + 2*ng);
  }
  PrintOutput(t, nx, ny, h_y0, h_y1, h_y2);

  /* Print final statistics */
  PrintFinalStats(cvode_mem, LS);


  /* Clean up */
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNNonlinSolFree(NLS);
  N_VDestroy(y);
  N_VDestroy(abstol);
  SUNContext_Free(&sunctx);

  cudaFree(udata.d_y0);
  cudaFree(udata.d_y1);
  cudaFree(udata.d_y2);
  cudaFree(udata.d_yd0);
  cudaFree(udata.d_yd1);
  cudaFree(udata.d_yd2);

  return 0;
}
