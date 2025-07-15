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

// constant memory

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
    sunrealtype *d_h;
    sunrealtype *d_mh;
} UserData;

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* Right hand side function evaluation kernel. */
__global__ static void f_kernel(
  const sunrealtype* y, 
  sunrealtype* yd, 
  sunrealtype* h,
  sunrealtype* mh,
  int nx, 
  int ny,
  bool is_red_phase)
{
    // compute 2D thread coordinates
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    if (ix >= nx || iy >= ny) return;

    // checkerboard partition: red if (i+j) % 2 == 0
    bool is_red = ((ix+iy)&1) == 0;
    if (is_red != is_red_phase) return;

    // linear group index and base pointer for 3 components
    int gid = iy * nx + ix;
    int base_idx = GROUPSIZE * gid;

    // neighbor group indices
    int ix_l = (ix > 0) ? ix - 1 : ix;
    int ix_r = (ix < nx-1) ? ix + 1 : ix;
    int iy_u = (iy > 0) ? iy - 1 : iy;
    int iy_d = (iy < ny-1) ? iy + 1 : iy;

    int base_l = GROUPSIZE * (iy * nx + ix_l);
    int base_r = GROUPSIZE * (iy * nx + ix_r);
    int base_u = GROUPSIZE * (iy_u * nx + ix);
    int base_d = GROUPSIZE * (iy_d * nx + ix);

    // compute h vector for each component
    for (int c = 0; c < GROUPSIZE; ++c) {
        h[base_idx + c] =
            che * (y[base_l + c] + y[base_r + c] + y[base_u + c] + y[base_d + c])
          + msk[c] * (chk * y[base_idx + 2] + cha)
          + nsk[c] * (y[base_r + c] + y[base_l + c]) * chb;
    }
    __syncthreads();
    // Dot product m*h for this group
    sunrealtype dot = 
        y[base_idx + 0] * h[base_idx + 0] +
        y[base_idx + 1] * h[base_idx + 1] +
        y[base_idx + 2] * h[base_idx + 2];
    mh[base_idx + 0] = dot;
    mh[base_idx + 1] = dot;
    mh[base_idx + 2] = dot;
    __syncthreads();
    yd[base_idx + 0] = chg * (y[base_idx + 2]*h[base_idx + 1] - y[base_idx + 1]*h[base_idx + 2])
                      + alpha * (h[base_idx + 0] - dot * y[base_idx + 0]);
    yd[base_idx + 1] = chg * (y[base_idx + 0]*h[base_idx + 2]
                             - y[base_idx + 2]*h[base_idx + 0])
                     + alpha * (h[base_idx + 1] - dot * y[base_idx + 1]);
    yd[base_idx + 2] = chg * (y[base_idx + 1]*h[base_idx + 0]
                             - y[base_idx + 0]*h[base_idx + 1])
                     + alpha * (h[base_idx + 2] - dot * y[base_idx + 2]);
}

/* Right hand side function. This just launches the CUDA kernel
   to do the actual computation. At the very least, doing this
   saves moving the vector data in y and ydot to/from the device
   every evaluation of f. */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    UserData* udata;
    sunrealtype *ydata, *ydotdata;

    udata    = (UserData*)user_data;
    ydata    = N_VGetDeviceArrayPointer_Cuda(y);
    ydotdata = N_VGetDeviceArrayPointer_Cuda(ydot);

    int nx = udata->nx, ny = udata->ny;
    dim3 block(16, 16);
    dim3 grid((nx + block.x - 1)/block.x, (ny + block.y - 1)/block.y);

    // Red phase
    f_kernel<<<grid, block>>>(ydata, ydotdata,
                                 udata->d_h, udata->d_mh,
                                 nx, ny, true);
    cudaDeviceSynchronize();

    // Blue phase
    f_kernel<<<grid, block>>>(ydata, ydotdata,
                                 udata->d_h, udata->d_mh,
                                 nx, ny, false);
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
static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3)
{
    #if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
    #elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
    #else
      printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
    #endif

  return;
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
    sunrealtype *ydata, *abstol_data; // Host-side pointers to solution and tolerance data
    sunrealtype t;
    sunrealtype tout;
    N_Vector y, abstol; // SUNDIALS vector structures for solution and absolute tolerance
    SUNLinearSolver LS; // Linear solver object (cuSolverSp QR)
    SUNNonlinearSolver NLS;
    void* cvode_mem; // CVODE integrator memory
    int retval, iout; // return status and output counter
    int neq, ngroups, groupj;// Problem size: number of equations, groups, and loop index
    UserData udata;

    /* Parse command-line to get number of groups */
    int nx = 128, ny = 128; 
    neq     = nx * ny * 3;

    /* Fill user data */
    udata.nx  = nx;
    udata.ny  = ny;
    udata.neq     = neq;
    cudaMalloc(&udata.d_h,  neq * sizeof(sunrealtype));
    cudaMalloc(&udata.d_mh, neq * sizeof(sunrealtype));

    /* Create SUNDIALS context */
    SUNContext_Create(SUN_COMM_NULL, &sunctx);

    /* Allocate CUDA vectors for solution and tolerances */
    y     = N_VNew_Cuda(neq, sunctx);
    abstol= N_VNew_Cuda(neq, sunctx);
    // get host pointers
    ydata       = N_VGetHostArrayPointer_Cuda(y);
    abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

    /* Initialize y and abstol on host then copy to device */
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int idx = 3 * (j * nx + i);
        if (j < ny / 2) {
          ydata[idx + 0] = 0.0;
          ydata[idx + 1] = 0.0175;
          ydata[idx + 2] = 0.998;
        } else {
          ydata[idx + 0] = 0.0;
          ydata[idx + 1] = 0.0175;
          ydata[idx + 2] = -0.998;
        }

        abstol_data[idx + 0] = ATOL1;
        abstol_data[idx + 1] = ATOL2;
        abstol_data[idx + 2] = ATOL3;
      }
    }
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(abstol);

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
    printf("number of groups = %d\n\n", nx, ny, nx * ny);

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
      // printf("%f\n",tout);
      N_VCopyFromDevice_Cuda(y);
      ydata = N_VGetHostArrayPointer_Cuda(y);
      int mid = (ny/2) * nx + (nx/2);
      printf("t = %.2f, m = [%g %g %g]\n", t,
              ydata[3*mid+0], ydata[3*mid+1], ydata[3*mid+2]);
    }

    /* Print final statistics */
    PrintFinalStats(cvode_mem, LS);

    /* Clean up */
    cudaFree(udata.d_h);
    cudaFree(udata.d_mh);
    N_VDestroy(y);
    N_VDestroy(abstol);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNNonlinSolFree(NLS);
    SUNContext_Free(&sunctx);

    return 0;
}
