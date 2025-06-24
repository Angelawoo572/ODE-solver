/** 
problem: three rate equations:
    dm1/dt = m3*f2 - m2*f3 + g1 - m*g*m1
    dm2/dt = m1*f3 - m3*f1 + g2 - m*g*m2
    dm3/dt = m2*f1 - m1*f2 + g3 - m*g*m3
on the interval from t = 0.0 to t = 4.e10, with initial
conditions: m1 = 1.0, m2 = 0.0, m3 = 0.0
This program solves the problem with the BDF method
*/

#include <cvode/cvode.h> /* prototypes for CVODE fcts., consts.           */
#include <nvector/nvector_cuda.h> /* access to cuda N_Vector                       */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h> /* defs. of sunrealtype, int                        */
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
/**
CUDA linear solver (if using cusolver)
access to cuSolverSp batch QR SUNLinearSolver */
#include <sunmatrix/sunmatrix_cusparse.h> /* access to cusparse SUNMatrix  */

__constant__ sunrealtype c_he;   // 前向/后向加权因子
__constant__ sunrealtype c_hk;   // 额外权重
__constant__ sunrealtype c_hap;  // 偏置
__constant__ sunrealtype c_ap;   // 放缩系数

/* Problem Constants */
#define GROUPSIZE 3               /* number of equations per group */
/* 我们每个 block 是 3×3，所以每组非零数 nnzper = 9 */
const int nnzper = GROUPSIZE * GROUPSIZE;
#define Y1        SUN_RCONST(1.0) /* initial y components */
#define Y2        SUN_RCONST(0.0)
#define Y3        SUN_RCONST(0.0)
#define RTOL      SUN_RCONST(1.0e-4) /* scalar relative tolerance            */
#define ATOL1     SUN_RCONST(1.0e-8) /* vector absolute tolerance components */
#define ATOL2     SUN_RCONST(1.0e-14)
#define ATOL3     SUN_RCONST(1.0e-6)
#define T0        SUN_RCONST(0.0)  /* initial time           */
#define T1        SUN_RCONST(0.4)  /* first output time      */
#define TMULT     SUN_RCONST(10.0) /* output time factor     */
#define NOUT      12               /* number of output times */

#define ZERO SUN_RCONST(0.0)

__constant__ float msk[3]=(0.0f,0.0f,1.0f);
__constant__ float chk=1.0f;
__constant__ float che=0.2f;
__constant__ float alpha=0.02;


/* Functions Called by the Solver */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to initialize the Jacobian sparsity pattern */
static int JacInit(SUNMatrix J);

/* Private function to output results */

static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3);

/* Private function to print final statistics */

static void PrintFinalStats(void* cvode_mem, SUNLinearSolver LS);


/* user data structure for parallel*/
typedef struct
{
  int ngroups; // number of groups
  int neq; // number of equations
  sunrealtype f1, f2, f3;
  sunrealtype g1, g2,g3;
  sunrealtype g,m;
} UserData;

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* Right hand side function evaluation kernel. */

__global__ static void f_kernel(const sunrealtype* __restrict__ y,
                                sunrealtype* __restrict__ yd,
                                int neq)
{
  sunrealtype gi, gj, gk;
  sunrealtype mg;
  sunindextype i, j, k, tid;

  tid = blockDim.x*blockIdx.x + threadIdx.x ;

  if ( tid > 2 && tid < blockDim - 3 )
  {
    iq=tid-3;
    ip=tid+3;
    ix=tid-(tid)%3;
    iy=ix+1;
    iz=iy+1;
    imsk=tid%3;
    h[tid] = che*(y[iq]+y[ip])+msk[imsk]*chk*y[iz];
  }
  __syncthreads()
  if ( tid > 2 && tid < blockDim - 3 )
  {
    i=tid-tid%3;
    j=i+1;
    k=j+1;
    mh[tid]=y[i]*h[i]+y[j]*h[j]+y[k]*h[k];

    j=tid+(tid+1)%3;
    k=tid+(tid+2)%3;
    yd[tid] = y[k]*h[j] - y[j]*h[k] + alpha*(h[tid] - mh[tid]*y[tid]);
  }
  else
  {
    yd[yid]=0;
  }
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

    unsigned block_size = 32;
    // total threads = grid_size * block_size
    // grid_size is ceil - (a+b-1)/b
    unsigned grid_size  = (udata->neq + block_size - 1) / block_size /3;
    f_kernel<<<grid_size, block_size>>>(ydata, ydotdata, udata->neq);

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
 * Jacobian initialization routine. This sets the sparisty pattern of
 * the blocks of the Jacobian J(t,y) = df/dy. This is performed on the CPU,
 * and only occurs at the beginning of the simulation.
 */

static int JacInit(SUNMatrix J)
{
  int    rowptrs[GROUPSIZE+1];
  int    colvals[nnzper];

  /* 全置零 */
  SUNMatZero(J);

  /* compressed sparse row 的 rowptrs */
  for (int i = 0; i <= GROUPSIZE; i++)
    rowptrs[i] = i * GROUPSIZE;

  /* 每行的列索引 0,1,2 */
  for (int i = 0; i < nnzper; i++)
    colvals[i] = i % GROUPSIZE;

  /* copy rowptrs, colvals to the device */
  SUNMatrix_cuSparse_CopyToDevice(J, NULL, rowptrs, colvals);
  cudaDeviceSynchronize();

  return (0);
}

/* Jacobian evaluation GPU kernel */
__global__ static void j_kernel(int ngroups,
                                sunrealtype f1, sunrealtype f2, sunrealtype f3,
                                sunrealtype g1, sunrealtype g2, sunrealtype g3,
                                sunrealtype m,  sunrealtype g,
                                sunrealtype* ydata,
                                sunrealtype* Jdata)
{
  int groupj;

  for (groupj = blockIdx.x * blockDim.x + threadIdx.x; groupj < ngroups;
       groupj += blockDim.x * gridDim.x)
  {

    /* first row of block: ∂f1/∂m1, ∂f1/∂m2, ∂f1/∂m3 */
    Jdata[nnzper * groupj + 0] = - m * g;
    Jdata[nnzper * groupj + 1] = - f3;
    Jdata[nnzper * groupj + 2] =   f2;

    /* second row of block: ∂f2/∂m1, ∂f2/∂m2, ∂f2/∂m3 */
    Jdata[nnzper * groupj + 3] =   f3;
    Jdata[nnzper * groupj + 4] = - m * g;
    Jdata[nnzper * groupj + 5] = - f1;

    /* third row of block: ∂f3/∂m1, ∂f3/∂m2, ∂f3/∂m3 */
    Jdata[nnzper * groupj + 6] = - f2;
    Jdata[nnzper * groupj + 7] =   f1;
    Jdata[nnzper * groupj + 8] = - m * g;
  }
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
 * This is done on the GPU.
 */

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata = (UserData*)user_data;
  sunrealtype *Jdata, *ydata;
  unsigned block_size, grid_size;

  Jdata  = SUNMatrix_cuSparse_Data(J);
  ydata  = N_VGetDeviceArrayPointer_Cuda(y);

  block_size = 32;
  grid_size  = (udata->neq + block_size - 1) / block_size;
  j_kernel<<<grid_size,block_size>>>(udata->ngroups,
                                     udata->f1, udata->f2, udata->f3,
                                     udata->g1, udata->g2, udata->g3,
                                     udata->m,  udata->g,
                                     ydata, Jdata);

  cudaDeviceSynchronize();
  cudaError_t cuerr = cudaGetLastError();
  if (cuerr != cudaSuccess)
  {
    fprintf(stderr, ">>> ERROR in Jac: cudaGetLastError returned %s\n",
            cudaGetErrorName(cuerr));
    return (-1);
  }

  return (0);
}


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main(int argc, char* argv[])
{
  SUNContext sunctx; // SUNDIALS context
  sunrealtype reltol, t, tout; // Solver tolerances and time variables
  sunrealtype *ydata, *abstol_data; // Host-side pointers to solution and tolerance data
  N_Vector y, abstol; // SUNDIALS vector structures for solution and absolute tolerance
  SUNMatrix A;
  SUNLinearSolver LS; // Linear solver object (cuSolverSp QR)
  void* cvode_mem; // CVODE integrator memory
  int retval, iout; // return status and output counter
  int neq, ngroups, groupj;// Problem size: number of equations, groups, and loop index
  UserData udata;
  cusparseHandle_t cusp_handle;
  cusolverSpHandle_t cusol_handle;

  sunrealtype h_he  = 1.0;
  sunrealtype h_hk  = 1.0;
  sunrealtype h_hap = 1.0;
  sunrealtype h_ap  = 1.0;

  cudaMemcpyToSymbol(c_he,  &h_he,  sizeof(sunrealtype));
  cudaMemcpyToSymbol(c_hk,  &h_hk,  sizeof(sunrealtype));
  cudaMemcpyToSymbol(c_hap, &h_hap, sizeof(sunrealtype));
  cudaMemcpyToSymbol(c_ap,  &h_ap,  sizeof(sunrealtype));

  y = abstol = NULL;// Initialize all pointers to NULL to ensure safe cleanup
  A = NULL;
  LS = NULL;  // Initialize linear solver pointer
  cvode_mem = NULL;  // Initialize CVODE memory


  /* Parse command line arguments */
  if (argc > 1) { ngroups = atoi(argv[1]); }
  else { ngroups = 100; }
  neq = ngroups * GROUPSIZE;

  udata.ngroups = ngroups;
  udata.neq     = neq;

  udata.f1 = 1.0;
  udata.f2 = 2.0;
  udata.f3 = 3.0;
  udata.g1 = 0.1;
  udata.g2 = 0.2;
  udata.g3 = 0.3;
  udata.g = 0.01;
  udata.m = 1.5;

  /* Initialize cuSOLVER and cuSPARSE handles */
  cusparseCreate(&cusp_handle);
  cusolverSpCreate(&cusol_handle);

  /* Create the SUNDIALS context */
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Create CUDA vector of length neq for I.C. and abstol */
  y = N_VNew_Cuda(neq, sunctx);
  abstol = N_VNew_Cuda(neq, sunctx);

  ydata       = N_VGetHostArrayPointer_Cuda(y);
  abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

  /* Initialize y */
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE)
  {
    ydata[groupj]     = Y1;
    ydata[groupj + 1] = Y2;
    ydata[groupj + 2] = Y3;
  }
  N_VCopyToDevice_Cuda(y);

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Set the vector absolute tolerance */
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE)
  {
    abstol_data[groupj]     = ATOL1;
    abstol_data[groupj + 1] = ATOL2;
    abstol_data[groupj + 2] = ATOL3;
  }
  N_VCopyToDevice_Cuda(abstol);

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  CVodeInit(cvode_mem, f, T0, y);

  /* Call CVodeSetUserData to attach the user data structure */
  CVodeSetUserData(cvode_mem, &udata);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  CVodeSVtolerances(cvode_mem, reltol, abstol);

  A = SUNMatrix_cuSparse_NewBlockCSR(ngroups, GROUPSIZE, GROUPSIZE,
                                     GROUPSIZE * GROUPSIZE, cusp_handle, sunctx);

  /* Set the sparsity pattern to be fixed so that the row pointers
   * and column indices are not zeroed out by SUNMatZero */
  SUNMatrix_cuSparse_SetFixedPattern(A, 1);
  /* Initialiize the Jacobian with its fixed sparsity pattern */
  JacInit(A);
  /* Create the SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_cuSolverSp_batchQR(y, A, cusol_handle, sunctx);

  CVodeSetLinearSolver(cvode_mem, LS, A);

  /* Set the user-supplied Jacobian routine Jac */
  CVodeSetJacFn(cvode_mem, Jac);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  printf(" \nGroup of independent 3-species kinetics problems\n\n");
  printf("number of groups = %d\n\n", ngroups);

  iout = 0;
  tout = T1;
  while (1)
  {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    N_VCopyFromDevice_Cuda(y);
    for (groupj = 0; groupj < ngroups; groupj += 10)
    {
      printf("group %d: ", groupj);
      PrintOutput(t, ydata[GROUPSIZE * groupj], ydata[1 + GROUPSIZE * groupj],
                  ydata[2 + GROUPSIZE * groupj]);
    }
    if (retval == CV_SUCCESS)
    {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) { break; }
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem, LS);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(A);

  SUNContext_Free(&sunctx);

  /* Destroy the cuSOLVER and cuSPARSE handles */
  cusparseDestroy(cusp_handle);
  cusolverSpDestroy(cusol_handle);

  return (0);
}
