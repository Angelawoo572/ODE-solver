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
#include <sunlinsol/sunlinsol_spgmr.h>

// constant memory
__constant__ sunrealtype c_he;
__constant__ sunrealtype c_hk;
__constant__ sunrealtype c_hap;
__constant__ sunrealtype c_ap;

/* Problem Constants */
#define GROUPSIZE 3               /* number of equations per group */
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

/* Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);    

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

__global__ static void f_kernel(const sunrealtype* y, sunrealtype* yd, int neq)
{
  sunrealtype fi, fj, fk;
  sunrealtype gi, gj, gk;
  sunrealtype mg;
  sunindextype i, j, k, tid;

  tid = blockDim.x*blockIdx.x + threadIdx.x ;

  i = tid % blockDim.x;
  j = blockDim.x + tid % blockDim.x;
  k = 2 * blockDim.x + tid % blockDim.x;
  /* only valid lanes do work; others write zero */
  if (tid>0 && tid<neq-1) {
    fi = c_he*(y[i+1] + y[i-1]);
    fj = c_he*(y[j+1] + y[j-1]);
    fk = c_he*(y[k+1] + y[k-1]) + c_hk*y[k] + c_hap;
    gi = c_ap*fi; gj = c_ap*fj; gk = c_ap*fk;
    mg = y[i]*gi + y[j]*gj + y[k]*gk;
    /* which equation depends on tid region */
    if (tid < blockDim.x) {
      yd[tid] = y[k]*fj - y[j]*fk + gi - mg*y[i];
    } else if (tid < 2*blockDim.x) {
      yd[tid] = y[i]*fk - y[k]*fi + gj - mg*y[j];
    } else {
      yd[tid] = y[j]*fi - y[i]*fj + gk - mg*y[k];
    }
  } else {
    yd[tid] = 0;
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
    unsigned grid_size  = (udata->neq + block_size - 1) / block_size;
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
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main(int argc, char* argv[])
{
    SUNContext      sunctx = NULL;
    sunrealtype     reltol, t, tout;
    sunrealtype    *ydata = NULL, *abstol_data = NULL;
    N_Vector        y = NULL, abstol = NULL;
    SUNLinearSolver LS = NULL;
    void           *cvode_mem = NULL;
    int             retval, iout;
    int             neq, ngroups, groupj;
    UserData        udata;

    /* Copy constants into device constant memory */
    sunrealtype h_he  = 1.0, h_hk = 1.0, h_hap = 1.0, h_ap = 1.0;
    cudaMemcpyToSymbol(c_he,  &h_he,  sizeof(sunrealtype));
    cudaMemcpyToSymbol(c_hk,  &h_hk,  sizeof(sunrealtype));
    cudaMemcpyToSymbol(c_hap, &h_hap, sizeof(sunrealtype));
    cudaMemcpyToSymbol(c_ap,  &h_ap,  sizeof(sunrealtype));

    /* Parse command-line to get number of groups */
    ngroups = (argc > 1 ? atoi(argv[1]) : 100);
    neq     = ngroups * GROUPSIZE;

    /* Fill user data */
    udata.ngroups = ngroups;
    udata.neq     = neq;
    udata.f1 = 1.0; udata.f2 = 2.0; udata.f3 = 3.0;
    udata.g1 = 0.1; udata.g2 = 0.2; udata.g3 = 0.3;
    udata.g  = 0.01; udata.m   = 1.5;

    /* Create SUNDIALS context */
    SUNContext_Create(SUN_COMM_NULL, &sunctx);

    /* Allocate CUDA vectors for solution and tolerances */
    y     = N_VNew_Cuda(neq, sunctx);
    abstol= N_VNew_Cuda(neq, sunctx);
    ydata       = N_VGetHostArrayPointer_Cuda(y);
    abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

    /* Initialize y and abstol on host then copy to device */
    for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
        ydata[groupj]     = Y1;
        ydata[groupj + 1] = Y2;
        ydata[groupj + 2] = Y3;
        abstol_data[groupj]     = ATOL1;
        abstol_data[groupj + 1] = ATOL2;
        abstol_data[groupj + 2] = ATOL3;
    }
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(abstol);

    /* Create and initialize CVODE solver memory */
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    CVodeInit(cvode_mem, f, T0, y);
    CVodeSetUserData(cvode_mem, &udata);
    CVodeSVtolerances(cvode_mem, RTOL, abstol);

    /* Matrix-free GMRES linear solver (no Jacobian needed) */
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
    CVodeSetLinearSolver(cvode_mem, LS, NULL);

    /* Print header */
    printf("\nGroup of independent 3-species kinetics problems\n\n");
    printf("number of groups = %d\n\n", ngroups);

    /* Time-stepping loop */
    iout = 0;
    tout = T1;
    while (iout < NOUT) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        N_VCopyFromDevice_Cuda(y);
        for (groupj = 0; groupj < ngroups; groupj += 10) {
            printf("group %d: ", groupj);
            PrintOutput(t,
                        ydata[GROUPSIZE * groupj],
                        ydata[1 + GROUPSIZE * groupj],
                        ydata[2 + GROUPSIZE * groupj]);
        }
        if (retval == CV_SUCCESS) {
            iout++;
            tout *= TMULT;
        }
    }

    /* Print final statistics */
    PrintFinalStats(cvode_mem, LS);

    /* Clean up */
    N_VDestroy(y);
    N_VDestroy(abstol);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNContext_Free(&sunctx);

    return 0;
}
