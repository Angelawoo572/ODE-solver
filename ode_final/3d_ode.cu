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

// constant memory
__constant__ sunrealtype c_he;
__constant__ sunrealtype c_hk;
__constant__ sunrealtype c_hap;
__constant__ sunrealtype c_ap;

/* Problem Constants */
#define GROUPSIZE 3               /* number of equations per group */
#define indexbound 2
#define RTOL      SUN_RCONST(1.0e-4) /* scalar relative tolerance            */
#define ATOL1     SUN_RCONST(1.0e-3) /* vector absolute tolerance components */
#define ATOL2     SUN_RCONST(1.0e-3)
#define ATOL3     SUN_RCONST(1.0e-3)
#define T0        SUN_RCONST(0.0)  /* initial time           */
#define T1        SUN_RCONST(0.4)  /* first output time      */
#define TMULT     SUN_RCONST(5.0) /* output time factor     */
#define NOUT      12               /* number of output times */

#define ZERO SUN_RCONST(0.0)

__constant__ float msk[3]={0.0f,0.0f,1.0f};
__constant__ float chk=1.0f;
__constant__ float che=0.2f;
__constant__ float alpha=0.02;   


/* user data structure for parallel*/
typedef struct
{
  int ngroups; // number of groups
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
  int neq)
{
  sunindextype i, j, k, tid,iq,ip,ix,iy,iz,imsk;
  // thread index
  tid = blockDim.x*blockIdx.x + threadIdx.x;
  if ( tid > indexbound && tid < blockDim.x - GROUPSIZE){
    iq=tid-3; // 前一组位置
    ip=tid+3; // 后一组位置
    ix=tid-(tid)%3; // ix = 3 * (tid / 3)
    iy=ix+1;
    iz=iy+1;
    imsk=tid%3; // tid在3个一组的thread的相对位置 x = 0, y = 1, z = 2
    /*
    normalize effective field, vector f
    che*(y[iq]+y[ip]); exchange interaction
    msk[imsk]*chk*y[iz]; AnisotropyTrem
     */
    h[tid] = che*(y[iq]+y[ip])+msk[imsk]*chk*y[iz];
  }
  __syncthreads();

 if (tid >= 3 && tid <= neq - 4) {
    int gid = tid / 3;         // group ID
    int base = 3 * gid;        // base index for this group's 3-vector
    int imsk = tid % 3;        // 0: x, 1: y, 2: z

    float mx = y[base];
    float my = y[base + 1];
    float mz = y[base + 2];

    float fx = h[base];
    float fy = h[base + 1];
    float fz = h[base + 2];

    // m × f
    float mxf_x = my * fz - mz * fy;
    float mxf_y = mz * fx - mx * fz;
    float mxf_z = mx * fy - my * fx;

    // m ⋅ f
    float mdotf = mx * fx + my * fy + mz * fz;

    // m × (m × f) = mdotf * m - f
    float mmxf_x = mdotf * mx - fx;
    float mmxf_y = mdotf * my - fy;
    float mmxf_z = mdotf * mz - fz;

    // Gilbert form: dm/dτ = - (m × f + α m × (m × f)) / (1 + α²)
    float denom = 1.0f + alpha * alpha;

    float dmx = (-mxf_x - alpha * mmxf_x) / denom;
    float dmy = (-mxf_y - alpha * mmxf_y) / denom;
    float dmz = (-mxf_z - alpha * mmxf_z) / denom;

    if (imsk == 0) yd[tid] = dmx;
    if (imsk == 1) yd[tid] = dmy;
    if (imsk == 2) yd[tid] = dmz;
  }
  else
  {
    yd[tid]=0;
    // printf("DEBUG: entering kernel, neq = %d\n", neq);
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

    unsigned block_size = 96;
    // total threads = grid_size * block_size
    // grid_size is ceil - (a+b-1)/b
    unsigned grid_size  = 1; // 1 (udata->neq + block_size - 1) / block_size
    f_kernel<<<grid_size, block_size>>>(ydata, ydotdata,udata->d_h,
      udata->d_mh, udata->neq);

    cudaDeviceSynchronize();

    //debug
    sunrealtype h_ydot[9];
    cudaMemcpy(h_ydot, ydotdata + 3, 3 * sizeof(sunrealtype), cudaMemcpyDeviceToHost);
    // printf("ydot sample (group 1): %f %f %f\n", h_ydot[0], h_ydot[1], h_ydot[2]);
    
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
    sunrealtype t, tout;
    N_Vector y, abstol; // SUNDIALS vector structures for solution and absolute tolerance
    SUNLinearSolver LS; // Linear solver object (cuSolverSp QR)
    SUNNonlinearSolver NLS;
    void* cvode_mem; // CVODE integrator memory
    int retval, iout; // return status and output counter
    int neq, ngroups, groupj;// Problem size: number of equations, groups, and loop index
    UserData udata;

    /* Parse command-line to get number of groups */
    ngroups = 32;
    neq     = ngroups * GROUPSIZE;

    /* Fill user data */
    udata.ngroups = ngroups;
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
    int nspin = ngroups; //32
    int ix, iy, iz;

    for(int i=0;i<nspin;i++)
    {
	    ix=3*i;
	    iy=ix+1;
	    iz=iy+1;

	    if(i==0)
	    {
		    ydata[ix]=0.0;
		    ydata[iy]=0.0;
		    ydata[iz]=1.0;
	    }
	    else if(i == nspin-1)
	    {
		    ydata[ix]=0.0;
		    ydata[iy]=0.0;
		    ydata[iz]=-1;;
	    }
	    else if(i < nspin/2)
	    {
		    ydata[ix]=0.0;
		    ydata[iy]=0.0175;
		    ydata[iz]=0.998;
	    }
	    else
	    {
		    ydata[ix]=0.0;
		    ydata[iy]=0.0175;
		    ydata[iz]=-0.998;
      }
    }

    for (int i = 0; i < neq; i += 3) {
        abstol_data[i]   = ATOL1;
        abstol_data[i+1] = ATOL2;
        abstol_data[i+2] = ATOL3;
    }
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(abstol);

    /* Create and initialize CVODE solver memory */
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    CVodeInit(cvode_mem, f, T0, y);
    CVodeSetUserData(cvode_mem, &udata);
    CVodeSVtolerances(cvode_mem, RTOL, abstol);

    CVodeSetMaxStep(cvode_mem, SUN_RCONST(1.0));  // Set max step size to 1.0 units of time
    CVodeSetMinStep(cvode_mem, 1e-12);

    /* Matrix-free GMRES linear solver (no Jacobian needed) */
    NLS = SUNNonlinSol_Newton(y, sunctx);
    CVodeSetNonlinearSolver(cvode_mem, NLS);
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
    CVodeSetLinearSolver(cvode_mem, LS, NULL);

    /* Print header */
    printf("\nGroup of independent 3-species kinetics problems\n\n");
    printf("number of groups = %d\n\n", ngroups);

    /* Time-stepping loop */
    iout = T0;
    tout = T1;
    while (iout < NOUT) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        // ********* FOR DEBUGGING ********** //
        long int nst, netf, ncfn;
        CVodeGetNumSteps(cvode_mem, &nst);
        CVodeGetNumErrTestFails(cvode_mem, &netf);
        CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
        printf("nst = %ld, netf = %ld, ncfn = %ld\n", nst, netf, ncfn);
        // ********* END DEBUGGING ********** //
        // copy solution back to host and print all groups
        N_VCopyFromDevice_Cuda(y);
        ydata = N_VGetHostArrayPointer_Cuda(y);
        for (groupj = 0; groupj < ngroups; groupj ++) {
            printf("group %d: ", groupj);
            PrintOutput(t,
                        ydata[GROUPSIZE * groupj],
                        ydata[1 + GROUPSIZE * groupj],
                        ydata[2 + GROUPSIZE * groupj]);
        }
        if (retval == CV_SUCCESS) {
            iout++;
            tout *= TMULT;
        }else {
            fprintf(stderr, "CVode error at output %d: retval = %d\n", iout, retval);
            break;
        }
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
