/* ode1.cu : Solve dy/dt = -y, y(0)=1 on GPU using SUNDIALS+CUDA (Krylov LS) */

#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>                      // CVODE main header
#include <nvector/nvector_cuda.h>             // CUDA N_Vector
#include <sunlinsol/sunlinsol_spgmr.h>        // Krylov (SPGMR) linear solver
#include <sundials/sundials_types.h>          // realtype definition

#define NEQ 1
#define RTOL  SUN_RCONST(1.0e-4)
#define ATOL  SUN_RCONST(1.0e-8)
#define T0    SUN_RCONST(0.0)
#define TEND  SUN_RCONST(1.0)
#define DT    SUN_RCONST(0.1)

__global__ void f_kernel(sunrealtype* y, sunrealtype* ydot, int n) {
    if (threadIdx.x == 0 && blockIdx.x == 0)
        printf("[kernel debug] kernel running: n=%d, y[0]=%g\n", n, y[0]);
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        ydot[i] = -y[i];
}

// ODE function: dy/dt = -y
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data) {
    sunrealtype* y_d = N_VGetDeviceArrayPointer_Cuda(y);
    sunrealtype* ydot_d = N_VGetDeviceArrayPointer_Cuda(ydot);
    int n = N_VGetLength(y);
    printf("[f debug] t=%g, y_d=%p, ydot_d=%p, n=%d\n", t, y_d, ydot_d, n);
    fflush(stdout);

    int block = 256;
    int grid = (n + block - 1) / block;
    f_kernel<<<grid, block>>>(y_d, ydot_d, n);
    cudaDeviceSynchronize();  // 保证kernel完成
    printf("hi");
    fflush(stdout);

    return 0;
}
// static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data) {
//     N_VCopyFromDevice_Cuda(y);
//     sunrealtype *y_h = N_VGetHostArrayPointer_Cuda(y);

//     N_VCopyFromDevice_Cuda(ydot);
//     sunrealtype *ydot_h = N_VGetHostArrayPointer_Cuda(ydot);

//     ydot_h[0] = -y_h[0];

//     N_VCopyToDevice_Cuda(ydot);

//     return 0;
// }

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
    int* retval;
    if (opt == 0 && returnvalue == NULL)
    {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }
    else if (opt == 1)
    {
        retval = (int*)returnvalue;
        if (*retval < 0)
        {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
            return 1;
        }
    }
    else if (opt == 2 && returnvalue == NULL)
    {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }
    return 0;
}

static void PrintOutput(sunrealtype t, N_Vector y)
{
    N_VCopyFromDevice_Cuda(y); // device -> host
    sunrealtype* y_h = N_VGetHostArrayPointer_Cuda(y); // get host array
    // printf("At t = %.4f    y = %14.8f\n", t, y_h[0]);
    printf("At t = %.4f    y[0] = %14.8f    y[%d] = %14.8f\n",
           t, y_h[0], NEQ-1, y_h[NEQ-1]); // 打印头尾两个分量
}

// static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
// {
//     N_VCopyFromDevice_Cuda(y);
//     sunrealtype* y_h = N_VGetHostArrayPointer_Cuda(y);
//     sunrealtype y_true = exp(-t);
//     sunrealtype y_num  = y_h[0];
//     sunrealtype err_weight = rtol * fabs(y_true) + atol;
//     sunrealtype err = fabs(y_num - y_true) / err_weight;
//     if (err < 1.0) {
//         printf("check_ans PASSED: err = %.6e < 1\n", err);
//         return 0;
//     } else {
//         printf("check_ans FAILED: err = %.6e >= 1\n", err);
//         return 1;
//     }
// }

static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
{
    N_VCopyFromDevice_Cuda(y);
    sunrealtype* y_h = N_VGetHostArrayPointer_Cuda(y);
    sunrealtype y_true = exp(-t);
    double maxerr = 0.0;
    for (int i = 0; i < NEQ; i++) {
        double err = fabs(y_h[i] - y_true) / (rtol * fabs(y_true) + atol);
        if (err > maxerr) maxerr = err;
    }
    if (maxerr < 1.0) {
        printf("check_ans PASSED: maxerr = %.6e < 1\n", maxerr);
        return 0;
    } else {
        printf("check_ans FAILED: maxerr = %.6e >= 1\n", maxerr);
        return 1;
    }
}
int main()
{
    SUNContext sunctx;
    N_Vector y = NULL;
    SUNLinearSolver LS = NULL;
    void* cvode_mem = NULL;
    sunrealtype t = T0, t_end = TEND, dt = DT;
    int retval;

    // SUNDIALS context
    retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
    if (retval) return 1;

    // 1. Create CUDA N_Vector
    y = N_VNew_Cuda(NEQ, sunctx);
    if (check_retval((void*)y, "N_VNew_Cuda", 0)) return 1;
    N_VConst_Cuda(SUN_RCONST(1.0), y);

    printf("[main debug] y->ops=%p, y->content=%p\n", (void*)y->ops, (void*)y->content);
    fflush(stdout);

    // 2. Create CVODE object
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) return 1;

    retval = CVodeInit(cvode_mem, f, t, y);
    if (retval < 0) return 1;

    retval = CVodeSStolerances(cvode_mem, RTOL, ATOL);
    if (retval < 0) return 1;

    retval = CVodeSetUserData(cvode_mem, NULL);
    if (retval < 0) return 1;

    // 3. Set up Krylov linear solver (matrix-free, no dense/banded/cuda matrix needed)
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0)) return 1;
    retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

    // 4. Time stepping
    while (t < t_end) {
        sunrealtype t_next = t + dt;
        retval = CVode(cvode_mem, t_next, y, &t, CV_NORMAL);
        if (retval < 0) break;
        PrintOutput(t, y);
    }

    check_ans(y, t, RTOL, ATOL);

    // 5. Clean up
    CVodeFree(&cvode_mem);
    N_VDestroy(y);
    SUNLinSolFree(LS);
    SUNContext_Free(&sunctx);

    return 0;
}

// int main()
// {
//     SUNContext sunctx;
//     N_Vector y = NULL;
//     SUNLinearSolver LS = NULL;
//     void* cvode_mem = NULL;
//     sunrealtype t = T0, t_end = TEND, dt = DT;
//     int retval;

//     // SUNDIALS context
//     retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
//     if (retval) return 1;

//     // 1. Create CUDA N_Vector
//     y = N_VNew_Cuda(NEQ, sunctx);
//     if (check_retval((void*)y, "N_VNew_Cuda", 0)) return 1;
//     N_VConst_Cuda(SUN_RCONST(1.0), y);

//     // 2. Create CVODE object
//     cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
//     if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) return 1;

//     retval = CVodeInit(cvode_mem, f, t, y);
//     if (retval < 0) return 1;

//     retval = CVodeSStolerances(cvode_mem, RTOL, ATOL);
//     if (retval < 0) return 1;

//     retval = CVodeSetUserData(cvode_mem, NULL);
//     if (retval < 0) return 1;

//     // 3. Set up Krylov linear solver (matrix-free, no dense/banded/cuda matrix needed)
//     LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
//     if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0)) return 1;
//     retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
//     if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

//     // 4. Time stepping
//     while (t < t_end) {
//         sunrealtype t_next = t + dt;
//         retval = CVode(cvode_mem, t_next, y, &t, CV_NORMAL);
//         if (retval < 0) break;
//         PrintOutput(t, y);
//     }

//     check_ans(y, t, RTOL, ATOL);

//     // 5. Clean up
//     CVodeFree(&cvode_mem);
//     N_VDestroy(y);
//     SUNLinSolFree(LS);
//     SUNContext_Free(&sunctx);

//     return 0;
// }