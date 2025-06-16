/* -----------------------------------------------------------------
 * Example problem:
 *
 * This is a minimal example of solving the ODE:
 *     dy/dt = -y,  with y(0) = 1
 * using CVODE with the Adams method (nonstiff).
 * -----------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>                  // CVODE main header
#include <nvector/nvector_serial.h>       // Serial N_Vector types, functions
#include <sundials/sundials_types.h>      // definition of realtype
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

#define Ith(v, i) NV_Ith_S(v, i - 1) /* i-th vector component i=1..NEQ 
                                        1-based indexing macro for readability
                                        NV_Ith_S - N_Vector第i个分量的值
                                        NV_Ith_S(v, i) v->data[i]
*/
#define IJth(A, i, j) \
  SM_ELEMENT_D(A, i - 1, j - 1) /* (i,j)-th matrix component i,j=1..NEQ */

#define NEQ   3               /* number of equations  */
#define Y0    SUN_RCONST(1.0)  // 初始条件 y(0) = 1
#define RTOL  SUN_RCONST(1.0e-4) // scalar relative tolerance
#define ATOL  SUN_RCONST(1.0e-8) // vector absolute tolerance components
#define T0    SUN_RCONST(0.0) // initial time
#define TEND  SUN_RCONST(1.0) // 

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static void PrintOutput(sunrealtype t, sunrealtype y);
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol);

// dy/dt = -y
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data){
    // compute derivative: dy/dt = -y
    NV_Ith_S(ydot,0) = -NV_Ith_S(y,0);
    return 0;
}

// opt == 0: 检查是否返回 NULL 指针
// opt == 1: 检查返回的 int 指针的值（负数是出错）
// opt == 2: 内存错误检查
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}

static void PrintOutput(sunrealtype t, sunrealtype y)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y = %14.6Le\n", t, y);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y = %14.6e\n", t, y);
#else
  printf("At t = %0.4e      y = %14.6e\n", t, y);
#endif
}

static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
{
  sunrealtype y_true, y_num, err_weight, err, ONE = SUN_RCONST(1.0);

  y_true = exp(-t);             // exact solution at time t
  y_num  = NV_Ith_S(y, 0);      // numerical solution from CVODE

  // Compute weighted error: |y_num - y_true| / (rtol * |y_true| + atol)
  err_weight = rtol * fabs(y_true) + atol;
  if (err_weight <= 0.0) {
    fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - err_weight <= 0\n\n");
    return -1;
  }

  err = fabs(y_num - y_true) / err_weight;

  if (err < ONE) {
    printf("check_ans PASSED: err = %.6e < 1\n", err);
    return 0; // pass
  } else {
    printf("check_ans FAILED: err = %.6e >= 1\n", err);
    return 1; // fail
  }
}


int main()
{
    SUNContext sunctx;           // SUNDIALS context object
    N_Vector y = NULL;           // Solution vector
    SUNMatrix A = NULL;
    SUNLinearSolver LS = NULL;
    void *cvode_mem = NULL;      // CVODE memory block
    sunrealtype t = 0.0;            // Initial time
    sunrealtype t_end = 1.0;        // Final time
    sunrealtype dt = 0.1;           // Time step for output
    int retval;                  // Return flag for SUNDIALS functions

    // create the Sundials context
    retval = SUNContext_Create(SUN_COMM_NULL,&sunctx);
    if (retval) return 1;

    // allocate and initialize the y vector with initial condition y(0) = 1
    y = N_VNew_Serial(1,sunctx); // y(0) = 1
    if (y == NULL) return 1;
    NV_Ith_S(y,0) = 1.0; // y's 0th - y0 initial value 1.0

    // Create CVODE memory and specify the Adams method for nonstiff problems
    cvode_mem = CVodeCreate(CV_ADAMS,sunctx);
    if (cvode_mem == NULL) return 1;

    // Initialize CVODE solver with the ODE function, initial time, and initial y
    retval = CVodeInit(cvode_mem, f, t, y);
    if (retval < 0) return 1;

    // Set solver tolerances: relative = 1e-4, absolute = 1e-8
    retval = CVodeSStolerances(cvode_mem, 1e-4, 1e-8);
    if (retval < 0) return 1;

    // Optionally set user data (not used here)
    retval = CVodeSetUserData(cvode_mem, NULL);
    if (retval < 0) return 1;

    // ADD: Set up dense linear solver
    A = SUNDenseMatrix(NEQ, NEQ, sunctx); // NEQ * NEQ dense matrix A
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return 1;

    LS = SUNLinSol_Dense(y, A, sunctx); // LU A * x = b Newton
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return 1;

    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;
    // SUNLinSol_Dense() 是创建一个合法 SUNLinearSolver 对象的标准入口


    // Time-stepping loop
    while (t < t_end) {
        sunrealtype t_next = t + dt;
        retval = CVode(cvode_mem, t_next, y, &t, CV_NORMAL); // Integrate to next time
        if (retval < 0) break;

        // Print current time and y(t)
        // printf("t = %.2f, y = %.6f\n", t, NV_Ith_S(y, 0));
        PrintOutput(t, NV_Ith_S(y, 0));
    }
    check_ans(y, t, 1e-4, 1e-8);  // check the final answer

    // Free memory resources
    CVodeFree(&cvode_mem);
    N_VDestroy(y);
    SUNContext_Free(&sunctx);

    return 0;

}