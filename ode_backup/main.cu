/*
 *  main.cpp
 */

//#define EMULATION_MODE
//#define use_export	// uncomment this if project is built using a compiler that
// supports the C++ keyword "export".  If not using such a 
// compiler, be sure not to add cuLsoda.cc to the target, or
// you will get redefinition errors.

#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include "../culsoda/trunk/Code/cuLsoda_kernel.cu"

__device__ __constant__ double msk[3]   = {0.0, 0.0, 1.0};
__device__ __constant__ double chk      = 1.0;
__device__ __constant__ double che      = 0.6;
__device__ __constant__ double alpha_ld = 0.02;


int main()
{
  // —————————————————————————————————————————
  // 一些常量／尺寸
  const int nspin = 32;
  const int neq   = 3 * nspin;
  // —————————————————————————————————————————
  // 1) 在主机上分配并初始化 y
  double *h_y = (double*)malloc(sizeof(double)*neq);
  for(int i=0;i<nspin;i++){
    int b = 3*i;
    if      (i==0)          { h_y[b+0]=0.0;   h_y[b+1]=0.0;   h_y[b+2]=1.0; }
    else if (i==nspin-1)    { h_y[b+0]=0.0;   h_y[b+1]=0.0;   h_y[b+2]=-1.0; }
    else if (i < nspin/2)   { h_y[b+0]=0.0;   h_y[b+1]=0.0175;h_y[b+2]=0.998; }
    else                    { h_y[b+0]=0.0;   h_y[b+1]=0.0175;h_y[b+2]=-0.998;}
  }
  double t0    = 0.0;
  double tout0 = 0.01;     // 初始步长
  // LSODA 公共参数
  int    jt     = 2;       // 数值 Jacobian
  int    itol   = 2;       // 向量公差
  int    itask  = 1;       // 单步
  int    istate = 1;       // 初始状态
  int    iopt   = 0;       // 不用可选
  double rtol   = 1e-4;
  double atol_v = 1e-6;    // 我们用标量 tol
  // LSODA workspace 大小
  int liw = 20 + neq;
  int lrw = 22 + 16*neq;
  // —————————————————————————————————————————
  // 2) 在设备上分配所有参数 + 拷贝
  double  *d_y,  *d_t,  *d_tout, *d_rtol, *d_atol;
  int     *d_neq, *d_jt, *d_itol, *d_itask, *d_istate, *d_iopt;
  int     *d_liw, *d_lrw, *d_iwork;
  double  *d_rwork;
  cuLsodaCommonBlock h_common;
  cuLsodaCommonBlockInit(&h_common);
  cuLsodaCommonBlock *d_common;

  #define CUDAMALLOC_AND_H2D(ptr, hostptr, sz) \
    cudaMalloc((void**)&ptr, sz);               \
    cudaMemcpy(ptr, hostptr, sz, cudaMemcpyHostToDevice);

  CUDAMALLOC_AND_H2D(d_neq,    &neq,    sizeof(int));
  CUDAMALLOC_AND_H2D(d_y,      h_y,     sizeof(double)*neq);
  CUDAMALLOC_AND_H2D(d_t,      &t0,     sizeof(double));
  CUDAMALLOC_AND_H2D(d_tout,   &tout0,  sizeof(double));
  CUDAMALLOC_AND_H2D(d_rtol,   &rtol,   sizeof(double));
  CUDAMALLOC_AND_H2D(d_atol,   &atol_v, sizeof(double));  // 标量版
  CUDAMALLOC_AND_H2D(d_jt,     &jt,     sizeof(int));
  CUDAMALLOC_AND_H2D(d_itol,   &itol,   sizeof(int));
  CUDAMALLOC_AND_H2D(d_itask,  &itask,  sizeof(int));
  CUDAMALLOC_AND_H2D(d_istate, &istate, sizeof(int));
  CUDAMALLOC_AND_H2D(d_iopt,   &iopt,   sizeof(int));
  CUDAMALLOC_AND_H2D(d_liw,    &liw,    sizeof(int));
  CUDAMALLOC_AND_H2D(d_lrw,    &lrw,    sizeof(int));
  cudaMalloc((void**)&d_iwork, sizeof(int)*liw);
  cudaMalloc((void**)&d_rwork, sizeof(double)*lrw);
  cudaMalloc((void**)&d_common, sizeof(h_common));
  cudaMemcpy(d_common, &h_common, sizeof(h_common), cudaMemcpyHostToDevice);

  #undef CUDAMALLOC_AND_H2D

  // 3) 构造 functor
  myFex fex;
  myJex jex;

  // 4) 启动 cuLsoda 求解器（单线程就行）
  cuLsoda<<<1,1>>>(
    fex, d_neq,
    d_y, d_t, d_tout,
    d_itol, d_rtol, d_atol,
    d_itask, d_istate, d_iopt,
    d_rwork, d_lrw, d_iwork, d_liw,
    jex, d_jt, d_common
  );
  cudaDeviceSynchronize();

  // 5) 拷回结果并打印
  cudaMemcpy(h_y, d_y, sizeof(double)*neq, cudaMemcpyDeviceToHost);
  cudaMemcpy(&t0, d_t, sizeof(double),   cudaMemcpyDeviceToHost);

  printf(">>> 归一化时间 τ = %e\n", t0);
  for(int i=0;i<nspin;i++){
    printf(" spin %2d : m = [%.6f, %.6f, %.6f]\n",
      i, h_y[3*i+0], h_y[3*i+1], h_y[3*i+2]
    );
  }

  return 0;
}



// int main(void)   /* Main program */ 
// {
	
//     /* Local variables */
//      double *t = (double*)malloc(sizeof(double));
// 	 double *y/*[3]*/ = (double*)malloc(sizeof(double)*3);
//      int *jt = (int*)malloc(sizeof(int));
//      int *neq = (int*)malloc(sizeof(int));
// 	 int *liw = (int*)malloc(sizeof(int));
// 	 int *lrw = (int*)malloc(sizeof(int));
//      double *atol/*[3]*/ = (double*)malloc(sizeof(double)*3);
//      int *itol =(int*) malloc(sizeof(int));
// 	 int *iopt =(int*) malloc(sizeof(int));
//      double *rtol = (double*)malloc(sizeof(double));
//      int *iout =(int*) malloc(sizeof(int));
//      double *tout =(double*) malloc(sizeof(double));
//      int *itask = (int*)malloc(sizeof(int));
// 	 int *iwork/*[23]*/ =(int*) malloc(sizeof(int)*23);
//      double *rwork/*[70]*/ = (double*)malloc(sizeof(double)*70);
// 	 int *istate = (int*)malloc(sizeof(int));
// 	struct cuLsodaCommonBlock common;
// 	struct cuLsodaCommonBlock *Hcommon = &common;
// 	/* End Local Block */

// /*
// 	cudaMallocHost((void**)&t,sizeof(double));
// 	cudaMallocHost((void**)&y,sizeof(double)*3);
// 	cudaMallocHost((void**)&jt,sizeof(int));
// 	cudaMallocHost((void**)&neq,sizeof(int));
// 	cudaMallocHost((void**)&liw,sizeof(int));		
// 	cudaMallocHost((void**)&lrw,sizeof(int));	
// 	cudaMallocHost((void**)&atol,sizeof(double)*13);
// 	cudaMallocHost((void**)&itol,sizeof(int));
// 	cudaMallocHost((void**)&iopt,sizeof(int));
// 	cudaMallocHost((void**)&rtol,sizeof(double));
// 	cudaMallocHost((void**)&iout,sizeof(int));
// 	cudaMallocHost((void**)&tout,sizeof(double));
// 	cudaMallocHost((void**)&itask,sizeof(int));
// 	cudaMallocHost((void**)&iwork,sizeof(int)*23);	
// 	cudaMallocHost((void**)&rwork,sizeof(double)*70);
// 	cudaMallocHost((void**)&istate,sizeof(int));

// */	
// 	/* Pointers to Device versions of Local variables */
// 	double	*_Dt;
// 	double	*_Dy;	// [3]
// 	int	*_Djt;
// 	int	*_Dneq;
// 	int	*_Dliw;
// 	int	*_Dlrw;
//     double	*_Datol;	//[3]
//     int	*_Ditol;
// 	int	*_Diopt;
//     double	*_Drtol;
//     int	*_Diout;
//     double	*_Dtout;
//     int	*_Ditask;
// 	int	*_Diwork;	// [23]
//     double	*_Drwork;	// [70]
// 	int	*_Distate;
// 	struct cuLsodaCommonBlock *_Dcommon;
// 	/* End Pointer Block */
	
	
	
// 	/* Method instantiations for Derivative and Jacobian functions to send to template */
// 	myFex fex;
// 	myJex jex;

	
// 	/* Assignment of initial values to locals */
//     *neq = 3;
// 	y[0] = (double)1.;
//     y[1] = (double)0.;
//     y[2] = (double)0.;
//     *t = (double)0.;
//     *tout = (double).4;
// 	*itol = 2;
//     *rtol = (double)1e-4;
//     atol[0] = (double)1e-6;
//     atol[1] = (double)1e-10;
//     atol[2] = (double) 1e-6;
//     *itask = 1;
//     *istate = 1;
//     *iopt = 0;
//     *lrw = 70;
//     *liw = 23;
//     *jt = 2;
// 	cuLsodaCommonBlockInit(Hcommon);
	
// 	/* Allocate device memory for each of the pointers, and copy the values from local to device */
// 	cudaMalloc((void**)&_Dt,sizeof(double));							cudaMemcpy(_Dt,t,sizeof(double),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Dy,sizeof(double)*3);							cudaMemcpy(_Dy,y,sizeof(double)*3,cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Djt,sizeof(int));								cudaMemcpy(_Djt,jt,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Dneq,sizeof(int));								cudaMemcpy(_Dneq,neq,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Dliw,sizeof(int));								cudaMemcpy(_Dliw,liw,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Dlrw,sizeof(int));								cudaMemcpy(_Dlrw,lrw,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Datol,sizeof(double)*13);						cudaMemcpy(_Datol,atol,sizeof(double)*13,cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Ditol,sizeof(int));							cudaMemcpy(_Ditol,itol,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Diopt,sizeof(int));							cudaMemcpy(_Diopt,iopt,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Drtol,sizeof(double));							cudaMemcpy(_Drtol,rtol,sizeof(double),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Diout,sizeof(int));							cudaMemcpy(_Diout,iout,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Dtout,sizeof(double));							cudaMemcpy(_Dtout,tout,sizeof(double),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Ditask,sizeof(int));							cudaMemcpy(_Ditask,itask,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Diwork,sizeof(int)*23);						cudaMemcpy(_Diwork,iwork,sizeof(int)*23,cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Drwork,sizeof(double)*70);						cudaMemcpy(_Drwork,rwork,sizeof(double)*70,cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Distate,sizeof(int));							cudaMemcpy(_Distate,istate,sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMalloc((void**)&_Dcommon,sizeof(struct cuLsodaCommonBlock));	cudaMemcpy(_Dcommon,Hcommon,sizeof(struct cuLsodaCommonBlock), cudaMemcpyHostToDevice);
// 	/* End Allocation and Copy Block */
	
// 	int error = -1;
// 	int error2 = -1;


	
//     for (*iout = 1; *iout <= 12; ++*iout) 
// 	{
	
// 		cuLsoda<<<1,1>>>(fex, _Dneq, _Dy, _Dt, _Dtout, _Ditol, _Drtol, _Datol, _Ditask, _Distate, _Diopt, _Drwork, _Dlrw, _Diwork, _Dliw, jex, _Djt, _Dcommon);

// 		error = cudaGetLastError();		
// 		error2 = cudaDeviceSynchronize();

// 		/* Copy memory back from Device to Host */
// 		cudaMemcpy(t,_Dt,sizeof(double),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(y,_Dy,sizeof(double)*3,cudaMemcpyDeviceToHost);
// 		cudaMemcpy(jt,_Djt,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(neq,_Dneq,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(liw,_Dliw,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(lrw,_Dlrw,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(atol,_Datol,sizeof(double)*13,cudaMemcpyDeviceToHost);
// 		cudaMemcpy(itol,_Ditol,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(iopt,_Diopt,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(rtol,_Drtol,sizeof(double),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(tout,_Dtout,sizeof(double),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(itask,_Ditask,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(iwork,_Diwork,sizeof(int)*23,cudaMemcpyDeviceToHost);
// 		cudaMemcpy(rwork,_Drwork,sizeof(double)*70,cudaMemcpyDeviceToHost);
// 		cudaMemcpy(istate,_Distate,sizeof(int),cudaMemcpyDeviceToHost);
// 		cudaMemcpy(Hcommon,_Dcommon,sizeof(struct cuLsodaCommonBlock), cudaMemcpyDeviceToHost);

// 		/* End Copy Block */


// 		printf("At t =\t%- #12.2ey = %- #16.12e\t%- #16.12e\t%- #16.12e\t%5d\n", *t, y[0], y[1], y[2],iwork[10]);
// 	/*	printf("\tCM_conit = %f\n",Hcommon->CM_conit);
// 		printf("\tCM_crate = %f\n",Hcommon->CM_crate);
// 		printf("\tCM_ccmax = %f\n",Hcommon->CM_ccmax);
// 		printf("\tCM_el0 = %f\n",Hcommon->CM_el0);
// 		printf("\tCM_h__ = %f\n",Hcommon->CM_h__);
// 		printf("\tCM_hmin = %f\n",Hcommon->CM_hmin);
// 		printf("\tCM_hmxi = %f\n",Hcommon->CM_hmxi);
// 		printf("\tCM_hu = %f\n",Hcommon->CM_hu);
// 		printf("\tCM_rc = %f\n",Hcommon->CM_rc);
// 		printf("\tCM_tn = %f\n",Hcommon->CM_tn);
// 		printf("\tCM_uround = %f\n",Hcommon->CM_uround);
// 		printf("\tCM_pdest = %f\n",Hcommon->CM_pdest);
// 		printf("\tCM_pdlast = %f\n",Hcommon->CM_pdlast);
// 		printf("\tCM_ratio = %f\n",Hcommon->CM_ratio);
// 		printf("\tCM_hold = %f\n",Hcommon->CM_hold);
// 		printf("\tCM_rmax = %f\n",Hcommon->CM_rmax);
// 		printf("\tCM_tsw = %f\n",Hcommon->CM_tsw);
// 		printf("\tCM_pdnorm = %f\n",Hcommon->CM_pdnorm);*/
		
// 		error = error2 = -1;
		
// 		if (istate < 0) {
// 			goto L80;
// 		}
		

// 		/* L40: */
// 		*tout *= 10.; 
// 		cudaMemcpy(_Dtout,tout,sizeof(double),cudaMemcpyHostToDevice);
//     }
//  	printf("Number of Steps:  %i\nNo. f-s: %i\nNo. J-s = %i\nMethod Last Used = %i\nLast switch was at t = %g\n",iwork[10],iwork[11],iwork[12],iwork[18],rwork[14]);
	
// L80:
//  	printf( "STOP istate is < 0\n");
//     return 0;
// } /* MAIN__ */