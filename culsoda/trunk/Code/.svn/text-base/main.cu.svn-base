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
#include "cuLsoda_kernel.cu"


int main(void)   /* Main program */ 
{
	
    /* Local variables */
     double *t = (double*)malloc(sizeof(double));
	 double *y/*[3]*/ = (double*)malloc(sizeof(double)*3);
     int *jt = (int*)malloc(sizeof(int));
     int *neq = (int*)malloc(sizeof(int));
	 int *liw = (int*)malloc(sizeof(int));
	 int *lrw = (int*)malloc(sizeof(int));
     double *atol/*[3]*/ = (double*)malloc(sizeof(double)*3);
     int *itol =(int*) malloc(sizeof(int));
	 int *iopt =(int*) malloc(sizeof(int));
     double *rtol = (double*)malloc(sizeof(double));
     int *iout =(int*) malloc(sizeof(int));
     double *tout =(double*) malloc(sizeof(double));
     int *itask = (int*)malloc(sizeof(int));
	 int *iwork/*[23]*/ =(int*) malloc(sizeof(int)*23);
     double *rwork/*[70]*/ = (double*)malloc(sizeof(double)*70);
	 int *istate = (int*)malloc(sizeof(int));
	struct cuLsodaCommonBlock common;
	struct cuLsodaCommonBlock *Hcommon = &common;
	/* End Local Block */

/*
	cudaMallocHost((void**)&t,sizeof(double));
	cudaMallocHost((void**)&y,sizeof(double)*3);
	cudaMallocHost((void**)&jt,sizeof(int));
	cudaMallocHost((void**)&neq,sizeof(int));
	cudaMallocHost((void**)&liw,sizeof(int));		
	cudaMallocHost((void**)&lrw,sizeof(int));	
	cudaMallocHost((void**)&atol,sizeof(double)*13);
	cudaMallocHost((void**)&itol,sizeof(int));
	cudaMallocHost((void**)&iopt,sizeof(int));
	cudaMallocHost((void**)&rtol,sizeof(double));
	cudaMallocHost((void**)&iout,sizeof(int));
	cudaMallocHost((void**)&tout,sizeof(double));
	cudaMallocHost((void**)&itask,sizeof(int));
	cudaMallocHost((void**)&iwork,sizeof(int)*23);	
	cudaMallocHost((void**)&rwork,sizeof(double)*70);
	cudaMallocHost((void**)&istate,sizeof(int));

*/	
	/* Pointers to Device versions of Local variables */
	double	*_Dt;
	double	*_Dy;	// [3]
	int	*_Djt;
	int	*_Dneq;
	int	*_Dliw;
	int	*_Dlrw;
    double	*_Datol;	//[3]
    int	*_Ditol;
	int	*_Diopt;
    double	*_Drtol;
    int	*_Diout;
    double	*_Dtout;
    int	*_Ditask;
	int	*_Diwork;	// [23]
    double	*_Drwork;	// [70]
	int	*_Distate;
	struct cuLsodaCommonBlock *_Dcommon;
	/* End Pointer Block */
	
	
	
	/* Method instantiations for Derivative and Jacobian functions to send to template */
	myFex fex;
	myJex jex;

	
	/* Assignment of initial values to locals */
    *neq = 3;
	y[0] = (double)1.;
    y[1] = (double)0.;
    y[2] = (double)0.;
    *t = (double)0.;
    *tout = (double).4;
	*itol = 2;
    *rtol = (double)1e-4;
    atol[0] = (double)1e-6;
    atol[1] = (double)1e-10;
    atol[2] = (double) 1e-6;
    *itask = 1;
    *istate = 1;
    *iopt = 0;
    *lrw = 70;
    *liw = 23;
    *jt = 2;
	cuLsodaCommonBlockInit(Hcommon);
	
	/* Allocate device memory for each of the pointers, and copy the values from local to device */
	cudaMalloc((void**)&_Dt,sizeof(double));							cudaMemcpy(_Dt,t,sizeof(double),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Dy,sizeof(double)*3);							cudaMemcpy(_Dy,y,sizeof(double)*3,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Djt,sizeof(int));								cudaMemcpy(_Djt,jt,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Dneq,sizeof(int));								cudaMemcpy(_Dneq,neq,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Dliw,sizeof(int));								cudaMemcpy(_Dliw,liw,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Dlrw,sizeof(int));								cudaMemcpy(_Dlrw,lrw,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Datol,sizeof(double)*13);						cudaMemcpy(_Datol,atol,sizeof(double)*13,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Ditol,sizeof(int));							cudaMemcpy(_Ditol,itol,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Diopt,sizeof(int));							cudaMemcpy(_Diopt,iopt,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Drtol,sizeof(double));							cudaMemcpy(_Drtol,rtol,sizeof(double),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Diout,sizeof(int));							cudaMemcpy(_Diout,iout,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Dtout,sizeof(double));							cudaMemcpy(_Dtout,tout,sizeof(double),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Ditask,sizeof(int));							cudaMemcpy(_Ditask,itask,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Diwork,sizeof(int)*23);						cudaMemcpy(_Diwork,iwork,sizeof(int)*23,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Drwork,sizeof(double)*70);						cudaMemcpy(_Drwork,rwork,sizeof(double)*70,cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Distate,sizeof(int));							cudaMemcpy(_Distate,istate,sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc((void**)&_Dcommon,sizeof(struct cuLsodaCommonBlock));	cudaMemcpy(_Dcommon,Hcommon,sizeof(struct cuLsodaCommonBlock), cudaMemcpyHostToDevice);
	/* End Allocation and Copy Block */
	
	int error = -1;
	int error2 = -1;


	
    for (*iout = 1; *iout <= 12; ++*iout) 
	{
	
		cuLsoda<<<1,1>>>(fex, _Dneq, _Dy, _Dt, _Dtout, _Ditol, _Drtol, _Datol, _Ditask, _Distate, _Diopt, _Drwork, _Dlrw, _Diwork, _Dliw, jex, _Djt, _Dcommon);

		error = cudaGetLastError();		
		error2 = cudaThreadSynchronize();

		/* Copy memory back from Device to Host */
		cudaMemcpy(t,_Dt,sizeof(double),cudaMemcpyDeviceToHost);
		cudaMemcpy(y,_Dy,sizeof(double)*3,cudaMemcpyDeviceToHost);
		cudaMemcpy(jt,_Djt,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(neq,_Dneq,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(liw,_Dliw,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(lrw,_Dlrw,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(atol,_Datol,sizeof(double)*13,cudaMemcpyDeviceToHost);
		cudaMemcpy(itol,_Ditol,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(iopt,_Diopt,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(rtol,_Drtol,sizeof(double),cudaMemcpyDeviceToHost);
		cudaMemcpy(tout,_Dtout,sizeof(double),cudaMemcpyDeviceToHost);
		cudaMemcpy(itask,_Ditask,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(iwork,_Diwork,sizeof(int)*23,cudaMemcpyDeviceToHost);
		cudaMemcpy(rwork,_Drwork,sizeof(double)*70,cudaMemcpyDeviceToHost);
		cudaMemcpy(istate,_Distate,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(Hcommon,_Dcommon,sizeof(struct cuLsodaCommonBlock), cudaMemcpyDeviceToHost);

		/* End Copy Block */


		printf("At t =\t%- #12.2ey = %- #16.12e\t%- #16.12e\t%- #16.12e\t%5d\n", *t, y[0], y[1], y[2],iwork[10]);
	/*	printf("\tCM_conit = %f\n",Hcommon->CM_conit);
		printf("\tCM_crate = %f\n",Hcommon->CM_crate);
		printf("\tCM_ccmax = %f\n",Hcommon->CM_ccmax);
		printf("\tCM_el0 = %f\n",Hcommon->CM_el0);
		printf("\tCM_h__ = %f\n",Hcommon->CM_h__);
		printf("\tCM_hmin = %f\n",Hcommon->CM_hmin);
		printf("\tCM_hmxi = %f\n",Hcommon->CM_hmxi);
		printf("\tCM_hu = %f\n",Hcommon->CM_hu);
		printf("\tCM_rc = %f\n",Hcommon->CM_rc);
		printf("\tCM_tn = %f\n",Hcommon->CM_tn);
		printf("\tCM_uround = %f\n",Hcommon->CM_uround);
		printf("\tCM_pdest = %f\n",Hcommon->CM_pdest);
		printf("\tCM_pdlast = %f\n",Hcommon->CM_pdlast);
		printf("\tCM_ratio = %f\n",Hcommon->CM_ratio);
		printf("\tCM_hold = %f\n",Hcommon->CM_hold);
		printf("\tCM_rmax = %f\n",Hcommon->CM_rmax);
		printf("\tCM_tsw = %f\n",Hcommon->CM_tsw);
		printf("\tCM_pdnorm = %f\n",Hcommon->CM_pdnorm);*/
		
		error = error2 = -1;

		if (istate < 0) {
			goto L80;
		}
		

		/* L40: */
		*tout *= 10.; 
		cudaMemcpy(_Dtout,tout,sizeof(double),cudaMemcpyHostToDevice);
    }
 	printf("Number of Steps:  %i\nNo. f-s: %i\nNo. J-s = %i\nMethod Last Used = %i\nLast switch was at t = %g\n",iwork[10],iwork[11],iwork[12],iwork[18],rwork[14]);
	
L80:
 	printf( "STOP istate is < 0\n");
    return 0;
} /* MAIN__ */


