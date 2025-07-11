/*
 *  cuLsoda.h
 *
 *	File Notes:	This file is a conversion of the double precision Livermore Solver for
 *	Ordinary Differential Equations with automatic switching for stiff and non-stiff
 *	problems (DLSODA)
 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
#ifndef CULSODA_CU_H
#define CULSODA_CU_H




/* Common Block Declarations */
struct cuLsodaCommonBlock
{
	double /*rowns[209],*/ CM_conit, CM_crate, CM_ccmax, CM_el0, CM_h__, CM_hmin, CM_hmxi, CM_hu, CM_rc, CM_tn, CM_uround, CM_pdest, CM_pdlast, CM_ratio, CM_hold, CM_rmax;
	double  CM_el[13], CM_elco[156]	/* was [13][12] */, CM_tesco[36]	/* was [3][12] */;
	double CM_rls[218];
	double CM_tsw, /*rowns2[20],*/ CM_pdnorm;
	double /*rownd2,*/ CM_cm1[12], CM_cm2[5];
	double CM_rlsa[22];
	double CM_sm1[12];
	int CM_init, CM_mxstep, CM_mxhnil, CM_nhnil, CM_nslast, CM_nyh, /*iowns[6],*/ CM_icf, 
	CM_ierpj, CM_iersl, CM_jcur, CM_jstart, CM_kflag, CM_l, CM_lyh, CM_lewt, CM_lacor, CM_lsavf,
	CM_lwm, CM_liwm, CM_meth, CM_miter, CM_maxord, CM_maxcor, CM_msbp, CM_mxncf, CM_n, CM_nq, 
	CM_nst, CM_nfe, CM_nje, CM_nqu;
	int /*iownd[6],*/ CM_ialth, CM_ipup, CM_lmax, /*meo,*/ CM_nqnyh, CM_nslp;
	int CM_ils[37];
	int CM_insufr, CM_insufi, CM_ixpr, /*iowns2[2],*/ CM_jtyp, CM_mused, CM_mxordn, CM_mxords; 
	int /*iownd2[3],*/ CM_icount, CM_irflag;
	int CM_ilsa[9];
};

/* End Common Block */ 


#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


#define Fex_and_Jex_definition
struct myFex
{
	__device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/)
	{
		ydot[0] = (double)1.0E4 * y[1] * y[2] - (double).04E0 * y[0];
		ydot[2] = (double)3.0E7 * y[1] * y[1];
		ydot[1] = (double)-1.0 * (ydot[0] + ydot[2]);
	}
};

struct myJex
{
	__device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/)
	{
		return;
	}
};




/* dlsoda.f -- translated by f2c (version 20090411).
 You must link the resulting object file with libf2c:
 on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm
 or, if you install libf2c.a in a standard place, with -lf2c -lm
 -- in that order, at the end of the command line, as in
 cc *.o -lf2c -lm
 Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
 
 http://www.netlib.org/f2c/libf2c.zip
 */




template<typename Fex, typename Jex>
__device__ int dlsoda_(Fex, int *, double *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, int *, int *, Jex, int *, struct cuLsodaCommonBlock *);

template<typename Fex, typename Jex> 
__device__ int dstoda_(int *neq, double *y, double *yh, int *NOT_nyh, double *yh1, double *ewt, double *savf, double *acor, double *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

template<typename Fex, typename Jex> 
__device__ int dprja_(int *neq, double *y, double *yh, int *NOT_nyh, double *ewt, double *ftem, double *savf, double *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

__device__ int dsolsy_(double *wm, int *iwm, double *x, double *tem, struct cuLsodaCommonBlock *common);
__device__ int dintdy_(double *t, int k, double *yh, int *NOT_nyh, double *dky, int *iflag, struct cuLsodaCommonBlock *common);
__device__ int dcfode_(int meth, double *DCFODE_elco, double *DCFODE_tesco, struct cuLsodaCommonBlock *common);
__device__ int dsolsy_(double *wm, int *iwm, double *x, double *tem, struct cuLsodaCommonBlock *common);
__device__ int dewset_(int *PARAM_n, int *itol, double *rtol, double *atol, double *ycur, double *ewt, struct cuLsodaCommonBlock *common);
__device__ double dmnorm_(int *PARAM_n, double *v, double *w, struct cuLsodaCommonBlock *common);
__device__ double dfnorm_(int *PARAM_n, double *a, double *w, struct cuLsodaCommonBlock *common);
__device__ double dbnorm_(int *PARAM_n, double *a, int *nra, int *ml, int *mu, double *w, struct cuLsodaCommonBlock *common);
__device__ int dsrcma_(double *rsav, int *isav, int *job, struct cuLsodaCommonBlock *common);
__device__ int dgefa_(double *a, int *lda, int *PARAM_n, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
__device__ int dgesl_(double *a, int *lda, int *PARAM_n, int *ipvt, double *b, int job, struct cuLsodaCommonBlock *common);
__device__ int dgbfa_(double *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
__device__ int dgbsl_(double *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, double *b, int job, struct cuLsodaCommonBlock *common);
__device__ double dumach_( struct cuLsodaCommonBlock *common);
//__device__ int xsetf_(int *mflag, struct CommonBlock *common);
//__device__ int xsetun_(int *lun, struct CommonBlock *common);
__device__ int ixsav_(int ipar, int *ivalue, int iset, struct cuLsodaCommonBlock *common);
__device__ int idamax_(int *PARAM_n, double *dx, int incx, struct cuLsodaCommonBlock *common);
__device__ int daxpy_(int *PARAM_n, double *da, double *dx, int incx, double *dy, int incy, struct cuLsodaCommonBlock *common);
__device__ int dumsum_(double a, double b, double *c__, struct cuLsodaCommonBlock *common);
__device__ int dscal_(int *PARAM_n, double *da, double *dx, int incx, struct cuLsodaCommonBlock *common);
__device__ double ddot_(int *PARAM_n, double *dx, int incx, double *dy, int incy, struct cuLsodaCommonBlock *common);
__device__ double d_sign(double *a, double *b);
__host__ __device__ void cuLsodaCommonBlockInit(struct cuLsodaCommonBlock *common);

#ifndef use_export
#include "cuLsoda.cu"
#endif

#endif



