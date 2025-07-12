#pragma once

#include "Define.h"

#include <cuda_runtime.h>

//#define MAX_Z			96
#define MAX_Z			64
//#define MAX_Z			32

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
__global__ void hamrcd(const int ZCount, const int nht, float* d_ymx, float* d_ymy, float* d_ymz,
					   float* d_gsx, float* d_gsy, const int kxofs, float* d_ttw,
					   const int ntpx, const int ntpy, float* d_tpf,
					   const int nhfx, const int nhfy, const float nps, const float MP_edge, const float tpf_edge,
					   float* d_tc, const float alftscaling, const float dtc, const float tambient, float* d_gv,
					   float* d_htho, float* d_hhw, float* d_hfx, float* d_hfy, float* d_hfz, float* d_hko, float* d_hks, float* d_ept, float* d_eptAex, float* d_heb,
					   const float tf, const float dtf, const float rumach,
					   float* d_alfg, int* d_rs,
					   const float fOffset,
					   float* d_rand, float* d_YHx, float* d_YHy, float* d_YHz,
					   const int animation, const int nmwskip, const int npwskip, float* d_amdx, float* d_amdy, float* d_amdz, float* d_ztmp, float* d_zhed,
					   float* d_ELCO, float* d_TESCO,
#ifdef GPU_LOG
					   const int nLog, float* d_asx, float* d_asy, float* d_bsy, float* d_alft, float* d_amk, float* d_hhwl, float* d_tmpt, float* d_htz, float* d_hax, float* d_hay, float* d_haz, float* d_huk, float* d_hea, float* d_bms, float* d_yy, float* d_yz, float* d_amy, float* d_amz, float* d_tmpc, float* d_hthi,
#endif
					   const int Grains);


