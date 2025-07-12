#pragma once

#pragma warning(disable : 4819)

#include "Define.h"
#include <cuda_runtime.h>


class CCUDAInterface
{
public:
	CCUDAInterface();
	virtual ~CCUDAInterface();

public:
	static cudaError_t GetProperty(int nDevice, cudaDeviceProp &Property);

	int CheckGPU();
	void ReleaseGPU();

	int GetProperty(cudaDeviceProp &Property);

	int IsGPUStart() {return (0 <= m_nDevice) ? 1: 0;}

	int Malloc(int ntpx, int ntpy, int ntpz, int ZCount, int nhfx, int nhfy, int nhfz, int Grains, int nht, int nmwskip, int npwskip, size_t &nSize);
	int Free();

	int Upload2(float* thickness, float* tpf);
	int Upload4(float* hfx, float* hfy, float* hfz);
	int Upload5(float* ept, float* eptAex, float* alfg);
	int Upload6(float* heb, float* hko);
	int Upload7(float* gv, float* htho);
	int Upload9(float* gsx, float* gsy);
	int Upload10(float* tc, float* hks);
	int Upload11(float*  ymx, float* ymy, float* ymz);
	int Upload12(float* ttw, float* hhw);
	int Upload13(int* rs);
	int Upload14(float* ELCO, float* TESCO);

	int Download1(float* ymx, float* ymy, float* ymz);
	int Download2(float* asx, float* asy, float* bsy, float* alft, float* amk, float* hhwl, float* tmpt, float* htz, float* hax, float* hay, float* haz, float* huk, float* hea, float* bms, float* yy, float* yz, float* amy, float* amz, float* tmpc, float* hthi);
	int Download3(float* amdz, float* ztmp, float* zhed);

	int Run(const int nLog, const int kxofs, const float dtc, const float tambient, const float tf, const float alftscaling, const float dtf,  const float rumach, const float fOffset, const float nps, const float MP_edge, const float tpf_edge, const int animation);

protected:
	int					m_nDevice;			// 使用デバイス

	cudaDeviceProp		m_Property;

	int					m_ZCount;			// 層数				nz
	int					m_Grains;			// セル数
	int					m_nht;				// 

	int					m_ntpx;
	int					m_ntpy;
	int					m_ntpz;

	int					m_nhfx;
	int					m_nhfy;
	int					m_nhfz;

	int					m_nmwskip;
	int					m_npwskip;

	// CUDAメモリ確保用

	float*				m_tpf;

	// Head Field Interpolation
	float*				m_hfx;
	float*				m_hfy;
	float*				m_hfz;

	// .med
	float*				m_ept;
	float*				m_eptAex;
	float*				m_alfg;

	// Exchange Field
	float*				m_heb;
	float*				m_hko;

	// 
	float*				m_gv;

	// Thermal Field
	float*				m_htho;

	// Grain Coordinate
	float*				m_gsx;
	float*				m_gsy;

	// TC distribution
	float*				m_tc;
	float*				m_hks;

	// 磁荷
	float*				m_ymx;
	float*				m_ymy;
	float*				m_ymz;

	// head field
	float*				m_ttw;
	float*				m_hhw;

	// register
	float*				m_rand;

	// shared
	float*				m_YHx;
	float*				m_YHy;
	float*				m_YHz;

	// 
	float*				m_amdx;
	float*				m_amdy;
	float*				m_amdz;
	float*				m_ztmp;
	float*				m_zhed;

#ifdef GPU_LOG
	// Start Simulation
	float*				m_asx;
	float*				m_asy;
	float*				m_bsy;
	float*				m_alft;
	float*				m_amk;
	float*				m_hhwl;
	float*				m_tmpt;
	float*				m_htz;
	float*				m_hax;
	float*				m_hay;
	float*				m_haz;
	float*				m_huk;
	float*				m_hea;
	float*				m_bms;
	float*				m_yy;
	float*				m_yz;
	float*				m_amy;
	float*				m_amz;
	float*				m_tmpc;
	float*				m_hthi;
#endif

	// Random Seed
	int*				m_rs;


	float*				m_ELCO;
	float*				m_TESCO;

};


// Beginning of GPU Architecture definitions
inline int _ConvertSMVer2Cores(int major, int minor)
{
	// Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
	typedef struct
	{
		int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] =
	{
		{ 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
		{ 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
		{ 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
		{ 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
		{ 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
		{ 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
		{ 0x53, 128}, // Maxwell Generation (SM 5.3) GM20x class
		{ 0x60, 64 }, // Pascal Generation (SM 6.0) GP100 class
		{ 0x61, 128}, // Pascal Generation (SM 6.1) GP10x class
		{ 0x62, 128}, // Pascal Generation (SM 6.2) GP10x class
		{ 0x70, 64 }, // Volta Generation (SM 7.0) GV100 class

		{   -1, -1 }
	};

	int index = 0;

	while (nGpuArchCoresPerSM[index].SM != -1)
	{
		if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
		{
			return nGpuArchCoresPerSM[index].Cores;
		}

		index++;
	}

	// If we don't find the values, we default use the previous one to run properly
//	printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[index-1].Cores);
	return nGpuArchCoresPerSM[index-1].Cores;
}

