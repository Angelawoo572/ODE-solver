#include "CUDAInterface.cuh"

#include "CUDAFunction.cuh"
#include "Define.h"

//#include <cuda_profiler_api.h>
//cudaProfilerStart();
//cudaProfilerStop();


//echo copy "$(CudaToolkitBinDir)\cudart*.dll" "$(OutDir)"
//copy "$(CudaToolkitBinDir)\cudart*.dll" "$(OutDir)"

template<class T>
cudaError_t CUDA_FREE(T *&p)
{
	cudaError_t error = cudaSuccess;
	if (nullptr != p)
		error = cudaFree(p);
	p = nullptr;
	return error;
}

template<class T>
cudaError_t CUDA_MALLOC(T *&p, const size_t nSize)
{
	CUDA_FREE<T>(p);
	cudaError_t error = cudaMalloc((void**)&p, (sizeof(T) * nSize));
	if (cudaSuccess != error)
		return error;

	error = cudaGetLastError();
	return cudaMemset(p, 0, (sizeof(T) * nSize));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
CCUDAInterface::CCUDAInterface()
{
	m_nDevice = -1;

	memset(&m_Property, 0, sizeof(cudaDeviceProp));

	m_ZCount = 0;
	m_Grains = 0;
	m_nht = 0;
	m_ntpx = 0;
	m_ntpy = 0;
	m_ntpz = 0;
	m_nmwskip = 0;
	m_npwskip = 0;

	m_tpf = nullptr;

	m_hfx = nullptr;
	m_hfy = nullptr;
	m_hfz = nullptr;

	m_ept = nullptr;
	m_eptAex = nullptr;
	m_alfg = nullptr;

	m_heb = nullptr;
	m_hko = nullptr;

	m_gv = nullptr;

	m_htho = nullptr;

	m_gsx = nullptr;
	m_gsy = nullptr;

	m_tc = nullptr;
	m_hks = nullptr;

	m_ymx = nullptr;
	m_ymy = nullptr;
	m_ymz = nullptr;

	m_ttw = nullptr;
	m_hhw = nullptr;

	m_rand = nullptr;

	m_YHx = nullptr;
	m_YHy = nullptr;
	m_YHz = nullptr;

	m_amdx = nullptr;
	m_amdy = nullptr;
	m_amdz = nullptr;
	m_ztmp = nullptr;
	m_zhed = nullptr;

#ifdef GPU_LOG
	m_asx = nullptr;
	m_asy = nullptr;
	m_bsy = nullptr;
	m_alft = nullptr;
	m_amk = nullptr;
	m_hhwl = nullptr;
	m_tmpt = nullptr;
	m_htz = nullptr;
	m_hax = nullptr;
	m_hay = nullptr;
	m_haz = nullptr;
	m_huk = nullptr;
	m_hea = nullptr;
	m_bms = nullptr;
	m_yy = nullptr;
	m_yz = nullptr;
	m_amy = nullptr;
	m_amz = nullptr;
	m_tmpc = nullptr;
	m_hthi = nullptr;
#endif

	m_rs = nullptr;

	m_ELCO = nullptr;
	m_TESCO = nullptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// デストラクタ
CCUDAInterface::~CCUDAInterface()
{
	ReleaseGPU();
	cudaDeviceReset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
cudaError_t CCUDAInterface::GetProperty(int nDevice, cudaDeviceProp &Property)
{
	cudaError_t error = cudaSuccess;

	error = cudaSetDevice(nDevice);
	if (cudaSuccess != error)
		return error;

	return cudaGetDeviceProperties(&Property, nDevice);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int CCUDAInterface::CheckGPU()
{
	cudaDeviceReset();

	m_nDevice = -1;
	memset(&m_Property, 0, sizeof(cudaDeviceProp));

	// 使用できる数を取得します。
	int nCount = 0;
	cudaError_t error = cudaGetDeviceCount(&nCount);
	if (cudaSuccess != error)
		return error;

	if (0 == nCount)
		return -1;

	// 使うのを決めます。
	cudaDeviceProp prop;
	for (int i = 0; i < nCount; i++)
	{
		error = cudaSetDevice(i);
		if (cudaSuccess != error)
			continue;

		error = cudaGetDeviceProperties(&prop, i);
		if (cudaSuccess != error)
			continue;

		// CCが7よりしたなら使用できなくします。
		if (6 > prop.major)
			continue;

		m_nDevice = i;
		m_Property = prop;
		break;
	}

	return cudaSuccess;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
void CCUDAInterface::ReleaseGPU()
{
	cudaDeviceSynchronize();

	Free();

	cudaDeviceReset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int CCUDAInterface::GetProperty(cudaDeviceProp &Property)
{
	memset(&Property, 0, sizeof(cudaDeviceProp));

	if (-1 == m_nDevice)
		return -1;

	Property = m_Property;

	return cudaSuccess;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int CCUDAInterface::Malloc(int ntpx, int ntpy, int ntpz, int ZCount, int nhfx, int nhfy, int nhfz, int Grains, int nht, int nmwskip, int npwskip, size_t &nSize)
{
	cudaError_t error = cudaSuccess;

	m_ZCount = ZCount;
	m_Grains = Grains;
	m_nht = nht;
	m_ntpx = ntpx;
	m_ntpy = ntpy;
	m_ntpz = ntpz;
	m_nhfx = nhfx;
	m_nhfy = nhfy;
	m_nmwskip = nmwskip;
	m_npwskip = npwskip;

	Free();

	// 使用するメモリサイズ
	nSize = 0;
	nSize += (sizeof(float) * m_ntpx * m_ntpy * m_ZCount);
	nSize += (sizeof(float) * m_ZCount * 9);
	nSize += (sizeof(float) * m_Grains * 3);
	nSize += (sizeof(float) * m_Grains * m_ZCount * 6);
	nSize += (sizeof(float) * m_nht * 2);
	nSize += (sizeof(float) * m_Grains * RAND_ARRAY);
	nSize += (sizeof(float) * m_Grains * MAX_Z * (MAXORD + 1) * 3);
#ifdef ALL_ANIMATION
	nSize += (sizeof(float) * m_Grains * m_nht * 3);
#else
	nSize += (sizeof(float) * m_Grains * m_nht / m_nmwskip + 1);
	nSize += (sizeof(float) * m_Grains * (m_nht / m_npwskip + 1) * 2);
#endif
	nSize += (sizeof(float) * m_Grains * SIZE_ELCO_D1 * SIZE_ELCO_D2);
	nSize += (sizeof(float) * m_Grains * SIZE_TESCO_D1 * SIZE_TESCO_D2);

#ifdef GPU_LOG
	nSize += (sizeof(float) * m_ZCount * m_nht * 17);
#endif

	// 基本
	if (cudaSuccess != (error = CUDA_MALLOC(m_tpf, (m_ntpx * m_ntpy * m_ZCount))))							return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hfx, (m_nhfx * m_nhfy * m_ZCount))))							return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hfy, (m_nhfx * m_nhfy * m_ZCount))))							return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hfz, (m_nhfx * m_nhfy * m_ZCount))))							return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ept, m_ZCount)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_eptAex, m_ZCount)))											return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_alfg, m_ZCount)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_heb, m_ZCount)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hko, m_ZCount)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_gv, m_Grains)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_htho, m_ZCount)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_gsx, m_Grains)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_gsy, m_Grains)))												return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_tc, (m_Grains * m_ZCount))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hks, (m_Grains * m_ZCount))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ymx, (m_Grains * m_ZCount))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ymy, (m_Grains * m_ZCount))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ymz, (m_Grains * m_ZCount))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ttw, m_nht)))													return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hhw, m_nht)))													return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_rand, (m_Grains * MAX_Z * RAND_ARRAY))))						return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_YHx, (m_Grains * MAX_Z * (MAXORD + 1)))))						return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_YHy, (m_Grains * MAX_Z * (MAXORD + 1)))))						return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_YHz, (m_Grains * MAX_Z * (MAXORD + 1)))))						return error;

	// アニメーション
#ifdef ALL_ANIMATION
	if (cudaSuccess != (error = CUDA_MALLOC(m_amdz, (m_Grains * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ztmp, (m_Grains * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_zhed, (m_Grains * m_nht))))									return error;
#else
	if (cudaSuccess != (error = CUDA_MALLOC(m_amdz, (m_Grains * (m_nht / m_nmwskip + 1)))))					return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_ztmp, (m_Grains * (m_nht / m_npwskip + 1)))))					return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_zhed, (m_Grains * (m_nht / m_npwskip + 1)))))					return error;
#endif

	// ログ出し時に必要
#ifdef GPU_LOG
	if (cudaSuccess != (error = CUDA_MALLOC(m_asx, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_asy, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_bsy, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_alft, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hhwl, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_tmpt, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_htz, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hax, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hay, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_haz, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_huk, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hea, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_amk, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_bms, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_yy, (m_ZCount * m_nht))))										return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_yz, (m_ZCount * m_nht))))										return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_amy, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_amz, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_tmpc, (m_ZCount * m_nht))))									return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_hthi, (m_ZCount * m_nht))))									return error;
#endif

	// 乱数種
	if (cudaSuccess != (error = CUDA_MALLOC(m_rs, (m_Grains * m_ZCount))))									return error;

	if (cudaSuccess != (error = CUDA_MALLOC(m_ELCO, (SIZE_ELCO_D1 * SIZE_ELCO_D2))))						return error;
	if (cudaSuccess != (error = CUDA_MALLOC(m_TESCO, (SIZE_TESCO_D1 * SIZE_TESCO_D2))))						return error;

	return cudaSuccess;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int CCUDAInterface::Free()
{
	CUDA_FREE(m_tpf);
	CUDA_FREE(m_hfx);
	CUDA_FREE(m_hfy);
	CUDA_FREE(m_hfz);
	CUDA_FREE(m_ept);
	CUDA_FREE(m_eptAex);
	CUDA_FREE(m_alfg);
	CUDA_FREE(m_heb);
	CUDA_FREE(m_hko);
	CUDA_FREE(m_gv);
	CUDA_FREE(m_htho);
	CUDA_FREE(m_gsx);
	CUDA_FREE(m_gsy);
	CUDA_FREE(m_tc);
	CUDA_FREE(m_hks);
	CUDA_FREE(m_ymx);
	CUDA_FREE(m_ymy);
	CUDA_FREE(m_ymz);
	CUDA_FREE(m_ttw);
	CUDA_FREE(m_hhw);
	CUDA_FREE(m_rand);
	CUDA_FREE(m_YHx);
	CUDA_FREE(m_YHy);
	CUDA_FREE(m_YHz);
	CUDA_FREE(m_amdx);
	CUDA_FREE(m_amdy);
	CUDA_FREE(m_amdz);
	CUDA_FREE(m_ztmp);
	CUDA_FREE(m_zhed);
#ifdef GPU_LOG
	CUDA_FREE(m_asx);
	CUDA_FREE(m_asy);
	CUDA_FREE(m_bsy);
	CUDA_FREE(m_alft);
	CUDA_FREE(m_amk);
	CUDA_FREE(m_hhwl);
	CUDA_FREE(m_tmpt);
	CUDA_FREE(m_htz);
	CUDA_FREE(m_hax);
	CUDA_FREE(m_hay);
	CUDA_FREE(m_haz);
	CUDA_FREE(m_huk);
	CUDA_FREE(m_hea);
	CUDA_FREE(m_bms);
	CUDA_FREE(m_yy);
	CUDA_FREE(m_yz);
	CUDA_FREE(m_amy);
	CUDA_FREE(m_amz);
	CUDA_FREE(m_tmpc);
	CUDA_FREE(m_hthi);
#endif
	CUDA_FREE(m_rs);

	CUDA_FREE(m_ELCO);
	CUDA_FREE(m_TESCO);

	return cudaGetLastError();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload2(float* thickness, float* tpf)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_tpf, tpf, (sizeof(float) * m_ntpx * m_ntpy * m_ZCount), cudaMemcpyHostToDevice)))		return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload4(float* hfx, float* hfy, float* hfz)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_hfx, hfx, (sizeof(float) * m_nhfx * m_nhfy * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_hfy, hfy, (sizeof(float) * m_nhfx * m_nhfy * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_hfz, hfz, (sizeof(float) * m_nhfx * m_nhfy * m_ZCount), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload5(float* ept, float* eptAex, float* alfg)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_ept, ept, (sizeof(float) * m_ZCount), cudaMemcpyHostToDevice)))				return error;
	if (cudaSuccess != (error = cudaMemcpy(m_eptAex, eptAex, (sizeof(float) * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_alfg, alfg, (sizeof(float) * m_ZCount), cudaMemcpyHostToDevice)))				return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload6(float* heb, float* hko)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_heb, heb, (sizeof(float) * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_hko, hko, (sizeof(float) * m_ZCount), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload7(float* gv, float* htho)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_gv, gv, (sizeof(float) * m_Grains), cudaMemcpyHostToDevice)))				return error;
	if (cudaSuccess != (error = cudaMemcpy(m_htho, htho, (sizeof(float) * m_ZCount), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload9(float* gsx, float* gsy)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_gsx, gsx, (sizeof(float) * m_Grains), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_gsy, gsy, (sizeof(float) * m_Grains), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload10(float* tc, float* hks)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_tc, tc, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_hks, hks, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload11(float* ymx, float* ymy, float* ymz)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_ymx, ymx, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_ymy, ymy, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_ymz, ymz, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload12(float* ttw, float* hhw)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_ttw, ttw, (sizeof(float) * m_nht), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_hhw, hhw, (sizeof(float) * m_nht), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload13(int* rs)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_rs, rs, (sizeof(int) * m_Grains * m_ZCount), cudaMemcpyHostToDevice)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算に必要なデータをコピー
int CCUDAInterface::Upload14(float* ELCO, float* TESCO)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(m_ELCO, ELCO, (sizeof(float) * SIZE_ELCO_D1 * SIZE_ELCO_D2), cudaMemcpyHostToDevice)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(m_TESCO, TESCO, (sizeof(float) * SIZE_TESCO_D1 * SIZE_TESCO_D2), cudaMemcpyHostToDevice)))		return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算結果の取得
int CCUDAInterface::Download1(float* ymx, float* ymy, float* ymz)
{
	cudaError_t error = cudaSuccess;

	if (cudaSuccess != (error = cudaMemcpy(ymx, m_ymx, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyDeviceToHost)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(ymy, m_ymy, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyDeviceToHost)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(ymz, m_ymz, (sizeof(float) * m_Grains * m_ZCount), cudaMemcpyDeviceToHost)))			return error;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算結果の取得
int CCUDAInterface::Download2(float* asx, float* asy, float* bsy, float* alft, float* amk, float* hhwl, float* tmpt, float* htz, float* hax, float* hay, float* haz, float* huk, float* hea, float* bms, float* yy, float* yz, float* amy, float* amz, float* tmpc, float* hthi)
{
	cudaError_t error = cudaSuccess;

#ifdef GPU_LOG
	if (cudaSuccess != (error = cudaMemcpy(asx, m_asx, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(asy, m_asy, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(bsy, m_bsy, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(alft, m_alft, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(amk, m_amk, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(hhwl, m_hhwl, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(tmpt, m_tmpt, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(htz, m_htz, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(hax, m_hax, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(hay, m_hay, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(haz, m_haz, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(huk, m_huk, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(hea, m_hea, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(bms, m_bms, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(yy, m_yy, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(yz, m_yz, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))			return error;
	if (cudaSuccess != (error = cudaMemcpy(amy, m_amy, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(amz, m_amz, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(tmpc, m_tmpc, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
	if (cudaSuccess != (error = cudaMemcpy(hthi, m_hthi, (sizeof(float) * m_ZCount * m_nht), cudaMemcpyDeviceToHost)))		return error;
#endif
	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算結果の取得
int CCUDAInterface::Download3(float* amdz, float* ztmp, float* zhed)
{
	cudaError_t error = cudaSuccess;

#ifdef ALL_ANIMATION
	if (cudaSuccess != (error = cudaMemcpy(amdz, m_amdz, (sizeof(float) * m_Grains * m_nht), cudaMemcpyDeviceToHost)))						return error;
	if (cudaSuccess != (error = cudaMemcpy(ztmp, m_ztmp, (sizeof(float) * m_Grains * m_nht), cudaMemcpyDeviceToHost)))						return error;
	if (cudaSuccess != (error = cudaMemcpy(zhed, m_zhed, (sizeof(float) * m_Grains * m_nht), cudaMemcpyDeviceToHost)))						return error;
#else
	if (cudaSuccess != (error = cudaMemcpy(amdz, m_amdz, (sizeof(float) * m_Grains * (m_nht / m_nmwskip + 1)), cudaMemcpyDeviceToHost)))	return error;
	if (cudaSuccess != (error = cudaMemcpy(ztmp, m_ztmp, (sizeof(float) * m_Grains * (m_nht / m_npwskip + 1)), cudaMemcpyDeviceToHost)))	return error;
	if (cudaSuccess != (error = cudaMemcpy(zhed, m_zhed, (sizeof(float) * m_Grains * (m_nht / m_npwskip + 1)), cudaMemcpyDeviceToHost)))	return error;
#endif
	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// GPUでの計算を実行
int CCUDAInterface::Run(const int nLog, const int kxofs, const float dtc, const float tambient, const float tf, const float alftscaling, const float dtf,  const float rumach, const float fOffset, const float nps, const float MP_edge, const float tpf_edge, const int animation)
{
	dim3 Grid(1, 1, 1);
	dim3 Block(1, 1, 1);

	Grid.x = m_Grains;
	Grid.y = 1;

	// Warp内での値まとめを使っているのでWarp単位で実行させる。
	if (32 == m_ZCount)			Block.x = 32;
	else if (64 == m_ZCount)	Block.x = 64;
	else						Block.x = (m_ZCount / 32 + 1) * 32;
	Block.y = 1;

#ifdef GPU_LOG
	Grid.x = 1;
	Grid.y = 1;
#endif

	hamrcd<<<Grid, Block>>>(m_ZCount, m_nht, m_ymx, m_ymy, m_ymz,
							m_gsx, m_gsy, kxofs, m_ttw,
							m_ntpx, m_ntpy, m_tpf,
							m_nhfx, m_nhfy, nps, MP_edge, tpf_edge,
							m_tc, alftscaling, dtc, tambient, m_gv,
							m_htho, m_hhw, m_hfx, m_hfy, m_hfz, m_hko, m_hks, m_ept, m_eptAex, m_heb,
							tf, dtf, rumach,
							m_alfg, m_rs,
							fOffset,
							m_rand, m_YHx, m_YHy, m_YHz,
							animation, m_nmwskip, m_npwskip, m_amdx, m_amdy, m_amdz, m_ztmp, m_zhed,
							m_ELCO, m_TESCO,
#ifdef GPU_LOG
							nLog, m_asx, m_asy, m_bsy, m_alft, m_amk, m_hhwl, m_tmpt, m_htz, m_hax, m_hay, m_haz, m_huk, m_hea, m_bms, m_yy, m_yz, m_amy, m_amz, m_tmpc, m_hthi,
#endif
							m_Grains);

	cudaDeviceSynchronize();


	return (int)cudaGetLastError();
}


