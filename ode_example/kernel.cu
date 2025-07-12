#pragma warning(disable:4819)

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <cuda_runtime.h>
//#include "device_launch_parameters.h"

#include "kernel.cuh"

#define SHARED_SIZE_GO	(32 * 8)


///////////////////////////////////////////////////////////////////////////////////////////////////
// 
__global__ void GetOutput(const int nTrackX, double *pY, double *pSenseX, double *pSenseY, double *pSenseZ, double *pMagX, double *pMagY, double *pMagZ, const int nLayer, const int nHeadXMax, const int nHeadYMax, const int nMagXMax, const int nMagYMax, const int nOffset)
{
	__shared__ double dY[SHARED_SIZE_GO + 1];

	if (threadIdx.x < nHeadXMax)
	{
		dY[threadIdx.x] = 0.0f;

		const int nHeadXCentre = nHeadXMax / 2;
		const int nHeadYCentre = nHeadYMax / 2;

		int nXStart = nHeadXCentre - nTrackX;
		if (nXStart < 0)
			nXStart = 0;

		int nXStop = nMagXMax + nHeadXCentre - nTrackX;
		if (nXStop > nHeadXMax)
			nXStop = nHeadXMax;

		int nYStart = nHeadYCentre - blockIdx.x;
		if (nYStart < 0)
			nYStart = 0;

		int nYStop = nMagYMax + nHeadYCentre - blockIdx.x;
		if (nYStop > nHeadYMax)
			nYStop = nHeadYMax;

		int i, j;
		const int nSXY = (nHeadXMax * nHeadYMax * nLayer) + nXStart + threadIdx.x;
		const int nMXY = (nMagXMax * nMagYMax * nLayer) + nTrackX + nXStart - nHeadXCentre + threadIdx.x;

		for (i = nYStart; i < nYStop; i++)
		{
			int nIndexS = nSXY + (nHeadXMax * i);
			int nIndexM = nMXY + (nMagXMax * (blockIdx.x + i - nHeadYCentre));

			for (j = (nXStart + threadIdx.x); j < nXStop; j += blockDim.x)
			{
				dY[threadIdx.x] += ((pSenseX[nIndexS] * pMagX[nIndexM]) + (pSenseY[nIndexS] * pMagY[nIndexM]) + (pSenseZ[nIndexS] * pMagZ[nIndexM]));

				nIndexM += blockDim.x;
				nIndexS += blockDim.x;
			}
		}

		// Warp内で値をまとめる
		for (i = (SHARED_SIZE_GO >> 1); i > 0; i >>= 1)
		{
			__syncthreads();

			if (i > threadIdx.x)
				dY[threadIdx.x] += dY[threadIdx.x + i];
		}

		if (0 == threadIdx.x)
			pY[nOffset + blockIdx.x] = dY[0];
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
void DeviceInit()
{
	cudaSetDevice(0);
	DeviceReset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
void DeviceReset()
{
	cudaDeviceReset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int GetMagMemory(const int nPoint, const int nCrossTrack, const int nLayer, const int nSizeMag, ReciprocityMemory &RM)
{
	RM.nPoint = nPoint;
	RM.nCrossTrack = nCrossTrack;
	RM.nLayerMax = nLayer;

	// 再設定は許すがメモリ確保は1回だけ。
	if (false != RM.bLoadMag)
		return 0;

	RM.nSize = (RM.nPoint + 50) * (RM.nCrossTrack + 20) * (RM.nLayerMax + 10);

	RM.h_Y = (double*)malloc(RM.nSize * sizeof(double));

	cudaError_t error = cudaSuccess;
	error = cudaMalloc((void**)&(RM.d_MagX), (nSizeMag * 2 * sizeof(double)));
	error = cudaMalloc((void**)&(RM.d_MagY), (nSizeMag * 2 * sizeof(double)));
	error = cudaMalloc((void**)&(RM.d_MagZ), (nSizeMag * 2 * sizeof(double)));
	error = cudaMalloc((void**)&(RM.d_Y), (RM.nSize * sizeof(double)));

	RM.bLoadMag = true;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int ReleaseMemoryMag(ReciprocityMemory &RM)
{
	cudaError_t error = cudaSuccess;
	error = cudaFree(RM.d_MagX);
	error = cudaFree(RM.d_MagY);
	error = cudaFree(RM.d_MagZ);
	error = cudaFree(RM.d_Y);

	free(RM.h_Y);

	RM.Clear();

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int ReleaseMemoryHead(HeadMemory &HM)
{
	cudaError_t error = cudaSuccess;
	error = cudaFree(HM.d_HeadX);
	error = cudaFree(HM.d_HeadY);
	error = cudaFree(HM.d_HeadZ);

	HM.Clear();

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int ReciprocitySetHead(HeadMemory &HM, int nSizeSense, sense *sSense)
{
	if (false != HM.bLoadHead)
		return 0;

	cudaError_t error = cudaSuccess;
	error = cudaMalloc((void**)&(HM.d_HeadX), (nSizeSense * sizeof(double)));
	error = cudaMalloc((void**)&(HM.d_HeadY), (nSizeSense * sizeof(double)));
	error = cudaMalloc((void**)&(HM.d_HeadZ), (nSizeSense * sizeof(double)));
	error = cudaMemcpy(HM.d_HeadX, sSense->dX, (nSizeSense * sizeof(double)), cudaMemcpyHostToDevice);
	error = cudaMemcpy(HM.d_HeadY, sSense->dY, (nSizeSense * sizeof(double)), cudaMemcpyHostToDevice);
	error = cudaMemcpy(HM.d_HeadZ, sSense->dZ, (nSizeSense * sizeof(double)), cudaMemcpyHostToDevice);

	HM.bLoadHead = true;

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
int ReciprocityGetOutput(ReciprocityMemory &RM, HeadMemory &HM, int nSizeMag, sense *sMag, int nLayerMin, const int nCrossTrackOffset, const int nLayerOffset)
{
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess)
		return error;

	error = cudaMemcpy(RM.d_MagX, sMag->dX, (nSizeMag * sizeof(double)), cudaMemcpyHostToDevice);
	error = cudaMemcpy(RM.d_MagY, sMag->dY, (nSizeMag * sizeof(double)), cudaMemcpyHostToDevice);
	error = cudaMemcpy(RM.d_MagZ, sMag->dZ, (nSizeMag * sizeof(double)), cudaMemcpyHostToDevice);
//	error = cudaMemset(RM.d_Y, 0, (RM.nPoint * RM.nCrossTrack * RM.nLayerMax * sizeof(double)));
//	error = cudaMemset(RM.d_Y, 0, (RM.nSize * sizeof(double)));

	dim3 Grid(1, 1, 1);
	dim3 Block(1, 1, 1);

	Grid.x = RM.nPoint;
	Grid.y = 1;
	Block.x = SHARED_SIZE_GO;
	Block.y = 1;

	for (int i = 0; i < RM.nCrossTrack; i++)
	{
		for (int j = nLayerMin; j <= RM.nLayerMax; j++)
			GetOutput<<<Grid, Block>>>(i, RM.d_Y, HM.d_HeadX, HM.d_HeadY, HM.d_HeadZ, RM.d_MagX, RM.d_MagY, RM.d_MagZ, j, HM.nHeadX, HM.nHeadY, sMag->nXMax, sMag->nYMax, ((nCrossTrackOffset * i) + (nLayerOffset * j)));
	}

//	error = cudaMemcpy(RM.h_Y, RM.d_Y, (RM.nPoint * RM.nCrossTrack * RM.nLayerMax * sizeof(double)), cudaMemcpyDeviceToHost);
	error = cudaMemcpy(RM.h_Y, RM.d_Y, (RM.nSize * sizeof(double)), cudaMemcpyDeviceToHost);

	return error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
bool CheckGPU()
{
	cudaDeviceProp Property;
	memset(&Property, 0, sizeof(cudaDeviceProp));

	// 使用できる数を取得します。
	int nCount = 0;
	cudaError_t error = cudaGetDeviceCount(&nCount);
	if (cudaSuccess != error)
		return false;

	if (0 == nCount)
		return false;

	// 使うのを決めます。
	for (int i = 0; i < nCount; i++)
	{
		error = cudaSetDevice(i);
		if (cudaSuccess != error)
			continue;

		error = cudaGetDeviceProperties(&Property, i);
		if (cudaSuccess != error)
			continue;

		// CCが7よりしたなら使用できなくします。
		if (6 > Property.major)
			continue;

		break;
	}

	return true;
}

