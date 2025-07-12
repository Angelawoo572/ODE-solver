#pragma once

#include <vector>
#include "reciprocity.h"
//#include "waveform.h"


struct HeadMemory
{
	HeadMemory() {Clear();}
	void Clear()
	{
		nHeadX = 0;
		nHeadY = 0;
		nHeadZ = 0;
		nHeadSize = 0;
		bLoadHead = false;
		vecHeadAttenuationZ.clear();
		d_HeadX = nullptr;
		d_HeadY = nullptr;
		d_HeadZ = nullptr;
	}

	int			nHeadX;
	int			nHeadY;
	int			nHeadZ;
	int			nHeadSize;

	bool		bLoadHead;

	std::vector<double>		vecHeadAttenuationZ;

	double*		d_HeadX;
	double*		d_HeadY;
	double*		d_HeadZ;
};


struct ReciprocityMemory
{
	ReciprocityMemory() {Clear();}
	void Clear()
	{
		nPoint = 0;
		nCrossTrack = 0;
		nLayerMax = 0;
		nSize = 0;
		bLoadMag = false;
		bFail = false;
		h_Y = nullptr;
		d_MagX = nullptr;
		d_MagY = nullptr;
		d_MagZ = nullptr;
		d_Y = nullptr;
		vecWholeWaveform.clear();
		vecWaveform.clear();
	}

	int			nPoint;
	int			nCrossTrack;
	int			nLayerMax;

	size_t		nSize;

	bool		bLoadMag;
	bool		bFail;

	double*		h_Y;

	double*		d_MagX;
	double*		d_MagY;
	double*		d_MagZ;
	double*		d_Y;

	std::vector<std::vector<double>> vecWholeWaveform;					// CrossTrack >> DownTrack
	std::vector<std::vector<std::vector<double>>> vecWaveform;			// Z >> CrossTrack >> DownTrack

};

void DeviceInit();
void DeviceReset();
int GetMagMemory(const int nPoint, const int nCrossTrack, const int nLayer, const int nSizeMag, ReciprocityMemory &RM);
int ReleaseMemoryMag(ReciprocityMemory &RM);
int ReleaseMemoryHead(HeadMemory &HM);
int ReciprocitySetHead(HeadMemory &HM, const int nSizeSense, sense *sSense);
int ReciprocityGetOutput(ReciprocityMemory &RM, HeadMemory &HM, int nSizeMag, sense *sMag, int nLayerMin, const int nCrossTrackOffset, const int nLayerOffset);
bool CheckGPU();

