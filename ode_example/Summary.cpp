// Summary.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include "Define.h"

//#include <crtdbg.h>
#include <stdio.h>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <locale>
#include <codecvt>
#include <regex>
#include <numeric>
#include <cmath>
#include <charconv>
#include <thread>

#include "Summary.h"
#include "polygon_io.h"
#include "reciprocity.h"
#include "kernel.cuh"

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
struct CELL
{
	CELL()
	{
		vertex.clear();
	}

	std::vector<point>		vertex;
};
using VectorCell = std::vector<CELL>;

struct PARAMETER
{
	PARAMETER() {Clear();}
	void Clear()
	{
		memset(fimed, 0, 256);
		memset(fihfd, 0, 256);
		memset(fitpf, 0, 256);
		nLayer = 30;
		nNWTemp = 0;
		nNWDynamic = 0;
		nNWTime = 0;
		nInitial = 0;
		dVhksigma = 0.0;
		dAlftscaling = 0.0;
		dTc = 0.0;
		dAmbient = 0.0;
		dLinerSpeed = 0.0;
		dFactionMotion = 0.0;
		dTPFOffset = 0.0;
		dOffset = 0.0;
		dGrainSize = 0.0;
		dThickness = 0.384;
		dThermalTimeStep = 0.0;
		dHeadField = 0.0;
		dRiseTime = 0.0;
		dAntf = 0.0;
		dBitLength = 18.0;
		dBitDelay = 0.0;
		dSigmaGS = 0.0;
		d0m = 0.0;
		dV = 0.0;
		vecThickness.clear();
	}

	char		fimed[256];
	char		fihfd[256];
	char		fitpf[256];
	int			nLayer;				// Z									m_ZCount
	int			nNWTemp;			// nw-temp								m_npwrite
	int			nNWDynamic;			// nw-dynamic							m_nmwrite
	int			nNWTime;			// nw-time								m_ntwrite
	int			nInitial;			// initial Mag							m_ichki
	double		dVhksigma;			//										m_vhksigma
	double		dAlftscaling;		//										m_alftscaling
	double		dTc;				// d_Tc									m_dtc
	double		dAmbient;			// T-ambient							m_tambient
	double		dLinerSpeed;		// Linear Speed(m/s)					m_spd
	double		dFactionMotion;		// fraction motion						m_pmove
	double		dTPFOffset;			// tpf-offset							m_pstart
	double		dOffset;			// offset(nm)							m_offset
	double		dGrainSize;			// grain size							m_gsz
	double		dThickness;			// c-lattice (nm)						m_asz
	double		dThermalTimeStep;	// thermal-time-step (ns) (dtsim)		m_dtsimi
	double		dAntf;				// antf									m_antf
	double		dHeadField;			// Head Field factor					m_hapi
	double		dRiseTime;			// rise time (ns)						m_trise
	double		dBitLength;			// bit-length (nm)						m_sbit
	double		dBitDelay;			// bit delay (nm)						m_sdelay
	double		dSigmaGS;			// sigma_GS_per							m_sigma_gs
	double		d0m;				// d0m									m_gstc_d
	double		dV;					// v									m_gstc_v

	std::vector<double>		vecThickness;	// 各層の厚さ

};

///////////////////////////////////////////////////////////////////////////////////////////////////
// 平均
template <template<class T, class Allocator = std::allocator<T>> class Container>
double Average(Container<double> &x)
{
	return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 分散
template <template<class T, class Allocator = std::allocator<T>> class Container>
double Variance(Container<double> &x)
{
	if (0 == x.size())
		return 0.0;

	const double size = (double)x.size();
	const double x_mean = Average(x);
	return (std::inner_product(x.begin(), x.end(), x.begin(), 0.0) - x_mean * x_mean * size)/ (size - 1.0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 標準偏差
template <template<class T, class Allocator = std::allocator<T>> class Container>
double StandardDeviation(Container<double> &x)
{
	if (0 == x.size())
		return 0.0;

	return std::sqrt(Variance(x));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// セル情報を作成する。
bool LoadParameter(const char* pszPARFile, PARAMETER &Parameter)
{
	FILE *pFile = fopen(pszPARFile, "rt");
	if (nullptr == pFile)
	{
		printf("**********************************************************************\n");
		printf("Failed to file open. File : %s\n", pszPARFile);
		printf("**********************************************************************\n");
		return false;
	}

	char szTemp[FILE_READ_BUFFER];

	// vertical_Hk_sigma
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf", &Parameter.dVhksigma, &Parameter.dAlftscaling);

	// d_Tc,T-ambient
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf", &Parameter.dTc, &Parameter.dAmbient);

	// Linear Speed(m/s),fraction motion, tpf-offset
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf %lf", &Parameter.dLinerSpeed, &Parameter.dFactionMotion, &Parameter.dTPFOffset);

	// offset(nm)
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf", &Parameter.dOffset);

	// grain size, c-lattice (nm), # of Fe monolayers
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf %d", &Parameter.dGrainSize, &Parameter.dThickness, &Parameter.nLayer);

	// thermal-time-step (ns) (dtsim)
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf", &Parameter.dThermalTimeStep, &Parameter.dAntf);

	// Head Field factor, rise time (ns)
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf", &Parameter.dHeadField, &Parameter.dRiseTime);

	// bit-length (nm), bit delay (nm)
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf", &Parameter.dBitLength, &Parameter.dBitDelay);

	// nw-temp, nw-dynamic, nw-time
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%d %d %d", &Parameter.nNWTemp, &Parameter.nNWDynamic, &Parameter.nNWTime);

	// initial Mag (1,-1,0)
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%d", &Parameter.nInitial);

	// sigma_GS_per, d0m, v
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%lf %lf %lf", &Parameter.dSigmaGS, &Parameter.d0m, &Parameter.dV);

	//'.med'
	char szFile[FILE_READ_BUFFER] = "";
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%s", szFile);
	sprintf(Parameter.fimed, "%s%s", szFile, EXTENSION_MED);

	//'.tpf'
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%s", szFile);
	sprintf(Parameter.fitpf, "%s%s", szFile, EXTENSION_TPF);

	//'.hfd'
	fgets(szTemp, FILE_READ_BUFFER, pFile);
	sscanf(szTemp, "%s", szFile);
	sprintf(Parameter.fihfd, "%s%s", szFile, EXTENSION_HFD);

	fclose(pFile);

	// 
	Parameter.vecThickness.clear();
	for (int i = 0; i < Parameter.nLayer; i++)
		Parameter.vecThickness.push_back(Parameter.dThickness);

	std::filesystem::path path(pszPARFile);
	std::string str = path.parent_path().string() + "/" + Parameter.fimed;
	pFile = fopen(str.c_str(), "rt");
	if (nullptr != pFile)
	{
		int nIndex = 0;
		while (NULL != fgets(szTemp, FILE_READ_BUFFER, pFile))
		{
			int iz;
			float dummy;
			float thickness = 0.0f;

			sscanf(szTemp, "%d %f %f %f %f %d %f %f %f %f %f", &iz, &dummy, &dummy, &dummy, &dummy, &iz, &dummy, &dummy, &dummy, &dummy, &thickness);

			if (0.0f != thickness)
				Parameter.vecThickness[nIndex] = thickness;

			nIndex++;
		}

		fclose(pFile);
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
bool LoadReadHead(const char* pszHeadFile, HeadMemory &HM)
{
	printf("Load Head.\n");

	sense sSense;
	HM.nHeadSize = Reciprocity_Load_Sense_Distribution(&sSense, pszHeadFile, 0);

	HM.nHeadX = sSense.nXMax;
	HM.nHeadY = sSense.nYMax;
	HM.nHeadZ = sSense.nZMax;

	const int nIndex = (sSense.nXMax / 2) + (sSense.nXMax * (sSense.nYMax / 2));
	const double dBase = sSense.dZ[nIndex];
	for (int i = 0; i < sSense.nZMax; i++)
		HM.vecHeadAttenuationZ.push_back(sSense.dZ[nIndex + (sSense.nXMax * sSense.nYMax * i)] / dBase);

	ReciprocitySetHead(HM, HM.nHeadSize, &sSense);

	Reciprocity_Sense_End(&sSense);

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// セル情報を作成する。
bool CreateBpol(const char* pszGPFile, const char* pszBpolFile, const int nLayer, int &nCrossTrack, int &nDownTrack)
{
	if (false != std::filesystem::exists(pszBpolFile))
	{
		POLYGONATTRIBUTES polyData;
		Polygon_IO_Load_Pol_Binary(&polyData, pszBpolFile);
		nCrossTrack = (int)(polyData.pSize.x / NM_2_A);
		nDownTrack = (int)(polyData.pSize.y / NM_2_A);
		Polygon_IO_End(&polyData);

		return true;
	}

	// gpファイルの読み込み。
	FILE *pFile = fopen(pszGPFile, "rt");
	if (nullptr == pFile)
	{
		printf("**********************************************************************\n");
		printf("Failed to file open. File : %s", pszGPFile);
		printf("**********************************************************************\n");
		return false;
	}

	VectorCell vecCell;
	int nVertexMax = 0;

	char szTemp[256];
	while (fgets(szTemp, 256, pFile))
	{
		// 番号と要素数。
		double d1, d2;
		sscanf(szTemp, " %lf %lf", &d1, &d2);

		CELL cell;
		for (int i = 0; i < (int)d2; i++)
		{
			fgets(szTemp, 256, pFile);

			// 頂点
			point pt;
			sscanf(szTemp, " %lf %lf", &(pt.x), &(pt.y));

			// そのままだと小さいので倍率を上げている。
			pt.x *= NM_2_A;
			pt.y *= NM_2_A;

			cell.vertex.push_back(pt);
		}

		vecCell.push_back(cell);

		// ついでに最大の頂点数を求める。
		if (nVertexMax < (int)cell.vertex.size())
			nVertexMax = (int)cell.vertex.size();
	}

	fclose(pFile);

	// セル数
	const int nPolygons = (int)vecCell.size();

	// 頂点の最大と最小を求める。
	double dXMin = 9999.0;
	double dXMax = -9999.0;
	double dYMin = 9999.0;
	double dYMax = -9999.0;

	for (size_t i = 0; i < vecCell.size(); i++)
	{
		for (size_t j = 0; j < vecCell[i].vertex.size(); j++)
		{
			double x = vecCell[i].vertex[j].x;
			if (dXMin > x)		dXMin = x;
			if (dXMax < x)		dXMax = x;

			double y = vecCell[i].vertex[j].y;
			if (dYMin > y)		dYMin = y;
			if (dYMax < y)		dYMax = y;
		}
	}

	const double dXSize = abs(dXMin) + abs(dXMax);;
	const double dYSize = abs(dYMin) + abs(dYMax);;

	// bpolは0始まりなのでマイナス値がなくなるようにオフセットする。
	for (size_t i = 0; i < vecCell.size(); i++)
	{
		double dSumX = 0.0;
		double dSumY = 0;
		for (size_t j = 0; j < vecCell[i].vertex.size(); j++)
		{
			vecCell[i].vertex[j].x += abs(dXMin);
			vecCell[i].vertex[j].y += abs(dYMin);
		}
	}

	// bpol用のバッファ。
	POLYGONATTRIBUTES* polyData = nullptr;
	polyData = new POLYGONATTRIBUTES;
	if (nullptr == polyData)
		return false;

	Polygon_IO_Reset_Polygon_Attributes(polyData);

	// 基本情報をぶっこむ。
	polyData->pSize.x = dXSize;
	polyData->pSize.y = dYSize;
	polyData->nPolygons = nPolygons * nLayer;
	polyData->nVertexMax = nVertexMax;
	Polygon_IO_Setup_Polygon_Attributes(&polyData, 0, 1, 0, 0);

	// gpから読み込んだ頂点を突っ込む。
	for (int i = 0; i < nLayer; i++)
	{
		for (size_t j = 0; j < vecCell.size(); j++)
		{
			const int nCell = i * nPolygons + (int)j;
			double dSumX = 0.0;
			double dSumY = 0.0;
			for (size_t k = 0; k < vecCell[j].vertex.size(); k++)
			{
				dSumX += vecCell[j].vertex[k].x;
				dSumY += vecCell[j].vertex[k].y;
			}

			// 母点は重心とする。
			polyData->pSeed[nCell].x = dSumX / (double)vecCell[j].vertex.size();
			polyData->pSeed[nCell].y = dSumY / (double)vecCell[j].vertex.size();

			polyData->nPolyLayer[nCell] = i;
			polyData->siPolyVertices[nCell] = (short)vecCell[j].vertex.size();

			for (size_t k = 0; k < vecCell[j].vertex.size(); k++)
				Polygon_IO_Set_Vertex(polyData, (int)nCell, (int)k, vecCell[j].vertex[k].x, vecCell[j].vertex[k].y);
		}
	}

	// bpolに書き出す。
	Polygon_IO_Save_Pol_Binary(pszBpolFile, polyData);

	nCrossTrack = (int)(polyData->pSize.x / NM_2_A);
	nDownTrack = (int)(polyData->pSize.y / NM_2_A);

	Polygon_IO_End(polyData);

	delete polyData;

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Waveformmを求める。
bool RunReciprocity(const char* pszEndFile, const char* pszBpolFile, const char* pszMsFile, int nCrossTrack, ReciprocityMemory &RM, HeadMemory &HM, PARAMETER &Parameter)
{
	printf("Reciprocity\n");

	POLYGONATTRIBUTES polyData;
	sense sMag;
	int nLayerMin,nLayerMax;
	int nSizeMag;

	// Load the magnetisation pattern
	Reciprocity_Load_Polygons(&polyData, pszBpolFile, pszEndFile, &nLayerMin, &nLayerMax);

	for (int i = 0; i < polyData.nPolygons; i++)
	{
		if (1 < abs(polyData.pM[i].z))
		{
			RM.bFail = true;
			break;
		}
	}

	const int nLayerSize = nLayerMax - nLayerMin + 1;
	const int nLayerCells = polyData.nPolygons / nLayerSize;

	// Load Ms data for each polygon and scale M in polyData
	Reciprocity_Load_Ms(&polyData, pszMsFile, nLayerCells, Parameter.vecThickness);

	nSizeMag = ((int)(polyData.pSize.x / NM_2_A) + 1 + 1) * ((int)(polyData.pSize.y / NM_2_A) + 1 + 1) * (nLayerMax + 1 + 1);

	// Set size for down-sampled magnetisation array
	Reciprocity_Get_Mag_Array(&sMag, ((int)(polyData.pSize.x / NM_2_A) + 1), ((int)(polyData.pSize.y / NM_2_A) + 1), (nLayerMax + 1));

	Reciprocity_Down_Sample(&polyData,&sMag,nLayerMin,nLayerMax);

	const int nPoint = ((int)polyData.pSize.y + 1) / 10;

	GetMagMemory(nPoint, nCrossTrack, nLayerMax, nSizeMag, RM);

	const int nCrossTrackOffset = RM.nPoint;
	const int nLayerOffset = RM.nCrossTrack * RM.nPoint;

	printf("GPU Start\n");
	int error = ReciprocityGetOutput(RM, HM, nSizeMag, &sMag, nLayerMin, nCrossTrackOffset, nLayerOffset);
	printf("GPU End %d\n", error);

	Reciprocity_Sense_End(&sMag);
	Polygon_IO_End(&polyData);

	if (0 != error)
		RM.bFail = true;

	double dMin = 99999.0;
	double dMax = -99999.0;
	for (int i = 0; i < RM.nSize; i++)
	{
		double d = RM.h_Y[i];
		if (d < dMin)		dMin = d;
		if (d > dMax)		dMax = d;
	}
	if (0.1 > (dMax - dMin))
	{
		RM.bFail = true;
		printf("Max:%f\nMin:%f\n", dMax, dMin);
	}

	// サイモン先生とジミー先生で書き込みの開始方向が違う。
	// ジミー先生に合わせるためWaveformのインデクスを逆転させている。

	// CrossTrack別のWaveform
	RM.vecWholeWaveform.clear();
	for(int j = 0; j < RM.nCrossTrack; j++)
	{
		std::vector<double> v;
		RM.vecWholeWaveform.push_back(v);
		std::vector<double> &wf = RM.vecWholeWaveform[RM.vecWholeWaveform.size() - 1];

		// Save waveform
//		for(int k = 0; k < RM.nPoint ; k++)
		for(int k = (RM.nPoint - 1); k >= 0; k--)
		{
			double d = 0.0;
			for (int i = nLayerMin; i <= RM.nLayerMax; i++)
				d += RM.h_Y[(nCrossTrackOffset * j) + (nLayerOffset * i) + k];

			wf.push_back(d);
		}
	}

	// 層毎のWaveform
	for (int i = nLayerMin; i <= RM.nLayerMax; i++)
	{
		std::vector<std::vector<double>> v1;
		RM.vecWaveform.push_back(v1);
		std::vector<std::vector<double>> &r1 = RM.vecWaveform[RM.vecWaveform.size() - 1];
		r1.reserve(RM.nCrossTrack + 5);

		std::vector<double> v2;
		for (int k = 0; k < RM.nCrossTrack; k++)
		{
			r1.push_back(v2);
			std::vector<double> &r2 = r1[r1.size() - 1];
			r2.reserve(RM.nPoint + 5);
		}

//		for(int k = 0; k < RM.nPoint ; k++)
		for(int k = (RM.nPoint - 1); k >= 0; k--)
		{
			for(int j = 0; j < RM.nCrossTrack; j++)
			{
				r1[j].push_back(RM.h_Y[(nCrossTrackOffset * j) + (nLayerOffset * i) + k]);
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
CSummary::CSummary()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
CSummary::~CSummary()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// stringのメモリリーク対策 (デストラクトを終了前に呼ばせる)
int CSummary::SubMain(const char* pszSimonPath, const char* pszResultPath, const char* pszFolder, const char* pszHeadFile, const char* pszCondition, const char* pszEndFile, const char* pszInputName)
{
	// ヘッドファイルが無ければ。
	if (false == std::filesystem::exists(pszHeadFile))
	{
		printf("Head file not found.\n");
		return 0;
	}

	printf("\n\n        Start tallying.\n\n");

	const std::filesystem::path pathEnd(pszEndFile);

	char szFile[256] = "";
	strcpy(szFile, pathEnd.stem().string().c_str());

	// 拡張子end、bpol、ms、kin、conを探す。
	char szMsFile[256] = "";
	sprintf(szMsFile, "%s/%s%s", pszSimonPath, pszInputName, EXTENSION_MS);

	char szConFile[256] = "";
	sprintf(szConFile, "%s/%s%s", pszSimonPath, pszInputName, EXTENSION_CON);

	char szBpolFile[256] = "";
	sprintf(szBpolFile, "%s/%s%s", pszSimonPath, pszInputName, EXTENSION_BPOL);

	char szLogFile[256] = "";
	sprintf(szLogFile, "%s/%s%s", pszResultPath, szFile, EXTENSION_LOG);

	if (false == std::filesystem::exists(pszEndFile))
	{
		printf("**********************************************************************\n");
		printf("End file not found. %s\n", pszEndFile);
		printf("**********************************************************************\n");
		return 0;
	}

	if (false == std::filesystem::exists(szMsFile))
	{
		printf("**********************************************************************\n");
		printf("Ms file not found. %s\n", szMsFile);
		printf("**********************************************************************\n");
		return 0;
	}

	if (false == std::filesystem::exists(szConFile))
	{
		printf("**********************************************************************\n");
		printf("Con file not found. %s\n", szConFile);
		printf("**********************************************************************\n");
		return 0;
	}

	// 
	DeviceInit();
	ReciprocityMemory RM;

	// 先に生成ヘッドを読み込む。
	HeadMemory HM;
	LoadReadHead(pszHeadFile, HM);

	RM.vecWaveform.clear();

	FILE *pFile = fopen(szConFile, "rt");
	if (nullptr == pFile)
	{
		printf("**********************************************************************\n");
		printf("Failed to file open. File : %s\n", szConFile);
		printf("**********************************************************************\n");
		return 0;
	}

	// 無駄なスペースがあるので
	std::regex r("\\s+");
	std::string strTemp;

	// gpファイル
	char szTemp[256] = "";
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");

	char szGpFile[256] = "";
	if (false == strTemp.empty())
		sprintf(szGpFile, "%s/%s", pszSimonPath, strTemp.c_str());

	// gaファイル
	szTemp[0] = '\0';
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");

	char szGaFile[256] = "";
	if (false == strTemp.empty())
		sprintf(szGaFile, "%s/%s", pszSimonPath, strTemp.c_str());

	// parファイル
	szTemp[0] = '\0';
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");

	char szParFile[256] = "";
	if (false == strTemp.empty())
		sprintf(szParFile, "%s/%s", pszSimonPath, strTemp.c_str());

	// 測定開始時間
	szTemp[0] = '\0';
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");

	char szStartTime[32] = "";
	if (false == strTemp.empty())
		sprintf(szStartTime, "%s/%s/%s %s:%s:%s", strTemp.substr(0, 4).c_str(), strTemp.substr(4, 2).c_str(), strTemp.substr(6, 2).c_str(), strTemp.substr(8, 2).c_str(), strTemp.substr(10, 2).c_str(), strTemp.substr(12, 2).c_str());

	// 測定終了時間
	szTemp[0] = '\0';
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");

	char szEndTime[32] = "";
	if (false == strTemp.empty())
		sprintf(szEndTime, "%s/%s/%s %s:%s:%s", strTemp.substr(0, 4).c_str(), strTemp.substr(4, 2).c_str(), strTemp.substr(6, 2).c_str(), strTemp.substr(8, 2).c_str(), strTemp.substr(10, 2).c_str(), strTemp.substr(12, 2).c_str());

	// nht
	szTemp[0] = '\0';
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");
	const int nht = strtol(strTemp.c_str(), NULL, 10);

	// ntpy
	szTemp[0] = '\0';
	fgets(szTemp, 256, pFile);
	strTemp = szTemp;
	strTemp = std::regex_replace(strTemp, r, "");
	const int ntpy = strtol(strTemp.c_str(), NULL, 10);

	fclose(pFile);
	pFile = nullptr;

	PARAMETER Parameter;
	LoadParameter(szParFile, Parameter);

	// bpol
	// 計算後のbpolは参照するファイルが書かれたファイルとなっている。
	// そのため、本来のバイナリファイルを参照ファイルから作成する。
	int nCrossTrack = 0, nDownTrack = 0;
	if (false == CreateBpol(szGpFile, szBpolFile, Parameter.nLayer, nCrossTrack, nDownTrack))
		return 0;

	// ギリギリだと微妙なので。
	nCrossTrack--;
	nDownTrack--;

	// reciprocity
	// bpol作成後にwaveformのファイルを作成する。
	RunReciprocity(pszEndFile, szBpolFile, szMsFile, nCrossTrack, RM, HM, Parameter);

	// 保存

	// Waveform
	char szWaveformFile[256] = "";
	pFile = nullptr;

	// offset0のwaveform
	// 計算が異常な時は集計されないようにするため作らない。
	sprintf(szWaveformFile, "%s/%s%s", pszResultPath, szFile, EXTENSION_OZ);
	if (false == RM.bFail)
	{
		pFile = fopen(szWaveformFile, "wt");
		if (nullptr != pFile)
		{
			printf("***   %s\n", szWaveformFile);

			std::vector<double> vecOffset0;
			const int nOffset0 = (int)RM.vecWaveform[0].size() / 2;
			for (size_t i = 0; i < RM.vecWaveform.size(); i++)
			{
				std::vector<std::vector<double>> &r1 = RM.vecWaveform[i];

				std::vector<double> &r2 = r1[nOffset0];
				for (size_t k = 0; k < r2.size(); k++)
				{
					// Layer、DownTrack
					fprintf(pFile, "%d,%d,%g\n", (int)(i + 1), (int)k, r2[k]);
					vecOffset0.push_back(r2[k]);
				}
			}

			fclose(pFile);
		}
	}
	pFile = nullptr;

	// 必要
	ReleaseMemoryMag(RM);
	ReleaseMemoryHead(HM);

	DeviceReset();

	printf("\n\n        The tally has been completed.\n\n");

	return 0;
}


