/* reciprocity.c - Calculate output waveforms from Voronoi simulations */
/*  6/ 9/2004 */
/* 30/ 5/2005 - Derived from reciprocity.c */
/* 29/ 8/2006 - Renamed to reciprocity.c, updated, much faster than before */
/*  4/ 2/2008 - Output routine transferred from CTK_R.c, main() routine separated to allow re-usability of this module */
/*  3/ 3/2009 - Can now use *.ms files to scale polygon Ms */

#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "polygon_io.h"
#include "waveform.h"

#include "reciprocity.h"

#include <filesystem>

void Reciprocity_Down_Sample(POLYGONATTRIBUTES *polyData,sense *sMag,int nLayerMin,int nLayerMax)
{			/* Convert M map from Angstroms to nm */
int nIndex,nX,nXY,nY,nZ;
double dX,dY,dZ;

#if 1
//	wprintf(L"  Downsizing polyData\n");
	nXY = sMag->nXMax * sMag->nYMax;
	for(nZ = nLayerMin; nZ <= nLayerMax; nZ++){
		Polygon_IO_Get_Bitmap(polyData,nZ);		/* Digitise layer */

		for(nX = 0; nX < sMag->nXMax; nX++){			/* Co-ordinates in averaged map are (nX,nY) */
			for(nY = 0; nY < sMag->nYMax; nY++){
				nIndex = nX + (nY * sMag->nXMax) + (nZ * nXY);
				Reciprocity_Get_Average_M(polyData,nX,nY,10,&dX,&dY,&dZ);
				sMag->dX[nIndex] = dX;
				sMag->dY[nIndex] = dY;
				sMag->dZ[nIndex] = dZ;
			}
		}
	}

#else

	nXY = sMag->nXMax * sMag->nYMax;
	for(nZ = nLayerMin; nZ <= nLayerMax; nZ++){
		Polygon_IO_Get_Bitmap(polyData,nZ);		/* Digitise layer */

		for(nX = 0; nX < sMag->nXMax; nX++){			/* Co-ordinates in averaged map are (nX,nY) */
			for(nY = 0; nY < sMag->nYMax; nY++){
				nIndex = nX + (nY * sMag->nXMax) + (nZ * nXY);
				Reciprocity_Get_Average_M(polyData,nX,nY,10,&dX,&dY,&dZ);
				sMag->dX[nIndex] = dX;
				sMag->dY[nIndex] = dY;
				sMag->dZ[nIndex] = dZ;
			}
		}
	}

#endif
}


void Reciprocity_Get_Average_M(POLYGONATTRIBUTES *polyData,int nGridX,int nGridY,int nDownSize,double *dMx,double *dMy,double *dMz)
{			/* Get average M in a 1nm square */
int nCount = 0,nX,nY,nXm,nYm;
double dX,dY,dZ;

	nXm = nGridX * nDownSize;		/* Converting from Angstroms to nm */
	nYm = nGridY * nDownSize;
	*dMx = 0.0;
	*dMy = 0.0;
	*dMz = 0.0;
	if(nXm >= 0){				/* Range checks */
		if(nXm + nDownSize < polyData->nPolyXMax){
			if(nYm >= 0){
				if(nYm + nDownSize < polyData->nPolyYMax){
					for(nX = 0; nX < nDownSize; nX++){
						for(nY = 0; nY < nDownSize; nY++){
							Polygon_IO_Get_M_X_Y_Z(polyData,nXm + nX,nYm + nY,&dX,&dY,&dZ);
							*dMx += dX;
							*dMy += dY;
							*dMz += dZ;
							nCount ++;
						}
					}
					*dMx /= (double)nCount;
					*dMy /= (double)nCount;
					*dMz /= (double)nCount;
				}
			}
		}
	}
}


void Reciprocity_Get_Dimensions(const char *cFilename,int *nXMax,int *nYMax, int *nZMax)
{			/* Find the dimensions of the sensitivity file */
FILE *fp;
int nX,nY,nZ;
double dX,dY,dZ;

	*nXMax = 0;
	*nYMax = 0;
	*nZMax = 0;

	// 全部読むのはアホらしいので最後だけ読みます。
	// 最後が一番大きい。
	fp = fopen(cFilename,"r");
	fseek(fp, -200, SEEK_END);

	char sz[100];
	fgets(sz, 100, fp);

	while(fscanf(fp,"%d,%d,%d,%lf,%lf,%lf",&nX,&nY,&nZ,&dX,&dY,&dZ)!=EOF){
		if(nX > *nXMax) *nXMax = nX;
		if(nY > *nYMax) *nYMax = nY;
		if(nZ > *nZMax) *nZMax = nZ;
	}
	fclose(fp);
}


void Reciprocity_Get_Mag_Array(sense *sMag,int nX,int nY,int nZ)
{			/* Get space for the downsized magnetisation distribution */
	Reciprocity_Sense_Setup(&sMag,nX,nY,nZ);
}


void Reciprocity_Load_Ms(POLYGONATTRIBUTES *polyData,const char *cMsFilename, int nLayerCells, std::vector<double> vecThickness)
{			/* Load Ms data for each polygon and scale M in polyData */
int n;
double dMs;
FILE *fp;

	std::uintmax_t size = std::filesystem::file_size(cMsFilename);

	char *psz = new char[size + 1];
	if (nullptr == psz)
		return;

	fp = fopen(cMsFilename, "rb");
	fread(psz, size, 1, fp);
	fclose(fp);

	int nIndex = 0;
	double d = 0.0;
	char *p = strtok(psz, "\r\n");
	for(n = 0; n < polyData->nPolygons; n++){
		sscanf(p,"%lf",&dMs);

		nIndex = (n % nLayerCells) + (n / nLayerCells * nLayerCells);
		d = dMs * vecThickness[n / nLayerCells] * 10e-4;
		polyData->pM[nIndex].x *= d;
		polyData->pM[nIndex].y *= d;
		polyData->pM[nIndex].z *= d;
	}

	delete [] psz;
	psz = nullptr;
}


void Reciprocity_Load_Polygons(POLYGONATTRIBUTES *polyData,const char *cPolyFilename,const char *cEndFilename,int *nLayerMin,int *nLayerMax)
{			/* Load and a set of polygons and magnetisation */
int n;

	Polygon_IO_Load_Pol_Binary(polyData,cPolyFilename);

	*nLayerMax = 0;
	*nLayerMin = 99999;
	for(n = 0; n < polyData->nPolygons; n++){
		if(polyData->nPolyLayer[n] > *nLayerMax) *nLayerMax = polyData->nPolyLayer[n];
		if(polyData->nPolyLayer[n] < *nLayerMin) *nLayerMin = polyData->nPolyLayer[n];
	}

	const int nLayerSize = *nLayerMax - *nLayerMin + 1;
	const int nLayerCells = polyData->nPolygons / nLayerSize;

	Polygon_IO_Load_Magnetisation(cEndFilename,polyData, nLayerCells);
}


int Reciprocity_Load_Sense_Distribution(sense *sSense,const char *cFilename, int nMagZMax)
{			/* Load a head sensitivity function */
int nIndex,nX,nY,nZ,nXMax,nYMax,nZMax;
int nSize;
double dX,dY,dZ;
FILE *fp;

	Reciprocity_Get_Dimensions(cFilename,&nXMax,&nYMax,&nZMax);
	if(nMagZMax > nZMax) nZMax = nMagZMax;		/* If the mag file has more layers we need extra space to avoid errors */
	nSize = Reciprocity_Sense_Setup(&sSense,nXMax,nYMax,nZMax);

	std::uintmax_t size = std::filesystem::file_size(cFilename);

	char *psz = new char[size + 1];
	if (nullptr == psz)
		return 0;

	fp = fopen(cFilename, "rb");
	fread(psz, size, 1, fp);
	fclose(fp);

	char *p = strtok(psz, "\r\n");
	while (NULL != p)
	{
		sscanf(p, "%d,%d,%d,%lf,%lf,%lf", &nX, &nY, &nZ, &dX, &dY, &dZ);
		nIndex = nX + (nXMax * nY) + (nXMax * nYMax * nZ);
		sSense->dX[nIndex] = dX;
		sSense->dY[nIndex] = dY;
		sSense->dZ[nIndex] = dZ;

		p = strtok(NULL, "\r\n");
	}

	delete [] psz;
	psz = nullptr;

	return nSize;
}


int Reciprocity_Sense_End(sense *sSense)
{			/* Free memory allocated to sensitivity function */
int nOK=0;

	if(sSense->nXMax > 0 && sSense->nYMax > 0 && sSense->nZMax > 0){
		delete [] sSense->dX;
		delete [] sSense->dY;
		delete [] sSense->dZ;
		nOK = 1;
	}

	return (nOK);
}


int Reciprocity_Sense_Setup(sense **sSense,int nXMax,int nYMax,int nZMax)
{			/* Allocate memory for a sensitivity function */
int nOK = 0, nSize = (nXMax + 1) * (nYMax + 1) * (nZMax + 1);		/* Add 1 because indices can run from 0 to nXMax etc */

	if(nXMax > 0 && nYMax > 0 && nZMax > 0){
		(*sSense)->nXMax = nXMax;
		(*sSense)->nYMax = nYMax;
		(*sSense)->nZMax = nZMax;
		(*sSense)->dX = new double[nSize];
		(*sSense)->dY = new double[nSize];
		(*sSense)->dZ = new double[nSize];
		memset((*sSense)->dX, 0, (sizeof(double) * nSize));
		memset((*sSense)->dY, 0, (sizeof(double) * nSize));
		memset((*sSense)->dZ, 0, (sizeof(double) * nSize));
		nOK = 1;
	}

	return (nSize);
}


