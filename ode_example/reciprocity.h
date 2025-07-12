/* reciprocity.h */

#ifndef __RECIPROCITY_H__
#define __RECIPROCITY_H__

#include "polygon_io.h"
#include "waveform.h"

#include <vector>

typedef struct{
	int nXMax;
	int nYMax;
	int nZMax;
	double *dX;
	double *dY;
	double *dZ;
	} sense;



void	Reciprocity_Down_Sample(POLYGONATTRIBUTES *polyData,sense *sMag,int nLayerMin,int nLayerMax);		/* Sample and average regions of polyData */
void	Reciprocity_Get_Average_M(POLYGONATTRIBUTES *polyData,int nGridX,int nGridY,int nDownSize,double *dMx,double *dMy,double *dMz);		/* Get average M over a 1nm square */
void	Reciprocity_Get_Dimensions(const char *cFilename,int *nXMax,int *nYMax,int *nZMax);			/* Find the dimensions of the sensitivity file */
void	Reciprocity_Get_Mag_Array(sense *sMag,int nX,int nY,int nZ);		/* Get space for the downsized magnetisation distribution */
void	Reciprocity_Load_Ms(POLYGONATTRIBUTES *polyData,const char *cMsFilename, int nLayerCells, std::vector<double> vecThickness);			/* Load Ms data for each polygon and scale M in polyData */
void	Reciprocity_Load_Polygons(POLYGONATTRIBUTES *polyData,const char *cPolyFilename,const char *cEndFilename,int *nLayerMin,int *nLayerMax);			/* Load and a set of polygons and magnetisation */
int		Reciprocity_Load_Sense_Distribution(sense *sSense,const char *cFilename,int nMagZMax);
int		Reciprocity_Sense_End(sense *sSense);
int		Reciprocity_Sense_Setup(sense **sSense,int nXMax,int nYMax,int nZMax);






#endif

