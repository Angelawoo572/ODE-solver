/* waveform.c - Calculate waveforms from Voronoi simulations */
/*  9/12/2003 */
/* 15/ 4/2005 - Added Get_DC_Noise() routine */
/* 26/ 5/2005 - No longer saves area under bit to Waveform.log file */
/* 16/ 1/2006 - Added Get_T50() */
/* 29/ 8/2006 - New version without dependence on digitise_3d.c */
/* 26/10/2006 - Fixed saving error for .bit and .noi files */

#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "polygon_io.h"

#include "waveform.h"


int Waveform_End(vector *v)
{			/* Free a vector */
int nOK=0;

	if(v->nPoints != 0){
		free(v->dY);
		nOK = 1;
	}
	return (nOK);
}


int Waveform_Init(vector *v,int nSize)
{			/* Instantiate a vector */
int n,nOK = 0;

	if(nSize > 0){
		v->nPoints = nSize;
		Waveform_Setup_Vector(&v);
		for(n = 0; n < v->nPoints; n++){			/* Initialise to zero */
			v->dY[n] = 0.0;
		}
		nOK = 1;
	}

	return (nOK);
}



void Waveform_Save_Log(const char *cLogFile,double *d,int nBitlength,int nYMin,int nYMax,int nX,const char *pszStartTime,const char *pszEndTime, const char *pszCondition, int nBaseBit)
{			/* Save the bit attributes to a log file */
			/* Signal, noise, saturation, jitter, T50, x min, x max, y min, y max, bitlength, layer, grain size, grain size sigma */
FILE *fp;

	fp = fopen(cLogFile,"a");
//	fwprintf(fp,L"%s,%s,%s,%d,%.4f,%.5f,%.2f,%.4f,%.3f,%d,%d,%d,%d,%d,%d,%.3f,%.3f,%.3f\n",pszStartTime,pszEndTime,pszCondition,nX,d[1],d[2],100.0 * d[0] / (double)(2 * nBitlength),d[3],d[4],nBitlength,nXMin,nXMax,nYMin,nYMax,nLayer,d[5],d[6],20.0 * log10(d[1] / d[2]));
	fprintf(fp,"%s,%s,%s,%d,%d,%d,%d,%f,%f,%f,%d,%f,%f,%f\n",pszStartTime,pszEndTime,pszCondition,nX,nYMin,nYMax,nBitlength,d[0],d[1],20.0 * log10(d[0] / d[1]),nBaseBit,d[2],d[3],20.0 * log10(d[2] / d[3]));
	fclose(fp);
}


void Waveform_Save_Vector(vector *v,const char *cFilename)
{			/* Save a vector */
FILE *fp;
int n;

	fp = fopen(cFilename,"w");
	for(n = 0; n < v->nPoints ; n++){
		fprintf(fp,"%d %lf\n",n,v->dY[n]);
	}
	fclose(fp);
}


void Waveform_Setup_Vector(vector **v)
{			/* Get memory for a vector */
	(*v)->dY = (double*)calloc((*v)->nPoints,sizeof(double));		/* y-axis data (double) */
}



