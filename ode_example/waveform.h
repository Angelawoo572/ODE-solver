/* waveform.h */


#ifndef __WAVEFORM_H__
#define __WAVEFORM_H__

#include "polygon_io.h"

typedef struct
{
	double *dY;
	int nPoints;
} vector;


int			Waveform_End(vector *v);							/* Free a vector */
int			Waveform_Init(vector *v,int nSize);		/* Setup a vector */
void		Waveform_Save_Log(const char *cLogFile,double *d,int nBitlength,int nYMin,int nYMax,int nX,const char *pszStartTime, const char *pszEndTime, const char *pszCondition, int nBaseBit);		/* Save the bit attributes to a log file */
void		Waveform_Save_Vector(vector *v,const char *cFilename);		/* Save a vector to a file */
void		Waveform_Setup_Vector(vector **v);		/* Get memory for a vector */


#endif

