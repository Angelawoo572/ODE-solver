/* Polygon_io.c */
/* Read / Write routines for polygon format files */
/* 18/ 4/2005 */
/* 12/ 7/2005 - Added BYTEORDER option to .h file and switches in the program to ease portability between machines */
/* 12/ 7/2005 - Extended POLYGONATTRIBUTES so that all polygon data can be stored in it */
/* 16/11/2005 - Added routine to copy layers of polygons */
/* 27/ 3/2006 - New Copy_Polygon() routine plus Process_Polygons() function added */
/*  6/ 4/2006 - Added functions to create, merge and translate polygon data */
/* 21/ 4/2006 - Added Is_Cubic() function to determine whether a cell is a cuboid */
/* 17/ 5/2006 - Added Get_Polygon_Array_Circle_Smooth() to make smooth circular objects */
/* 11/ 8/2006 - Added Bitmap function to replace digitise_3d file */
/*  5/10/2006 - Removed BYTEORDER and added Get_Byte_Order() function */
/*  6/10/2006 - Added fix to Convert_Pol_To_P3D() to remove zero-area polygons */
/* 15/12/2006 - Added functions to merge polygon sets, crop polygon sets and display statistics of polygon sets */
/* 31/ 1/2007 - Removed Polygon_IO_Add_Layer() as it doesn't work when compiled with sxcc */
/*  9/ 5/2007 - Updated functions */
/*  1/ 2/2008 - Added ellipse objects with sub-cell discretisation */
/*  2/ 7/2008 - DTM mesh update for improved flexibility */
/* 16/ 3/2010 - Improved sub-cell discretisation of ellipses with smoother edges */
/*  7/12/2010 - Added boundary cell allocation routine Add_Boundary() */
/* 16/ 3/2012 - Added Spin Ice mesh construction routines */


#define _CRT_SECURE_NO_WARNINGS

#include "point.h"

#include "polygon_io.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <filesystem>

#define SAFE_DELETE(x)			{if (NULL != x) {delete x;} x = NULL;}
#define SAFE_DELETES(x)			{if (NULL != x) {delete [] x;} x = NULL;}
#define SAFE_DELETE_CLEAR(x)	{if (NULL != x) {x->Clear(); delete x;} x = NULL;}
#define SAFE_RELEASE(x)			{if (NULL != x) {x->Release();} x = NULL;}

#define NEW_NONE(a, b, c)		{int FWE = c; if (0 == FWE) FWE++; a = new b[FWE];}
#define NEW_ZERO(a, b, c)		{int FWE = c; if (0 == FWE) FWE++; a = new b[FWE]; memset(a, 0, (sizeof(b) * FWE));}


void Polygon_IO_Add_Boundary(POLYGONATTRIBUTES *polyData,int nLayerFrom,int nLayerTo,double dBoundary)
{			/* Add boundary cells to polyData */
POLYGONATTRIBUTES polyTemp,polyTemp2;
int n,nLayer,nLayers;
double dLongest,dX,dY;

	Polygon_IO_Reset_Polygon_Attributes(&polyTemp);
	Polygon_IO_Reset_Polygon_Attributes(&polyTemp2);

	if(nLayerTo < nLayerFrom){
		n = nLayerFrom;
		nLayerFrom = nLayerTo;
		nLayerTo = n;
	}
	nLayers = nLayerTo - nLayerFrom + 1;
	polyTemp.nPolygons = (4 * nLayers);			/* Set the header details */
	polyTemp.nVertexMax = 4;
	Polygon_IO_Setup_Polygon_Attributes_From_P(&polyTemp,0,1,0,0);
	dX = polyData->pSize.x;
	dY = polyData->pSize.y;

	polyTemp.pSize.x = dX + (2.0 * dBoundary);
	polyTemp.pSize.y = dY + (2.0 * dBoundary);

	if(polyData->pSize.x > polyData->pSize.y) dLongest = polyData->pSize.x;
	else dLongest = polyData->pSize.y;
	polyTemp.dSeedToVertexMax2 = 0.25 * ((dBoundary * dBoundary) + (dLongest * dLongest));

	n = 0;
	for(nLayer = nLayerFrom; nLayer <= nLayerTo; nLayer++){
		polyTemp.nPolyLayer[n] = nLayer;
		polyTemp.siPolyVertices[n] = 4;
		polyTemp.fBoundary[n] = 0.0;
		polyTemp.pSeed[n].x = (dX + dBoundary) * 0.5;		/* Seed point is centre of polygon */
		polyTemp.pSeed[n].y = dBoundary * 0.5;

		polyTemp.pVertex[n].x = 0.0;				/* Set vertices (lower left cell) */
		polyTemp.pVertex[n].y = 0.0;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].x = 0.0;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].y = dBoundary;
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].x = dX + dBoundary;
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].y = dBoundary;
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].x = dX + dBoundary;
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].y = 0.0;

		n++;
		polyTemp.nPolyLayer[n] = nLayer;
		polyTemp.siPolyVertices[n] = 4;
		polyTemp.fBoundary[n] = 0.0;
		polyTemp.pSeed[n].x = dBoundary * 0.5;		/* Seed point is centre of polygon */
		polyTemp.pSeed[n].y = dBoundary + ((dY + dBoundary) * 0.5);

		polyTemp.pVertex[n].x = 0.0;				/* Set vertices (upper left cell) */
		polyTemp.pVertex[n].y = dBoundary;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].x = 0.0;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].y = dY + (2.0 * dBoundary);
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].x = dBoundary;
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].y = dY + (2.0 * dBoundary);
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].x = dBoundary;
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].y = dBoundary;

		n++;
		polyTemp.nPolyLayer[n] = nLayer;
		polyTemp.siPolyVertices[n] = 4;
		polyTemp.fBoundary[n] = 0.0;
		polyTemp.pSeed[n].x = dBoundary + ((dX + dBoundary) * 0.5);		/* Seed point is centre of polygon */
		polyTemp.pSeed[n].y = dY + (1.5 * dBoundary);

		polyTemp.pVertex[n].x = dBoundary;				/* Set vertices (upper right cell) */
		polyTemp.pVertex[n].y = dY + dBoundary;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].x = dBoundary;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].y = dY + (2.0 * dBoundary);
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].x = dX + (2.0 * dBoundary);
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].y = dY + (2.0 * dBoundary);
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].x = dX + (2.0 * dBoundary);
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].y = dY + dBoundary;

		n++;
		polyTemp.nPolyLayer[n] = nLayer;
		polyTemp.siPolyVertices[n] = 4;
		polyTemp.fBoundary[n] = 0.0;
		polyTemp.pSeed[n].x = dX + (1.5 * dBoundary);		/* Seed point is centre of polygon */
		polyTemp.pSeed[n].y = (dY + dBoundary) * 0.5;

		polyTemp.pVertex[n].x = dX + dBoundary;				/* Set vertices (lower right cell) */
		polyTemp.pVertex[n].y = 0.0;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].x = dX + dBoundary;
		polyTemp.pVertex[n + (polyTemp.nPolygons)].y = dY + dBoundary;
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].x = dX + (2.0 * dBoundary);
		polyTemp.pVertex[n + (2 * polyTemp.nPolygons)].y = dY + dBoundary;
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].x = dX + (2.0 * dBoundary);
		polyTemp.pVertex[n + (3 * polyTemp.nPolygons)].y = 0.0;
	}

	Polygon_IO_Copy_All_Polygons(polyData,&polyTemp2);		/* Copy polyData to polyTemp2 */
	Polygon_IO_Translate_Polygons(&polyTemp2,dBoundary,dBoundary);		/* Offset existing polygons */

	Polygon_IO_Merge_Polygons(&polyTemp2,&polyTemp,polyData);		/* Merge polyTemp2 and polyTemp */
	Polygon_IO_End(&polyTemp);						/* Free polyTemp */
	Polygon_IO_End(&polyTemp2);
	Polygon_IO_Update_Header(polyData);		/* Update the header */
}


void Polygon_IO_Check_Vertex_Ordering(POLYGONATTRIBUTES *polyData)
{			/* Make vertex chirality the same for all polygons */
int nPolygon,nBaseVertex;
point p1,p2,pZ;

	p1.z = 0.0;		/* Vertices are all in the same plane */
	p2.z = 0.0;
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		nBaseVertex = 0;
		do{
			Polygon_IO_Get_Vectors_From_Vertices(polyData,nPolygon,nBaseVertex,&p1,&p2);
			Point_Cross_Product(&pZ,&p1,&p2);
			nBaseVertex++;
		}
		while(pZ.z == 0.0);
		if(pZ.z < 0.0) Polygon_IO_Reverse_Vertex_Order(polyData,nPolygon);
	}
}


void Polygon_IO_Clone_Cells(POLYGONATTRIBUTES *polyData,int nSourceLayer,int nTargetLayer)
{			/* Clone cells from one layer to another */
POLYGONATTRIBUTES polyTemp,polyTemp2;

	Polygon_IO_Reset_Polygon_Attributes(&polyTemp);
	Polygon_IO_Reset_Polygon_Attributes(&polyTemp2);
	Polygon_IO_Copy_Polygon_Layer(polyData,&polyTemp,nSourceLayer,nTargetLayer);		/* Copy polygons from nSourceLayer into polyTemp */
	Polygon_IO_Copy_All_Polygons(polyData,&polyTemp2);
	Polygon_IO_Merge_Polygons(&polyTemp2,&polyTemp,polyData);
	Polygon_IO_End(&polyTemp);						/* Free polyTemp */
	Polygon_IO_End(&polyTemp2);
}


void Polygon_IO_Convert_Pol_To_P3D(char *cPolFilename,char *cP3DFilename,POLYGONATTRIBUTES *polyData,int nL)
{			/* Convert a .pol file to the new .p3d file format */
int n;
POLYGONATTRIBUTES polyNew;

	Polygon_IO_Reset_Polygon_Attributes(&polyNew);
	Polygon_IO_Get_File_Attributes(cPolFilename,polyData);
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	Polygon_IO_Load_Pol_Polygons_Old(cPolFilename,polyData->nPolygons,polyData->pSeed,polyData->pVertex);
	polyData->nVertexMax--;		/* Old format repeated the first vertex */
	for(n = 0; n < polyData->nPolygons; n++){
		polyData->siPolyVertices[n] = Polygon_IO_Get_Vertex_Count(n,polyData->nPolygons,polyData->pVertex);		/* Fill siVertices array */
		polyData->nPolyLayer[n] = nL;		/* Set layer to nL */
	}
	Polygon_IO_Remove_Zero_Area_Polygons(polyData,&polyNew);
	Polygon_IO_Save_P3D(cP3DFilename,&polyNew);
}


int Polygon_IO_Copy_All_Polygons(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyCopy)
{			/* Copy polygons from polySource to polyCopy */
int nOK = 0;

	Polygon_IO_End(polyCopy);
	Polygon_IO_Copy_Header(polySource,polyCopy);			/* First copy header */
	Polygon_IO_Setup_Polygon_Attributes(&polyCopy,polySource->nGroupSet,1,0,0);	/* Then get memory for copy */
	Polygon_IO_Copy_Polygons(polySource,polyCopy,0);	/* Copy all the polygons */
	Polygon_IO_Update_Header(polyCopy);		/* Update the header */

	return (nOK);
}


void Polygon_IO_Copy_Header(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew)
{			/* Copy a header from one data set to another */
	polyNew->nPolygons = polySource->nPolygons;			/* Total number of polygons */
	polyNew->nVertexMax = polySource->nVertexMax;		/* Maximum number of vertices for a single polygon */
	polyNew->pSize.x = polySource->pSize.x;			/* Maximum extent of polygon region (A) */
	polyNew->pSize.y = polySource->pSize.y;			/* Maximum extent of polygon region (A) */
	polyNew->dSeedToVertexMax2 = polySource->dSeedToVertexMax2;		/* Maximum distance from seed point to vertex (squared) */
}


void Polygon_IO_Copy_Polygon(POLYGONATTRIBUTES *polyFrom,int nPolyFrom,POLYGONATTRIBUTES *polyTo,int nPolyTo)
{			/* Copy a polygon from polyFrom to polyTo */
int nVertex;

	if(nPolyFrom < polyFrom->nPolygons && nPolyTo < polyTo->nPolygons){
		polyTo->nPolyLayer[nPolyTo] = polyFrom->nPolyLayer[nPolyFrom];		/* Copy layer assignment */
		polyTo->siPolyVertices[nPolyTo] = polyFrom->siPolyVertices[nPolyFrom];	/* Copy number of vertices */
		polyTo->pSeed[nPolyTo] = polyFrom->pSeed[nPolyFrom];		/* Copy seed point */
		polyTo->fBoundary[nPolyTo] = polyFrom->fBoundary[nPolyFrom];		/* Non-magnetic boundary */
		if(polyFrom->nGroupSet == 1) polyTo->nGroup[nPolyTo] = polyFrom->nGroup[nPolyFrom];		/* Copy group index (if available) */
		for(nVertex = 0; nVertex < polyTo->siPolyVertices[nPolyTo]; nVertex++){		/* Copy vertices */
			polyTo->pVertex[nPolyTo + (nVertex * polyTo->nPolygons)] = polyFrom->pVertex[nPolyFrom + (nVertex * polyFrom->nPolygons)];
		}
	}
}


void Polygon_IO_Copy_Polygon_Layer(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nLayerFrom,int nLayerTo)
{			/* Copy polygons from one layer in polySource to a new polygon object */
int n,nT,nCopyPolygons = 0,nVertexMax = 0;

	for(n = 0; n < polySource->nPolygons; n++){
		if(polySource->nPolyLayer[n] == nLayerFrom){
			nCopyPolygons++;
			if(polySource->siPolyVertices[n] > nVertexMax) nVertexMax = polySource->siPolyVertices[n];
		}
	}
	printf("    %d polygons to be copied from layer %d to layer %d\n",nCopyPolygons,nLayerFrom,nLayerTo);

	Polygon_IO_End(polyNew);
	polyNew->nPolygons = nCopyPolygons;		/* Number of polygons in new object */
	polyNew->nVertexMax = nVertexMax;
	Polygon_IO_Setup_Polygon_Attributes(&polyNew,polySource->nGroupSet,1,0,0);		/* Space for vertex data */

	nT=0;
	for(n=0;n<polySource->nPolygons;n++){
		if(polySource->nPolyLayer[n] == nLayerFrom){
			Polygon_IO_Copy_Polygon(polySource,n,polyNew,nT);		/* Copy polygon */
			polyNew->nPolyLayer[nT] = nLayerTo;				/* Store new layer */
			nT++;
		}
	}

	Polygon_IO_Update_Header(polyNew);		/* Update size, seed to vertex max etc */
}


int Polygon_IO_Copy_Polygon_Layer_As_Group_Core(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nLayerFrom,int nLayerTo)
{			/* Copy polygons from nLayerFrom in polySource into nLayerTo of polyNew */
int nLastGroup = -1,nPolygon,nSource,nV;

	if(polySource->nGroupSet == 0){		/* No groups were set, copy cells directly */
		Polygon_IO_Copy_Polygon_Layer(polySource,polyNew,nLayerFrom,nLayerTo);
	}
	else{			/* Need to find boundaries of groups */
		Polygon_IO_End(polyNew);
		polyNew->nPolygons = Polygon_IO_Get_Group_Count_In_Layer(polySource,nLayerFrom);
		polyNew->nVertexMax = Polygon_IO_Get_Max_Vertex_Count_In_Layer(polySource,nLayerFrom);
		printf("    Group layer %d: %d group cells with a maximum of %d vertices\n",nLayerTo,polyNew->nPolygons,polyNew->nVertexMax);
		Polygon_IO_Setup_Polygon_Attributes(&polyNew,1,1,0,0);		/* Space for vertex data */

		nPolygon = -1;
		for(nSource = 0; nSource < polySource->nPolygons; nSource++){
			if(polySource->nPolyLayer[nSource] == nLayerFrom){
				if(polySource->nGroup[nSource] != nLastGroup){		/* True for first cell in group */
					nLastGroup = polySource->nGroup[nSource];
					nV = 0;
					nPolygon++;
					polyNew->siPolyVertices[nPolygon] = polySource->siPolyVertices[nSource];
					polyNew->nPolyLayer[nPolygon] = nLayerTo;
					polyNew->pSeed[nPolygon] = polySource->pSeed[nSource];
					polyNew->nGroup[nPolygon] = polySource->nGroup[nSource];
				}
				else{		/* True for second and subsequent cells in a group */
					polyNew->pVertex[nPolygon + (nV * polyNew->nPolygons)] = polySource->pVertex[nSource];
					nV++;
				}
			}
		}
		Polygon_IO_Update_Header(polyNew);		/* Update size, seed to vertex max etc */
	}

	return (nPolygon);
}



int Polygon_IO_Copy_Polygon_Layer_As_Group_XY(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nLayerFrom,int nLayerTo)
{			/* Copy polygons from nLayerFrom in polySource into nLayerTo of polyNew */
int nLastGroup = -1,nPolygon,nSource;
double dX1,dX2,dXMin,dXMax,dY1,dY2,dYMin,dYMax;

	if(polySource->nGroupSet == 0){		/* No groups were set, copy cells directly */
		Polygon_IO_Copy_Polygon_Layer(polySource,polyNew,nLayerFrom,nLayerTo);
	}
	else{			/* Need to find boundaries of groups */
		Polygon_IO_End(polyNew);
		polyNew->nPolygons = Polygon_IO_Get_Group_Count_In_Layer(polySource,nLayerFrom);
		polyNew->nVertexMax = Polygon_IO_Get_Max_Vertex_Count_In_Layer(polySource,nLayerFrom);
		printf("    Group layer %d: %d group cells with a maximum of %d vertices\n",nLayerTo,polyNew->nPolygons,polyNew->nVertexMax);
		Polygon_IO_Setup_Polygon_Attributes(&polyNew,1,1,0,0);		/* Space for vertex data */

		nPolygon = -1;
		for(nSource = 0; nSource < polySource->nPolygons; nSource++){
			if(polySource->nPolyLayer[nSource] == nLayerFrom){
				if(polySource->nGroup[nSource] != nLastGroup){		/* True for first cell in group */
					nLastGroup = polySource->nGroup[nSource];
					nPolygon++;
					polyNew->siPolyVertices[nPolygon] = 4;
					polyNew->nPolyLayer[nPolygon] = nLayerTo;
					polyNew->nGroup[nPolygon] = polySource->nGroup[nSource];
					polyNew->fBoundary[nPolygon] = 0.0;			/* No boundary for cells */
					Polygon_IO_Get_Polygon_Limits(polySource,nSource,&dXMin,&dXMax,&dYMin,&dYMax);

					polyNew->pSeed[nPolygon].x = 0.5 * (dXMin + dXMax);
					polyNew->pSeed[nPolygon].y = 0.5 * (dYMin + dYMax);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,0,dXMin,dYMin);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,1,dXMin,dYMax);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,2,dXMax,dYMax);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,3,dXMax,dYMin);
				}
				else{		/* True for second and subsequent cells in a group */
					Polygon_IO_Get_Polygon_Limits(polySource,nSource,&dX1,&dX2,&dY1,&dY2);
					if(dX1 < dXMin) dXMin = dX1;
					if(dX2 > dXMax) dXMax = dX2;
					if(dY1 < dYMin) dYMin = dY1;
					if(dY2 > dYMax) dYMax = dY2;

					polyNew->pSeed[nPolygon].x = 0.5 * (dXMin + dXMax);
					polyNew->pSeed[nPolygon].y = 0.5 * (dYMin + dYMax);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,0,dXMin,dYMin);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,1,dXMin,dYMax);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,2,dXMax,dYMax);
					Polygon_IO_Set_Vertex(polyNew,nPolygon,3,dXMax,dYMin);
				}
			}
		}
		Polygon_IO_Update_Header(polyNew);		/* Update size, seed to vertex max etc */
	}

	return (nPolygon);
}


void Polygon_IO_Copy_Polygons(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nOffset)
{			/* Copy polygons from one data set to another. nOffset changes the polygon indices */
int nPolygon;

	for(nPolygon = 0; nPolygon < polySource->nPolygons; nPolygon++){
		Polygon_IO_Copy_Polygon(polySource,nPolygon,polyNew,nPolygon+nOffset);
	}
}


int Polygon_IO_Crop_Polygon_Region(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut,double dXStart,double dXStop,double dYStart,double dYStop)
{			/* Crop a polygon region */
int n,nCount = 0,nNewPolygon = 0,nVertexMax = 0;
double dTolerance = 1.0;		/* 1A tolerance on in/out decision */

	for(n = 0; n < polyIn->nPolygons; n++){
		if(Polygon_IO_Is_Polygon_In_Region(polyIn,n,dXStart,dXStop,dYStart,dYStop,dTolerance) == 1){
			nCount++;
			if(polyIn->siPolyVertices[n] > nVertexMax) nVertexMax = polyIn->siPolyVertices[n];
		}
	}

	printf("    Crop polygons : Region = (%.3f,%.3f) - (%.3f,%.3f)\n",dXStart,dYStart,dXStop,dYStop);
	printf("    Crop polygons : Old count %d, New count = %d\n",polyIn->nPolygons,nCount);
	polyOut->nPolygons = nCount;
	polyOut->nVertexMax = nVertexMax;
	Polygon_IO_Setup_Polygon_Attributes(&polyOut,0,1,0,0);

	for(n = 0; n < polyIn->nPolygons; n++){		/* Copy polygons from polyIn */
		if(Polygon_IO_Is_Polygon_In_Region(polyIn,n,dXStart,dXStop,dYStart,dYStop,dTolerance) == 1){
			Polygon_IO_Copy_Polygon(polyIn,n,polyOut,nNewPolygon);
			nNewPolygon++;
		}
	}
	Polygon_IO_Update_Header(polyOut);

	return (nNewPolygon);
}


void Polygon_IO_End(POLYGONATTRIBUTES *polyData)
{			/* Free memory in polyData */
	if(polyData->nPolygons != 0){
		SAFE_DELETES(polyData->pCentre);
		SAFE_DELETES(polyData->pSeed);
		SAFE_DELETES(polyData->pVertex);
		SAFE_DELETES(polyData->siPolyVertices);
		SAFE_DELETES(polyData->nPolyLayer);
		if(polyData->nBitmapSet == 1){
			SAFE_DELETES(polyData->nBitmap);
			polyData->nBitmapSet = 0;
		}

		if(polyData->nBoundarySet == 1){
			SAFE_DELETES(polyData->fBoundary);
			polyData->nBoundarySet = 0;
		}

		Polygon_IO_End_Group(polyData);

		if(polyData->nMagnetisationSet == 1){
			SAFE_DELETES(polyData->pM);
			polyData->nMagnetisationSet = 0;
		}
		Polygon_IO_Reset_Polygon_Attributes(polyData);
	}
}


int Polygon_IO_End_Group(POLYGONATTRIBUTES *polyData)
{			/* Delete support for group indexing in polyData */
int nOK = 0;

	if(polyData->nGroupSet == 1){
		SAFE_DELETES(polyData->nGroup);
		polyData->nGroupSet = 0;
		nOK = 1;
	}

	return (nOK);
}


void Polygon_IO_Get_Aspect_Ratio(POLYGONATTRIBUTES *polyData)
{			/* Calculate the aspect ratio of a set of polygons */
int nPolygon,nV;
double dX,dXMax,dXMin,dY,dYMax,dYMin;
point p;

	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		p = Polygon_IO_Get_Vertex(polyData,nPolygon,0);
		dXMin = p.x;
		dXMax = p.x;
		dYMin = p.y;
		dYMax = p.y;
		for(nV = 1; nV < polyData->siPolyVertices[nPolygon]; nV++){
			p = Polygon_IO_Get_Vertex(polyData,nPolygon,nV);
			if(p.x < dXMin) dXMin = p.x;
			if(p.x > dXMax) dXMax = p.x;
			if(p.y < dYMin) dYMin = p.y;
			if(p.y > dYMax) dYMax = p.y;
		}
		dX = dXMax - dXMin;
		dY = dYMax - dYMin;

		printf("%d %lf %lf\n",nPolygon,dX/dY,dY/dX);
	}
}


double Polygon_IO_Get_Average_Polygon_Area(POLYGONATTRIBUTES *polyData)
{			/* Return the average area of polygons in polyData */
	return (Polygon_IO_Get_Average_Polygon_Area_Region(polyData,0.0,polyData->pSize.x,0.0,polyData->pSize.y));
}


double	Polygon_IO_Get_Average_Polygon_Area_Region(POLYGONATTRIBUTES *polyData,double dXStart,double dXStop,double dYStart,double dYStop)
{			/* Return the average area of polygons in a region */
int nPolygon,nPolygons=0;
double dArea,dTotalArea=0.0;

	if(polyData->nPolygons == 0) return (0.0);

	for(nPolygon=0;nPolygon < polyData->nPolygons; nPolygon++){
		if(polyData->pSeed[nPolygon].x > dXStart){
			if(polyData->pSeed[nPolygon].x < dXStop){
				if(polyData->pSeed[nPolygon].y > dYStart){
					if(polyData->pSeed[nPolygon].y < dYStop){
						dArea = Polygon_IO_Get_Polygon_Area(polyData,nPolygon);
						if(dArea < 0.0) dArea *= -1.0;
						dTotalArea += dArea;
						nPolygons++;
					}
				}
			}
		}
	}

	if(nPolygons > 0){
		dArea = dTotalArea/ (double)nPolygons;
	}
	else{
		dArea = 0.0;
	}

	return (dArea);
}


void Polygon_IO_Get_Bitmap(POLYGONATTRIBUTES *polyData, int nLayer)
{			/* Make a bitmap from the polygon data */
int n2,nPolygon,nPolygons,nVertex,nY;
int nXMin,nXMax,nYMin,nYMax;
double dA,dB,dC,dX,dY,dXMin,dXMax,dYMin,dYMax;
point p0,p1,p2;
long lIndex,lPointCount = 0,lX;

	Polygon_IO_Setup_Polygon_Attributes_Bitmap(&polyData);		/* Create space for the bitmap */
	Polygon_IO_Reset_Bitmap(polyData);

//	printf("    Digitising -");
	nPolygons = polyData->nPolygons;
	for(nPolygon = 0; nPolygon < nPolygons; nPolygon++){
		if(polyData->nPolyLayer[nPolygon] == nLayer){		/* If polygon is in the desired layer ... */
			p0 = polyData->pSeed[nPolygon];
			n2 = polyData->siPolyVertices[nPolygon] - 1;			/* Last vertex in this polygon */
			for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
				p1 = Polygon_IO_Get_Vertex(polyData,nPolygon,nVertex);
				p2 = Polygon_IO_Get_Vertex(polyData,nPolygon,n2);
//				p1=polyData->pVertex[nPolygon+(nVertex*nPolygons)];
//				p2=polyData->pVertex[nPolygon+(n2*nPolygons)];
				dXMin = p0.x;
				if(p1.x < dXMin) dXMin = p1.x;
				if(p2.x < dXMin) dXMin = p2.x;
				dYMin = p0.y;
				if(p1.y < dYMin) dYMin = p1.y;
				if(p2.y < dYMin) dYMin = p2.y;
				dXMax = p0.x;
				if(p1.x > dXMax) dXMax = p1.x;
				if(p2.x > dXMax) dXMax = p2.x;
				dYMax = p0.y;
				if(p1.y > dYMax) dYMax = p1.y;
				if(p2.y > dYMax) dYMax = p2.y;
				if(dXMin < 0.0) dXMin = 0.0;
				if(dYMin < 0.0) dYMin = 0.0;
				if(dXMax >= polyData->pSize.x) dXMax = polyData->pSize.x;
				if(dYMax >= polyData->pSize.y) dYMax = polyData->pSize.y;
				nXMin = (int)rint(dXMin);
				nXMax = (int)rint(dXMax);
				nYMin = (int)rint(dYMin);
				nYMax = (int)rint(dYMax);
				for(nY = nYMin; nY <= nYMax; nY++){
					dY = (double)nY;
					lIndex = (nY * polyData->nPolyXMax);
					for(lX = nXMin; lX <= nXMax; lX++){
						dX = (double)lX;
						dA = (p2.x - p1.x) * (dY - p1.y) - (p2.y - p1.y) * (dX - p1.x);
						dB = (p0.x - p2.x) * (dY - p2.y) - (p0.y - p2.y) * (dX - p2.x);
						dC = (p1.x - p0.x) * (dY - p0.y) - (p1.y - p0.y) * (dX - p0.x);
						if((dA * dB >= 0.0) && (dB * dC >= 0.0) && (dC * dA >= 0.0)){
							polyData->nBitmap[lX + lIndex] = nPolygon;
							lPointCount++;
						}
					}
				}
				n2 = nVertex;
			}
		}
	}
//	printf(" %ld points assigned\n",lPointCount);
}


int Polygon_IO_Get_Byte_Order(char *cFilename)
{			/* Determine the byte order of a binary (.bp3d) file */
			/* Isn't guaranteed to work */
FILE *fp;
int nByteOrder = 0,nPolygons,nVertexMax;
double dPolyXMax,dPolyYMax,dSeedToVertexMax2;

	fp = fopen(cFilename,"r");
	nPolygons = Polygon_IO_Get_Int(fp);		/* Read the first integer (nPolygons) */
	nVertexMax = Polygon_IO_Get_Int(fp);
	dPolyXMax = Polygon_IO_Get_Float(fp);
	dPolyYMax = Polygon_IO_Get_Float(fp);
	dSeedToVertexMax2 = Polygon_IO_Get_Float(fp);
	fclose(fp);
	if(nPolygons < 0) nByteOrder = 1;		/* Tests to determine the byte order */
	if(nVertexMax > 1000) nByteOrder = 1;
	if(dPolyXMax < 1.0) nByteOrder = 1;
	if(dPolyYMax < 1.0) nByteOrder = 1;

	return (nByteOrder);
}


double Polygon_IO_Get_Double(FILE *fp)
{			/* Return a double extracted from a binary file */
double *d;
char s[sizeof(double)];
int c,n;

	for(n = 0; n < sizeof(double); n++){
		c = fgetc(fp);
		s[n] = c;
	}
	d = (double*)&s[0];

	return (*d);
}


double Polygon_IO_Get_Double_Reverse(FILE *fp)
{			/* Return a double extracted from a binary file using reverse byte ordering */
double *d;
char s[sizeof(double)];
int c,n;

	for(n = sizeof(double) - 1; n >= 0; n--){
		c = fgetc(fp);
		s[n] = c;
	}
	d = (double*)&s[0];

	return (*d);
}


void Polygon_IO_Get_File_Attributes(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Return the number of polygons in a file */
FILE *fp;
int nLoadCentre,nVertex = 0;
double d,dX,dY,dSeedX,dSeedY;

	fp = fopen(cFilename,"r");
	nLoadCentre = 1;
	Polygon_IO_Reset_Polygon_Attributes(polyData);
	while(fscanf(fp,"%lf %lf",&dX,&dY) != EOF){
		switch(nLoadCentre){
		case 0:
			if(dX == -1000.0 && dY == -1000.0){		/* This is the polygon terminator, (-1000,-1000) */
				nLoadCentre = 1;
				polyData->nPolygons++;
				if(nVertex > polyData->nVertexMax) polyData->nVertexMax = nVertex;
			}
			else{				/* This is a vertex */
				nVertex++;
				d = ((dSeedX - dX) * (dSeedX - dX)) + ((dSeedY - dY) * (dSeedY - dY));		/* Centre-vertex distance squared */
				if(d > polyData->dSeedToVertexMax2) polyData->dSeedToVertexMax2 = d;
				if(dX > polyData->pSize.x) polyData->pSize.x = dX;
				if(dY > polyData->pSize.y) polyData->pSize.y = dY;
			}
			break;

		case 1:
			if(dX != -1000.0 && dY != -1000.0){			/* This is the polygon seed point */
				nLoadCentre = 0;
				dSeedX = dX;			/* Store seed point position */
				dSeedY = dY;
			}
			nVertex = 0;
			break;
		}
	}
	fclose(fp);
}


void Polygon_IO_Get_File_Attributes_Binary(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Return the number of polygons in a binary file */
FILE *fp;
long lOld,lNew;
int nLoadCentre,nVertex = 0;
double d,dSeedX,dSeedY;
float fX,fY;

	fp = fopen(cFilename,"r");
	nLoadCentre = 1;
	Polygon_IO_Reset_Polygon_Attributes(polyData);
	do{
		lOld = ftell(fp);
		fX = Polygon_IO_Get_Float(fp);
		fY = Polygon_IO_Get_Float(fp);
		switch(nLoadCentre){
		case 0:
			if(fX == -1000.0 && fY == -1000.0){		/* This is the polygon terminator, (-1000,-1000) */
				nLoadCentre = 1;
				polyData->nPolygons++;
			}
			else{				/* This is a vertex */
				nVertex++;
				d = ((dSeedX - fX) * (dSeedX - fX)) + ((dSeedY - fY) * (dSeedY - fY));		/* Centre-vertex distance squared */
				if(d > polyData->dSeedToVertexMax2) polyData->dSeedToVertexMax2 = d;
				if(fX > polyData->pSize.x) polyData->pSize.x = fX;
				if(fY > polyData->pSize.y) polyData->pSize.y = fY;
			}
			break;

		case 1:
			if(fX != -1000.0 && fY != -1000.0){			/* This is the polygon seed point */
				nLoadCentre = 0;
				dSeedX = fX;			/* Store seed point position */
				dSeedY = fY;
			}
			if(nVertex > polyData->nVertexMax) polyData->nVertexMax = nVertex;
			nVertex = 0;
			break;
		}
		lNew = ftell(fp);
	}
	while(lOld != lNew);
	
	fclose(fp);
}


float Polygon_IO_Get_Float(FILE *fp)
{			/* Return a float extracted from a binary file */
float *f;
char s[sizeof(float)];
int c,n;

//	for(n = 0; n < sizeof(float); n++){
	for(n = 0; n < 4; n++){
		c = fgetc(fp);
		s[n] = c;
	}
	f = (float*)&s[0];

	return (*f);
}


float Polygon_IO_Get_Float_Reverse(FILE *fp)
{			/* Return a float extracted from a binary file using reverse byte ordering */
float *f;
char s[sizeof(float)];
int c,n;

	for(n = sizeof(float) - 1; n >= 0; n--){
		c = fgetc(fp);
		s[n] = c;
	}
	f = (float*)&s[0];

	return (*f);
}


int Polygon_IO_Get_Group_Cell_Count_Max(POLYGONATTRIBUTES *polyData, int nLayer)
{			/* Get the maximum number of cells in a group in nLayer */
int n,nCount = 0,nCountMax = 1,nLastGroup = -1;

	if(polyData->nGroupSet == 0){			/* Groups were not set, one cell per group per layer */
		return (nCountMax);
	}
	else{
		for(n = 0; n < polyData->nPolygons; n++){		/* Count unique groups in nLayer */
			if(polyData->nPolyLayer[n] == nLayer){
				if(polyData->nGroup[n] != nLastGroup){
					nLastGroup = polyData->nGroup[n];
					nCount = 1;
				}
				else{
					nCount++;
					if(nCount > nCountMax) nCountMax = nCount;
				}
			}
		}
	}

	return (nCountMax);
}


int Polygon_IO_Get_Group_Count(POLYGONATTRIBUTES *polyData)
{			/* How many groups are in polyData ? */
int n,nGroups = 0;

	if(polyData->nGroupSet == 0){			/* Groups were not set */
		return (0);
	}
	else{
		for(n = 0; n < polyData->nPolygons; n++){		/* Assumes that groups are numbered consecutively, without gaps */
			if(polyData->nGroup[n] > nGroups) nGroups = polyData->nGroup[n];
		}
	}
	nGroups++;		/* First group is group 0 */

	return (nGroups);
}


int Polygon_IO_Get_Group_Count_In_Layer(POLYGONATTRIBUTES *polyData,int nLayer)
{			/* Get the number of groups in nLayer */
int n,nCount = 0,nLastGroup = -1;

	if(polyData->nGroupSet == 0){			/* Groups were not set, use cell count instead */
		nCount = Polygon_IO_Get_Polygon_Count_In_Layer(polyData,nLayer);
	}
	else{
		for(n = 0; n < polyData->nPolygons; n++){		/* Count unique groups in nLayer */
			if(polyData->nPolyLayer[n] == nLayer){
				if(polyData->nGroup[n] != nLastGroup){
					nCount++;
					nLastGroup = polyData->nGroup[n];
				}
			}
		}
	}

	return (nCount);
}


int Polygon_IO_Get_Int(FILE *fp)
{			/* Return an integer extracted from a binary file */
char s[sizeof(int)];
int *i,c,n;

//	for(n = 0; n < sizeof(int); n++){
	for(n = 0; n < 4; n++){
		c = fgetc(fp);
		s[n] = c;
	}
	i = (int*)&s[0];

	return (*i);
}


int Polygon_IO_Get_Int_Reverse(FILE *fp)
{			/* Return an integer extracted from a binary file (reverse byte ordering) */
char s[sizeof(int)];
int *i,c,n;

	for(n = sizeof(int) - 1; n >= 0; n--){
		c = fgetc(fp);
		s[n] = c;
	}
	i = (int*)&s[0];

	return (*i);
}


double Polygon_IO_Get_M(POLYGONATTRIBUTES *polyData,int nX,int nY,int nMComponent)
{			/* Get M at (nX,nY) (bitmap must be defined first and M data loaded) */
int nPolygon = polyData->nBitmap[nX + (nY * polyData->nPolyXMax)];
double dM = 0.0;

	if(nPolygon != -1){
		switch (nMComponent){
		case 0:
			dM = polyData->pM[nPolygon].x;
			break;
		case 1:
			dM = polyData->pM[nPolygon].y;
			break;
		case 2:
			dM = polyData->pM[nPolygon].z;
			break;
		}
	}

	return (dM);
}


void Polygon_IO_Get_M_X_Y_Z(POLYGONATTRIBUTES *polyData,int nX,int nY,double *dMx,double *dMy,double *dMz)
{			/* Get Mx, My and Mz at (nX,nY) */
int nPolygon = polyData->nBitmap[nX + (nY * polyData->nPolyXMax)];

	if(nPolygon != -1){
		*dMx = polyData->pM[nPolygon].x;
		*dMy = polyData->pM[nPolygon].y;
		*dMz = polyData->pM[nPolygon].z;
	}
	else{
		*dMx = 0.0;
		*dMy = 0.0;
		*dMz = 0.0;
	}
}


int Polygon_IO_Get_Max_Vertex_Count_In_Layer(POLYGONATTRIBUTES *polyData, int nLayer)
{			/* Return the maximum number of vertices for the cells in nLayer */
int n,nVMax = 0;

	for(n = 0; n < polyData->nPolygons; n++){
		if(polyData->nPolyLayer[n] == nLayer){
			if(polyData->siPolyVertices[n] > nVMax) nVMax = polyData->siPolyVertices[n];		/* Core cell has the same number of vertices as the perimeter */
		}
	}

	return (nVMax);
}


double Polygon_IO_Get_Polygon_Area(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Calculate the area of a polygon */
double dArea = 0.0;
int nVertex,nV2;
point p1,p2;

	if(nPolygon < polyData->nPolygons){
		nV2 = polyData->siPolyVertices[nPolygon] - 1;			/* Last vertex */
		for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
			p1 = Polygon_IO_Get_Vertex(polyData,nPolygon,nVertex);
			p2 = Polygon_IO_Get_Vertex(polyData,nPolygon,nV2);
			dArea += ((p1.x * p2.y) - (p2.x * p1.y));
			nV2 = nVertex;
		}
		dArea *= 0.5;
	}

	return (dArea);
}


void Polygon_IO_Get_Polygon_Array_Circle(POLYGONATTRIBUTES *polyData, int nZLayer,double dXCentre,double dYCentre,double dRadius,double dXSize,double dYSize)
{			/* Create a circular object consiting of cubic (cuboid) polygons in nLayer */
			/* dXSize and dYSize are cell sizes, in Angstroms */
int nCellCount = 0;
double dX,dXC,dY,dYC,dXMax = 0.0,dYMax = 0.0;

	dX = -dRadius;
	while(dX <= dRadius){			/* Find how many cells there will be */
		dY = -dRadius;
		dXC = dX + (0.5 * dXSize);
		while(dY <= dRadius){
			dYC = dY + (0.5 * dYSize);
			if((dXC * dXC) + (dYC * dYC) < (dRadius * dRadius)){
				if(dX + dXSize + dXCentre > dXMax) dXMax = dX + dXSize + dXCentre;
				if(dY + dYSize + dYCentre > dYMax) dYMax = dY + dYSize + dYCentre;
				nCellCount++;
			}
			dY += dYSize;
		}
		dX += dXSize;
	}

	polyData->nPolygons = nCellCount;			/* Set the header details */
	polyData->nVertexMax = 4;
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = dXMax;
	polyData->pSize.y = dYMax;
	polyData->dSeedToVertexMax2 = 0.25 * ((dXSize * dXSize) + (dYSize * dYSize));

	nCellCount = 0;
	dX = -dRadius;
	while(dX <= dRadius){			/* Put the cells into the structure */
		dY = -dRadius;
		dXC = dX + (0.5 * dXSize);
		while(dY <= dRadius){
			dYC = dY + (0.5 * dYSize);
			if((dXC * dXC) + (dYC * dYC) < (dRadius * dRadius)){
				polyData->nPolyLayer[nCellCount] = nZLayer;
				polyData->siPolyVertices[nCellCount] = 4;
				polyData->fBoundary[nCellCount] = 0.0;
				polyData->pSeed[nCellCount].x = dXC + dXCentre;		/* Seed point is centre of polygon */
				polyData->pSeed[nCellCount].y = dYC + dYCentre;
				polyData->pVertex[nCellCount].x = dX + dXCentre;				/* Set vertices */
				polyData->pVertex[nCellCount].y = dY + dYCentre;
				polyData->pVertex[nCellCount + (polyData->nPolygons)].x = dX + dXCentre;
				polyData->pVertex[nCellCount + (polyData->nPolygons)].y = dY + dYSize + dYCentre;
				polyData->pVertex[nCellCount + (2 * polyData->nPolygons)].x = dX + dXSize + dXCentre;
				polyData->pVertex[nCellCount + (2 * polyData->nPolygons)].y = dY + dYSize + dYCentre;
				polyData->pVertex[nCellCount + (3 * polyData->nPolygons)].x = dX + dXSize + dXCentre;
				polyData->pVertex[nCellCount + (3 * polyData->nPolygons)].y = dY + dYCentre;
				nCellCount++;
			}
			dY += dYSize;
		}
		dX += dXSize;
	}
}


void Polygon_IO_Get_Polygon_Array_Circle_Outline(POLYGONATTRIBUTES *polyData,int nZLayer,int nDotsX,int nDotsY,double dPitchX,double dPitchY,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle)
{			/* Create an outline of an array of circles, or ellipses */
			/* dPitchX and dPitchY are unit cell sizes, in Angstroms */
int nCellCount = 0,nV,nVMax,nX,nY;
double dQ,dX,dXCentre,dXCorner,dY,dYCentre,dYCorner;
double dPi = 3.1415926535;
double dPhi = dEllipseAngle * dPi / 180.0;	/* Angle of ellipse from y-axis */

	nCellCount = (nDotsX * nDotsY);	/* There may be multiple dots */
	nVMax = (int)(360.0 / dTheta);
	polyData->nPolygons = nCellCount;			/* Set the header details */
	polyData->nVertexMax = nVMax;					/* Number of vertices */
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = (dPitchX * (double)nDotsX);
	polyData->pSize.y = (dPitchY * (double)nDotsY);
	polyData->dSeedToVertexMax2 = (dRadius * dRadius);
	dXCorner = (0.5 * dPitchX);		/* Centre dot in bit cell */
	dYCorner = (0.5 * dPitchY);

	nCellCount = 0;
	for(nX = 0; nX < nDotsX; nX++){
		dXCentre = ((double)nX * dPitchX) + dXCorner;
		for(nY = 0; nY < nDotsY; nY++){
			dYCentre = ((double)nY * dPitchY) + dYCorner;

			polyData->nPolyLayer[nCellCount] = nZLayer;			/* First quadrant */
			polyData->siPolyVertices[nCellCount] = nVMax;
			polyData->fBoundary[nCellCount] = 0.0;
			polyData->pSeed[nCellCount].x = dXCentre;		/* Seed point is centre of polygon */
			polyData->pSeed[nCellCount].y = dYCentre;
			for(nV = 0; nV < nVMax; nV++){
				dQ = ((double)nV) * dPi * 2.0 / ((double)nVMax);
				dX = (dRadius * cos (dQ) * cos (dPhi)) - (dRadius * dEllipticity * sin (dQ) * sin (dPhi));
				dY = (dRadius * cos (dQ) * sin (dPhi)) + (dRadius * dEllipticity * sin (dQ) * cos (dPhi));
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x = dX + dXCentre;
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y = dY + dYCentre;
			}
			nCellCount++;
		}
	}
}


void Polygon_IO_Get_Polygon_Array_Cubic(POLYGONATTRIBUTES *polyData,int nCells,int *nX,int *nY,int *nZ,int nXCount,int nYCount,double dXSize,double dYSize)
{			/* Create an array (nXCount x nYCount) of cubic (cuboid) polygons in nLayer */
			/* dXSize and dYSize are cell sizes, in Angstroms */
int n;

	polyData->nPolygons = nCells;			/* Set the header details */
	polyData->nVertexMax = 4;
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = dXSize * (double)nXCount;
	polyData->pSize.y = dYSize * (double)nYCount;
	polyData->dSeedToVertexMax2 = 0.25 * ((dXSize * dXSize) + (dYSize * dYSize));

	for(n = 0; n < polyData->nPolygons; n++){		/* Add each polygon */
		polyData->nPolyLayer[n] = nZ[n];
		polyData->siPolyVertices[n] = 4;
		polyData->fBoundary[n] = 0.0;
		polyData->pSeed[n].x = dXSize * (0.5 + (double)nX[n]);		/* Seed point is centre of polygon */
		polyData->pSeed[n].y = dYSize * (0.5 + (double)nY[n]);

		polyData->pVertex[n].x = dXSize * (double)nX[n];				/* Set vertices */
		polyData->pVertex[n].y = dYSize * (double)nY[n];
		polyData->pVertex[n + (polyData->nPolygons)].x = dXSize * (double)nX[n];
		polyData->pVertex[n + (polyData->nPolygons)].y = dYSize * (1.0 + (double)nY[n]);
		polyData->pVertex[n + (2 * polyData->nPolygons)].x = dXSize * (1.0 + (double)nX[n]);
		polyData->pVertex[n + (2 * polyData->nPolygons)].y = dYSize * (1.0 + (double)nY[n]);
		polyData->pVertex[n + (3 * polyData->nPolygons)].x = dXSize * (1.0 + (double)nX[n]);
		polyData->pVertex[n + (3 * polyData->nPolygons)].y = dYSize * (double)nY[n];
	}
}


int Polygon_IO_Get_Polygon_Array_Cubic_Dot(POLYGONATTRIBUTES *polyData, int nDotsX, int nDotsY, double dPitchX, double dPitchY, double dDotX, double dDotY, int nSubCellsX, int nSubCellsY, int nLayerFrom, int nLayerTo)
{			/* Generate an array of cubic dots with sub-dot discretisation */
int n,nCells,nGroup,nLayer,nSX,nSY,nX,nY;
double dX,dXCorner,dXOffset,dXSize;
double dY,dYCorner,dYOffset,dYSize;

	nCells = nDotsX * nDotsY * nSubCellsX * nSubCellsY * (1 + nLayerTo - nLayerFrom);
	dXSize = dDotX / (double)nSubCellsX;		/* Size of cells within each dot */
	dYSize = dDotY / (double)nSubCellsY;
	polyData->nPolygons = nCells;						/* Set the header details */
	polyData->nVertexMax = 4;
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = (double)nDotsX * dPitchX;
	polyData->pSize.y = (double)nDotsY * dPitchY;
	polyData->dSeedToVertexMax2 = 0.25 * ((dXSize * dXSize) + (dYSize * dYSize));

	n = 0;
	dXCorner = (0.5 * (dPitchX - dDotX));		/* Centre dot in bit cell */
	dYCorner = (0.5 * (dPitchY - dDotY));
	for(nLayer = nLayerFrom; nLayer <= nLayerTo; nLayer++){
		nGroup = 0;
		for(nX = 0; nX < nDotsX; nX++){
			dXOffset = ((double)nX * dPitchX) + dXCorner;
			for(nY = 0; nY < nDotsY; nY++){
				dYOffset = ((double)nY * dPitchY) + dYCorner;
				for(nSX = 0; nSX < nSubCellsX; nSX++){
					dX = (double)nSX * dXSize;
					for(nSY = 0; nSY < nSubCellsY; nSY++){
						dY = (double)nSY * dYSize;
						polyData->nPolyLayer[n] = nLayer;
						polyData->siPolyVertices[n] = 4;
						polyData->fBoundary[n] = 0.0;
						polyData->pSeed[n].x = dX + (0.5 * dXSize) + dXOffset;		/* Seed point is centre of polygon */
						polyData->pSeed[n].y = dY + (0.5 * dYSize) + dYOffset;

						polyData->pVertex[n].x = dX + dXOffset;				/* Set vertices */
						polyData->pVertex[n].y = dY + dYOffset;
						polyData->pVertex[n + (polyData->nPolygons)].x = dX + dXOffset;
						polyData->pVertex[n + (polyData->nPolygons)].y = dY + dYSize + dYOffset;
						polyData->pVertex[n + (2 * polyData->nPolygons)].x = dX + dXOffset + dXSize;
						polyData->pVertex[n + (2 * polyData->nPolygons)].y = dY + dYSize + dYOffset;
						polyData->pVertex[n + (3 * polyData->nPolygons)].x = dX + dXOffset + dXSize;
						polyData->pVertex[n + (3 * polyData->nPolygons)].y = dY + dYOffset;
						n++;
					}
				}
				nGroup++;
			}
		}
	}

	return (nCells);
}


void Polygon_IO_Get_Polygon_Array_Hexagon(POLYGONATTRIBUTES *polyData,int nYCount,int *nZLayer,int nLayers,double dXSize,double dYSize)
{			/* Create an array hexagonal polygons in nLayer */
			/* dXSize and dYSize are cell sizes, in Angstroms */
int n,nSeedX,nSeedY,nXCount,nZ = 0;
double dX,dXMin = 0.0,dXStep,dY,dYMin = 0.0,dYOffset,dYStep,dOneOverRoot3 = 1.0 / sqrt(3.0);

	dYStep = dYSize / ((double)nYCount);
	dXStep = pow(dYStep * dYStep * 0.75,0.5);
	nXCount = (int)(dXSize / dXStep) + 1;
	polyData->nPolygons = (nXCount * nYCount * nLayers);			/* Set the header details */
	polyData->nVertexMax = 6;
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = dXSize;
	polyData->pSize.y = dYSize;
	polyData->dSeedToVertexMax2 = (dYStep * dOneOverRoot3) * (dYStep * dOneOverRoot3);

	n = 0;
	for(nSeedX = 0; nSeedX < nXCount; nSeedX++){
		dX = (0.25 + (double)nSeedX) * dXStep;
		dYOffset = (dYStep * 0.5) * (double)(nSeedX % 2);
		for(nSeedY = 0; nSeedY < nYCount; nSeedY++){
			dY = ((0.25 + (double)nSeedY) * dYStep) + dYOffset;
			if(n < polyData->nPolygons){
				polyData->nPolyLayer[n] = nZLayer[nZ];
				polyData->siPolyVertices[n] = 6;
				polyData->fBoundary[n] = 0.0;
				polyData->pSeed[n].x = dX;		/* Seed point is centre of polygon */
				polyData->pSeed[n].y = dY;
				polyData->pVertex[n].x = dX	+ (dYStep * dOneOverRoot3);			/* Set vertices */
				polyData->pVertex[n].y = dY;
				polyData->pVertex[n + (polyData->nPolygons)].x = dX + (0.5 * dYStep * dOneOverRoot3);
				polyData->pVertex[n + (polyData->nPolygons)].y = dY - (0.5 * dYStep);
				polyData->pVertex[n + (2 * polyData->nPolygons)].x = dX - (0.5 * dYStep * dOneOverRoot3);
				polyData->pVertex[n + (2 * polyData->nPolygons)].y = dY - (0.5 * dYStep);
				polyData->pVertex[n + (3 * polyData->nPolygons)].x = dX - (dYStep * dOneOverRoot3);
				polyData->pVertex[n + (3 * polyData->nPolygons)].y = dY;
				polyData->pVertex[n + (4 * polyData->nPolygons)].x = dX - (0.5 * dYStep * dOneOverRoot3);
				polyData->pVertex[n + (4 * polyData->nPolygons)].y = dY + (0.5 * dYStep);
				polyData->pVertex[n + (5 * polyData->nPolygons)].x = dX + (0.5 * dYStep * dOneOverRoot3);
				polyData->pVertex[n + (5 * polyData->nPolygons)].y = dY + (0.5 * dYStep);
				if((dX - (dYStep * dOneOverRoot3)) < dXMin) dXMin = (dX - (dYStep * dOneOverRoot3));		/* All the vertices should be >= 0, if they aren't they must be shifted later */
				if((dY - (0.5 * dYStep)) < dYMin) dYMin = (dY - (0.5 * dYStep));
			}
			n++;
		}
	}

	Polygon_IO_Translate_Polygons(polyData,-dXMin,-dYMin);		/* Shift polygons so vertices are all >= 0 */
}


void Polygon_IO_Get_Polygon_Array_Pattern(POLYGONATTRIBUTES *polyData,int nCells,n_point *pPattern,int nXMax,int nYMax,double dXSize,double dYSize)
{			/* Create an array of cubic polygons from a pattern */
			/* dXSize and dYSize are cell sizes, in Angstroms */
int n,nX,nY;

	polyData->nPolygons = nCells;			/* Set the header details */
	polyData->nVertexMax = 4;
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = dXSize * (double)(nXMax + 1);
	polyData->pSize.y = dYSize * (double)(nYMax + 1);
	polyData->dSeedToVertexMax2 = 0.25 * ((dXSize * dXSize) + (dYSize * dYSize));

	for(n = 0; n < polyData->nPolygons; n++){		/* Add each polygon */
		nX = pPattern[n].x;
		nY = pPattern[n].y;
		polyData->nPolyLayer[n] = pPattern[n].z;
		polyData->siPolyVertices[n] = 4;
		polyData->fBoundary[n] = 0.0;
		polyData->pSeed[n].x = dXSize * (0.5 + (double)nX);		/* Seed point is centre of polygon */
		polyData->pSeed[n].y = dYSize * (0.5 + (double)nY);

		polyData->pVertex[n].x = dXSize * (double)nX;				/* Set vertices */
		polyData->pVertex[n].y = dYSize * (double)nY;
		polyData->pVertex[n + (polyData->nPolygons)].x = dXSize * (double)nX;
		polyData->pVertex[n + (polyData->nPolygons)].y = dYSize * (1.0 + (double)nY);
		polyData->pVertex[n + (2 * polyData->nPolygons)].x = dXSize * (1.0 + (double)nX);
		polyData->pVertex[n + (2 * polyData->nPolygons)].y = dYSize * (1.0 + (double)nY);
		polyData->pVertex[n + (3 * polyData->nPolygons)].x = dXSize * (1.0 + (double)nX);
		polyData->pVertex[n + (3 * polyData->nPolygons)].y = dYSize * (double)nY;
	}
}


void Polygon_IO_Get_Polygon_Array_Spin_Ice_Hexagon(POLYGONATTRIBUTES *polyData,int nZLayer,int nCellsX,int nCellsY,double dPitchX,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle)
{			/* Create an outline of an array of circles, or ellipses in a spin-ice-like pattern with hexagonal cells */
			/* dPitchX is the spacing between cell centres, in Angstroms */
int nCellCount = 0,nDotsX,nDotsY,nV,nVMax,nX,nY;
double dQ,dX,dXCentre,dY,dYCentre,dYOffset,dYStep;
double dPitchY = 0.25 * dPitchX * sqrt(3.0);
double dPhi[2],dPi = 3.1415926535;
double dPhi0 = 0.0;									/* Angles of major axes */
double dPhiP60 = 60.0 * dPi / 180.0;
double dPhiM60 = -60.0 * dPi / 180.0;

	nCellCount = 1 + (5 * nCellsX);
	nCellCount += (nCellsY - 1) * (2 + (3 * nCellsX));
	nVMax = (int)(360.0 / dTheta);
	polyData->nPolygons = nCellCount;			/* Set the header details */
	polyData->nVertexMax = nVMax;					/* Number of vertices */
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = (2.0 * dPitchX) + (1.5 * dPitchX * (double)nCellsX);
	polyData->pSize.y = (2.0 * dPitchY) + (4.0 * dPitchY * (double)nCellsY);
	polyData->dSeedToVertexMax2 = (dRadius * dRadius);

	nCellCount = 0;
	nDotsX = (2 * nCellsX) + 1;
	for(nX = 0; nX < nDotsX; nX++){
		dXCentre = (0.5 * dPitchX) + (0.75 * dPitchX * (double)nX);
		switch(nX % 4){			/* Set y starting point */
		case 0:
			dYOffset = 2.0 * dPitchY;
			dYStep = 2.0 * dPitchY;
			nDotsY = 1 + (2 * nCellsY);
			dPhi[0] = dPhiM60;
			dPhi[1] = dPhiP60;
			if(nX == 0) nDotsY--;			/* One less element in the first column */
			if(nX == nDotsX - 1){
				dYOffset = 4.0 * dPitchY;		/* Start higher up and add one less dot if this is the last column */
				nDotsY--;
				dPhi[0] = dPhiP60;
				dPhi[1] = dPhiM60;
			}
			break;
		case 1:
			dYOffset = dPitchY;
			dYStep = 4.0 * dPitchY;
			nDotsY = 1 + nCellsY;
			dPhi[0] = dPhi0;
			dPhi[1] = dPhi0;
			break;
		case 2:
			dYOffset = 2.0 * dPitchY;
			dYStep = 2.0 * dPitchY;
			nDotsY = 1 + (2 * nCellsY);
			dPhi[0] = dPhiP60;
			dPhi[1] = dPhiM60;
			if(nX == nDotsX - 1) nDotsY--;		/* One less element in the last column */
			break;
		case 3:
			dYOffset = 3.0 * dPitchY;
			dYStep = 4.0 * dPitchY;
			nDotsY = 1 + nCellsY;
			dPhi[0] = dPhi0;
			dPhi[1] = dPhi0;
			break;
		}

		for(nY = 0; nY < nDotsY; nY++){
			dYCentre = dYOffset + (dYStep * (double)nY);
			polyData->nPolyLayer[nCellCount] = nZLayer;
			polyData->siPolyVertices[nCellCount] = nVMax;
			polyData->fBoundary[nCellCount] = 0.0;
			polyData->pSeed[nCellCount].x = dXCentre;		/* Seed point is centre of polygon */
			polyData->pSeed[nCellCount].y = dYCentre;
			for(nV = 0; nV < nVMax; nV++){
				dQ = ((double)nV) * dPi * 2.0 / ((double)nVMax);
				dX = (dRadius * cos (dQ) * cos (dPhi[nY % 2])) - (dRadius * dEllipticity * sin (dQ) * sin (dPhi[nY % 2]));
				dY = (dRadius * cos (dQ) * sin (dPhi[nY % 2])) + (dRadius * dEllipticity * sin (dQ) * cos (dPhi[nY % 2]));
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x = dX + dXCentre;
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y = dY + dYCentre;
			}
			nCellCount++;
		}
	}
}


void Polygon_IO_Get_Polygon_Array_Spin_Ice_Square(POLYGONATTRIBUTES *polyData,int nZLayer,int nCellsX,int nCellsY,double dPitchX,double dPitchY,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle)
{			/* Create an outline of an array of circles, or ellipses in a spin-ice-like pattern */
			/* dPitchX and dPitchY are unit cell sizes, in Angstroms */
int nCellCount = 0,nDotsX,nDotsY,nV,nVMax,nX,nY;
double dQ,dX,dXCentre,dXCorner,dY,dYCentre,dYCorner;
double dPi = 3.1415926535;
double dPhi = dEllipseAngle * dPi / 180.0;	/* Angle of ellipse from y-axis */

	nDotsX = (2 * nCellsX) + 1;
	nDotsY = nCellsY;
	nCellCount = (nDotsX * nDotsY) + nCellsX;	/* There may be multiple dots */
	nVMax = (int)(360.0 / dTheta);
	polyData->nPolygons = nCellCount;			/* Set the header details */
	polyData->nVertexMax = nVMax;					/* Number of vertices */
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = (dPitchX * (double)nDotsX);
	polyData->pSize.y = (dPitchY * (double)nDotsY) + dPitchX;
	polyData->dSeedToVertexMax2 = (dRadius * dRadius);
	dXCorner = (0.5 * dPitchX);		/* Centre dot in bit cell */
	dYCorner = (0.5 * dPitchY);

	nCellCount = 0;
	for(nX = 0; nX < nDotsX; nX++){
		dXCentre = ((double)nX * dPitchX) + dXCorner;
		for(nY = 0; nY < nDotsY; nY++){
			dYCentre = ((double)nY * dPitchY) + dYCorner;

			polyData->nPolyLayer[nCellCount] = nZLayer;			/* First quadrant */
			polyData->siPolyVertices[nCellCount] = nVMax;
			polyData->fBoundary[nCellCount] = 0.0;
			polyData->pSeed[nCellCount].x = dXCentre;		/* Seed point is centre of polygon */
			polyData->pSeed[nCellCount].y = dYCentre;
			for(nV = 0; nV < nVMax; nV++){
				dQ = ((double)nV) * dPi * 2.0 / ((double)nVMax);
				dX = (dRadius * cos (dQ) * cos (dPhi)) - (dRadius * dEllipticity * sin (dQ) * sin (dPhi));
				dY = (dRadius * cos (dQ) * sin (dPhi)) + (dRadius * dEllipticity * sin (dQ) * cos (dPhi));
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x = dX + dXCentre;
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y = dY + dYCentre;
			}

			if(nX % 2 == 0){			/* Rotate dots in alternate columns by 90 degrees */
				for(nV = 0; nV < polyData->siPolyVertices[nCellCount]; nV++){
					dX = polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x - polyData->pSeed[nCellCount].x;
					dY = polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y - polyData->pSeed[nCellCount].y;
					polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x = polyData->pSeed[nCellCount].x + dY;
					polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y = polyData->pSeed[nCellCount].y - dX;
				}
			}

			nCellCount++;
		}
	}

	for(nX = 0; nX < nCellsX; nX++){			/* Add top row of dots to close top cells */
		dXCentre = ((double)((2 * nX) + 1) * dPitchX) + dXCorner;
		nY = nDotsY;
		dYCentre = ((double)nY * dPitchY) + dYCorner;

		polyData->nPolyLayer[nCellCount] = nZLayer;			/* First quadrant */
		polyData->siPolyVertices[nCellCount] = nVMax;
		polyData->fBoundary[nCellCount] = 0.0;
		polyData->pSeed[nCellCount].x = dXCentre;		/* Seed point is centre of polygon */
		polyData->pSeed[nCellCount].y = dYCentre;
		for(nV = 0; nV < nVMax; nV++){
			dQ = ((double)nV) * dPi * 2.0 / ((double)nVMax);
			dX = (dRadius * cos (dQ) * cos (dPhi)) - (dRadius * dEllipticity * sin (dQ) * sin (dPhi));
			dY = (dRadius * cos (dQ) * sin (dPhi)) + (dRadius * dEllipticity * sin (dQ) * cos (dPhi));
			polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x = dX + dXCentre;
			polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y = dY + dYCentre;
		}
		nCellCount++;
	}
}


void Polygon_IO_Get_Polygon_Array_Spin_Ice_Triangle(POLYGONATTRIBUTES *polyData,int nZLayer,int nCellsX,int nCellsY,double dPitchX,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle)
{			/* Create an outline of an array of circles, or ellipses in a spin-ice-like pattern */
			/* dPitchX is the spacing between cell centres, in Angstroms */
int nCellCount = 0,nDotsX,nDotsY,nV,nVMax,nX,nY;
double dQ,dX,dXCentre,dXOffset,dXPitch,dY,dYCentre;
double dPitchY = 0.25 * dPitchX * sqrt(3.0);
double dPhi[2],dPi = 3.1415926535;
double dPhi0 = 0.0;									/* Angles of major axes */
double dPhiP60 = 60.0 * dPi / 180.0;
double dPhiM60 = -60.0 * dPi / 180.0;

	nCellCount = 3 + (2 * (nCellsX - 1));		/* Triangular cell count */
	if(nCellsX %2 == 0){
		nCellCount += ((((3 * nCellsX) / 2) + 1) * (nCellsY - 1));
	}
	else{
		nCellCount += ((3 * ((nCellsX + 1) / 2)) * (nCellsY - 1));
		if(nCellsY % 2 == 0){
			nCellCount -= ((nCellsY - 2) / 2);
		}
		else{
			nCellCount -= ((nCellsY - 1) / 2);
		}
	}

	nVMax = (int)(360.0 / dTheta);
	polyData->nPolygons = nCellCount;			/* Set the header details */
	polyData->nVertexMax = nVMax;					/* Number of vertices */
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
	polyData->pSize.x = dPitchX + (0.5 * dPitchX * (double)nCellsX);
	polyData->pSize.y = dPitchY + (2.0 * dPitchY * (double)nCellsY);
	polyData->dSeedToVertexMax2 = (dRadius * dRadius);

	nCellCount = 0;
	nDotsY = (2 * nCellsY) + 1;
	for(nY = 0; nY < nDotsY; nY++){
		dYCentre = (0.5 * dPitchY) + (dPitchY * (double)nY);
		switch(nY % 4){
		case 0:		/* First row */
			nDotsX = ((nCellsX - 1) / 2) + 1;
			dXOffset = (0.5 * dPitchX);
			dXPitch = dPitchX;
			dPhi[0] = dPhi0;
			dPhi[1] = dPhi0;
			break;
		case 1:
			nDotsX = nCellsX + 1;
			dXOffset = (0.25 * dPitchX);
			dXPitch = (0.5 * dPitchX);
			dPhi[0] = dPhiP60;
			dPhi[1] = dPhiM60;
			break;
		case 2:
			nDotsX = (nCellsX / 2);
			dXOffset = dPitchX;
			dXPitch = dPitchX;
			dPhi[0] = dPhi0;
			dPhi[1] = dPhi0;
			break;
		case 3:
			nDotsX = nCellsX + 1;
			dXOffset = (0.25 * dPitchX);
			dXPitch = (0.5 * dPitchX);
			dPhi[0] = dPhiM60;
			dPhi[1] = dPhiP60;
			break;
		}
		for(nX = 0; nX < nDotsX; nX++){
			dXCentre = dXOffset + (dXPitch * (double)nX);
			polyData->nPolyLayer[nCellCount] = nZLayer;
			polyData->siPolyVertices[nCellCount] = nVMax;
			polyData->fBoundary[nCellCount] = 0.0;
			polyData->pSeed[nCellCount].x = dXCentre;		/* Seed point is centre of polygon */
			polyData->pSeed[nCellCount].y = dYCentre;
			for(nV = 0; nV < nVMax; nV++){
				dQ = ((double)nV) * dPi * 2.0 / ((double)nVMax);
				dX = (dRadius * cos (dQ) * cos (dPhi[nX % 2])) - (dRadius * dEllipticity * sin (dQ) * sin (dPhi[nX % 2]));
				dY = (dRadius * cos (dQ) * sin (dPhi[nX % 2])) + (dRadius * dEllipticity * sin (dQ) * cos (dPhi[nX % 2]));
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].x = dX + dXCentre;
				polyData->pVertex[nCellCount + (nV * polyData->nPolygons)].y = dY + dYCentre;
			}
			nCellCount++;
		}
	}
}


int Polygon_IO_Get_Polygon_At_XY(POLYGONATTRIBUTES *polyData,int nX,int nY)
{			/* Get the polygon at (nX,nY) (bitmap must be defined first) */
	return (polyData->nBitmap[nX+(nY * polyData->nPolyXMax)]);
}


point Polygon_IO_Get_Polygon_Centre(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Calculate the centre (centroid) of a polygon */
double dArea,dX = 0.0,dY = 0.0;
int nVertex,nV2;
point p1,p2,pCentre;

	if(nPolygon < polyData->nPolygons){
		dArea = Polygon_IO_Get_Polygon_Area(polyData,nPolygon);
		if(dArea != 0.0){
			nV2 = polyData->siPolyVertices[nPolygon] - 1;			/* Last vertex */
			for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
				p1 = Polygon_IO_Get_Vertex(polyData,nPolygon,nVertex);
				p2 = Polygon_IO_Get_Vertex(polyData,nPolygon,nV2);
				dX += (p1.x + p2.x) * ((p1.x * p2.y) - (p2.x * p1.y));
				dY += (p1.y + p2.y) * ((p1.x * p2.y) - (p2.x * p1.y));
				nV2 = nVertex;
			}
			dX /= (6.0 * dArea);
			dY /= (6.0 * dArea);
		}
	}
	pCentre.x = dX;
	pCentre.y = dY;

	return (pCentre);
}


int Polygon_IO_Get_Polygon_Count_In_Layer(POLYGONATTRIBUTES *polyData,int nLayer)
{			/* Get the number of polygons in nLayer */
int n,nCount=0;

	for(n=0;n<polyData->nPolygons;n++){
		if(polyData->nPolyLayer[n] == nLayer) nCount++;
	}

	return (nCount);
}


int Polygon_IO_Get_Polygon_Count_In_Region(POLYGONATTRIBUTES *polyData,int *nVertexMax,double dXStart,double dXStop,double dYStart,double dYStop)
{			/* Count the number of polygons in a region (based on seed point location) */
int n,nCount = 0;

	for(n = 0; n < polyData->nPolygons; n++){
		if(polyData->pSeed[n].x > dXStart){
			if(polyData->pSeed[n].x < dXStop){
				if(polyData->pSeed[n].y > dYStart){
					if(polyData->pSeed[n].y < dYStop){
						if(polyData->siPolyVertices[n] > *nVertexMax) *nVertexMax = polyData->siPolyVertices[n];
						nCount++;		/* Polygon is in the defined region */
					}
				}
			}
		}
	}

	return (nCount);
}


void Polygon_IO_Get_Polygon_Extent(POLYGONATTRIBUTES *polyData,int nPolygon,double *dSeedToVertex,double *dXMax,double *dYMax)
{			/* Get the maximum extent of a polygon */
int n;
double d,dX,dY;

	*dSeedToVertex = 0.0;
	*dXMax = 0.0;
	*dYMax = 0.0;
	for(n = 0; n < polyData->siPolyVertices[nPolygon]; n++){		/* Loop around vertices */
		dX = polyData->pVertex[nPolygon + (n * polyData->nPolygons)].x;		/* Vertex */
		dY = polyData->pVertex[nPolygon + (n * polyData->nPolygons)].y;
		if(dX > *dXMax) *dXMax = dX;
		if(dY > *dYMax) *dYMax = dY;
		d = ((dX - polyData->pSeed[nPolygon].x) * (dX - polyData->pSeed[nPolygon].x)) + ((dY - polyData->pSeed[nPolygon].y) * (dY - polyData->pSeed[nPolygon].y));
		if(d > *dSeedToVertex) *dSeedToVertex = d;
	}
}


void Polygon_IO_Get_Polygon_Limits(POLYGONATTRIBUTES *polyData,int nPolygon,double *dXMin,double *dXMax,double *dYMin,double *dYMax)
{			/* Get the vertices at the limits of a polygon */
int nV;
point p = Polygon_IO_Get_Vertex(polyData,nPolygon,0);

	*dXMin = p.x;
	*dXMax = p.x;
	*dYMin = p.y;
	*dYMax = p.y;
	for(nV = 1; nV < polyData->siPolyVertices[nPolygon]; nV++){
		p = Polygon_IO_Get_Vertex(polyData,nPolygon,nV);
		if(p.x < *dXMin) *dXMin = p.x;
		if(p.x > *dXMax) *dXMax = p.x;
		if(p.y < *dYMin) *dYMin = p.y;
		if(p.y > *dYMax) *dYMax = p.y;
	}
}


double Polygon_IO_Get_Polygon_Perimeter(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Get the perimeter length of a polygon */
int n2,nVertex;
double dBoundary = 0.0,dSegment;
point p1,p2;

	n2 = polyData->siPolyVertices[nPolygon] - 1;			/* Last vertex */
	for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
		p1 = Polygon_IO_Get_Vertex(polyData,nPolygon,nVertex);
		p2 = Polygon_IO_Get_Vertex(polyData,nPolygon,n2);
		dSegment = sqrt(((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)));		/* Length of common boundary */
		dBoundary += dSegment;
		n2 = nVertex;
	}

	return (dBoundary);
}


double Polygon_IO_Get_Polygon_Statistics(POLYGONATTRIBUTES *polyData)
{			/* Get average area and standard deviation of a set of polygons */
	return (Polygon_IO_Get_Polygon_Statistics_Region(polyData,0.0,polyData->pSize.x,0.0,polyData->pSize.y));
}


double Polygon_IO_Get_Polygon_Statistics_Region(POLYGONATTRIBUTES *polyData,double dXStart,double dXStop,double dYStart,double dYStop)
{			/* Get average area and standard deviation of a set of polygons within a defined region */
int nBiggest = 0,nPolygon,nPolygons = 0,nSmallest = 0;
double dArea,dAverageArea = 0.0,dBiggest = 0.0,dSigmaArea = 0.0,dSmallest = 1.0e100;

	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		if(polyData->pSeed[nPolygon].x > dXStart){
			if(polyData->pSeed[nPolygon].x < dXStop){
				if(polyData->pSeed[nPolygon].y > dYStart){
					if(polyData->pSeed[nPolygon].y < dYStop){
						dArea = Polygon_IO_Get_Polygon_Area(polyData,nPolygon);
						if(dArea < 0.0) dArea *= -1.0;
						if(dArea > dBiggest){
							dBiggest = dArea;
							nBiggest = nPolygon;
						}
						if(dArea < dSmallest){
							dSmallest = dArea;
							nSmallest = nPolygon;
						}
						dAverageArea += dArea;
						nPolygons++;
					}
				}
			}
		}
	}
	dAverageArea /= (double)nPolygons;

	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		if(polyData->pSeed[nPolygon].x > dXStart){
			if(polyData->pSeed[nPolygon].x < dXStop){
				if(polyData->pSeed[nPolygon].y > dYStart){
					if(polyData->pSeed[nPolygon].y < dYStop){
						dArea = Polygon_IO_Get_Polygon_Area(polyData,nPolygon);
						if(dArea < 0.0) dArea *= -1.0;
						dSigmaArea += ((dAverageArea - dArea) * (dAverageArea - dArea));
					}
				}
			}
		}
	}
	dSigmaArea /= (double)nPolygons;
	dSigmaArea = sqrt(dSigmaArea);

	printf("  Polygon statistics for region\n");
	printf("    ( %.3f , %.3f ) to ( %.3f , %.3f )\n",dXStart,dYStart,dXStop,dYStop);
	printf("    Number of polygons : %d\n",nPolygons);
	printf("    Total area : %.3f\n",dAverageArea * nPolygons);
	printf("    Average area : %.3f\n",dAverageArea);
	printf("    Sigma area : %.3f\n",dSigmaArea);
	printf("    <D> : %.3f\n",sqrt(dAverageArea));
	printf("    Sigma D : %.3f ( %.3f%% )\n",dSigmaArea / (2.0 * sqrt(dAverageArea)),50.0 * dSigmaArea / dAverageArea);
	printf("    Biggest grain : %.3f, (%d)\n",sqrt(dBiggest),nBiggest);
	printf("    Smallest grain : %.3f, (%d)\n",sqrt(dSmallest),nSmallest);

	return (50.0 * dSigmaArea / dAverageArea);			/* Return the size distribution (%) */
}


short int Polygon_IO_Get_Short_Int(FILE *fp)
{			/* Return a short integer extracted from a binary file */
char s[sizeof(short int)];
int c,n;
short int *i;

	for(n = 0; n < sizeof(short int); n++){
		c = fgetc(fp);
		s[n] = c;
	}
	i = (short int*)&s[0];

	return (*i);
}


short int Polygon_IO_Get_Short_Int_Reverse(FILE *fp)
{			/* Return a short integer extracted from a binary file (reverse byte ordering) */
char s[sizeof(short int)];
int c,n;
short int *i;

	for(n = sizeof(short int) - 1; n >= 0; n--){
		c = fgetc(fp);
		s[n] = c;
	}
	i = (short int*)&s[0];

	return (*i);
}


unsigned int Polygon_IO_Get_State_Bytes(POLYGONATTRIBUTES *polyData)
{			/* Return an integer representing the state of the polygon object */
unsigned int uiState = 0;

	if(polyData->nGroupSet == 1) uiState += 1;
	if(polyData->nMagnetisationSet == 1) uiState += 2;
	if(polyData->nBoundarySet == 1) uiState += 4;
	if(polyData->nThicknessSet == 1) uiState += 8;

	return (uiState);
}


point Polygon_IO_Get_Tetrahedron_Centre(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Return the centroid of a tetrahedron */
point p0,pA,pB,pC,pCentre;

	p0 = Polygon_IO_Get_Vertex_No_Boundary(polyData,nPolygon,0);
	pA = Polygon_IO_Get_Vertex_No_Boundary(polyData,nPolygon,1);
	pB = Polygon_IO_Get_Vertex_No_Boundary(polyData,nPolygon,2);
	pC = Polygon_IO_Get_Vertex_No_Boundary(polyData,nPolygon,3);

	pA.x -= p0.x;
	pA.y -= p0.y;
	pA.z -= p0.z;

	pB.x -= p0.x;
	pB.y -= p0.y;
	pB.z -= p0.z;

	pC.x -= p0.x;
	pC.y -= p0.y;
	pC.z -= p0.z;

	pCentre.x = p0.x + (0.25 * (pA.x + pB.x + pC.x));
	pCentre.y = p0.y + (0.25 * (pA.y + pB.y + pC.y));
	pCentre.z = p0.z + (0.25 * (pA.z + pB.z + pC.z));

	return (pCentre);
}


double Polygon_IO_Get_Total_Polygon_Area(POLYGONATTRIBUTES *polyData)
{			/* Get the total area of all polygons */
int n;
double dA,dTotal = 0.0;

	for(n = 0; n < polyData->nPolygons; n++){
		dA = Polygon_IO_Get_Polygon_Area(polyData,n);
		if(dA < 0.0) dA *= -1.0;
		dTotal += dA;
	}

	return (dTotal);
}


unsigned int Polygon_IO_Get_Unsigned_Int(FILE *fp)
{			/* Return an unsigned integer extracted from a binary file */
char s[sizeof(unsigned int)];
int c,n;
unsigned int *ui;

	for(n = 0; n < sizeof(int); n++){
		c = fgetc(fp);
		s[n] = c;
	}
	ui = (unsigned int*)&s[0];

	return (*ui);
}


unsigned int Polygon_IO_Get_Unsigned_Int_Reverse(FILE *fp)
{			/* Return an unsigned integer extracted from a binary file (reverse byte ordering) */
char s[sizeof(unsigned int)];
int c,n;
unsigned int *ui;

	for(n = sizeof(int) - 1; n >= 0; n--){
		c = fgetc(fp);
		s[n] = c;
	}
	ui = (unsigned int*)&s[0];

	return (*ui);
}


void Polygon_IO_Get_Vectors_From_Vertices(POLYGONATTRIBUTES *polyData,int nPolygon,int nBaseVertex,point *p1,point *p2)
{			/* Get two vectors from p0->p1 and p1->p2 (sides of a polygon) */
		p1->x = polyData->pVertex[nPolygon + ((nBaseVertex + 1) * polyData->nPolygons)].x - polyData->pVertex[nPolygon + (nBaseVertex * polyData->nPolygons)].x;
		p1->y = polyData->pVertex[nPolygon + ((nBaseVertex + 1) * polyData->nPolygons)].y - polyData->pVertex[nPolygon + (nBaseVertex * polyData->nPolygons)].y;
		p2->x = polyData->pVertex[nPolygon + ((nBaseVertex + 2) * polyData->nPolygons)].x - polyData->pVertex[nPolygon + ((nBaseVertex + 1)*polyData->nPolygons)].x;
		p2->y = polyData->pVertex[nPolygon + ((nBaseVertex + 2) * polyData->nPolygons)].y - polyData->pVertex[nPolygon + ((nBaseVertex + 1)*polyData->nPolygons)].y;
}


point Polygon_IO_Get_Vertex(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex)
{			/* Return a vertex of a polygon */
point pV;

	pV.x =  ((double)polyData->fBoundary[nPolygon] * (polyData->pSeed[nPolygon].x - polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x)) + polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x;		/* Vertex after allowing for boundary */
	pV.y =  ((double)polyData->fBoundary[nPolygon] * (polyData->pSeed[nPolygon].y - polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y)) + polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y;
	pV.z =  ((double)polyData->fBoundary[nPolygon] * (polyData->pSeed[nPolygon].z - polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].z)) + polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].z;

//	pV.x =  polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x;		/* Vertex after allowing for boundary */
//	pV.y =  polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y;
//	pV.z =  polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].z;

	return (pV);
}


int Polygon_IO_Get_Vertex_Count(int nPolygon,int nPolygons,point *pVertex)
{			/* Return the number of vertices in nPolygon */
int nVertices = 0;
double dX0,dY0;

	dX0 = pVertex[nPolygon].x;		/* Initial vertex */
	dY0 = pVertex[nPolygon].y;
	do{
		nVertices++;
	}
	while((pVertex[nPolygon + (nVertices * nPolygons)].x != dX0) || (pVertex[nPolygon + (nVertices * nPolygons)].y != dY0));

	return (nVertices);
}


point Polygon_IO_Get_Vertex_No_Boundary(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex)
{			/* Return a vertex of a polygon without applying the boundary factor */
point pV;

	pV.x =  polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x;
	pV.y =  polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y;
	pV.z =  polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].z;

	return (pV);
}


int Polygon_IO_Get_Zero_Area_Polygon_Count(POLYGONATTRIBUTES *polyData)
{			/* Count the number of polygons with zero area - these should be removed */
int nCount = 0,nPolygon;

	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		if(polyData->siPolyVertices[nPolygon] < 3) nCount++;
	}
	printf("    There are %d polygons with zero area\n",nCount);

	return (nCount);
}


void Polygon_IO_Init(POLYGONATTRIBUTES *polyData)
{			/* Initialise a polygon object */
	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);
}


int Polygon_IO_Init_Group(POLYGONATTRIBUTES **polyData)
{			/* Add support for group indexing to polyData */
int nOK = 0;

	if((*polyData)->nGroupSet == 0){
		if((*polyData)->nPolygons > 0){
			NEW_ZERO((*polyData)->nGroup, int, (*polyData)->nPolygons);
			(*polyData)->nGroupSet = 1;
			nOK = 1;
		}
	}

	return (nOK);
}


void Polygon_IO_Init_Group_Default(POLYGONATTRIBUTES *polyData)
{			/* Initialise polyData with default groups (nGroup = nPolygon) */
int n;

	Polygon_IO_Init_Group(&polyData);
	for(n = 0; n < polyData->nPolygons; n++){
		polyData->nGroup[n] = n;
	}
}


int Polygon_IO_Is_Cubic(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Return 1 if a polygon is cubic/cuboid, 0 otherwise */
int nCubic = 0;
double dXDiff = -1.0,dYDiff = -1.0,dXSame,dYSame;

	if(nPolygon < polyData->nPolygons){
		if(polyData->siPolyVertices[nPolygon] == 4){
			dXSame = polyData->pVertex[nPolygon].x;
			dYSame = polyData->pVertex[nPolygon].y;

			if(polyData->pVertex[nPolygon + (polyData->nPolygons)].x == dXSame){		/* Point 1 has same x -> different y */
				dYDiff = polyData->pVertex[nPolygon + (polyData->nPolygons)].y;
			}

			if(polyData->pVertex[nPolygon + (polyData->nPolygons)].y == dYSame){		/* Point 1 has same y -> different x */
				dXDiff = polyData->pVertex[nPolygon + (polyData->nPolygons)].x;
			}

			if(polyData->pVertex[nPolygon + (3 * polyData->nPolygons)].x == dXSame){		/* Point 3 has same x -> different y */
				dYDiff = polyData->pVertex[nPolygon + (3 * polyData->nPolygons)].y;
			}

			if(polyData->pVertex[nPolygon + (3 * polyData->nPolygons)].y == dYSame){		/* Point 3 has same y -> different x */
				dXDiff = polyData->pVertex[nPolygon + (3 * polyData->nPolygons)].x;
			}

			if(polyData->pVertex[nPolygon + (2 * polyData->nPolygons)].x == dXDiff){
				if(polyData->pVertex[nPolygon + (2 * polyData->nPolygons)].y == dYDiff){
					nCubic = 1;
				}
			}
		}
	}

	return (nCubic);
}


int Polygon_IO_Is_Point_In_Polygon(POLYGONATTRIBUTES *polyData,int nPolygon,double dX,double dY)
{			/* Return 1 if (dX,dY) is inside polygon nPolygon */
int nI, nJ, nIn = 0, nV = polyData->siPolyVertices[nPolygon];

	for(nI = 0, nJ = nV - 1; nI < nV; nJ = nI++){
		if((((polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y <= dY) &&
		(dY < polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y)) ||
		((polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y <= dY) &&
		(dY < polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y))) &&
		(dX < (polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].x - polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].x) * (dY - polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y) / (polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y - polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y) + polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].x)) nIn = !nIn;
	}

	return (nIn);
}


int Polygon_IO_Is_Point_In_Polygon_With_Boundary(POLYGONATTRIBUTES *polyData,int nPolygon,double dX,double dY)
{			/* Return 1 if (dX,dY) is inside polygon nPolygon with boundary */
int nI, nJ, nIn = 0, nV = polyData->siPolyVertices[nPolygon];
point nVI, nVJ;
//FILE *fp;

// point Polygon_IO_Get_Vertex(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex)
//	printf("Polygon_IO_Is_Point_In_Polygon_With_Boundary\n");
	for(nI = 0, nJ = nV - 1; nI < nV; nJ = nI++){

		nVI = Polygon_IO_Get_Vertex(polyData, nPolygon, nI);
		nVJ = Polygon_IO_Get_Vertex(polyData, nPolygon, nJ);
		
//		fp = fopen("Polygon_IO_Is_Point_In_Polygon_With_Boundary.txt","at");
//		fprintf(fp, "nI %d %lf %lf %lf w/b %lf %lf %lf nJ %d %lf %lf %lf w/b %lf %lf %lf\n",nI, polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].x, polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y, polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].z, nVI.x, nVI.y, nVI.z, nJ, polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].x, polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y, polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].z, nVJ.x, nVJ.y, nVJ.z);
//		fprintf(fp, "nI %d nJ %d\n",nI, nJ);
//		fclose(fp);
		
//		if((((polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y <= dY) &&
//		(dY < polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y)) ||
//		((polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y <= dY) &&
//		(dY < polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y))) &&
//		(dX < (polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].x - polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].x) * (dY - polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y) / (polyData->pVertex[nPolygon + (nJ * polyData->nPolygons)].y - polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].y) + polyData->pVertex[nPolygon + (nI * polyData->nPolygons)].x)) nIn = !nIn;

		if((((nVI.y <= dY) &&
		(dY < nVJ.y)) ||
		((nVJ.y <= dY) &&
		(dY < nVI.y))) &&
		(dX < (nVJ.x - nVI.x) * (dY - nVI.y) / (nVJ.y - nVI.y) + nVI.x)) nIn = !nIn;

	
	
	}

	return (nIn);
}


int Polygon_IO_Is_Polygon_In_Region(POLYGONATTRIBUTES *polyData,int nPolygon,double dXStart,double dXStop,double dYStart,double dYStop,double dTolerance)
{			/* Determine if a polygon (vertices) is within a region (dXStart,dYStart) to (dXStop, dYStop) */
int nIn = 0;
double dXMin,dXMax,dYMin,dYMax;
point p = Polygon_IO_Get_Vertex(polyData,nPolygon,0);

	if(polyData->pSeed[nPolygon].x > dXStart){
		if(polyData->pSeed[nPolygon].x < dXStop){
			if(polyData->pSeed[nPolygon].y > dYStart){
				if(polyData->pSeed[nPolygon].y < dYStop){
					Polygon_IO_Get_Polygon_Limits(polyData,nPolygon,&dXMin,&dXMax,&dYMin,&dYMax);		/* Get limits of polygon */
					if(dXMin > (dXStart - dTolerance)){
						if(dXMax < (dXStop + dTolerance)){
							if(dYMin > (dYStart - dTolerance)){
								if(dYMax < (dYStop + dTolerance)){
									nIn = 1;
								}
							}
						}
					}
				}
			}
		}
	}

	return (nIn);
}


void Polygon_IO_Load_Magnetisation(const char *cFilename,POLYGONATTRIBUTES *polyData, int nLayerCells)
{			/* Load the magnetisation data into the map */
FILE *fp;
int nPolygon,nZ;
double dMx,dMy,dMz,dDuff;

	Polygon_IO_Setup_Polygon_Attributes_M(&polyData);		/* Allocate memory for M */

	std::uintmax_t size = std::filesystem::file_size(cFilename);

	char *psz = new char[size + 1];
	if (nullptr == psz)
		return;

	fp = fopen(cFilename, "rb");
	fread(psz, size, 1, fp);
	fclose(fp);

	int nIndex = 0;
	nPolygon=0;
	char *p = strtok(psz, "\r\n");
	while (NULL != p)
	{
		sscanf(p, "%lf %lf %d %lf %lf %lf",&dDuff,&dDuff,&nZ,&dMx,&dMy,&dMz);
		if(nPolygon<polyData->nPolygons)
		{
			nIndex = (nPolygon % nLayerCells) + (nPolygon / nLayerCells * nLayerCells);
			polyData->pM[nIndex].x = dMx;
			polyData->pM[nIndex].y = dMy;
			polyData->pM[nIndex].z = dMz;
			nPolygon++;
		}

		p = strtok(NULL, "\r\n");
	}

	delete [] psz;
	psz = nullptr;
}


void Polygon_IO_Load_Mesh(POLYGONATTRIBUTES *polyData,char *cFilename)
{			/* Load a neutral format tetrahedronal mesh */
FILE *fp;
int n,n0,n1,n2,n3,n4,nPoints,nTets,nVertex,nRank = 0;
int *nTet;
double dX,dY,dZ;
point pSize,*p;

	fp = fopen(cFilename,"r");
	fscanf(fp,"%d",&nPoints);
	printf("    %d points\n",nPoints);
	NEW_ZERO(p, point, nPoints);				/* Array to hold vertices of tetrahedrons */

	for(n = 0; n < nPoints; n++){			/* Load points */
		fscanf(fp,"%lf %lf %lf",&dX,&dY,&dZ);
		p[n].x = dX;
		p[n].y = dY;
		p[n].z = dZ;
	}

	fscanf(fp,"%d",&nTets);
	printf("    %d tetrahedrons\n",nTets);
	NEW_ZERO(nTet, int, (nTets * 4));			/* Array to hold list of vertices in each tetrahedron */

	for(n = 0; n < nTets; n++){				/* Load list of tetrahedrons */
		fscanf(fp,"%d %d %d %d %d\n",&n0,&n1,&n2,&n3,&n4);
		nTet[(n * 4)] = n1 - 1;				/* Subtract 1 because points are numbered from 1 in neutral file */
		nTet[(n * 4) + 1] = n2 - 1;
		nTet[(n * 4) + 2] = n3 - 1;
		nTet[(n * 4) + 3] = n4 - 1;
	}
	fclose(fp);

	Polygon_IO_Reset_Polygon_Attributes(polyData);
	polyData->nPolygons = nTets;
	polyData->nVertexMax = 4;
	polyData->nBoundarySet = 1;
	polyData->nTetrahedronSet = 1;
	pSize = p[0];
	for(n = 1; n < nPoints; n++){			/* Get maximum size of object */
		if(p[n].x > pSize.x) pSize.x = p[n].x;
		if(p[n].x > pSize.y) pSize.y = p[n].x;
		if(p[n].x > pSize.z) pSize.z = p[n].x;
	}
	polyData->pSize = pSize;

//	polyData->dSeedToVertexMax2 = 1.0;
	Polygon_IO_Setup_Polygon_Attributes(&polyData,polyData->nGroupSet,polyData->nBoundarySet,polyData->nThicknessSet,polyData->nMagnetisationSet);	/* Dimension the arrays in polyData */
	for(n = 0; n < polyData->nPolygons; n++){
		polyData->siPolyVertices[n] = 4;		/* Set no. of vertices */
		polyData->nPolyLayer[n] = 3;							/* Set layer */
		for(nVertex = 0; nVertex < 4; nVertex++){		/* Set vertices */
			polyData->pVertex[n + (nVertex * polyData->nPolygons)].x = p[ nTet[(n * 4) + nVertex] ].x;
			polyData->pVertex[n + (nVertex * polyData->nPolygons)].y = p[ nTet[(n * 4) + nVertex] ].y;
			polyData->pVertex[n + (nVertex * polyData->nPolygons)].z = p[ nTet[(n * 4) + nVertex] ].z;
		}
		polyData->pSeed[n] = Polygon_IO_Get_Tetrahedron_Centre(polyData,n);
	}

	SAFE_DELETES(p);
	SAFE_DELETES(nTet);
}


void Polygon_IO_Load_P3D(POLYGONATTRIBUTES *polyData,char *cFilename)
{			/* Load a .p3d file */
	Polygon_IO_Load_P3D_Header(cFilename,polyData);
	Polygon_IO_Load_P3D_Polygons(cFilename,polyData);
}


void Polygon_IO_Load_P3D_Header(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Load the header information from a .p3d file */
FILE *fp;
int nPolygons,nVertexMax;
double dPolyXMax,dPolyYMax,dSeedToVertexMax2;

	Polygon_IO_Reset_Polygon_Attributes(polyData);
	fp = fopen(cFilename,"r");
	fscanf(fp,"%d %d %lf %lf %lf",&nPolygons,&nVertexMax,&dPolyXMax,&dPolyYMax,&dSeedToVertexMax2);
	fclose(fp);

	polyData->nPolygons = nPolygons;
	polyData->nVertexMax = nVertexMax;
	polyData->pSize.x = dPolyXMax;
	polyData->pSize.y = dPolyYMax;
	polyData->dSeedToVertexMax2 = dSeedToVertexMax2;
}


void Polygon_IO_Load_P3D_Header_Binary(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Load the header from a .bp3d file */
FILE *fp;
int nByteOrder,nPolygons,nVertexMax;
double dPolyXMax,dPolyYMax,dSeedToVertexMax2;

	Polygon_IO_Reset_Polygon_Attributes(polyData);
	nByteOrder = Polygon_IO_Get_Byte_Order(cFilename);
	fp = fopen(cFilename,"r");
	switch(nByteOrder){			/* Byte order depends on the machine that saved the file */
	case 0:
		nPolygons = Polygon_IO_Get_Int(fp);
		nVertexMax = Polygon_IO_Get_Int(fp);
		dPolyXMax = Polygon_IO_Get_Float(fp);
		dPolyYMax = Polygon_IO_Get_Float(fp);
		dSeedToVertexMax2 = Polygon_IO_Get_Float(fp);
		break;
	case 1:
		nPolygons = Polygon_IO_Get_Int_Reverse(fp);
		nVertexMax = Polygon_IO_Get_Int_Reverse(fp);
		dPolyXMax = Polygon_IO_Get_Float_Reverse(fp);
		dPolyYMax = Polygon_IO_Get_Float_Reverse(fp);
		dSeedToVertexMax2 = Polygon_IO_Get_Float_Reverse(fp);
		break;
	}
	fclose(fp);

	polyData->nPolygons = nPolygons;
	polyData->nVertexMax = nVertexMax;
	polyData->pSize.x = dPolyXMax;
	polyData->pSize.y = dPolyYMax;
	polyData->dSeedToVertexMax2 = dSeedToVertexMax2;
}


void Polygon_IO_Load_P3D_Polygons(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Load the polygon data from a .p3d file */
FILE *fp;
int n,nL,nPolygon,nVertex,nVertices;
double d,dX,dY;

	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);	/* Dimension the arrays in polyData */
	fp = fopen(cFilename,"r");
	fscanf(fp,"%d %d %lf %lf %lf",&n,&n,&d,&d,&d);		/* Skip the header (it should have been read separately) */
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		polyData->fBoundary[nPolygon] = 0.0;		/* This isn't saved */
		fscanf(fp,"%d %d",&nL,&nVertices);	/* Load layer number and no. of vertices */
		polyData->siPolyVertices[nPolygon] = nVertices;		/* Set no. of vertices */
		polyData->nPolyLayer[nPolygon] = nL;							/* Set layer */
		fscanf(fp,"%lf %lf",&dX,&dY);				/* Load seed point */
		polyData->pSeed[nPolygon].x = dX;
		polyData->pSeed[nPolygon].y = dY;
		for(nVertex = 0; nVertex < nVertices; nVertex++){
			fscanf(fp,"%lf %lf",&dX,&dY);
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x = dX;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y = dY;
		}
	}
	fclose(fp);
}


void Polygon_IO_Load_P3D_Polygons_Binary(char *cFilename, POLYGONATTRIBUTES *polyData)
{			/* Load the polygon data from a .bp3d file */
FILE *fp;
int n,nByteOrder,nL,nPolygon,nVertex,nVertices;
double dX,dY;

	Polygon_IO_Setup_Polygon_Attributes(&polyData,0,1,0,0);			/* Dimension the arrays in polyData */
	nByteOrder = Polygon_IO_Get_Byte_Order(cFilename);
	if(nByteOrder) Polygon_IO_Load_P3D_Polygons_Binary_Reverse(cFilename,polyData);
	else{
		fp = fopen(cFilename,"r");
		for(n = 0; n < 20; n++) fgetc(fp);			/* Skip the header (5x4 bytes) */
		for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
			polyData->fBoundary[nPolygon] = 0.0;			/* This isn't saved */
			nL = Polygon_IO_Get_Int(fp);
			nVertices = Polygon_IO_Get_Int(fp);
			polyData->siPolyVertices[nPolygon] = nVertices;	/* Set no. of vertices */
			polyData->nPolyLayer[nPolygon] = nL;						/* Set layer */
			dX = Polygon_IO_Get_Float(fp);
			dY = Polygon_IO_Get_Float(fp);
			polyData->pSeed[nPolygon].x = dX;
			polyData->pSeed[nPolygon].y = dY;
			for(nVertex = 0; nVertex < nVertices; nVertex++){
				dX = Polygon_IO_Get_Float(fp);
				dY = Polygon_IO_Get_Float(fp);
				polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x = dX;
				polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y = dY;
			}
		}
		fclose(fp);
	}
}


void Polygon_IO_Load_P3D_Polygons_Binary_Reverse(char *cFilename, POLYGONATTRIBUTES *polyData)
{			/* Load the polygon data from a .bp3d file */
FILE *fp;
int n,nL,nPolygon,nVertex,nVertices;
double dX,dY;

	fp = fopen(cFilename,"r");
	for(n = 0; n < 20; n++) fgetc(fp);			/* Skip the header (5x4 bytes) */
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		polyData->fBoundary[nPolygon] = 0.0;		/* This isn't saved */
		nL = Polygon_IO_Get_Int_Reverse(fp);
		nVertices = Polygon_IO_Get_Int_Reverse(fp);
		polyData->siPolyVertices[nPolygon] = nVertices;	/* Set no. of vertices */
		polyData->nPolyLayer[nPolygon] = nL;						/* Set layer */
		dX = Polygon_IO_Get_Float_Reverse(fp);
		dY = Polygon_IO_Get_Float_Reverse(fp);
		polyData->pSeed[nPolygon].x = dX;
		polyData->pSeed[nPolygon].y = dY;
		for(nVertex = 0; nVertex < nVertices; nVertex++){
			dX = Polygon_IO_Get_Float_Reverse(fp);
			dY = Polygon_IO_Get_Float_Reverse(fp);
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x = dX;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y = dY;
		}
	}
	fclose(fp);
}


void Polygon_IO_Load_Pol(POLYGONATTRIBUTES *polyData,char *cFilename)
{			/* Load a .pol file */
	Polygon_IO_Load_Pol_Header(cFilename,polyData);
	Polygon_IO_Load_Pol_Polygons(cFilename,polyData);
}


void Polygon_IO_Load_Pol_Binary(POLYGONATTRIBUTES *polyData,const char *cFilename)
{			/* Load a .bpol file */

	Polygon_IO_Load_Pol_Header_Binary(cFilename,polyData);
	Polygon_IO_Load_Pol_Polygons_Binary(cFilename,polyData);
}


void Polygon_IO_Load_Pol_Header(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Load the header information from a .pol file */
FILE *fp;
int nPolygons,nVertexMax;
unsigned int uiState;
double dPolyXMax,dPolyYMax,dSeedToVertexMax2;

	Polygon_IO_Reset_Polygon_Attributes(polyData);
	fp = fopen(cFilename,"r");
	fscanf(fp,"%d %d %u %lf %lf %lf",&nPolygons,&nVertexMax,&uiState,&dPolyXMax,&dPolyYMax,&dSeedToVertexMax2);
	fclose(fp);

	polyData->nPolygons = nPolygons;
	polyData->nVertexMax = nVertexMax;
	if(uiState & 1) polyData->nGroupSet = 1;
	if(uiState & 2) polyData->nMagnetisationSet = 1;
	if(uiState & 4) polyData->nBoundarySet = 1;
	if(uiState & 8) polyData->nThicknessSet = 1;
	if(uiState & 16) polyData->nTetrahedronSet = 1;

	polyData->pSize.x = dPolyXMax;
	polyData->pSize.y = dPolyYMax;
	polyData->dSeedToVertexMax2 = dSeedToVertexMax2;
}


int Polygon_IO_Load_Pol_Header_Binary(const char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Load the header from a .bpol file */
FILE *fp;
int n123,nByteOrder,nPolygons,nVertexMax;
unsigned int uiState;
double dPolyXMax,dPolyYMax,dSeedToVertexMax2;

	Polygon_IO_Reset_Polygon_Attributes(polyData);
	fp = fopen(cFilename,"rb");

	n123 = Polygon_IO_Get_Int(fp);
	if(n123 == 123) nByteOrder = 0;		/* Select byte order based on value of n123 */
	else nByteOrder = 1;

	switch(nByteOrder){			/* Byte order depends on the machine that saved the file */
	case 0:
		nPolygons = Polygon_IO_Get_Int(fp);
		nVertexMax = Polygon_IO_Get_Int(fp);
		uiState = Polygon_IO_Get_Unsigned_Int(fp);
		dPolyXMax = Polygon_IO_Get_Float(fp);
		dPolyYMax = Polygon_IO_Get_Float(fp);
		dSeedToVertexMax2 = Polygon_IO_Get_Float(fp);
		break;
	case 1:
		nPolygons = Polygon_IO_Get_Int_Reverse(fp);
		nVertexMax = Polygon_IO_Get_Int_Reverse(fp);
		uiState = Polygon_IO_Get_Unsigned_Int_Reverse(fp);
		dPolyXMax = Polygon_IO_Get_Float_Reverse(fp);
		dPolyYMax = Polygon_IO_Get_Float_Reverse(fp);
		dSeedToVertexMax2 = Polygon_IO_Get_Float_Reverse(fp);
		break;
	}
	fclose(fp);

	polyData->nPolygons = nPolygons;
	polyData->nVertexMax = nVertexMax;
	if(uiState & 1) polyData->nGroupSet = 1;
	if(uiState & 2) polyData->nMagnetisationSet = 1;
	if(uiState & 4) polyData->nBoundarySet = 1;
	if(uiState & 8) polyData->nThicknessSet = 1;
	if(uiState & 16) polyData->nTetrahedronSet = 1;
	polyData->pSize.x = dPolyXMax;
	polyData->pSize.y = dPolyYMax;
	polyData->dSeedToVertexMax2 = dSeedToVertexMax2;

	return (nByteOrder);
}


void Polygon_IO_Load_Pol_Polygons(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Load the polygon data from a .pol file */
FILE *fp;
int n,nGroup = 0,nL,nPolygon,nVertex,nVertices;
unsigned int uiState;
float fBoundary = 0.0;
double d,dX,dY,dZ;

	Polygon_IO_Setup_Polygon_Attributes(&polyData,polyData->nGroupSet,polyData->nBoundarySet,polyData->nThicknessSet,polyData->nMagnetisationSet);	/* Dimension the arrays in polyData */
	fp = fopen(cFilename,"r");
	fscanf(fp,"%d %d %u %lf %lf %lf",&n,&n,&uiState,&d,&d,&d);		/* Skip the header (it should have been read separately) */
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		fscanf(fp,"%d %d",&nL,&nVertices);	/* Load layer number, no. of vertices */
		polyData->siPolyVertices[nPolygon] = (short int)nVertices;		/* Set no. of vertices */
		polyData->nPolyLayer[nPolygon] = nL;							/* Set layer */

		if(polyData->nGroupSet == 1){
			fscanf(fp,"%d",&nGroup);
			polyData->nGroup[nPolygon] = nGroup;		/* Set group */
		}

		if(polyData->nBoundarySet == 1){
			fscanf(fp,"%f",&fBoundary);
			polyData->fBoundary[nPolygon] = fBoundary;		/* Set boundary */
		}

		if(polyData->nThicknessSet == 1){
			fscanf(fp,"%lf",&d);
			polyData->dThickness[nPolygon] = d;		/* Set thickness */
		}

		if(polyData->nMagnetisationSet == 1){
			fscanf(fp,"%lf %lf %lf",&dX,&dY,&dZ);
			polyData->pM[nPolygon].x = dX;		/* Set magnetisation */
			polyData->pM[nPolygon].y = dY;
			polyData->pM[nPolygon].z = dZ;
		}

		fscanf(fp,"%lf %lf %lf",&dX,&dY,&dZ);				/* Load seed point */
		polyData->pSeed[nPolygon].x = dX;
		polyData->pSeed[nPolygon].y = dY;
		polyData->pSeed[nPolygon].z = dZ;
		for(nVertex = 0; nVertex < nVertices; nVertex++){		/* Load vertices */
			if(polyData->nTetrahedronSet == 0){
				fscanf(fp,"%lf %lf",&dX,&dY);
				dZ = 0.0;
			}
			else{
				fscanf(fp,"%lf %lf %lf",&dX,&dY,&dZ);
			}
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x = dX;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y = dY;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].z = dZ;
		}
	}
	fclose(fp);
}


void Polygon_IO_Load_Pol_Polygons_Binary(const char *cFilename, POLYGONATTRIBUTES *polyData)
{			/* Load the polygon data from a .bpol file */
FILE *fp;
int n,nGroup,nL,nPolygon,nSkipBytes,nVertex,nVertices;
float fBoundary;
double dX,dY,dZ;

	Polygon_IO_Setup_Polygon_Attributes(&polyData,polyData->nGroupSet,polyData->nBoundarySet,polyData->nThicknessSet,polyData->nMagnetisationSet);			/* Dimension the arrays in polyData */
//	nSkipBytes = (3 * sizeof(int)) + sizeof(unsigned int) + (3 * sizeof(float));
	nSkipBytes = (3 * 4) + 4 + (3 * 4);
	fp = fopen(cFilename,"rb");

	for(n = 0; n < nSkipBytes; n++) fgetc(fp);			/* Skip the header (4*4 + 3*8 bytes) */
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		nL = Polygon_IO_Get_Int(fp);
		nVertices = Polygon_IO_Get_Int(fp);

		polyData->nPolyLayer[nPolygon] = nL;						/* Set layer */
		polyData->siPolyVertices[nPolygon] = (short int)nVertices;	/* Set no. of vertices */

		if(polyData->nGroupSet == 1){
			nGroup = Polygon_IO_Get_Int(fp);
			polyData->nGroup[nPolygon] = nGroup;		/* Set group */
		}

		if(polyData->nBoundarySet == 1){
			fBoundary = Polygon_IO_Get_Float(fp);
			polyData->fBoundary[nPolygon] = fBoundary;		/* Set boundary */
		}

		if(polyData->nThicknessSet == 1){
			dZ = (double)Polygon_IO_Get_Float(fp);
			polyData->dThickness[nPolygon] = dZ;		/* Set thickness */
		}

		if(polyData->nMagnetisationSet == 1){
			dX = (double)Polygon_IO_Get_Float(fp);
			dY = (double)Polygon_IO_Get_Float(fp);
			dZ = (double)Polygon_IO_Get_Float(fp);
			polyData->pM[nPolygon].x = dX;		/* Set magnetisation */
			polyData->pM[nPolygon].y = dY;
			polyData->pM[nPolygon].z = dZ;
		}

		dX = (double)Polygon_IO_Get_Float(fp);
		dY = (double)Polygon_IO_Get_Float(fp);
		dZ = (double)Polygon_IO_Get_Float(fp);
		polyData->pSeed[nPolygon].x = dX;
		polyData->pSeed[nPolygon].y = dY;
		polyData->pSeed[nPolygon].z = dZ;
		for(nVertex = 0; nVertex < nVertices; nVertex++){
			if(polyData->nTetrahedronSet == 0){
				dX = (double)Polygon_IO_Get_Float(fp);
				dY = (double)Polygon_IO_Get_Float(fp);
				dZ = 0.0;
			}
			else{
				dX = (double)Polygon_IO_Get_Float(fp);
				dY = (double)Polygon_IO_Get_Float(fp);
				dZ = (double)Polygon_IO_Get_Float(fp);
			}
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x = dX;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y = dY;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].z = dZ;
		}
	}
	fclose(fp);
}



void Polygon_IO_Load_Pol_Polygons_Binary_Reverse(char *cFilename, POLYGONATTRIBUTES *polyData)
{			/* Load the polygon data from a .bpol file */
FILE *fp;
int n,nGroup,nL,nPolygon,nSkipBytes,nVertex,nVertices;
float fBoundary;
double dX,dY,dZ;

	Polygon_IO_Setup_Polygon_Attributes(&polyData,polyData->nGroupSet,polyData->nBoundarySet,polyData->nThicknessSet,polyData->nMagnetisationSet);			/* Dimension the arrays in polyData */
	nSkipBytes = (3 * sizeof(int)) + sizeof(unsigned int) + (3 * sizeof(float));
	fp = fopen(cFilename,"r");

	for(n = 0; n < nSkipBytes; n++) fgetc(fp);			/* Skip the header (4*4 + 3*8 bytes) */
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		nL = Polygon_IO_Get_Int_Reverse(fp);
		nVertices = Polygon_IO_Get_Int_Reverse(fp);
		polyData->nPolyLayer[nPolygon] = nL;						/* Set layer */
		polyData->siPolyVertices[nPolygon] = (short int)nVertices;	/* Set no. of vertices */

		if(polyData->nGroupSet == 1){
			nGroup = Polygon_IO_Get_Int_Reverse(fp);
			polyData->nGroup[nPolygon] = nGroup;		/* Set group */
		}

		if(polyData->nBoundarySet == 1){
			fBoundary = Polygon_IO_Get_Float_Reverse(fp);
			polyData->fBoundary[nPolygon] = fBoundary;		/* Set boundary */
		}

		if(polyData->nThicknessSet == 1){
			dZ = (double)Polygon_IO_Get_Float_Reverse(fp);
			polyData->dThickness[nPolygon] = dZ;		/* Set thickness */
		}

		if(polyData->nMagnetisationSet == 1){
			dX = (double)Polygon_IO_Get_Float_Reverse(fp);
			dY = (double)Polygon_IO_Get_Float_Reverse(fp);
			dZ = (double)Polygon_IO_Get_Float_Reverse(fp);
			polyData->pM[nPolygon].x = dX;		/* Set magnetisation */
			polyData->pM[nPolygon].y = dY;
			polyData->pM[nPolygon].z = dZ;
		}

		dX = (double)Polygon_IO_Get_Float_Reverse(fp);
		dY = (double)Polygon_IO_Get_Float_Reverse(fp);
		dZ = (double)Polygon_IO_Get_Float_Reverse(fp);
		polyData->pSeed[nPolygon].x = dX;
		polyData->pSeed[nPolygon].y = dY;
		polyData->pSeed[nPolygon].z = dZ;
		for(nVertex = 0; nVertex < nVertices; nVertex++){
			dX = (double)Polygon_IO_Get_Float_Reverse(fp);
			dY = (double)Polygon_IO_Get_Float_Reverse(fp);
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].x = dX;
			polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)].y = dY;
		}
	}
	fclose(fp);
}


void Polygon_IO_Load_Pol_Polygons_Old(char *cFilename,int nPolygons,point *PolySeed,point *PolyVertex)
{			/* Load a *.pol polygon data file (old format) */
FILE *fp;
int nPolygon,nVertex,nLoadCentre;
double dX,dY;

	nPolygon = 0;
	fp = fopen(cFilename,"r");
	nVertex = 0;
	nLoadCentre = 1;
  
	while(fscanf(fp,"%lf %lf",&dX,&dY)!=EOF){
		switch(nLoadCentre){
		case 0:
			if(dX == -1000.0 && dY == -1000.0){		/* This is the polygon terminator, (-1000,-1000) */
				nLoadCentre = 1;
				nPolygon++;
			}
			else{				/* This is a vertex */
				if(nVertex != 0){		/* Stop consecutive occurences of the same vertex */
					if(dX == PolyVertex[nPolygon + ((nVertex - 1) * nPolygons)].x && dY == PolyVertex[nPolygon + ((nVertex - 1) * nPolygons)].y) nVertex--;
				}
				PolyVertex[nPolygon + (nVertex * nPolygons)].x = dX;
				PolyVertex[nPolygon + (nVertex * nPolygons)].y = dY;
				nVertex++;
			}
			break;

		case 1:
			if(dX != -1000.0 && dY != -1000.0){			/* This is the polygon seed point */
				PolySeed[nPolygon].x = dX;
				PolySeed[nPolygon].y = dY;
				nLoadCentre = 0;
			}
			nVertex = 0;
			break;
		}
	}
	fclose(fp);
}


void Polygon_IO_Merge_Polygons(POLYGONATTRIBUTES *polyA,POLYGONATTRIBUTES *polyB,POLYGONATTRIBUTES *polyNew)
{			/* Merge polyA and polyB together */
int nGroupSet = 0;

	Polygon_IO_End(polyNew);
	polyNew->nPolygons = polyA->nPolygons + polyB->nPolygons;		/* Total number of polygons */
	if(polyA->nVertexMax > polyB->nVertexMax) polyNew->nVertexMax = polyA->nVertexMax;		/* Maximum number of vertices for a single polygon */
	else polyNew->nVertexMax = polyB->nVertexMax;
	if(polyA->pSize.x > polyB->pSize.x) polyNew->pSize.x = polyA->pSize.x;				/* Maximum extent of polygon region (A) */
	else polyNew->pSize.x = polyB->pSize.x;
	if(polyA->pSize.y > polyB->pSize.y) polyNew->pSize.y = polyA->pSize.y;				/* Maximum extent of polygon region (A) */
	else polyNew->pSize.y = polyB->pSize.y;
	if(polyA->dSeedToVertexMax2 > polyB->dSeedToVertexMax2) polyNew->dSeedToVertexMax2 = polyA->dSeedToVertexMax2;		/* Maximum distance from seed point to vertex (squared) */
	else polyNew->dSeedToVertexMax2 = polyB->dSeedToVertexMax2;
	if(polyA->nGroupSet == 1 || polyB->nGroupSet == 1) nGroupSet = 1;
	Polygon_IO_Setup_Polygon_Attributes(&polyNew,nGroupSet,1,0,0);
	Polygon_IO_Copy_Polygons(polyA,polyNew,0);			/* Copy polygons from polyA to polyNew */
	Polygon_IO_Copy_Polygons(polyB,polyNew,polyA->nPolygons);	/* Add polyB to polyNew; */
	Polygon_IO_Update_Header(polyNew);			/* Update the header for the new polygon set */
}


void Polygon_IO_Print_Polygon_Attributes(POLYGONATTRIBUTES *polyData)
{			/* Print the details of a POLYGONATTRIBUTES structure */
	printf("    Number of polygons : %d\n",polyData->nPolygons);
	printf("    Maximum number of vertices : %d\n",polyData->nVertexMax);
	printf("    Maximum seed to vertex distance : %lf\n",sqrt(polyData->dSeedToVertexMax2));
	printf("    Polygon region dimensions ( %.3f , %.3f )\n",polyData->pSize.x,polyData->pSize.y);
	printf("    Group set %d\n",polyData->nGroupSet);
	printf("    Boundary set %d\n",polyData->nBoundarySet);
	printf("    Thickness set %d\n",polyData->nThicknessSet);
	printf("    Magnetisation set %d\n",polyData->nMagnetisationSet);
}


void Polygon_IO_Print_Polygon_Details(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Display the details of a polygon */
int n;

	if(nPolygon < polyData->nPolygons){
		printf("    Details of polygon %d\n\n",nPolygon);
		if(polyData->nTetrahedronSet == 0){
			printf("    Seed point ( %.3f , %.3f )\n",polyData->pSeed[nPolygon].x,polyData->pSeed[nPolygon].y);
		}
		else{
			printf("    Seed point ( %.3f , %.3f , %.3f )\n",polyData->pSeed[nPolygon].x,polyData->pSeed[nPolygon].y,polyData->pSeed[nPolygon].z);
		}
		printf("    Layer %d\n",polyData->nPolyLayer[nPolygon]);
		if(polyData->nGroupSet == 1) printf("    Group %d\n",polyData->nGroup[nPolygon]);
		printf("    Area %.1f\n",Polygon_IO_Get_Polygon_Area(polyData,nPolygon));
		printf("    Boundary %.3f\n",polyData->fBoundary[nPolygon]);
		printf("    Number of vertices %d\n",polyData->siPolyVertices[nPolygon]);
		printf("    Vertex  (x,y)\n");
		for(n = 0; n < polyData->siPolyVertices[nPolygon]; n++){
			if(polyData->nTetrahedronSet == 0){
				printf("    %d  ( %.3f , %.3f )\n",n,polyData->pVertex[nPolygon + (n * polyData->nPolygons)].x,polyData->pVertex[nPolygon + (n * polyData->nPolygons)].y);
			}
			else{
				printf("    %d  ( %.3f , %.3f , %.3f )\n",n,polyData->pVertex[nPolygon + (n * polyData->nPolygons)].x,polyData->pVertex[nPolygon + (n * polyData->nPolygons)].y,polyData->pVertex[nPolygon + (n * polyData->nPolygons)].z);
			}
		}
		printf("\n");
	}
}


void Polygon_IO_Process_Polygons(POLYGONATTRIBUTES *polyOld,POLYGONATTRIBUTES *polyNew)
{			/* Do some processing on polyOld, save the result in polyNew */
int nCount = 0,nPolygon;
double dSeedToVertex,dXMax,dYMax;

	Polygon_IO_Reset_Polygon_Attributes(polyNew);
	for(nPolygon = 0; nPolygon < polyOld->nPolygons; nPolygon++){		/* First set up the new header */
		if(Polygon_IO_Process_Polygons_Rule(polyOld,nPolygon) == 1){
			Polygon_IO_Get_Polygon_Extent(polyOld,nPolygon,&dSeedToVertex,&dXMax,&dYMax);
			if(polyOld->siPolyVertices[nPolygon] > polyNew->nVertexMax) polyNew->nVertexMax = polyOld->siPolyVertices[nPolygon];	/* Set maximum number of vertices */
			if(dXMax > polyNew->pSize.x) polyNew->pSize.x = dXMax;
			if(dYMax > polyNew->pSize.y) polyNew->pSize.y = dYMax;
			if(dSeedToVertex > polyNew->dSeedToVertexMax2) polyNew->dSeedToVertexMax2 = dSeedToVertex;
			nCount++;
		}
	}
	polyNew->nPolygons = nCount;		/* Set number of polygons */
	Polygon_IO_Setup_Polygon_Attributes(&polyNew,0,1,0,0);

	nCount = 0;
	for(nPolygon = 0; nPolygon < polyOld->nPolygons; nPolygon++){		/* Now copy the polygons */
		if(Polygon_IO_Process_Polygons_Rule(polyOld,nPolygon) == 1){
			Polygon_IO_Copy_Polygon(polyOld,nPolygon,polyNew,nCount);
			nCount++;
		}
	}
}


int Polygon_IO_Process_Polygons_Rule(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Rules to use when processing polygons */
int nOK = 0;

	/* Rules for inclusion */
	if(polyData->pSeed[nPolygon].x > 250.0 && polyData->pSeed[nPolygon].x < 750.0) nOK = 1;
	if(polyData->pSeed[nPolygon].x > 1000.0 && polyData->pSeed[nPolygon].x < 1500.0) nOK = 1;
	if(polyData->pSeed[nPolygon].x > 1750.0 && polyData->pSeed[nPolygon].x < 2250.0) nOK = 1;
	if(polyData->pSeed[nPolygon].x > 2500.0 && polyData->pSeed[nPolygon].x < 3000.0) nOK = 1;
	if(polyData->pSeed[nPolygon].x > 3250.0 && polyData->pSeed[nPolygon].x < 3750.0) nOK = 1;

//	if(polyData->nPolyLayer[nPolygon] == 4) nOK=1;

	return (nOK);
}


int Polygon_IO_Remove_Zero_Area_Polygons(POLYGONATTRIBUTES *polyData,POLYGONATTRIBUTES *polyNew)
{			/* Strip any zero area polygons from polyData, store result in polyNew */
			/* A zero area polygon is one with 2 or less vertices */
int n,nNewPolygon,nNewPolygons,nZeroCount;

	Polygon_IO_End(polyNew);
	nZeroCount = Polygon_IO_Get_Zero_Area_Polygon_Count(polyData);
	nNewPolygons = polyData->nPolygons - nZeroCount;
	Polygon_IO_Copy_Header(polyData,polyNew);
	polyNew->nPolygons = nNewPolygons;
	Polygon_IO_Setup_Polygon_Attributes(&polyNew,polyData->nGroupSet,1,0,0);

	nNewPolygon = 0;
	for(n = 0; n < polyData->nPolygons; n++){		/* Copy polygons from polyData */
		if(polyData->siPolyVertices[n] > 2){
			Polygon_IO_Copy_Polygon(polyData,n,polyNew,nNewPolygon);
			nNewPolygon++;
		}
	}
	Polygon_IO_Update_Header(polyNew);

	return (nZeroCount);
}


int Polygon_IO_Renumber_Groups(POLYGONATTRIBUTES *polyData)
{			/* Renumber groups in polyData to start from zero */
int n,nGroups = 0,nGroupsMax = 0,*nNewGroup;

	for(n = 0; n < polyData->nPolygons;n++){		/* Find max. no. of groups */
		if(polyData->nGroup[n] > nGroupsMax) nGroupsMax = polyData->nGroup[n];
	}
	nGroupsMax++;
	NEW_ZERO(nNewGroup, int, nGroupsMax);
	for(n = 0; n < nGroupsMax; n++) nNewGroup[n] = -1;

	for(n = 0; n < polyData->nPolygons; n++){
		if(nNewGroup[ polyData->nGroup[n] ] == -1){
			nNewGroup[ polyData->nGroup[n] ] = nGroups;
			nGroups++;
		}
	}

	for(n = 0; n < polyData->nPolygons; n++){
		polyData->nGroup[n] = nNewGroup[ polyData->nGroup[n] ];
	}
	SAFE_DELETES(nNewGroup);

	return (nGroups);
}


void Polygon_IO_Reset_Bitmap(POLYGONATTRIBUTES *polyData)
{			/* Set the bitmap to an initial state */
int n,nSize;

	if(polyData->nBitmapSet == 1){
		nSize = (polyData->nPolyXMax * polyData->nPolyYMax);
		for(n = 0; n < nSize; n++){
			polyData->nBitmap[n] = -1;		/* Set to -1 as 0 is a polygon */
		}
	}
}


void Polygon_IO_Reset_Polygon_Attributes(POLYGONATTRIBUTES *polyData)
{			/* Reset the attributes of a set of polygons */
			/* Call before first use of structure */
	polyData->nPolygons = 0;
	polyData->nPolyXMax = 0;
	polyData->nPolyYMax = 0;
	polyData->nVertexMax = 0;
	polyData->pSize.x = 0.0;
	polyData->pSize.y = 0.0;
	polyData->dSeedToVertexMax2 = 0.0;
	polyData->nBitmapSet = 0;
	polyData->nBoundarySet = 0;
	polyData->nGroupSet = 0;
	polyData->nMagnetisationSet = 0;
	polyData->nThicknessSet = 0;
	polyData->nTetrahedronSet = 0;
	strcpy(polyData->cName,"");
}


void Polygon_IO_Reverse_Vertex_Order(POLYGONATTRIBUTES *polyData,int nPolygon)
{			/* Reverse the order of the vertices in a polygon (change chirality) */
point p;
int n,nSwap;

	nSwap=polyData->siPolyVertices[nPolygon]-1;
	for(n=1;n<polyData->siPolyVertices[nPolygon];n++){		/* Swap point n with point nSwap */
		if(n<nSwap){
			p.x=polyData->pVertex[nPolygon+(n*polyData->nPolygons)].x;
			p.y=polyData->pVertex[nPolygon+(n*polyData->nPolygons)].y;
			polyData->pVertex[nPolygon+(n*polyData->nPolygons)].x=polyData->pVertex[nPolygon+(nSwap*polyData->nPolygons)].x;
			polyData->pVertex[nPolygon+(n*polyData->nPolygons)].y=polyData->pVertex[nPolygon+(nSwap*polyData->nPolygons)].y;
			polyData->pVertex[nPolygon+(nSwap*polyData->nPolygons)].x=p.x;
			polyData->pVertex[nPolygon+(nSwap*polyData->nPolygons)].y=p.y;
			nSwap--;
		}
	}
}


void Polygon_IO_Save_P3D(char *cSavefile,POLYGONATTRIBUTES *polyData)
{			/* Save a .p3d file */
	Polygon_IO_Save_P3D_Header(cSavefile,polyData);
	Polygon_IO_Save_P3D_Data(cSavefile,polyData);
}


void Polygon_IO_Save_P3D_Binary(char *cSavefile,POLYGONATTRIBUTES *polyData)
{			/* Save a .bp3d file in binary format */
	Polygon_IO_Save_P3D_Header_Binary(cSavefile,polyData);
	Polygon_IO_Save_P3D_Data_Binary(cSavefile,polyData);
}


void Polygon_IO_Save_P3D_Data(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Save the data points in .p3d file format */
FILE *fp;
int nPolygon,nVertex;
point p;

	fp = fopen(cFilename,"a");
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		fprintf(fp,"%d %d ",polyData->nPolyLayer[nPolygon],polyData->siPolyVertices[nPolygon]);
		fprintf(fp,"%lf %lf ",polyData->pSeed[nPolygon].x,polyData->pSeed[nPolygon].y);
		for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
			p = Polygon_IO_Get_Vertex(polyData,nPolygon,nVertex);
			fprintf(fp,"%.3f %.3f ",p.x,p.y);
		}
	}
	fclose(fp);
}


void Polygon_IO_Save_P3D_Data_Binary(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Save data in .bp3d format */
FILE *fp;
int n[2],nByte,nPolygon,nVertex,nSizeInt,nSizeFloat;
float f[2];
point p;
char *c;

	nSizeInt = sizeof(int) * 2;
	nSizeFloat = sizeof(float) * 2;
	fp = fopen(cFilename,"a");
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		n[0] = polyData->nPolyLayer[nPolygon];			/* Save layer and no. of vertices */
		n[1] = polyData->siPolyVertices[nPolygon];
		c = (char*)n;
		for(nByte = 0; nByte < nSizeInt; nByte++){
			fprintf(fp,"%c",*c);
			c++;
		}

		f[0] = (float)polyData->pSeed[nPolygon].x;		/* Save seed point */
		f[1] = (float)polyData->pSeed[nPolygon].y;
		c = (char*)f;
		for(nByte = 0; nByte < nSizeFloat; nByte++){
			fprintf(fp,"%c",*c);
			c++;
		}

		for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
			p = Polygon_IO_Get_Vertex(polyData,nPolygon,nVertex);
			f[0] = (float)p.x;		/* Save vertices */
			f[1] = (float)p.y;
			c = (char*)f;
			for(nByte = 0; nByte < nSizeFloat; nByte++){
				fprintf(fp,"%c",*c);
				c++;
			} 
		}
	}
	fclose(fp);
}


void Polygon_IO_Save_P3D_Header(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Save the header information for a .p3d file */
FILE *fp;

	fp=fopen(cFilename,"w");
	fprintf(fp,"%d\n",polyData->nPolygons);
	fprintf(fp,"%d\n",polyData->nVertexMax);
	fprintf(fp,"%lf\n",polyData->pSize.x);
	fprintf(fp,"%lf\n",polyData->pSize.y);
	fprintf(fp,"%lf\n",polyData->dSeedToVertexMax2);
	fclose(fp);
}


void Polygon_IO_Save_P3D_Header_Binary(char *cSavefile,POLYGONATTRIBUTES *polyData)
{			/* Save a .bp3d header */
FILE *fp;
int n[2],nByte,nSizeInt,nSizeFloat;
float f[3];
char *c;

	nSizeInt = sizeof(int) * 2;
	nSizeFloat = sizeof(float) * 3;
	fp = fopen(cSavefile,"w");
	n[0] = polyData->nPolygons;
	n[1] = polyData->nVertexMax;
	c = (char*)n;
	for(nByte = 0; nByte < nSizeInt; nByte++){
		fprintf(fp,"%c",*c);
		c++;
	}

	f[0] = (float)polyData->pSize.x;
	f[1] = (float)polyData->pSize.y;
	f[2] = (float)polyData->dSeedToVertexMax2;
	c = (char*)f;
	for(nByte = 0; nByte < nSizeFloat; nByte++){
		fprintf(fp,"%c",*c);
		c++;
	}

	fclose(fp);
}


void Polygon_IO_Save_Pol(char *cSavefile,POLYGONATTRIBUTES *polyData)
{			/* Save a .pol file */
	Polygon_IO_Save_Pol_Header(cSavefile,polyData);
	Polygon_IO_Save_Pol_Data(cSavefile,polyData);
}


void Polygon_IO_Save_Pol_Binary(const char *cSavefile,POLYGONATTRIBUTES *polyData)
{			/* Save a .bpol file in binary format */
	Polygon_IO_Save_Pol_Header_Binary(cSavefile,polyData);
	Polygon_IO_Save_Pol_Data_Binary(cSavefile,polyData);
}


void Polygon_IO_Save_Pol_Data(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Save the data points in .pol file format */
FILE *fp;
int nPolygon,nVertex;
point p;

	fp = fopen(cFilename,"a");
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		fprintf(fp,"%d %d ",polyData->nPolyLayer[nPolygon],(int)polyData->siPolyVertices[nPolygon]);
		if(polyData->nGroupSet == 1) fprintf(fp,"%d ",polyData->nGroup[nPolygon]);		/* Save group */
		if(polyData->nBoundarySet == 1) fprintf(fp,"%f ",polyData->fBoundary[nPolygon]);		/* Save boundary */
		if(polyData->nThicknessSet == 1) fprintf(fp,"%lf ",polyData->dThickness[nPolygon]);		/* Save thickness */
		if(polyData->nMagnetisationSet == 1) fprintf(fp,"%lf %lf %lf ",polyData->pM[nPolygon].x,polyData->pM[nPolygon].y,polyData->pM[nPolygon].z);		/* Save magnetisation */
		fprintf(fp,"%lf %lf %lf ",polyData->pSeed[nPolygon].x,polyData->pSeed[nPolygon].y,polyData->pSeed[nPolygon].z);
		for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
			p = polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)];		/* Save vertices without boundary factor */
			if(polyData->nTetrahedronSet == 0){
				fprintf(fp,"%.3f %.3f ",p.x,p.y);		/* Save vertices */
			}
			else{
				fprintf(fp,"%.3f %.3f %.3f ",p.x,p.y,p.z);		/* Save Z data for tetrahedrons */
			}
		}
	}
	fclose(fp);
}


void Polygon_IO_Save_Pol_Data_Binary(const char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Save data in .bpol format */
FILE *fp;
int n[3],nByte,nPolygon,nVertex,nSizeInt,nSizeFloat,nSizeVertex;
int nIndex,nBoundaryIndex,nThicknessIndex,nMagnetisationIndex;
float f[8];
point p;
char *c;

	nSizeInt = sizeof(int) * 2;			/* Integer part */
	if(polyData->nGroupSet == 1) nSizeInt += sizeof(int);

	nSizeFloat = sizeof(float) * 3;		/* Float part */
	nIndex = 0;
	if(polyData->nBoundarySet == 1){
		nSizeFloat += sizeof(float);
		nBoundaryIndex = nIndex;
		nIndex++;
	}

	if(polyData->nThicknessSet == 1){
		nSizeFloat += sizeof(float);
		nThicknessIndex = nIndex;
		nIndex++;
	}

	if(polyData->nMagnetisationSet == 1){
		nSizeFloat += (sizeof(float) * 3);
		nMagnetisationIndex = nIndex;
		nIndex += 3;
	}

	if(polyData->nTetrahedronSet == 0){
		nSizeVertex = sizeof(float) * 2;		/* For vertices */
	}
	else{
		nSizeVertex = sizeof(float) * 3;		/* Save Z data for tetrahedrons */
	}
	fp = fopen(cFilename,"ab");
	for(nPolygon = 0; nPolygon < polyData->nPolygons; nPolygon++){
		n[0] = polyData->nPolyLayer[nPolygon];			/* Save layer and no. of vertices */
		n[1] = (int)polyData->siPolyVertices[nPolygon];
		if(polyData->nGroupSet == 1) n[2] = polyData->nGroup[nPolygon];					/* Will save group if nGroupSet = 1 */
		c = (char*)n;
		for(nByte = 0; nByte < nSizeInt; nByte++){
			fprintf(fp,"%c",*c);
			c++;
		}

		if(polyData->nBoundarySet == 1) f[nBoundaryIndex] = polyData->fBoundary[nPolygon];
		if(polyData->nThicknessSet == 1) f[nThicknessIndex] = (float)polyData->dThickness[nPolygon];
		if(polyData->nMagnetisationSet == 1){
			f[nMagnetisationIndex] = (float)polyData->pM[nPolygon].x;
			f[nMagnetisationIndex + 1] = (float)polyData->pM[nPolygon].y;
			f[nMagnetisationIndex + 2] = (float)polyData->pM[nPolygon].z;
		}
		f[nIndex] = (float)polyData->pSeed[nPolygon].x;		/* Save seed point */
		f[nIndex + 1] = (float)polyData->pSeed[nPolygon].y;
		f[nIndex + 2] = (float)polyData->pSeed[nPolygon].z;

		c = (char*)f;
		for(nByte = 0; nByte < nSizeFloat; nByte++){
			fprintf(fp,"%c",*c);
			c++;
		}

		for(nVertex = 0; nVertex < polyData->siPolyVertices[nPolygon]; nVertex++){
			p = polyData->pVertex[nPolygon + (nVertex * polyData->nPolygons)];		/* Save vertices without boundary factor */
			f[0] = (float)p.x;		/* Save vertices */
			f[1] = (float)p.y;
			if(polyData->nTetrahedronSet == 1) f[2] = (float)p.z;	/* Save Z data for tetrahedrons */
			c = (char*)f;
			for(nByte = 0; nByte < nSizeVertex; nByte++){
				fprintf(fp,"%c",*c);
				c++;
			} 
		}
	}
	fclose(fp);
}


void Polygon_IO_Save_Pol_Header(char *cFilename,POLYGONATTRIBUTES *polyData)
{			/* Save the header information for a .pol file */
FILE *fp;
unsigned int uiState = Polygon_IO_Get_State_Bytes(polyData);		/* Get state of polyData object */

	fp = fopen(cFilename,"w");
	fprintf(fp,"%d\n",polyData->nPolygons);			/* No. of polygons */
	fprintf(fp,"%d\n",polyData->nVertexMax);		/* Max. no. of vertices */
	fprintf(fp,"%u\n",uiState);									/* State (are other items saved?) */
	fprintf(fp,"%lf\n",polyData->pSize.x);		/* X extent */
	fprintf(fp,"%lf\n",polyData->pSize.y);		/* Y extent */
	fprintf(fp,"%lf\n",polyData->dSeedToVertexMax2);
	fclose(fp);
}


void Polygon_IO_Save_Pol_Header_Binary(const char *cSavefile,POLYGONATTRIBUTES *polyData)
{			/* Save a .bpol header */
FILE *fp;
int n[4],nByte,nSizeInt,nSizeFloat;
unsigned int uiState[1];
float f[3];
char *c;

	nSizeInt = sizeof(int) * 3;
	nSizeFloat = sizeof(float) * 3;
	uiState[0] = Polygon_IO_Get_State_Bytes(polyData);		/* Get state of polyData object */
	fp = fopen(cSavefile,"wb");
	n[0] = 123;			/* For checking endianness */
	n[1] = polyData->nPolygons;
	n[2] = polyData->nVertexMax;
	c = (char*)n;
	for(nByte = 0; nByte < nSizeInt; nByte++){
		fprintf(fp,"%c",*c);
		c++;
	}

	c = (char*)uiState;
	for(nByte = 0; nByte < sizeof(unsigned int); nByte++){
		fprintf(fp,"%c",*c);
		c++;
	}

	f[0] = (float)polyData->pSize.x;
	f[1] = (float)polyData->pSize.y;
	f[2] = (float)polyData->dSeedToVertexMax2;
	c = (char*)f;
	for(nByte = 0; nByte < nSizeFloat; nByte++){
		fprintf(fp,"%c",*c);
		c++;
	}
	fclose(fp);
}


void Polygon_IO_Save_Polygon_Areas(char *cFilename, POLYGONATTRIBUTES *polyData)
{			/* Save the areas of the polygons */
int n;
FILE *fp;

	fp = fopen(cFilename,"w");
	for(n = 0; n < polyData->nPolygons; n++){
		fprintf(fp,"%lf\n",fabs(Polygon_IO_Get_Polygon_Area(polyData,n)));
	}
	fclose(fp);
}


void Polygon_IO_Save_Seed_Points(char *cFilename, POLYGONATTRIBUTES *polyData)
{			/* Save the seed points */
int n;
FILE *fp;

	fp = fopen(cFilename,"w");
	for(n = 0; n < polyData->nPolygons; n++){
		fprintf(fp,"%f %f\n",(float)polyData->pSeed[n].x,(float)polyData->pSeed[n].y);
	}
	fclose(fp);
}


void Polygon_IO_Set_Group_Indexing_Layer(POLYGONATTRIBUTES *polyData, int nLayer, int nStart)
{			/* Set group index for polyData cells */
int n,nGroup = nStart - 1,nLastGroup = -1;

	for(n = 0; n < polyData->nPolygons; n++){
		if(polyData->nPolyLayer[n] == nLayer){
			if(polyData->nGroup[n] != nLastGroup){		/* Cell is in a new group, increment nGroup */
				nLastGroup = polyData->nGroup[n];
				nGroup++;
				polyData->nGroup[n] = nGroup;		/* Set new group */
			}
			else{		/* Cells are in the same group, keep same nGroup */
				polyData->nGroup[n] = nGroup;
			}
		}
	}
}


int Polygon_IO_Set_Layer(POLYGONATTRIBUTES *polyData,int nPolygon,int nLayer)
{			/* Set the layer for a polygon */
int nOK = 0;

	if(nPolygon < polyData->nPolygons){
		polyData->nPolyLayer[nPolygon] = nLayer;
		nOK = 1;
	}

	return (nOK);
}


int Polygon_IO_Set_Name(POLYGONATTRIBUTES *polyData,char *c)
{			/* Set the name of the polygon structure */
int nOK = 0;

	if(strlen(c) < 255){
		strcpy(polyData->cName,c);
		nOK = 1;
	}

	return (nOK);
}

int Polygon_IO_Set_Seed_Point(POLYGONATTRIBUTES *polyData,int nPolygon,double dX,double dY)
{			/* Set a seed point */
int nOK=0;

	if(nPolygon < polyData->nPolygons){
		polyData->pSeed[nPolygon].x=dX;
		polyData->pSeed[nPolygon].y=dY;
		nOK=1;
	}

	return (nOK);
}


int Polygon_IO_Set_Vertex(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex,double dX,double dY)
{			/* Set a vertex for a polygon */
int nOK=0;

	if(nPolygon < polyData->nPolygons){
		if(nVertex < polyData->nVertexMax){
			polyData->pVertex[nPolygon+(nVertex*polyData->nPolygons)].x=dX;
			polyData->pVertex[nPolygon+(nVertex*polyData->nPolygons)].y=dY;
			nOK=1;
		}
	}

	return (nOK);
}


int Polygon_IO_Set_Vertex_Count(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertices)
{			/* Set the number of vertices for a polygon */
int nOK=0;

	if(nPolygon < polyData->nPolygons){
		polyData->siPolyVertices[nPolygon]=(short int)nVertices;
		nOK=1;
	}

	return (nOK);
}


int Polygon_IO_Set_Virtual(POLYGONATTRIBUTES *polyData,double dXFrom,double dXTo,double dYFrom,double dYTo)
{			/* Mark which polygons are virtual */
int n,nIn = 0;

	for(n = 0; n < polyData->nPolygons; n++){
		polyData->nGroup[n] = 1;		/* 1 indicates polygon is virtual */
		if(polyData->pSeed[n].x > dXFrom){
			if(polyData->pSeed[n].x < dXTo){
				if(polyData->pSeed[n].y > dYFrom){
					if(polyData->pSeed[n].y < dYTo){
						polyData->nGroup[n] = 0;		/* 0 = polygons are included in full LLG calculation */
						nIn++;
					}
				}
			}
		}
	}

	return (polyData->nPolygons - nIn);
}


void Polygon_IO_Setup_Polygon_Attributes(POLYGONATTRIBUTES **polyData, int nGroupSet, int nBoundarySet, int nThicknessSet, int nMagnetisationSet)
{			/* Dimension the arrays in polyData before loading polygons */
	if((*polyData)->nPolygons > 0){
		NEW_ZERO((*polyData)->pCentre, point, (*polyData)->nPolygons);	/* Polygon centres */
		NEW_ZERO((*polyData)->pSeed, point, (*polyData)->nPolygons);		/* Seed points */
		NEW_ZERO((*polyData)->pVertex, point, ((*polyData)->nPolygons * (*polyData)->nVertexMax));		/* Vertices */
		NEW_ZERO((*polyData)->siPolyVertices, short int, (*polyData)->nPolygons);					/* Number of vertices per polygon */
		NEW_ZERO((*polyData)->nPolyLayer, int, (*polyData)->nPolygons);		/* Layer polygon is in */

		if(nBoundarySet == 1){
			NEW_ZERO((*polyData)->fBoundary, float, (*polyData)->nPolygons);		/* Boundary thicknesses */
			(*polyData)->nBoundarySet = 1;		/* Indicate boundary array is set */
		}

		if(nGroupSet == 1){
			NEW_ZERO((*polyData)->nGroup, int, (*polyData)->nPolygons);		/* Group numbering */
			(*polyData)->nGroupSet = 1;		/* Indicate group array is set */
		}
	}
}


void Polygon_IO_Setup_Polygon_Attributes_Bitmap(POLYGONATTRIBUTES **polyData)
{			/* Make space for a bitmap of grain number vs. (x,y) */
int nX,nY;

	if((*polyData)->nBitmapSet != 1){
		nX = (int)((*polyData)->pSize.x) + 1;		/* Add 1 to create a buffer */
		nY = (int)((*polyData)->pSize.y) + 1;
		printf("    Bitmap size ( %d , %d )\n",nX,nY);
		NEW_ZERO((*polyData)->nBitmap, int, (nX * nY));
		(*polyData)->nBitmapSet = 1;
		(*polyData)->nPolyXMax = nX;						/* Store size for faster access to bitmap */
		(*polyData)->nPolyYMax = nY;
	}
}


void Polygon_IO_Setup_Polygon_Attributes_From_P(POLYGONATTRIBUTES *polyData, int nGroupSet, int nBoundarySet, int nThicknessSet, int nMagnetisationSet)
{			/* Dimension the arrays in polyData before loading polygons */
	Polygon_IO_Setup_Polygon_Attributes(&polyData,nGroupSet,nBoundarySet,nThicknessSet,nMagnetisationSet);
}


void Polygon_IO_Setup_Polygon_Attributes_M(POLYGONATTRIBUTES **polyData)
{			/* Dimension the arrays in polyData for magnetisation */
	NEW_ZERO((*polyData)->pM, point, (*polyData)->nPolygons);	/* Polygon magnetisation */
	(*polyData)->nMagnetisationSet = 1;
}


void Polygon_IO_Shift_Alternate_Columns(POLYGONATTRIBUTES *polyData,double dWidth,double dOffset,int nOddEven,int nRotate)
{			/* Shift alternate columns of dots up or down track by dOffset */
int nColumn,nPolygon,nPolygons = polyData->nPolygons,nV;
double dX,dY;

	for(nPolygon = 0; nPolygon < nPolygons; nPolygon++){
		nColumn = (int)(polyData->pSeed[nPolygon].x / dWidth);
		if((nColumn % 2 == 0 && nOddEven == 0) || (nColumn % 2 == 1 && nOddEven == 1)){		/* Decide if cell should be moved */
			polyData->pSeed[nPolygon].y += dOffset;
			polyData->pCentre[nPolygon].y += dOffset;
			for(nV = 0; nV < polyData->siPolyVertices[nPolygon]; nV++){
				polyData->pVertex[nPolygon + (nV * nPolygons)].y += dOffset;
			}
			if(nRotate == 1){
				for(nV = 0; nV < polyData->siPolyVertices[nPolygon]; nV++){
					dX = polyData->pVertex[nPolygon + (nV * nPolygons)].x - polyData->pSeed[nPolygon].x;
					dY = polyData->pVertex[nPolygon + (nV * nPolygons)].y - polyData->pSeed[nPolygon].y;
					polyData->pVertex[nPolygon + (nV * nPolygons)].x = polyData->pSeed[nPolygon].x + dY;
					polyData->pVertex[nPolygon + (nV * nPolygons)].y = polyData->pSeed[nPolygon].y - dX;
				}
			}
		}
	}
}


int Polygon_IO_Subdivide_Polygons(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut, int nBoundary)
{			/* Subdivide each polygon in a set of polygons into smaller polygons */
			/* nBoundary = 1 creates new polygons with the intended boundary */
			/* nBoundary = 0 ignores boundary settings */
int nCore,nNew = 0,nNewPolygons = 0,nPolygon,nV;
point p0,p1,p2,p3;
double dArea,dB,dCoreArea;

	Polygon_IO_End(polyOut);
	for(nPolygon = 0; nPolygon < polyIn->nPolygons; nPolygon++){
		if(polyIn->siPolyVertices[nPolygon] > 2) nNewPolygons += (polyIn->siPolyVertices[nPolygon] + 1);
	}
	printf("    %d polygons in subdivided mesh\n",nNewPolygons);
	polyOut->nPolygons = nNewPolygons;
	if(polyIn->nVertexMax == 3){		/* If all cells are triangles, some quadrilaterals will be created */
		polyOut->nVertexMax = 4;
	}
	else{
		polyOut->nVertexMax = polyIn->nVertexMax;
	}
	Polygon_IO_Setup_Polygon_Attributes(&polyOut,1,1,0,0);

	for(nPolygon = 0; nPolygon < polyIn->nPolygons; nPolygon++){
		if(polyIn->siPolyVertices[nPolygon] > 2){
			dArea = Polygon_IO_Get_Polygon_Area(polyIn,nPolygon);
			dCoreArea = dArea / (double)(polyIn->siPolyVertices[nPolygon] + 1);		/* Area of new core polygon */
			dB = 1.0 - sqrt(dCoreArea / dArea);			/* Boundary for new vertices */

			polyOut->siPolyVertices[nNew] = polyIn->siPolyVertices[nPolygon];				/* First create the core polygon */
			polyOut->nPolyLayer[nNew] = polyIn->nPolyLayer[nPolygon];
			polyOut->pSeed[nNew] = Polygon_IO_Get_Polygon_Centre(polyIn,nPolygon);
			polyOut->fBoundary[nNew] = 0.0;			/* No boundary for core polygon */
//			polyOut->nGroup[nNew] = nPolygon;
			polyOut->nGroup[nNew] = polyIn->nGroup[nPolygon];		/* Copy group from old polygon */
			nCore = nNew;		/* Index of current core polygon */
			for(nV = 0; nV < polyIn->siPolyVertices[nPolygon]; nV++){
				p0.x = (dB * ( polyOut->pSeed[nNew].x - polyIn->pVertex[nPolygon + (nV * polyIn->nPolygons)].x)) + polyIn->pVertex[nPolygon + (nV * polyIn->nPolygons)].x;
				p0.y = (dB * ( polyOut->pSeed[nNew].y - polyIn->pVertex[nPolygon + (nV * polyIn->nPolygons)].y)) + polyIn->pVertex[nPolygon + (nV * polyIn->nPolygons)].y;
				polyOut->pVertex[nNew + (nV * nNewPolygons)] = p0;
			}
			nNew++;

			p1 = polyIn->pVertex[nPolygon + ((polyIn->siPolyVertices[nPolygon] - 1) * polyIn->nPolygons)];		/* Get last vertex */
			if(nBoundary == 1) p1 = Polygon_IO_Get_Vertex(polyIn,nPolygon,polyIn->siPolyVertices[nPolygon] - 1);
			p3 = polyOut->pVertex[nCore + ((polyOut->siPolyVertices[nCore] - 1) * polyOut->nPolygons)];				/* Get core vertex */
			for(nV = 0; nV < polyIn->siPolyVertices[nPolygon]; nV++){
				p0 = polyIn->pVertex[nPolygon + (nV * polyIn->nPolygons)];			/* p0 and p1 are the outside vertices */
				if(nBoundary == 1) p0 = Polygon_IO_Get_Vertex(polyIn,nPolygon,nV);
				p2 = polyOut->pVertex[nCore + (nV * polyOut->nPolygons)];				/* p2 and p3 are the inside vertices */
				polyOut->siPolyVertices[nNew] = 4;
				polyOut->nPolyLayer[nNew] = polyIn->nPolyLayer[nPolygon];
				polyOut->pVertex[nNew] = p0;
				polyOut->pVertex[nNew + polyOut->nPolygons] = p1;
				polyOut->pVertex[nNew + (2 * polyOut->nPolygons)] = p3;
				polyOut->pVertex[nNew + (3 * polyOut->nPolygons)] = p2;
				polyOut->pSeed[nNew] = Polygon_IO_Get_Polygon_Centre(polyOut,nNew);
				polyOut->fBoundary[nNew] = 0.0;			/* No boundary */
//				polyOut->nGroup[nNew] = nPolygon;
				polyOut->nGroup[nNew] = polyIn->nGroup[nPolygon];		/* Copy group from old polygon */

				nNew++;
				p1 = p0;
				p3 = p2;
			}
		}
	}
	Polygon_IO_Update_Header(polyOut);

	return (polyOut->nPolygons);
}


int Polygon_IO_Subdivide_Polygons_Circles(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut, double dSubCellSize,double dXFrom,double dXTo)
{			/* Subdivide circular polygons into smaller sub-cells */
int n,nIn,nMax = 10,nNew = 0,nPolygon;
double dX,dScale,dY,dXMax,dXMin,dYMax,dYMin;
double dDelta = dSubCellSize / pow(2.0,2 * nMax);
double dPi = 3.1415926535;
double dPx1,dPx2,dPx3,dPx4,dPy1,dPy2,dPy3,dPy4;
double dPx1p,dPx2p,dPx3p,dPx4p,dPy1p,dPy2p,dPy3p,dPy4p;
double dCos30 = cos(dPi / 6.0),dCos60 = cos(dPi / 3.0);
double dSin30 = sin(dPi / 6.0),dSin60 = sin(dPi / 3.0);

	Polygon_IO_End(polyOut);

	for(nPolygon = 0; nPolygon < polyIn->nPolygons; nPolygon++){
		Polygon_IO_Get_Polygon_Limits(polyIn,nPolygon,&dXMin,&dXMax,&dYMin,&dYMax);
		if(dXMin < dXFrom || dXMax > dXTo){			/* Don't subdivide polygon */
			nNew++;			/* Just add a single copy of the polygon */
		}
		else{
			for(dX = dXMin - dSubCellSize; dX < dXMax + dSubCellSize; dX += dSubCellSize){
				for(dY = dYMin - dSubCellSize; dY < dYMax + dSubCellSize; dY += dSubCellSize){
					nIn = 0;
					nIn += Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX,dY);		/* Check if any of the four vertices of a cell are inside the polygon */
					nIn += Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX + dSubCellSize,dY);
					nIn += Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX + dSubCellSize,dY + dSubCellSize);
					nIn += Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX,dY + dSubCellSize);
					if(nIn != 0){		/* Create a new polygon */
						nNew++;
					}
				}
			}
		}
	}
	printf("    %d polygons in subdivided mesh\n",nNew);
	polyOut->nPolygons = nNew;
	polyOut->nVertexMax = polyIn->nVertexMax;			/* Keep same nVertexMax in case not all polygons are subdivided */
	if(polyOut->nVertexMax < 6) polyOut->nVertexMax = 6;		/* Some edge cells will have 6 vertices */
	Polygon_IO_Setup_Polygon_Attributes(&polyOut,1,1,0,0);

	nNew = 0;
	for(nPolygon = 0; nPolygon < polyIn->nPolygons; nPolygon++){
		Polygon_IO_Get_Polygon_Limits(polyIn,nPolygon,&dXMin,&dXMax,&dYMin,&dYMax);

		if(dXMin < dXFrom || dXMax > dXTo){			/* Don't subdivide polygon */
			Polygon_IO_Copy_Polygon(polyIn, nPolygon, polyOut, nNew);		/* Copy the polygon as-is */
			nNew++;
		}
		else{
			for(dX = dXMin - dSubCellSize; dX < dXMax + dSubCellSize; dX += dSubCellSize){
				for(dY = dYMin - dSubCellSize; dY < dYMax + dSubCellSize; dY += dSubCellSize){
					nIn = 0;
					nIn += Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX,dY);		/* Check if any of the four vertices of a cell are inside the polygon */
					nIn += 2 * Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX + dSubCellSize,dY);
					nIn += 4 * Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX + dSubCellSize,dY + dSubCellSize);
					nIn += 8 * Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dX,dY + dSubCellSize);
					if(nIn != 0){		/* Create a new polygon */
						polyOut->siPolyVertices[nNew] = 4;				/* Sub-cells have 4 vertices (usually) */
						polyOut->nPolyLayer[nNew] = polyIn->nPolyLayer[nPolygon];
						polyOut->nGroup[nNew] = polyIn->nGroup[nPolygon];		/* Copy group from old polygon */
						polyOut->fBoundary[nNew] = 0.0;			/* No boundary for sub-cells */
						polyOut->pSeed[nNew].x = dX + (0.5 * dSubCellSize);
						polyOut->pSeed[nNew].y = dY + (0.5 * dSubCellSize);

						Polygon_IO_Set_Vertex(polyOut,nNew,0,dX,dY);			/* Default vertex assignment */
						Polygon_IO_Set_Vertex(polyOut,nNew,1,dX + dSubCellSize,dY);
						Polygon_IO_Set_Vertex(polyOut,nNew,2,dX + dSubCellSize,dY + dSubCellSize);
						Polygon_IO_Set_Vertex(polyOut,nNew,3,dX,dY + dSubCellSize);

						if(nIn == 1){		/* Only point at (dX,dY) is inside the ellipse */
							polyOut->siPolyVertices[nNew] = 5;				/* This cell will have 5 vertices */
							dPx1 = dX + dSubCellSize;
							dPy1 = dY;
							dPx2 = dX + (dSubCellSize * dCos30);
							dPy2 = dY + (dSubCellSize * dSin30);
							dPx3 = dX + (dSubCellSize * dCos60);
							dPy3 = dY + (dSubCellSize * dSin60);
							dPx4 = dX;
							dPy4 = dY + dSubCellSize;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPx1p = dPx1 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1p,dPy1) == 0){
									dPx1 = dPx1p;
								}

								dPx2p = dPx2 - (dScale * dCos30);
								dPy2p = dPy2 - (dScale * dSin30);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPx3p = dPx3 - (dScale * dCos60);
								dPy3p = dPy3 - (dScale * dSin60);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3p) == 0){
									dPx3 = dPx3p;
									dPy3 = dPy3p;
								}

								dPy4p = dPy4 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4,dPy4p) == 0){
									dPy4 = dPy4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + dDelta;
							polyOut->pSeed[nNew].y = dY + dDelta;
						}

						if(nIn == 2){		/* Only point at (dX+d,dY) is inside the ellipse */
							polyOut->siPolyVertices[nNew] = 5;				/* This cell will have 5 vertices */
							dPx1 = dX + dSubCellSize;
							dPy1 = dY + dSubCellSize;
							dPx2 = dX + dSubCellSize - (dSubCellSize * dSin30);
							dPy2 = dY + (dSubCellSize * dCos30);
							dPx3 = dX + dSubCellSize - (dSubCellSize * dSin60);
							dPy3 = dY + (dSubCellSize * dCos60);
							dPx4 = dX;
							dPy4 = dY;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPy1p = dPy1 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1,dPy1p) == 0){
									dPy1 = dPy1p;
								}

								dPx2p = dPx2 + (dScale * dSin30);
								dPy2p = dPy2 - (dScale * dCos30);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPx3p = dPx3 + (dScale * dSin60);
								dPy3p = dPy3 - (dScale * dCos60);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3p) == 0){
									dPx3 = dPx3p;
									dPy3 = dPy3p;
								}

								dPx4p = dPx4 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4p,dPy4) == 0){
									dPx4 = dPx4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX + dSubCellSize,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + dSubCellSize - dDelta;
							polyOut->pSeed[nNew].y = dY + dDelta;
						}

						if(nIn == 3){		/* Points at (dX,dY+d) and (dX+d,dY+d) are outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX + dSubCellSize;
							dPy1 = dY + dSubCellSize;
							dPx2 = dX + (0.666 * dSubCellSize);
							dPy2 = dY + dSubCellSize;
							dPx3 = dX + (0.333 * dSubCellSize);
							dPy3 = dY + dSubCellSize;
							dPx4 = dX;
							dPy4 = dY + dSubCellSize;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPy1p = dPy1 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1,dPy1p) == 0){
									dPy1 = dPy1p;
								}

								dPy2p = dPy2 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2,dPy2p) == 0){
									dPy2 = dPy2p;
								}

								dPy3p = dPy3 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3,dPy3p) == 0){
									dPy3 = dPy3p;
								}

								dPy4p = dPy4 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4,dPy4p) == 0){
									dPy4 = dPy4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX + dSubCellSize,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + 0.5 * dSubCellSize;
							polyOut->pSeed[nNew].y = dY + dDelta;
						}

						if(nIn == 4){		/* Only point at (dX+d,dY+d) is inside the ellipse */
							polyOut->siPolyVertices[nNew] = 5;				/* This cell will have 5 vertices */
							dPx1 = dX;
							dPy1 = dY + dSubCellSize;
							dPx2 = dX + dSubCellSize - (dSubCellSize * dCos30);
							dPy2 = dY + dSubCellSize - (dSubCellSize * dSin30);
							dPx3 = dX + dSubCellSize - (dSubCellSize * dCos60);
							dPy3 = dY + dSubCellSize - (dSubCellSize * dSin60);
							dPx4 = dX + dSubCellSize;
							dPy4 = dY;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPx1p = dPx1 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1p,dPy1) == 0){
									dPx1 = dPx1p;
								}

								dPx2p = dPx2 + (dScale * dCos30);
								dPy2p = dPy2 + (dScale * dSin30);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPx3p = dPx3 + (dScale * dCos60);
								dPy3p = dPy3 + (dScale * dSin60);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3p) == 0){
									dPx3 = dPx3p;
									dPy3 = dPy3p;
								}

								dPy4p = dPy4 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4,dPy4p) == 0){
									dPy4 = dPy4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX + dSubCellSize,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + dSubCellSize - dDelta;
							polyOut->pSeed[nNew].y = dY + dSubCellSize - dDelta;
						}

						if(nIn == 6){		/* Points at (dX,dY) and (dX,dY+d) are outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX;
							dPy1 = dY + dSubCellSize;
							dPx2 = dX;
							dPy2 = dY + (0.666 * dSubCellSize);
							dPx3 = dX;
							dPy3 = dY + (0.333 * dSubCellSize);
							dPx4 = dX;
							dPy4 = dY;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPx1p = dPx1 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1p,dPy1) == 0){
									dPx1 = dPx1p;
								}

								dPx2p = dPx2 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2) == 0){
									dPx2 = dPx2p;
								}

								dPx3p = dPx3 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3) == 0){
									dPx3 = dPx3p;
								}

								dPx4p = dPx4 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4p,dPy4) == 0){
									dPx4 = dPx4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX + dSubCellSize,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX + dSubCellSize,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + dSubCellSize - dDelta;
							polyOut->pSeed[nNew].y = dY + 0.5 * dSubCellSize;
						}

						if(nIn == 7){		/* Only point at (dX,dY+d) is outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX;
							dPy1 = dY + dSubCellSize;
							dPx2 = dX;
							dPy2 = dY + dSubCellSize;
							dPx3 = dX;
							dPy3 = dY + dSubCellSize;

							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPx1p = dPx1 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1p,dPy1) == 0){
									dPx1 = dPx1p;
								}

								dPx2p = dPx2 + dScale;
								dPy2p = dPy2 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPy3p = dPy3 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3,dPy3p) == 0){
									dPy3 = dPy3p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX + dSubCellSize,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dX + dSubCellSize,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx3,dPy3);
							polyOut->pSeed[nNew].x = dX + dSubCellSize - dDelta;
							polyOut->pSeed[nNew].y = dY + dSubCellSize + dDelta;
						}

						if(nIn == 8){		/* Only point at (dX,dY+d) is inside the ellipse */
							polyOut->siPolyVertices[nNew] = 5;				/* This cell will have 5 vertices */
							dPx1 = dX;
							dPy1 = dY;
							dPx2 = dX + (dSubCellSize * dSin30);
							dPy2 = dY + dSubCellSize - (dSubCellSize * dCos30);
							dPx3 = dX + (dSubCellSize * dSin60);
							dPy3 = dY + dSubCellSize - (dSubCellSize * dCos60);
							dPx4 = dX + dSubCellSize;
							dPy4 = dY + dSubCellSize;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPy1p = dPy1 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1,dPy1p) == 0){
									dPy1 = dPy1p;
								}

								dPx2p = dPx2 - (dScale * dSin30);
								dPy2p = dPy2 + (dScale * dCos30);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPx3p = dPx3 - (dScale * dSin60);
								dPy3p = dPy3 + (dScale * dCos60);
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3p) == 0){
									dPx3 = dPx3p;
									dPy3 = dPy3p;
								}

								dPx4p = dPx4 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4p,dPy4) == 0){
									dPx4 = dPx4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + dDelta;
							polyOut->pSeed[nNew].y = dY + dSubCellSize - dDelta;
						}

						if(nIn == 9){		/* Points at (dX+d,dY) and (dX+d,dY+d) are outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX + dSubCellSize;
							dPy1 = dY;
							dPx2 = dX + dSubCellSize;
							dPy2 = dY + (0.333 * dSubCellSize);
							dPx3 = dX + dSubCellSize;
							dPy3 = dY + (0.666 * dSubCellSize);
							dPx4 = dX + dSubCellSize;
							dPy4 = dY + dSubCellSize;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPx1p = dPx1 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1p,dPy1) == 0){
									dPx1 = dPx1p;
								}

								dPx2p = dPx2 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2) == 0){
									dPx2 = dPx2p;
								}

								dPx3p = dPx3 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3) == 0){
									dPx3 = dPx3p;
								}

								dPx4p = dPx4 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4p,dPy4) == 0){
									dPx4 = dPx4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX, dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + dDelta;
							polyOut->pSeed[nNew].y = dY + 0.5 * dSubCellSize;
						}

						if(nIn == 11){		/* Only point at (dX+d,dY+d) is outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX + dSubCellSize;
							dPy1 = dY + dSubCellSize;
							dPx2 = dX + dSubCellSize;
							dPy2 = dY + dSubCellSize;
							dPx3 = dX + dSubCellSize;
							dPy3 = dY + dSubCellSize;

							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPy1p = dPy1 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1,dPy1p) == 0){
									dPy1 = dPy1p;
								}

								dPx2p = dPx2 - dScale;
								dPy2p = dPy2 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPx3p = dPx3 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3) == 0){
									dPx3 = dPx3p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dX + dSubCellSize,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx3,dPy3);
							polyOut->pSeed[nNew].x = dX + dDelta;
							polyOut->pSeed[nNew].y = dY + dDelta;
						}

						if(nIn == 12){		/* Points at (dX,dY) and (dX+d,dY) are outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX;
							dPy1 = dY;
							dPx2 = dX + (0.333 * dSubCellSize);
							dPy2 = dY;
							dPx3 = dX + (0.666 * dSubCellSize);
							dPy3 = dY;
							dPx4 = dX + dSubCellSize;
							dPy4 = dY;
							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPy1p = dPy1 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1,dPy1p) == 0){
									dPy1 = dPy1p;
								}

								dPy2p = dPy2 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2,dPy2p) == 0){
									dPy2 = dPy2p;
								}

								dPy3p = dPy3 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3,dPy3p) == 0){
									dPy3 = dPy3p;
								}

								dPy4p = dPy4 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx4,dPy4p) == 0){
									dPy4 = dPy4p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX + dSubCellSize,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx3,dPy3);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx4,dPy4);
							polyOut->pSeed[nNew].x = dX + 0.5 * dSubCellSize;
							polyOut->pSeed[nNew].y = dY + dSubCellSize - dDelta;
						}

						if(nIn == 13){		/* Only point at (dX+d,dY) is outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX + dSubCellSize;
							dPy1 = dY;
							dPx2 = dX + dSubCellSize;
							dPy2 = dY;
							dPx3 = dX + dSubCellSize;
							dPy3 = dY;

							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPx1p = dPx1 - dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1p,dPy1) == 0){
									dPx1 = dPx1p;
								}

								dPx2p = dPx2 - dScale;
								dPy2p = dPy2 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPy3p = dPy3 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3,dPy3p) == 0){
									dPy3 = dPy3p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX + dSubCellSize,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dX,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx3,dPy3);
							polyOut->pSeed[nNew].x = dX + dDelta;
							polyOut->pSeed[nNew].y = dY + dSubCellSize - dDelta;
						}

						if(nIn == 14){		/* Only point at (dX,dY) is outside the ellipse */
							polyOut->siPolyVertices[nNew] = 6;				/* This cell will have 6 vertices */
							dPx1 = dX;
							dPy1 = dY;
							dPx2 = dX;
							dPy2 = dY;
							dPx3 = dX;
							dPy3 = dY;

							dScale = 0.5 * dSubCellSize;
							n = 0;
							do{
								dPy1p = dPy1 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx1,dPy1p) == 0){
									dPy1 = dPy1p;
								}

								dPx2p = dPx2 + dScale;
								dPy2p = dPy2 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx2p,dPy2p) == 0){
									dPx2 = dPx2p;
									dPy2 = dPy2p;
								}

								dPx3p = dPx3 + dScale;
								if(Polygon_IO_Is_Point_In_Polygon(polyIn,nPolygon,dPx3p,dPy3) == 0){
									dPx3 = dPx3p;
								}

								dScale *= 0.5;
								n++;
							}
							while(n < nMax);

							Polygon_IO_Set_Vertex(polyOut,nNew,0,dX + dSubCellSize,dY);
							Polygon_IO_Set_Vertex(polyOut,nNew,1,dX + dSubCellSize,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,2,dX,dY + dSubCellSize);
							Polygon_IO_Set_Vertex(polyOut,nNew,3,dPx1,dPy1);
							Polygon_IO_Set_Vertex(polyOut,nNew,4,dPx2,dPy2);
							Polygon_IO_Set_Vertex(polyOut,nNew,5,dPx3,dPy3);
							polyOut->pSeed[nNew].x = dX + dSubCellSize - dDelta;
							polyOut->pSeed[nNew].y = dY + dSubCellSize - dDelta;
						}
						nNew++;
					}
				}
			}
		}
	}
	Polygon_IO_Update_Header(polyOut);

	return (polyOut->nPolygons);
}


int Polygon_IO_Subdivide_Polygons_XY(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut, int nSubCellsX, int nSubCellsY,double dXFrom,double dXTo,double dYFrom,double dYTo,int nVirtual)
{			/* Subdivide rectangular polygons into nSubCellsX * nSubCellsY sub-cells */
int nNew,nNewPolygons = 0,nPolygon,nV,nX,nY;
double dNX,dNY,dX,dY,dXMin,dXMax,dYMin,dYMax;
point p0;

	Polygon_IO_End(polyOut);

	nNew = 0;
	for(nPolygon = 0; nPolygon < polyIn->nPolygons; nPolygon++){		/* Find polygons which can be sub-divided */
		if(polyIn->pSeed[nPolygon].x > dXFrom && polyIn->pSeed[nPolygon].x < dXTo){
			if(polyIn->pSeed[nPolygon].y > dYFrom && polyIn->pSeed[nPolygon].y < dYTo){
				nNew++;
			}
		}
	}

	nNewPolygons = (nNew * nSubCellsX * nSubCellsY) + (polyIn->nPolygons - nNew);
	if(nVirtual == 1) nNewPolygons = (nNew * nSubCellsX * nSubCellsY);	/* Only copy subdivided polygons if nVirtual = 1 */
	printf("    %d polygons in subdivided mesh (limits : x = %.1f to x = %.1f, y = %.1f to y = %.1f)\n",nNewPolygons,dXFrom,dXTo,dYFrom,dYTo);
	polyOut->nPolygons = nNewPolygons;
	polyOut->nVertexMax = polyIn->nVertexMax;		/* This should be 4 */
	Polygon_IO_Setup_Polygon_Attributes(&polyOut,1,1,0,0);

	nNew = 0;
	for(nPolygon = 0; nPolygon < polyIn->nPolygons; nPolygon++){
		if(polyIn->pSeed[nPolygon].x > dXFrom && polyIn->pSeed[nPolygon].x < dXTo && polyIn->pSeed[nPolygon].y > dYFrom && polyIn->pSeed[nPolygon].y < dYTo){		/* Cell is within range for sub-division */
			p0 = Polygon_IO_Get_Vertex(polyIn,nPolygon,0);
			dXMin = p0.x;
			dXMax = p0.x;
			dYMin = p0.y;
			dYMax = p0.y;
				for(nV = 1; nV < 4; nV++){
				p0 = Polygon_IO_Get_Vertex(polyIn,nPolygon,nV);
				if(p0.x < dXMin) dXMin = p0.x;
				if(p0.x > dXMax) dXMax = p0.x;
				if(p0.y < dYMin) dYMin = p0.y;
				if(p0.y > dYMax) dYMax = p0.y;
			}

			dX = (dXMax - dXMin) / (double)nSubCellsX;		/* Size of new cells */
			dY = (dYMax - dYMin) / (double)nSubCellsY;
			for(nX = 0; nX < nSubCellsX; nX++){
				dNX = (double)nX;
				for(nY = 0; nY < nSubCellsY; nY++){
					dNY = (double)nY;
					polyOut->siPolyVertices[nNew] = 4;				/* All sub-cells have 4 vertices */
					polyOut->nPolyLayer[nNew] = polyIn->nPolyLayer[nPolygon];
					polyOut->nGroup[nNew] = polyIn->nGroup[nPolygon];		/* Copy group from old polygon */
					polyOut->fBoundary[nNew] = 0.0;			/* No boundary for sub-cells */
					polyOut->pSeed[nNew].x = dXMin + (dX * (dNX + 0.5));
					polyOut->pSeed[nNew].y = dYMin + (dY * (dNY + 0.5));

					Polygon_IO_Set_Vertex(polyOut,nNew,0,dXMin + (dX * dNX),dYMin + (dY * dNY));
					Polygon_IO_Set_Vertex(polyOut,nNew,1,dXMin + (dX * dNX),dYMin + (dY * (dNY + 1.0)));
					Polygon_IO_Set_Vertex(polyOut,nNew,2,dXMin + (dX * (dNX + 1.0)),dYMin + (dY * (dNY + 1.0)));
					Polygon_IO_Set_Vertex(polyOut,nNew,3,dXMin + (dX * (dNX + 1.0)),dYMin + (dY * dNY));

					nNew++;
				}
			}
		}
		else{		/* Just copy existing cell without sub-division */
			if(nVirtual == 0){			/* Only add undivided cells if nVirtual = 0 */
				Polygon_IO_Copy_Polygon(polyIn,nPolygon,polyOut,nNew);
				nNew++;
			}
		}
	}

	if(nVirtual == 1) Polygon_IO_Renumber_Groups(polyOut);		/* Some groups may be missing, renumber the existing groups */
	Polygon_IO_Update_Header(polyOut);

	return (polyOut->nPolygons);
}


void Polygon_IO_Swap_Polygons(POLYGONATTRIBUTES *polyData,int nA,int nB)
{			/* Swap polygons nA and nB in polyData */
point pCentre,pSeed;
int n,nGroup,nPolyLayer,nV;
short int siPolyVertices;
float fBoundary;
//double dThickness;

point *pVertex = new point[polyData->nVertexMax];

	siPolyVertices = polyData->siPolyVertices[nA];			/* Copy nA to temporary memory */
	for(n = 0; n < siPolyVertices; n++){
		pVertex[n] = polyData->pVertex[nA + (n * polyData->nPolygons)];
	}
	nPolyLayer = polyData->nPolyLayer[nA];
	nGroup = polyData->nGroup[nA];
	fBoundary = polyData->fBoundary[nA];
//	dThickness = polyData->dThickness[nA];
	pCentre = polyData->pCentre[nA];
	pSeed = polyData->pSeed[nA];

	nV = polyData->siPolyVertices[nB];									/* Copy nB to nA */
	for(n = 0; n < nV; n++){
		polyData->pVertex[nA + (n * polyData->nPolygons)] = polyData->pVertex[nB + (n * polyData->nPolygons)];
	}
	polyData->nPolyLayer[nA] = polyData->nPolyLayer[nB];
	polyData->nGroup[nA] = polyData->nGroup[nB];
	polyData->siPolyVertices[nA] = polyData->siPolyVertices[nB];
	polyData->fBoundary[nA] = polyData->fBoundary[nB];
//	polyData->dThickness[nA] = polyData->dThickness[nB];
	polyData->pCentre[nA] = polyData->pCentre[nB];
	polyData->pSeed[nA] = polyData->pSeed[nB];

	for(n = 0; n < siPolyVertices; n++){								/* Copy temporary memory to nB */
		polyData->pVertex[nB + (n * polyData->nPolygons)] = pVertex[n];
	}
	polyData->nPolyLayer[nB] = nPolyLayer;
	polyData->nGroup[nB] = nGroup;
	polyData->siPolyVertices[nB] = siPolyVertices;
	polyData->fBoundary[nB] = fBoundary;
//	polyData->dThickness[nB] = dThickness;
	polyData->pCentre[nB] = pCentre;
	polyData->pSeed[nB] = pSeed;

delete [] pVertex;
}


void Polygon_IO_Translate_Polygon(POLYGONATTRIBUTES *polyData,int nPolygon,double dXShift,double dYShift)
{			/* Move nPolygon by a vector (dXShift,dYShift) */
int nV;
double dXMax = 0.0,dYMax = 0.0;

	polyData->pSeed[nPolygon].x += dXShift;		/* Move the seed point */
	polyData->pSeed[nPolygon].y += dYShift;
	for(nV = 0; nV < polyData->siPolyVertices[nPolygon]; nV++){		/* Move the vertices */
		polyData->pVertex[nPolygon + (nV * polyData->nPolygons)].x += dXShift;
		polyData->pVertex[nPolygon + (nV * polyData->nPolygons)].y += dYShift;
		if(polyData->pVertex[nPolygon + (nV * polyData->nPolygons)].x > dXMax) dXMax = polyData->pVertex[nPolygon + (nV * polyData->nPolygons)].x;
		if(polyData->pVertex[nPolygon + (nV * polyData->nPolygons)].y > dYMax) dYMax = polyData->pVertex[nPolygon + (nV * polyData->nPolygons)].y;
	}

	if(dXMax > polyData->pSize.x) polyData->pSize.x = dXMax;
	if(dYMax > polyData->pSize.y) polyData->pSize.y = dYMax;
}


void Polygon_IO_Translate_Polygons(POLYGONATTRIBUTES *polyData,double dXShift,double dYShift)
{			/* Move the polygons by a vector (dXShift,dYShift) */
int n,nV;

	polyData->pSize.x += dXShift;
	polyData->pSize.y += dYShift;
	for(n = 0; n < polyData->nPolygons; n++){
		polyData->pSeed[n].x += dXShift;		/* Move the seed point */
		polyData->pSeed[n].y += dYShift;
		for(nV = 0; nV < polyData->siPolyVertices[n]; nV++){		/* Move the vertices */
			polyData->pVertex[n + (nV * polyData->nPolygons)].x += dXShift;
			polyData->pVertex[n + (nV * polyData->nPolygons)].y += dYShift;
		}
	}
}


void Polygon_IO_Translate_Polygon_Layer(POLYGONATTRIBUTES *polyData,int nZOld,int nZNew)
{			/* Move polygons from nZOld to nZNew */
int n;

	for(n = 0; n < polyData->nPolygons; n++){
		if(polyData->nPolyLayer[n] == nZOld) polyData->nPolyLayer[n] = nZNew;
	}
}


int Polygon_IO_Update_Header(POLYGONATTRIBUTES *polyData)
{			/* Update the header */
int n,nChanged = 0;
double dXMax = 0.0,dYMax = 0.0,dSeedToVertexMax = 0.0;
double dS,dX,dY;

	for(n = 0; n < polyData->nPolygons; n++){			/* Find size limits */
		Polygon_IO_Get_Polygon_Extent(polyData,n,&dS,&dX,&dY);
		if(dS > dSeedToVertexMax) dSeedToVertexMax = dS;
		if(dX > dXMax) dXMax = dX;
		if(dY > dYMax) dYMax = dY;
	}

	if(dXMax != polyData->pSize.x){
		polyData->pSize.x = dXMax;
		nChanged = 1;
	}

	if(dYMax != polyData->pSize.y){
		polyData->pSize.y = dYMax;
		nChanged = 1;
	}

	if(dSeedToVertexMax != polyData->dSeedToVertexMax2){
		polyData->dSeedToVertexMax2 = dSeedToVertexMax;
		nChanged = 1;
	}

	return (nChanged);
}


