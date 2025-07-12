/* Polygon_io.h */


#ifndef __POLYGON_IO_H__
#define __POLYGON_IO_H__

#include <stdio.h>
#include "point.h"

typedef struct{								/* A structure containing information about the polygons */
	int nBitmapSet;							/* A bitmap exists */
	int nBoundarySet;						/* Is the boundary set or not */
	int nGroupSet;							/* Array indicating group membership exists */
	int nMagnetisationSet;			/* Magnetisation vectors exist */
	int nPolygons;							/* Total number of polygons */
	int nPolyXMax;							/* Integer extent of polygon region (x) or bitmap */
	int nPolyYMax;							/* Integer extent of polygon region (y) */
	int nTetrahedronSet;						/* Is the file describing tetrahedrons, or flat objects */
	int nThicknessSet;					/* Is thickness data available ? */
	int nVertexMax;							/* Maximum number of vertices for a single polygon */
	double dPolyXMax;						/* Maximum extent of polygon region (A) */
	double dPolyYMax;						/* Maximum extent of polygon region (A) */
	double dPolyZMax;						/* Maximum extent of polygon region (A) */
	double dSeedToVertexMax2;		/* Maximum distance from seed point to vertex (squared) */
	point pSize;								/* Maximum extent of polygon region (A) */
	int *nBitmap;								/* Used for making a bitmap of which points are covered by which polygon */
	int *nGroup;								/* Indicates to which group a cell belongs */
	int *nPolyLayer;						/* Layer in which each polygon lies */
	short int *siPolyVertices;	/* Number of vertices per polygon */
	float *fBoundary;						/* Non-magnetic boundary thickness */
	double *dThickness;					/* Thickness of polygon */
	point *pSeed;								/* Polygon seed points */
	point *pVertex;							/* Polygon vertices */
	point *pCentre;							/* Polygon centres (centroids) */
	point *pM;									/* Magnetisation of polygons */
	char cName[255];						/* A name for the structure (to help with debugging) */
} POLYGONATTRIBUTES;



void					Polygon_IO_Add_Boundary(POLYGONATTRIBUTES *polyData,int nLayerFrom,int nLayerTo,double dBoundary);			/* Add boundary cells to polyData */
void					Polygon_IO_Check_Vertex_Ordering(POLYGONATTRIBUTES *polyData);		/* Make the vertex chirality uniform for a set of polygons */
void					Polygon_IO_Clone_Cells(POLYGONATTRIBUTES *polyData,int nSourceLayer,int nTargetLayer);		/* Clone cells from one layer to another */
void					Polygon_IO_Convert_Pol_To_P3D(char *cPolFilename,char *cP3DFilename,POLYGONATTRIBUTES *polyData,int nL);		/* Convert a .pol file to a .p3d file */
int						Polygon_IO_Copy_All_Polygons(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyCopy);		/* Copy polygons from polySource to polyCopy */
void					Polygon_IO_Copy_Header(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew);		/* Copy a header from one data set to another */
void					Polygon_IO_Copy_Polygon(POLYGONATTRIBUTES *polyFrom,int nPolyFrom,POLYGONATTRIBUTES *polyTo,int nPolyTo);		/* Copy a polygon from polyFrom to polyTo */
void					Polygon_IO_Copy_Polygon_Layer(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nLayerFrom,int nLayerTo);		/* Copy polygons from one layer in polySource to a new polygon object */
int						Polygon_IO_Copy_Polygon_Layer_As_Group_Core(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nLayerFrom,int nLayerTo);		/* Copy polygons from nLayerFrom in polySource into nLayerTo of polyNew */
int						Polygon_IO_Copy_Polygon_Layer_As_Group_XY(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nLayerFrom,int nLayerTo);			/* Copy polygons from nLayerFrom in polySource into nLayerTo of polyNew */
void					Polygon_IO_Copy_Polygons(POLYGONATTRIBUTES *polySource,POLYGONATTRIBUTES *polyNew,int nOffset);			/* Copy polygons from one data set to another. nOffset changes the polygon indices */
int						Polygon_IO_Crop_Polygon_Region(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut,double dXStart,double dXStop,double dYStart,double dYStop);		/* Crop a polygon region */
void					Polygon_IO_End(POLYGONATTRIBUTES *polyData);
int						Polygon_IO_End_Group(POLYGONATTRIBUTES *polyData);			/* Delete support for group indexing in polyData */
void					Polygon_IO_Get_Aspect_Ratio(POLYGONATTRIBUTES *polyData);		/* Calculate the aspect ratio of a set of polygons */
double				Polygon_IO_Get_Average_Polygon_Area(POLYGONATTRIBUTES *polyData);		/* Return the average area of polygons in polyData */
double				Polygon_IO_Get_Average_Polygon_Area_Region(POLYGONATTRIBUTES *polyData,double dXStart,double dXStop,double dYStart,double dYStop);		/* Return the average area of polygons in a region */
void					Polygon_IO_Get_Bitmap(POLYGONATTRIBUTES *polyData, int nLayer);			/* Make a bitmap from the polygon data */
int						Polygon_IO_Get_Byte_Order(char *cFilename);		/* Determine the byte order of a binary file */
double				Polygon_IO_Get_Double(FILE *fp);							/* Return a double extracted from a binary file */
double				Polygon_IO_Get_Double_Reverse(FILE *fp);			/* Return a double extracted from a binary file using reverse byte ordering */
void					Polygon_IO_Get_File_Attributes(char *cFilename,POLYGONATTRIBUTES *polyData);	/* Get attributes of a .pol file */
void					Polygon_IO_Get_File_Attributes_Binary(char *cFilename,POLYGONATTRIBUTES *polyData);		/* Get attributes of a .bpol binary file */
float					Polygon_IO_Get_Float(FILE *fp);				/* Get a float from a binary file */
float					Polygon_IO_Get_Float_Reverse(FILE *fp);			/* Get a float from a binary file, using reverse byte ordering */
int						Polygon_IO_Get_Group_Cell_Count_Max(POLYGONATTRIBUTES *polyData, int nLayer);			/* Get the maximum number of cells in a group in nLayer */
int						Polygon_IO_Get_Group_Count(POLYGONATTRIBUTES *polyData);			/* How many groups are in polyData ? */
int						Polygon_IO_Get_Group_Count_In_Layer(POLYGONATTRIBUTES *polyData,int nLayer);		/* Get the number of groups in nLayer */
int						Polygon_IO_Get_Int(FILE *fp);			/* Get an integer from a binary file */
int						Polygon_IO_Get_Int_Reverse(FILE *fp);		/* Get an integer from a binary file, using reverse byte ordering */
double				Polygon_IO_Get_M(POLYGONATTRIBUTES *polyData,int nX,int nY,int nMComponent);		/* Get M at (nX,nY) */
void					Polygon_IO_Get_M_X_Y_Z(POLYGONATTRIBUTES *polyData,int nX,int nY,double *dMx,double *dMy,double *dMz);		/* Get Mx, My and Mz at (nX,nY) */
int						Polygon_IO_Get_Max_Vertex_Count_In_Layer(POLYGONATTRIBUTES *polyData, int nLayer);			/* Return the maximum number of vertices for the cells in nLayer */
double				Polygon_IO_Get_Polygon_Area(POLYGONATTRIBUTES *polyData,int nPolygon);
void					Polygon_IO_Get_Polygon_Array_Circle(POLYGONATTRIBUTES *polyData, int nZLayer,double dXCentre,double dYCentre,double dRadius,double dXSize,double dYSize);		/* Create a circular object */
void					Polygon_IO_Get_Polygon_Array_Circle_Outline(POLYGONATTRIBUTES *polyData, int nZLayer,int nDotsX,int nDotsY,double dPitchX,double dPitchY,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle);			/* Create an outline of an array of circles */
void					Polygon_IO_Get_Polygon_Array_Cubic(POLYGONATTRIBUTES *polyData,int nCells,int *nX,int *nY,int *nZ,int nXCount,int nYCount,double dXSize,double dYSize);		/* Get a cubic array of polygons */
int						Polygon_IO_Get_Polygon_Array_Cubic_Dot(POLYGONATTRIBUTES *polyData, int nDotsX, int nDotsY, double dPitchX, double dPitchY, double dDotX, double dDotY, int nSubCellsX, int nSubCellsY, int nLayerFrom, int nLayerTo);			/* Generate an array of cubic dots with sub-dot discretisation */
void					Polygon_IO_Get_Polygon_Array_Hexagon(POLYGONATTRIBUTES *polyData,int nYCount,int *nZLayer,int nLayers,double dXSize,double dYSize);			/* Create an array hexagonal polygons in nLayer */
void					Polygon_IO_Get_Polygon_Array_Pattern(POLYGONATTRIBUTES *polyData,int nCells,n_point *pPattern,int nXMax,int nYMax,double dXSize,double dYSize);			/* Create an array of cubic polygons from a pattern */
void					Polygon_IO_Get_Polygon_Array_Spin_Ice_Hexagon(POLYGONATTRIBUTES *polyData,int nZLayer,int nCellsX,int nCellsY,double dPitchX,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle);					/* Create an outline of an array of circles, or ellipses in a spin-ice-like pattern with hexagonal cells */
void					Polygon_IO_Get_Polygon_Array_Spin_Ice_Square(POLYGONATTRIBUTES *polyData,int nZLayer,int nCellsX,int nCellsY,double dPitchX,double dPitchY,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle);			/* Create an outline of an array of circles, or ellipses in a square spin-ice-like pattern */
void					Polygon_IO_Get_Polygon_Array_Spin_Ice_Triangle(POLYGONATTRIBUTES *polyData,int nZLayer,int nCellsX,int nCellsY,double dPitchX,double dRadius,double dTheta,double dEllipticity,double dEllipseAngle);			/* Create an outline of an array of circles, or ellipses in a triangular spin-ice-like pattern */
int						Polygon_IO_Get_Polygon_At_XY(POLYGONATTRIBUTES *polyData,int nX,int nY);		/* Get the polygon at (nX,nY) */
point					Polygon_IO_Get_Polygon_Centre(POLYGONATTRIBUTES *polyData,int nPolygon);
int						Polygon_IO_Get_Polygon_Count_In_Layer(POLYGONATTRIBUTES *polyData,int nLayer);		/* Get the number of polygons in nLayer */
int 					Polygon_IO_Get_Polygon_Count_In_Region(POLYGONATTRIBUTES *polyData,int *nVertexMax,double dXStart,double dXStop,double dYStart,double dYStop);		/* Return the number of polygons within the defined region */
void					Polygon_IO_Get_Polygon_Extent(POLYGONATTRIBUTES *polyData,int nPolygon,double *dSeedToVertex,double *dXMax,double *dYMax);		/* Get maximum extent of a polygon */
void					Polygon_IO_Get_Polygon_Limits(POLYGONATTRIBUTES *polyData,int nPolygon,double *dXMin,double *dXMax,double *dYMin,double *dYMax);		/* Get the vertices at the limits of a polygon */
double				Polygon_IO_Get_Polygon_Perimeter(POLYGONATTRIBUTES *polyData,int nPolygon);		/* Get the perimeter length of a polygon */
double				Polygon_IO_Get_Polygon_Statistics(POLYGONATTRIBUTES *polyData);		/* Get average area and standard deviation of a set of polygons */
double				Polygon_IO_Get_Polygon_Statistics_Region(POLYGONATTRIBUTES *polyData,double dXStart,double dXStop,double dYStart,double dYStop);		/* Get average area and standard deviation of a set of polygons within a defined region */
short int			Polygon_IO_Get_Short_Int(FILE *fp);							/* Return a short integer extracted from a binary file */
short int			Polygon_IO_Get_Short_Int_Reverse(FILE *fp);			/* Return a short integer extracted from a binary file (reverse byte ordering) */
unsigned int	Polygon_IO_Get_State_Bytes(POLYGONATTRIBUTES *polyData);		/* Return an integer representing the state of the polygon object */
point					Polygon_IO_Get_Tetrahedron_Centre(POLYGONATTRIBUTES *polyData,int nPolygon);			/* Return the centroid of a tetrahedron */
double				Polygon_IO_Get_Total_Polygon_Area(POLYGONATTRIBUTES *polyData);		/* Get the total area of all polygons */
unsigned int	Polygon_IO_Get_Unsigned_Int(FILE *fp);			/* Get an unsigned integer from a binary file */
unsigned int	Polygon_IO_Get_Unsigned_Int_Reverse(FILE *fp);		/* Get an unsigned integer from a binary file, using reverse byte ordering */
void					Polygon_IO_Get_Vectors_From_Vertices(POLYGONATTRIBUTES *polyData,int nPolygon,int nBaseVertex,point *p1,point *p2);		/* Get two vectors formed from the sides of a polygon */
point					Polygon_IO_Get_Vertex(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex);			/* Return a vertex of a polygon */
int						Polygon_IO_Get_Vertex_Count(int nPolygon,int nPolygons,point *pVertex);	/* Return the number of vertices in nPolygon */
point					Polygon_IO_Get_Vertex_No_Boundary(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex);			/* Return a vertex of a polygon without applying the boundary factor */
int						Polygon_IO_Get_Zero_Area_Polygon_Count(POLYGONATTRIBUTES *polyData);		/* How many polygons have zero area ? */
void					Polygon_IO_Init(POLYGONATTRIBUTES *polyData);			/* Initialise a polygon object */
int						Polygon_IO_Init_Group(POLYGONATTRIBUTES **polyData);			/* Add support for group indexing to polyData */
void					Polygon_IO_Init_Group_Default(POLYGONATTRIBUTES *polyData);			/* Initialise polyData with default groups (nGroup = nPolygon) */
int						Polygon_IO_Is_Cubic(POLYGONATTRIBUTES *polyData,int nPolygon);
int						Polygon_IO_Is_Point_In_Polygon(POLYGONATTRIBUTES *polyData,int nPolygon,double dX,double dY);		/* Return 1 if (dX,dY) is inside polygon nPolygon */
int						Polygon_IO_Is_Point_In_Polygon_With_Boundary(POLYGONATTRIBUTES *polyData,int nPolygon,double dX,double dY);
int						Polygon_IO_Is_Polygon_In_Region(POLYGONATTRIBUTES *polyData,int nPolygon,double dXStart,double dXStop,double dYStart,double dYStop,double dTolerance);			/* Determine if a polygon (vertices) is within a region (dXStart,dYStart) to (dXStop, dYStop) */
void					Polygon_IO_Load_Magnetisation(const char *cFilename,POLYGONATTRIBUTES *polyData, int nLayerCells);	/* Load magnetisation data */
void					Polygon_IO_Load_Mesh(POLYGONATTRIBUTES *polyData,char *cFilename);					/* Load a neutral format tetrahedronal mesh */
void					Polygon_IO_Load_P3D(POLYGONATTRIBUTES *polyData,char *cFilename);						/* Load a .p3d file */
void					Polygon_IO_Load_P3D_Header(char *cFilename, POLYGONATTRIBUTES *polyData);		/* Load the header information from a .p3d file */
void					Polygon_IO_Load_P3D_Header_Binary(char *cFilename,POLYGONATTRIBUTES *polyData);		/* Load the header from a .bp3d file */
void					Polygon_IO_Load_P3D_Polygons(char *cFilename, POLYGONATTRIBUTES *polyData);			/* Load the polygon data from a .p3d file */
void					Polygon_IO_Load_P3D_Polygons_Binary(char *cFilename, POLYGONATTRIBUTES *polyData);		/* Load polygons from a .bp3d file */
void					Polygon_IO_Load_P3D_Polygons_Binary_Reverse(char *cFilename, POLYGONATTRIBUTES *polyData);		/* Load binary polygons with reversy byte ordering */
void					Polygon_IO_Load_Pol(POLYGONATTRIBUTES *polyData,char *cFilename);						/* Load a .pol file */
void					Polygon_IO_Load_Pol_Binary(POLYGONATTRIBUTES *polyData,const char *cFilename);		/* Load a .bpol file */
void					Polygon_IO_Load_Pol_Header(char *cFilename, POLYGONATTRIBUTES *polyData);		/* Load the header information from a .pol file */
int						Polygon_IO_Load_Pol_Header_Binary(const char *cFilename,POLYGONATTRIBUTES *polyData);		/* Load the header from a .bpol file */
void					Polygon_IO_Load_Pol_Polygons(char *cFilename,POLYGONATTRIBUTES *polyData);			/* Load the polygon data from a .pol file */
void					Polygon_IO_Load_Pol_Polygons_Binary(const char *cFilename, POLYGONATTRIBUTES *polyData);			/* Load the polygon data from a .bpol file */
void					Polygon_IO_Load_Pol_Polygons_Binary_Reverse(char *cFilename, POLYGONATTRIBUTES *polyData);			/* Load the polygon data from a .bpol file */
void					Polygon_IO_Load_Pol_Polygons_Old(char *cFilename,int nPolygons,point *PolySeed,point *PolyVertex);				/* Load polygons from a .pol file */
void					Polygon_IO_Merge_Polygons(POLYGONATTRIBUTES *polyA,POLYGONATTRIBUTES *polyB,POLYGONATTRIBUTES *polyNew);		/* Merge polyA and polyB together */
void					Polygon_IO_Print_Polygon_Attributes(POLYGONATTRIBUTES *polyData);	/* Print the details of a POLYGONATTRIBUTES structure */
void					Polygon_IO_Print_Polygon_Details(POLYGONATTRIBUTES *polyData,int nPolygon);			/* Print details of a specific polygon */
void					Polygon_IO_Process_Polygons(POLYGONATTRIBUTES *polyOld,POLYGONATTRIBUTES *polyNew);		/* Process polyOld and put the result in polyNew */
int						Polygon_IO_Process_Polygons_Rule(POLYGONATTRIBUTES *polyData,int nPolygon);	/* Rules to use when processing polygons */
int						Polygon_IO_Remove_Zero_Area_Polygons(POLYGONATTRIBUTES *polyData,POLYGONATTRIBUTES *polyNew);	/* Strip any zero area polygons from polyData, store result in polyNew */
int						Polygon_IO_Renumber_Groups(POLYGONATTRIBUTES *polyData);			/* Renumber groups in polyData to start from zero */
void					Polygon_IO_Reset_Bitmap(POLYGONATTRIBUTES *polyData);
void					Polygon_IO_Reset_Polygon_Attributes(POLYGONATTRIBUTES *polyData);
void					Polygon_IO_Reverse_Vertex_Order(POLYGONATTRIBUTES *polyData,int nPolygon);	/* Swap the vertex points around */
void					Polygon_IO_Save_P3D(char *cSavefile, POLYGONATTRIBUTES *polyData);		/* Save a .p3d file */
void					Polygon_IO_Save_P3D_Binary(char *cSavefile, POLYGONATTRIBUTES *polyData);		/* Save a .bp3d file */
void					Polygon_IO_Save_P3D_Data(char *cFilename, POLYGONATTRIBUTES *polyData);			/* Save data in .p3d format */
void					Polygon_IO_Save_P3D_Data_Binary(char *cFilename, POLYGONATTRIBUTES *polyData);		/* Save a .bp3d data file */
void					Polygon_IO_Save_P3D_Header(char *cFilename, POLYGONATTRIBUTES *polyData);		/* Write .p3d format header */
void					Polygon_IO_Save_P3D_Header_Binary(char *cSavefile,POLYGONATTRIBUTES *polyData);		/* Write .bp3d format header */
void					Polygon_IO_Save_Pol(char *cFilename, POLYGONATTRIBUTES *polyData);								/* Save polygons in .pol format */
void					Polygon_IO_Save_Pol_Binary(const char *cFilename, POLYGONATTRIBUTES *polyData);					/* Save polygons in .bpol format */
void					Polygon_IO_Save_Pol_Data(char *cFilename, POLYGONATTRIBUTES *polyData);						/* Save data in .pol format */
void					Polygon_IO_Save_Pol_Data_Binary(const char *cFilename, POLYGONATTRIBUTES *polyData);		/* Save data in .bpol format */
void					Polygon_IO_Save_Pol_Header(char *cFilename, POLYGONATTRIBUTES *polyData);					/* Save header in .pol format */
void					Polygon_IO_Save_Pol_Header_Binary(const char *cFilename, POLYGONATTRIBUTES *polyData);	/* Save header in .bpol format */
void					Polygon_IO_Save_Polygon_Areas(char *cFilename, POLYGONATTRIBUTES *polyData);			/* Save the areas of the polygons */
void					Polygon_IO_Save_Seed_Points(char *cFilename, POLYGONATTRIBUTES *polyData);				/* Save the seed points */
void					Polygon_IO_Set_Group_Indexing_Layer(POLYGONATTRIBUTES *polyData, int nLayer, int nStart);			/* Set group index for polyData cells */
int 					Polygon_IO_Set_Layer(POLYGONATTRIBUTES *polyData,int nPolygon,int nLayer);	/* Set the layer for a polygon */
int						Polygon_IO_Set_Name(POLYGONATTRIBUTES *polyData,char *c);			/* Set the name of the polygon structure */
int						Polygon_IO_Set_Seed_Point(POLYGONATTRIBUTES *polyData,int nPolygon,double dX,double dY);		/* Set a seed point */
int 					Polygon_IO_Set_Vertex(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertex,double dX,double dY);		/* Set a vertex of a polygon */
int 					Polygon_IO_Set_Vertex_Count(POLYGONATTRIBUTES *polyData,int nPolygon,int nVertices);		/* Set the number of vertices for a polygon */
int						Polygon_IO_Set_Virtual(POLYGONATTRIBUTES *polyData,double dXFrom,double dXTo,double dYFrom,double dYTo);		/* Mark which polygons are virtual */
void					Polygon_IO_Setup_Polygon_Attributes(POLYGONATTRIBUTES **polyData, int nGroupSet, int nBoundarySet, int nThicknessSet, int nMagnetisationSet);			/* Dimension the arrays in polyData before loading polygons */
void					Polygon_IO_Setup_Polygon_Attributes_Bitmap(POLYGONATTRIBUTES **polyData);	/* Make space for bitmap data */
void					Polygon_IO_Setup_Polygon_Attributes_From_P(POLYGONATTRIBUTES *polyData, int nGroupSet, int nBoundarySet, int nThicknessSet, int nMagnetisationSet);			/* Dimension the arrays in polyData before loading polygons */
void					Polygon_IO_Setup_Polygon_Attributes_M(POLYGONATTRIBUTES **polyData);		/* Dimension space for M data */
void					Polygon_IO_Shift_Alternate_Columns(POLYGONATTRIBUTES *polyData,double dWidth,double dOffset,int nOddEven,int nRotate);	/* Shift alternate columns of dots up or down track by dOffset */
int						Polygon_IO_Subdivide_Polygons(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut, int nBoundary);					/* Subdivide each polygon in a set of polygons into four smaller polygons */
int						Polygon_IO_Subdivide_Polygons_Circles(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut, double dSubCellSize,double dXFrom,double dXTo);			/* Subdivide circular polygons into smaller sub-cells */
int						Polygon_IO_Subdivide_Polygons_XY(POLYGONATTRIBUTES *polyIn,POLYGONATTRIBUTES *polyOut, int nSubCellsX, int nSubCellsY,double dXFrom,double dXTo,double dYFrom,double dYTo,int nVirtual);		/* Subdivide rectangular polygons into nSubCellsX * nSubCellsY sub-cells */
void					Polygon_IO_Swap_Polygons(POLYGONATTRIBUTES *polyData,int nA,int nB);		/* Swap polygons nA and nB in polyData */
void					Polygon_IO_Translate_Polygon(POLYGONATTRIBUTES *polyData,int nPolygon,double dXShift,double dYShift);			/* Move nPolygon by a vector (dXShift,dYShift) */
void					Polygon_IO_Translate_Polygons(POLYGONATTRIBUTES *polyData,double dXShift,double dYShift);		/* Move polygons */
void					Polygon_IO_Translate_Polygon_Layer(POLYGONATTRIBUTES *polyData,int nZOld,int nZNew);	/* Move polygons from nZOld to nZNew */
int						Polygon_IO_Update_Header(POLYGONATTRIBUTES *polyData);		/* Update the header */


#endif

