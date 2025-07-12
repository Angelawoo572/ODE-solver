/* Point.h */


#ifndef __POINT_H__
#define __POINT_H__

typedef struct{
	double x;
	double y;
	double z;
	} point;


typedef struct{
	float x;
	float y;
	float z;
	} f_point;


typedef struct{
	int x;
	int y;
	int z;
	} n_point;


typedef struct{
	long x;
	long y;
	long z;
	} l_point;


typedef struct{
	point pMin;
	point pMax;
} box;


typedef struct{
	point p0;
	point p1;
	point p2;
} triangle;



void		Point_Cross_Product(point *p0,point *p1,point *p2);			/* Calculate p0 = p1 x p2 */
void		Point_Cross_Product_N(n_point *np0,n_point *np1,n_point *np2);		/* Calculate p0 = p1 x p2 */
double	Point_Dot_Product(point *p1,point *p2);									/* Calculate p1.p2 */
void		Point_Get_Line_ABC_2_Points(point p1,point p2,double *dA,double *dB,double *dC);		/* Get equation of line passing through two points Ax + By = C */
void		Point_Get_Line_ABC_Point_And_Gradient(point p1,double dInA,double dInB,double *dA,double *dB,double *dC);			/* Get equation of line passing through a point p1 with gradient A/B (Ax + By = C) */
double	Point_Get_Point_To_Line_Distance(point p0,point p1,point p2);		/* Distance between p0 and line p1 - p2 */
double	Point_Get_Point_To_Point_Distance(point p0,point p1);		/* Distance between two points */
double	Point_Get_Triangle_Area(point p0,point p1,point p2);		/* Return the (x-y) area of a triangle given three vertices */
double	Point_Get_Triangle_Area_3D(point *p);										/* Area of a 3D triangle */
point		Point_Get_Triangle_Centre(point p0,point p1,point p2);	/* Centre of a triangle */
int			Point_Is_Same(point *p1,point *p2,double dTolerance);		/* Check if points are the same */
void		Point_Normalise(point *p);															/* Normalise the components of a point */
void		Point_Print_Properties(point p);												/* Show details of a point */
void		Point_Set(point *p,double dX,double dY,double dZ);			/* Set p(x,y,z) */
void		Point_Set_F(f_point *p,float fX,float fY,float fZ);			/* Set p(x,y,z) */
void		Point_Set_N(n_point *np,int nX,int nY,int nZ);					/* Set p(x,y,z) */


#endif

