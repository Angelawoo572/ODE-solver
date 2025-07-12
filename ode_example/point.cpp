/* Point.c - routines associated with point objects (could also be vectors) */
/* 19/10/2004 */
/*  7/ 7/2005 - Added function for 3D triangle area */
/* 11/ 7/2005 - Added dot and cross products for vectors */


#include "point.h"

#include <math.h>



void Point_Cross_Product(point *p0,point *p1,point *p2)
{			/* Calculate p0 = p1 x p2 */
	p0->x = (p1->y * p2->z) - (p1->z * p2->y);		/* y1z2 - z1y2 */
	p0->y = (p1->z * p2->x) - (p1->x * p2->z);		/* - ( x1z2-z1x2 ) */
	p0->z = (p1->x * p2->y) - (p1->y * p2->x);		/* x1y2-y1x2 */
}


void Point_Cross_Product_N(n_point *np0, n_point *np1, n_point *np2)
{			/* Calculate p0 = p1 x p2 */
	np0->x = (np1->y * np2->z) - (np1->z * np2->y);		/* y1z2 - z1y2 */
	np0->y = (np1->z * np2->x) - (np1->x * np2->z);		/* - ( x1z2-z1x2 ) */
	np0->z = (np1->x * np2->y) - (np1->y * np2->x);		/* x1y2-y1x2 */
}


double Point_Dot_Product(point *p1,point *p2)
{			/* Calculate p1.p2 */
	return ((p1->x * p2->x) + (p1->y * p2->y) + (p1->z * p2->z));
}


void Point_Get_Line_ABC_2_Points(point p1,point p2,double *dA,double *dB,double *dC)
{			/* Get equation of line passing through two points (Ax + By = C) */
	*dA = p1.y - p2.y;
	*dB = p2.x - p1.x;
	*dC = (p2.x * p1.y) - (p1.x * p2.y);
}


void Point_Get_Line_ABC_Point_And_Gradient(point p1,double dInA,double dInB,double *dA,double *dB,double *dC)
{			/* Get equation of line passing through a point p1 with gradient A/B (Ax + By = C) */
	*dA = -dInA;
	*dB = dInB;
	*dC = (p1.y * dInB) - (p1.x * dInA);
}


double Point_Get_Point_To_Line_Distance(point p0,point p1,point p2)
{			/* Distance between p0 and line p1 - p2 */
double d,dX12 = p2.x - p1.x,dY12 = p2.y - p1.y;

	d = fabs((dX12 * (p1.y - p0.y)) - (dY12 * (p1.x - p0.x)));
	d /= sqrt((dX12 * dX12) + (dY12 * dY12));

	return (d);
}


double Point_Get_Point_To_Point_Distance(point p0,point p1)
{			/* Return the distance between two points */
	return (sqrt(((p1.x - p0.x) * (p1.x - p0.x)) + ((p1.y - p0.y) * (p1.y - p0.y)) + ((p1.z - p0.z) * (p1.z - p0.z))));
}


double Point_Get_Triangle_Area(point p0,point p1,point p2)
{			/* Return the (x-y) area of a triangle given three vertices */
double dA,dB,dC,dS,dArea = 0.0;

	dA = sqrt(((p1.x - p0.x) * (p1.x - p0.x)) + ((p1.y - p0.y) * (p1.y - p0.y)));
	dB = sqrt(((p2.x - p0.x) * (p2.x - p0.x)) + ((p2.y - p0.y) * (p2.y - p0.y)));
	dC = sqrt(((p2.x - p1.x) * (p2.x - p1.x)) + ((p2.y - p1.y) * (p2.y - p1.y)));
	dS = ((dA + dB + dC) * 0.5);
	dArea = sqrt(dS * (dS - dA) * (dS - dB) * (dS - dC));

	return (dArea);
}


double Point_Get_Triangle_Area_3D(point *p)
{			/* Return the area of a triangle given three vertices */
double dA,dB,dC,dS,dArea = 0.0;

	dA = Point_Get_Point_To_Point_Distance(p[0],p[1]);
	dB = Point_Get_Point_To_Point_Distance(p[1],p[2]);
	dC = Point_Get_Point_To_Point_Distance(p[2],p[0]);
	dS = ((dA + dB + dC) * 0.5);
	dArea = sqrt(dS * (dS - dA) * (dS - dB) * (dS - dC));

	return (dArea);
}


point Point_Get_Triangle_Centre(point p0,point p1,point p2)
{			/* Return the centre of gravity of a triangle given its vertices */
point p01,p02,pCentre;
double dM1,dM2,dC1,dC2;

	p01.x = 0.5 * (p0.x + p1.x);					/* Mid point p0-p1 */
	p01.y = 0.5 * (p0.y + p1.y);
	p02.x = 0.5 * (p0.x + p2.x);					/* Mid point p0-p2 */
	p02.y = 0.5 * (p0.y + p2.y);
	dM1 = (p01.y - p2.y) / (p01.x - p2.x);	/* Slopes */
	dM2 = (p02.y - p1.y) / (p02.x - p1.x);
	dC1 = p01.y - (dM1 * p01.x);						/* Intercepts */
	dC2 = p02.y - (dM2 * p02.x);
	pCentre.x = (dC2 - dC1) / (dM1 - dM2);		/* Centre of triangle */
	pCentre.y = (dM1 * pCentre.x) + dC1;

	return (pCentre);
}


int Point_Is_Same(point *p1,point *p2,double dTolerance)
{			/* Check if points are the same, within tolerance */
int nSame = 0;
//double dTolerance = 0.001;

	if(fabs(p1->x - p2->x) < dTolerance){
		if(fabs(p1->y - p2->y) < dTolerance){
			nSame = 1;
		}
	}

	return (nSame);
}


void Point_Normalise(point *p)
{			/* Normalise the components of p */
double dR = (p->x * p->x) + (p->y * p->y) + (p->z * p->z);

	if(dR > 0.0){
		dR = sqrt(dR);
		p->x /= dR;
		p->y /= dR;
		p->z /= dR;
	}
}


void Point_Print_Properties(point p)
{			/* Display the properties of a point */
//	printf("Point co-ordinates : ( %.3f , %.3f , %.3f )\n",p.x,p.y,p.z);
}


void Point_Set(point *p,double dX,double dY,double dZ)
{			/* Set a point variable to (x,y,z) */
	p->x = dX;
	p->y = dY;
	p->z = dZ;
}


void Point_Set_F(f_point *p,float fX,float fY,float fZ)
{			/* Set a point variable to (x,y,z) */
	p->x = fX;
	p->y = fY;
	p->z = fZ;
}


void Point_Set_N(n_point *np,int nX,int nY,int nZ)
{			/* Set a point variable to (x,y,z) */
	np->x = nX;
	np->y = nY;
	np->z = nZ;
}

