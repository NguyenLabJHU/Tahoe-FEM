/* $Id: QuadL4FaceT.cpp,v 1.4 2001-04-10 00:13:57 rjones Exp $ */

#include "QuadL4FaceT.h"
#include "FaceT.h"

#include "dArrayT.h"
#include "dMatrixT.h"

/* vector functions */
inline static void Cross(const double* A, const double* B, double* AxB)
{
        AxB[0] = A[1]*B[2] - A[2]*B[1];
        AxB[1] = A[2]*B[0] - A[0]*B[2];
        AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B)
{       return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

inline static void Diff(const double* start, const double* end, double* v)
{
        v[0] = end[0] - start[0]; 
        v[1] = end[1] - start[1];
        v[2] = end[2] - start[2];
};

inline static void Add(const double* start, const double* end, double* v)
{
        v[0] = end[0] + start[0]; 
        v[1] = end[1] + start[1];
        v[2] = end[2] + start[2];
};

inline static void Ave(double* v1,double* v2, double* v3,double* v4, double* v)
{
	v[0] = 0.25*(v1[0] + v2[0] + v3[0] + v4[0]); 
	v[1] = 0.25*(v1[1] + v2[1] + v3[1] + v4[1]);
	v[2] = 0.25*(v1[2] + v2[2] + v3[2] + v4[2]);
};

inline static double Mag(const double* v)
{
	return  sqrt (Dot(v,v)) ;
};


inline static void Normalize(double* v)
{
	double scale = 1.0/ Mag(v) ;
	v[0] *= scale ;
	v[1] *= scale ;
	v[2] *= scale ;
};



/* use vector functions----------------------------------------*/


QuadL4FaceT::QuadL4FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
        FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	for (int i = 0; i < fNumNodes; i++) {
		fx[i] = fSurfaceCoordinates(fConnectivity[i]);
	}
}

QuadL4FaceT::~QuadL4FaceT (void)
{
}

void
QuadL4FaceT::ComputeCentroid(double& centroid)
{
	Ave(fx[0],fx[1],fx[2],fx[3],&centroid);
}

double
QuadL4FaceT::ComputeRadius(void)
{
	double diagonal[3];
	Diff (fx[0],fx[2],diagonal);
	double radius = 0.5*Mag(diagonal);
	return radius;
}

#if 0
void
QuadL4FaceT::ComputeNormal(double& normal)
{ /* compute face average normal */
	double e1[3], e2[3], e3[3], e4[3], v1[3], v2[3];
	Add(fx[0],fx[1],e1);
	Add(fx[1],fx[2],e2);
	Add(fx[2],fx[3],e3);
	Add(fx[3],fx[0],e4);
	Diff(e4,e2,v1);
	Diff(e1,e3,v2);
	Cross(v1,v2,&normal);
	Normlize(&normal);
}
#endif

#if 0
void
QuadL4FaceT::ComputeNormal(int local_node, double& normal)
{ /* compute face node normal */
}
#endif

void
QuadL4FaceT::ComputeNormal(dArrayT& local_coordinates,double& normal)
{
}

void
QuadL4FaceT::ComputeShapeFunctions
(dArrayT& local_coordinates, dArrayT& shape_functions)
{
	double xi  = local_coordinates[0];
	double eta = local_coordinates[1];
	shape_functions[0] = 0.25 * (1.0 - xi ) * (1.0 - eta) ;
	shape_functions[1] = 0.25 * (1.0 + xi ) * (1.0 - eta) ;
	shape_functions[2] = 0.25 * (1.0 + xi ) * (1.0 + eta) ;
	shape_functions[3] = 0.25 * (1.0 - xi ) * (1.0 + eta) ;
}

void
QuadL4FaceT::ComputeShapeFunctions
(dArrayT& local_coordinates, dMatrixT& shape_functions)
{
	dArrayT shape_f;
	ComputeShapeFunctions(local_coordinates, shape_f);
//MORE
}

double
QuadL4FaceT::ComputeJacobian (dArrayT& local_coordinates)
{
	//HACK
	return 1.0;
}

bool
QuadL4FaceT::Projection 
(double& point, double& normal, dArrayT& local_coordinates, double gap)
{
        //HACK
        return 0;
}

