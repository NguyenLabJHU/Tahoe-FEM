/* $Id: LineL2FaceT.cpp,v 1.4 2001-04-11 14:48:58 rjones Exp $ */

#include "LineL2FaceT.h"
#include "FaceT.h"

#include "dArrayT.h"
#include "dMatrixT.h"

/* vector functions */
inline static void Cross(const double* v,  double* vXe3)
{
	vXe3[0] = -v[1];
	vXe3[1] =  v[0];
};

inline static double Dot(const double* v1, const double* v2)
{ 	
	return v1[0]*v2[0] + v1[1]*v2[1]; 
};

inline static void Diff(const double* start, const double* end, double* v)
{
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
};

inline static void Add(const double* v1, const double* v2, double* v)
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
};

inline static void Ave(const double* v1, const double* v2, double* v)
{
	v[0] = 0.5 * ( v1[0] - v2[0]);
	v[1] = 0.5 * ( v1[1] - v2[1]);
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
};


/* use vector functions -----------------------------------------*/

LineL2FaceT::LineL2FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
	FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
}

LineL2FaceT::~LineL2FaceT (void)
{
}

void
LineL2FaceT::Initialize(void)
{
        for (int i = 0; i < fConnectivity.Length(); i++) {
                fx[i] = fSurfaceCoordinates(fConnectivity[i]);
        }
}


void
LineL2FaceT::ComputeCentroid(double& centroid)
{
	Ave(fx[0],fx[1],&centroid); 
}

double
LineL2FaceT::ComputeRadius(void)
{
	double diagonal[2];
	Diff (fx[0],fx[1],diagonal);
	double radius = 0.5* Mag(diagonal);
	return radius;
}

void
LineL2FaceT::ComputeNormal(dArrayT& local_coordinates, double& normal)
{
	double t1[2];
	Diff(fx[0],fx[1],t1);
	Cross(t1,&normal);
	Normalize(&normal);
}

void
LineL2FaceT::ComputeShapeFunctions
(dArrayT& local_coordinates, dArrayT& shape_functions)
{
	double xi  = local_coordinates[0];
	shape_functions[0] = 0.5 * (1.0 - xi );
	shape_functions[1] = 0.5 * (1.0 + xi );
}

void
LineL2FaceT::ComputeShapeFunctions
(dArrayT& local_coordinates, dMatrixT& shape_functions)
{
	dArrayT shape_f;
	ComputeShapeFunctions(local_coordinates, shape_f);
// MORE
}

#if 0
void
LineL2FaceT::InterpolateVector
(dArrayT& local_coordinates, dArray2DT& nodal_vectors, double& vector);
{
	ComputeShapeFunctions
		(dArrayT& local_coordinates, dArrayT& shape_functions);
	vector[0] = shape_function[0]*nodal_vector[0][0];
	          + shape_function[1]*nodal_vector[1][0];
	vector[1] = shape_function[0]*nodal_vector[0][1];
	          + shape_function[1]*nodal_vector[1][1];
}
#endif


double
LineL2FaceT::ComputeJacobian (dArrayT& local_coordinates)
{
	//HACK
	return 1.0;
}

bool
LineL2FaceT::Projection 
(double& point, double& normal, dArrayT& local_coordinates, double gap)
{
	//HACK
	return 0;
}
