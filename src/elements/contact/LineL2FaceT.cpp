/* LineL2FaceT.cpp */

#include "FaceT.h"
#include "Vector.h"

/* vector functions */
inline static void Cross(const double* v,  double* vXe3)
{
	vXe3[0] = -v[1];
	vXe3[1] =  v[0];
};

inline static double Dot(const double* v1, const double* v2)
{ 	return v1[0]*v2[0] + v1[1]*v2[1]; };

inline static void Diff(const double* start, const double* end, double* v)
{
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
};

inline static double Mag(const double* v)
{
	return  sqrt (Dot(v,v)) ;
};

inline static void Cross(const double* v,  double* nv)
{
	scale = 1.0/ Mag(v) ;
	nv[0] *= scale ;
	nv[1] *= scale ;
};





LineL2FaceT::LineL2FaceT
(SurfaceT& surface,iArrayT& connectivity, dArrayT& coordinates);
{
	fSurface = surface ;
	fConnectivity.Allocate(4);
	for (i = 0; i<fNumNodes ; i++) {
		fConnectivity[i] = connectivity[i];
	}
	fCoordinates = coordinates;
}

LineL2FaceT::LineL2FaceT (void)
{
	delete [] fSurface;
	delete [] fConnectivity;
	delete [] fCoordinates;
	delete [] fNeigbors;
}

void
LineL2FaceT::ComputeCentroid(Vector& centroid)
{
	centroid[0] = 0.5 * (fCoordinates[0][0] 
			   + fCoordinates[1][0] );
	centroid[1] = 0.5 * (fCoordinates[0][1] 
			   + fCoordinates[1][1] );
}

double
LineL2FaceT::ComputeRadius(void)
{
	Diff (fCoordinates[0],fCoordinates[1],diagonal);
	radius = 0.5* Mag(diagonal);
	return radius;
}

void
LineL2FaceT::ComputeNormal(double& local_coordinates, Vector& normal)
{
	t1[0] = fCoordinates[1][0] - fCoordinates[0][0];
	t1[1] = fCoordinates[1][1] - fCoordinates[0][1];
	Cross(t1,normal);
	Norm(normal,normal);
}

void
LineL2FaceT::ComputeShapeFunctions
(double& local_coordinates, double& shape_functions);
{
	xi  = local_coordinate[0];
	shape_function[0] = 0.5 * (1.0 - xi );
	shape_function[1] = 0.5 * (1.0 + xi );
}

void
LineL2FaceT::InterpolateVector
(double& local_coordinates, double& nodal_vectors, double& vector);
{
	ComputeShapeFunctions
		(double& local_coordinates, double& shape_functions);
	vector[0] = shape_function[0]*nodal_vector[0][0];
	          + shape_function[1]*nodal_vector[1][0];
	vector[1] = shape_function[0]*nodal_vector[0][1];
	          + shape_function[1]*nodal_vector[1][1];
}


void
LineL2FaceT::ComputeShapeFunctions
(double& local_coordinates, Matrix& shape_functions);
{
	ComputeShapeFunctions
		(double& local_coordinates, double& shape_functions);
	// MORE
}

double
LineL2FaceT::ComputeJacobian (double& local_coordinates);
{
	
