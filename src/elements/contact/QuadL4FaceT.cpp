/* QuadL4FaceT.cpp */

#include "FaceT.h"
#include "Vector.h"

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


QuadL4FaceT::QuadL4FaceT
(SurfaceT& surface,iArrayT& connectivity, dArrayT& coordinates);
{
	FaceT::FaceT
	   (SurfaceT& surface,iArrayT& connectivity, dArrayT& coordinates);
	fx1 =  fConnectivity[0];
	fx2 =  fConnectivity[1];
	fx3 =  fConnectivity[2];
	fx4 =  fConnectivity[3];
}

QuadL4FaceT::QuadL4FaceT (void)
{
	delete [] fSurface;
	delete [] fConnectivity;
	delete [] fCoordinates;
	delete [] fNeigbors;
}

void
QuadL4FaceT::ComputeCentroid(Vector& centroid)
{
	centroid[0] = 0.25 * (fx1[0] + fx2[0] + fx3[0] + fx4[0] );
	centroid[1] = 0.25 * (fx1[1] + fx2[1] + fx3[1] + fx4[1] );
	centroid[2] = 0.25 * (fx1[2] + fx2[2] + fx3[2] + fx4[2] );
}

double
QuadL4FaceT::ComputeRadius(void)
{
	Diff (fx1,fx3,diagonal);
	return radius = 0.5*Mag(diagonal);
}

void
QuadL4FaceT::ComputeNormal(Vector& normal)
{ /* compute face average normal */
	Add(fx1,fx2,e1);
	Add(fx2,fx3,e2);
	Add(fx3,fx4,e3);
	Add(fx4,fx1,e4);
	Diff(e4,e2,v1);
	Diff(e1,e3,v2);
	Cross(v1,v2,normal);
	Norm(normal,normal);
}

void
QuadL4FaceT::ComputeNormal(int local_node, Vector& normal)
{ /* compute face node normal */
}
void
QuadL4FaceT::ComputeNormal(double& local_coordinates,Vector& normal)
{
}

void
QuadL4FaceT::ComputeShapeFunctions
(double& local_coordinates, double& shape_functions);
{
	xi  = local_coordinate[0];
	eta = local_coordinate[1];
	shape_function[0] = 0.25 * (1.0 - xi ) * (1.0 - eta) ;
	shape_function[1] = 0.25 * (1.0 + xi ) * (1.0 - eta) ;
	shape_function[2] = 0.25 * (1.0 + xi ) * (1.0 + eta) ;
	shape_function[3] = 0.25 * (1.0 - xi ) * (1.0 + eta) ;
}

void
QuadL4FaceT::ComputeShapeFunctions
(double& local_coordinates, Matrix& shape_functions);
{
	ComputeShapeFunctions
		(double& local_coordinates, double& shape_functions);
	// MORE
}

double
QuadL4FaceT::ComputeJacobian (double& local_coordinates);
{
	

}
