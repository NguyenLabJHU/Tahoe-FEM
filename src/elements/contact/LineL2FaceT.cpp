/* $Id: LineL2FaceT.cpp,v 1.8 2001-04-24 00:33:22 rjones Exp $ */

#include "LineL2FaceT.h"
#include "FaceT.h"

#include "dArrayT.h"
#include "dMatrixT.h"

/* vector functions */
#include "vector2D.h"

LineL2FaceT::LineL2FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
	FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	fNumVertexNodes = 2;
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
LineL2FaceT::NodeNormal(int local_node_number, double& normal)
{
	int curr = local_node_number;
	int next = Next(curr);
	/* get sense of left and right */
	if (next > curr) {
		Diff(fx[next],fx[curr],t1);
	}
	else {
		Diff(fx[curr],fx[next],t1);
	}
	RCross(t1,&normal);
}

void
LineL2FaceT::FaceNormal(void)
{
        /* right to left */
        Diff(fx[0],fx[1],t1);
        RCross(t1,fnormal);
}


void
LineL2FaceT::ComputeNormal(dArrayT& local_coordinates, double& normal)
{
	Diff(fx[0],fx[1],t1);
	/* this assumes a CW parameterization of the boundary */
	RCross(t1,&normal);
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
(ContactNodeT* node, dArrayT& parameters) 
{
	//HACK
	return 0;
}

void
LineL2FaceT::LocalBasis
(double* normal, double* tangent1, double* tangent2)
{
}
