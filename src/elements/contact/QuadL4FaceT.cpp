/* $Id: QuadL4FaceT.cpp,v 1.7 2001-04-19 23:47:01 rjones Exp $ */

#include "QuadL4FaceT.h"
#include "FaceT.h"

#include "ContactElementT.h"
#include "dArrayT.h"
#include "dMatrixT.h"

/* vector functions */
#include "vector3D.h"


QuadL4FaceT::QuadL4FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
        FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	fNumVertexNodes = 4;
}

QuadL4FaceT::~QuadL4FaceT (void)
{
}

void
QuadL4FaceT::Initialize(void)
{
        for (int i = 0; i < fConnectivity.Length(); i++) {
                fx[i] = fSurfaceCoordinates(fConnectivity[i]);
        }
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

void
QuadL4FaceT::NodeNormal(int local_node_number,double& normal)
{ /* computes (unnormalized) outward normal at vertex node */
	int curr = local_node_number;
	int prev = Prev(local_node_number);	
	int next = Next(local_node_number);	
	Diff(fx[next],fx[curr],t1);
	Diff(fx[prev],fx[curr],t2);
	Cross(t1,t2,&normal); 
}

void
QuadL4FaceT::FaceNormal(void)
{ /* compute face average normal */
	double e1[3], e2[3], e3[3], e4[3], v1[3], v2[3];
	Add(fx[0],fx[1],e1);
	Add(fx[1],fx[2],e2);
	Add(fx[2],fx[3],e3);
	Add(fx[3],fx[0],e4);
	Diff(e4,e2,v1);
	Diff(e1,e3,v2);
	Cross(v1,v2,fnormal);
//Normalize(&fnormal);
}

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
(ContactNodeT* node,dArrayT& parameters)
{
	double tol_g  = parameters[ContactElementT::kGapTol];
	double tol_xi = parameters[ContactElementT::kXiTol];

	/* check normal opposition */
	if ( Dot(node->Normal(),fnormal) > 0.0 ) {
#if 0
	  /* compute local coordinates */

	  if( CheckLocalCoordinates(xi,tolxi) ) { 
	    /* compute gap */
	    Interpolate (xi, fx, x_proj );
	    g = Gap(x_proj, x0, nm); // inline for FaceT
	    if (CheckGap(g,tolg) ) {
		//assign opposite (chooses closest)
		isbetter = node->AssignOpposing(fSurface,this,xi,g);
		return 1;
	    }
	  }
#endif
	}
        return 0;
}

