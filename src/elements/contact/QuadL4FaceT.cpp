/* $Id: QuadL4FaceT.cpp,v 1.6 2001-04-16 17:30:52 rjones Exp $ */

#include "QuadL4FaceT.h"
#include "FaceT.h"

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
#if 0
// declare these in header
	/* check normal opposition */
	FaceNormal(fnm) ; // recompute??

	if ( CheckOpposition(nm,fnm) ) {

	  /* compute (approximate) local coordinates */
	  Diff(x0, fx1,v1);
	  Diff(x0, fx2,v2);
	  Diff(x0, fx3,v3);
	  Diff(x0, fx4,v4);
	  a1 = TripleProduct(v1,v2,n) ;
	  a2 = TripleProduct(v2,v3,n) ;
	  a3 = TripleProduct(v3,v4,n) ;
	  a4 = TripleProduct(v4,v1,n) ;
	
	  xi(0) = ( a4 - a2)/ (a4 - a2) ;
	  xi(1) = ( a1 - a3)/ (a1 - a3) ;
	  if( CheckLocalCoordinates(xi,tolxi) ) { // inline for FaceT
	    /* compute gap */
	    Interpolate (xi, fx, x_check );
	    g = Gap(x_check, x0, nm); // inline for FaceT
	    if (CheckGap(g,tolg) ) {
		//assign opposite (chooses closest)
		isbetter = node->AssignOpposing(fSurface,this,xi,g);
		return 1;
	    }
	  }
	}

#endif
        return 0;
}

