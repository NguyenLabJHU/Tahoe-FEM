/* $Id: QuadL4FaceT.cpp,v 1.13 2001-04-30 21:27:17 rjones Exp $ */

#include "QuadL4FaceT.h"
#include "FaceT.h"

#include "ContactElementT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include <math.h>

/* vector functions */
#include "vector3D.h"

/* parameters */
static const double kTol_Quad = 0.00000001;
static const double kTol_One  = 1.00000001;

dArray2DT QuadL4FaceT::fIntegrationPoints;

QuadL4FaceT::QuadL4FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
        FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	fNumVertexNodes = 4;
	if (!fIntegrationPoints.IsAllocated()) {
		fIntegrationPoints.Allocate(4,2);
		double* ip;
		ip = fIntegrationPoints(0);	
		ip[0] = -1.0 ; ip[1] = -1.0;
		ip = fIntegrationPoints(1);	
		ip[0] =  1.0 ; ip[1] = -1.0;
		ip = fIntegrationPoints(2);	
		ip[0] =  1.0 ; ip[1] =  1.0;
		ip = fIntegrationPoints(3);	
		ip[0] = -1.0 ; ip[1] =  1.0;
	}
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
QuadL4FaceT::ComputeCentroid(double& centroid) const
{
	Ave(fx[0],fx[1],fx[2],fx[3],&centroid);
}

double
QuadL4FaceT::ComputeRadius(void) const
{
	double diagonal[3];
	Diff (fx[0],fx[2],diagonal);
	double radius = 0.5*Mag(diagonal);
	return radius;
}

void
QuadL4FaceT::NodeNormal(int local_node_number,double& normal) const
{ /* computes (unnormalized) outward normal at vertex node */
	int curr = local_node_number;
	int prev = Prev(local_node_number);	
	int next = Next(local_node_number);	
	double t1[3], t2[3];
	Diff(fx[next],fx[curr],t1);
	Diff(fx[prev],fx[curr],t2);
	Cross(t1,t2,&normal); 
}

void
QuadL4FaceT::CalcFaceNormal(void)
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
QuadL4FaceT::ComputeNormal(double* local_coordinates,double& normal) const
{
}

void
QuadL4FaceT::ComputeShapeFunctions 
(double* local_coordinates, dArrayT& shape_functions) const
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
(double* local_coordinates, dMatrixT& shape_functions) const
{
	dArrayT shape_f;
	ComputeShapeFunctions(local_coordinates, shape_f);
//MORE
}

double
QuadL4FaceT::ComputeJacobian (double* local_coordinates) const
{
	//HACK
	return 1.0;
}

bool
QuadL4FaceT::Projection 
(ContactNodeT* node,dArrayT& parameters)  const
{
	double tol_g  = parameters[ContactElementT::kGapTol];
	double tol_xi = parameters[ContactElementT::kXiTol];

	const double* nm = node->Normal();
	/* check normal opposition */
	if ( Dot(nm,fnormal) < 0.0 ) {
	  const double* x0 = node->Position();
	  /* compute local coordinates */
	  double a[3], b[3], c[3], d[3];
	  Polynomial(a,b,c,d);
	  /* components */
	  double a1,b1,c1,d1,a2,b2,c2,d2,x1,x2;
	  const double* t1 = node->Tangent1();
	  const double* t2 = node->Tangent2();
	  x1 = Dot(x0,t1); 
	  x2 = Dot(x0,t2);
	  a1 = Dot(a,t1); b1 = Dot(b,t1); c1 = Dot(c,t1); d1 = Dot(d,t1);
	  a2 = Dot(a,t2); b2 = Dot(b,t2); c2 = Dot(c,t2); d2 = Dot(d,t2);
	  double p0,p1,p2,p3,m0,m1,m2,m3;
	  /*difference*/
	  m0 = a1 - a2 - x1 + x2;
	  m1 = b1 - b2; m2 = c1 - c2; m3 = d1 - d2;
	  /*average*/
	  p0 = a1 + a2 - x1 - x2;
	  p1 = b1 + b2; p2 = c1 + c2; p3 = d1 + d2;
	  /* reduced equation for xi, valid for p0 - p2*eta != 0 */
	  double qua = p3*m1 - p1*m3;
	  double lin = p2*m1 - p1*m2 + p3*m0 - p0*m3;
	  double con = p2*m0 - p0*m2;
	  double xi[2];
	  if (fabs(qua) < kTol_Quad) { xi[0] = -con/lin; }
	  else {
		double b2a = 0.5*lin/qua;
		double discrim = lin*lin - 4.0*con*qua;
		if (discrim < 0.0) { return 0; }
		else if (b2a > kTol_One ) 
		  { xi[0] = -b2a + sqrt(b2a*b2a - qua/con); }
		else if (b2a <-kTol_One ) 
		  { xi[0] = -b2a - sqrt(b2a*b2a - qua/con); }
		else {
		  double xi1 = 0.5*(-lin + sqrt(discrim))/qua; 
		  double xi2 = 0.5*(-lin - sqrt(discrim))/qua; 
		  fabs(xi1) < kTol_One ? xi[0] = xi1 : xi[0] = xi2;}
	  }
	  if (p2 + p3*xi[0] != 0.0) 
		{xi[1] = -(p0 + p1*xi[0])/(p2 + p3*xi[0]);}
	  else
		{xi[1] = -(m0 + m1*xi[0])/(m2 + m3*xi[0]);}
	  if( CheckLocalCoordinates(xi,tol_xi) ) { 
	    double a3,b3,c3,d3,x3;
	    x3 = Dot(x0,nm);
	    a3 = Dot(a,nm); b3 = Dot(b,nm); c3 = Dot(c,nm); d3 = Dot(d,nm);
	    /* compute gap */
	    double g =  a3 + b3*xi[0] + c3*xi[1]+ d3*xi[0]*xi[1] - x3;
	    if (CheckGap(g,tol_g) ) {
		/*assign opposite (chooses closest)*/
		bool isbetter = node->AssignOpposing(fSurface,*this,xi,g);
		return isbetter;
	    }
	  }
	}
        return 0;
}


void
QuadL4FaceT::LocalBasis  
(double* normal, double* tangent1, double* tangent2) const
{
	double t2[3];
	/* calculate (approx) face tangent */
	Diff(fx[0],fx[3],t2); 	
	/* calculate tangents */
	Cross(normal,t2,tangent1);
	Normalize(tangent1);
	Cross(normal,tangent1,tangent2);
	Normalize(tangent2);
}

void
QuadL4FaceT::Quadrature
(dArray2DT& points, dArrayT& weights) const
{
	points = fIntegrationPoints;
	for (int i = 0; i < fIntegrationPoints.Length(); i++) {
		weights[i] = ComputeJacobian(points(i));
	}
}

