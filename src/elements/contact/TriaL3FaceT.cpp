/* $Id: TriaL3FaceT.cpp,v 1.14 2011-12-01 20:38:01 beichuan Exp $ */
#include "TriaL3FaceT.h"

/* suppress CW warning messages */
#pragma warn_unusedarg off

#include "ContactElementT.h"
#include "ContactNodeT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include <cmath>

/* vector functions */
#include "vector3D.h"

using namespace Tahoe;

/* parameters */
dArray2DT TriaL3FaceT::fIntegrationPoints;

TriaL3FaceT::TriaL3FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
        FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	fNumVertexNodes = 3;
	if (!fIntegrationPoints.IsAllocated()) {
		fIntegrationPoints.Dimension(3,2);
		double* ip;
		ip = fIntegrationPoints(0);	
		ip[0] = 1.0 ; ip[1] = 0.0;
		ip = fIntegrationPoints(1);	
		ip[0] = 0.0 ; ip[1] = 1.0;
		ip = fIntegrationPoints(2);	
		ip[0] = 0.0 ; ip[1] = 0.0;
	}
	fGeometryType = GeometryT::kQuadrilateral;
}

TriaL3FaceT::~TriaL3FaceT (void)
{
}

void
TriaL3FaceT::Initialize(void)
{
        for (int i = 0; i < fConnectivity.Length(); i++) {
                fx[i] = fSurfaceCoordinates(fConnectivity[i]);
        }
}

void
TriaL3FaceT::ComputeCentroid(double* centroid) const
{
	Ave(fx[0],fx[1],fx[2],centroid);
}

double
TriaL3FaceT::ComputeRadius(void) const
{
	double centroid[3];
	ComputeCentroid(centroid);
	Diff (fx[0],centroid,centroid);
	double radius = Magnitude(centroid);
	return radius;
}

void
TriaL3FaceT::NodeNormal(int local_node_number,double* normal) const
{ /* computes (unnormalized) outward normal at vertex node */
	double t1[3], t2[3];
	Diff(fx[1],fx[0],t1);
	Diff(fx[2],fx[0],t2);
	Cross(t1,t2,normal); 
}

void
TriaL3FaceT::CalcFaceNormal(void)
{ /* compute face average normal */
    double t1[3], t2[3];
    Diff(fx[1],fx[0],t1);
    Diff(fx[2],fx[0],t2);
    Cross(t1,t2,fnormal);

}

void
TriaL3FaceT::ComputeNormal
(const double* local_coordinates,double* normal) const
{
cout << "not implemented";
throw;
}

void
TriaL3FaceT::ComputeTangent1
(const double* local_coordinates,double* tangent1) const
{
cout << "not implemented";
throw;
}

void
TriaL3FaceT::ComputeTangent2
(const double* local_coordinates,double* tangent2) const
{
cout << "not implemented";
throw;
}

void
TriaL3FaceT::InterpolatePosition
(const double* local_coordinates, double* x)
const
{
        dArrayT shape_f(3);
        ComputeShapeFunctions (local_coordinates, shape_f);
        x[0] = shape_f[0]*fx[0][0]
             + shape_f[1]*fx[1][0]
             + shape_f[2]*fx[2][0];
        x[1] = shape_f[0]*fx[0][1]
             + shape_f[1]*fx[1][1]
             + shape_f[2]*fx[2][1];
        x[2] = shape_f[0]*fx[0][2]
             + shape_f[1]*fx[1][2]
             + shape_f[2]*fx[2][2];
}

double
TriaL3FaceT::Interpolate
(const double* local_coordinates, dArrayT& nodal_values)
const
{
        dArrayT shape_f(3);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*nodal_values[0]
                     + shape_f[1]*nodal_values[1] 
                     + shape_f[2]*nodal_values[2] ;
        return value;
}

double
TriaL3FaceT::Interpolate
(const double* local_coordinates, ArrayT<double*>& nodal_values)
const
{
        dArrayT shape_f(3);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*(*nodal_values[0])
                     + shape_f[1]*(*nodal_values[1])
                     + shape_f[2]*(*nodal_values[2]) ;
        return value;
}


void
TriaL3FaceT::InterpolateVector
(const double* local_coordinates, dArray2DT& nodal_vectors, double* vector)
const
{
        dArrayT shape_f(4);
        ComputeShapeFunctions (local_coordinates, shape_f);
        vector[0] = shape_f[0]*nodal_vectors(0)[0]
                  + shape_f[1]*nodal_vectors(1)[0] 
                  + shape_f[2]*nodal_vectors(2)[0] ;
        vector[1] = shape_f[0]*nodal_vectors(0)[1]
                  + shape_f[1]*nodal_vectors(1)[1] 
                  + shape_f[2]*nodal_vectors(2)[1] ;
        vector[2] = shape_f[0]*nodal_vectors(0)[2]
                  + shape_f[1]*nodal_vectors(1)[2] 
                  + shape_f[2]*nodal_vectors(2)[2] ;
}

void
TriaL3FaceT::ComputeShapeFunctions 
(const double* local_coordinates, dArrayT& shape_functions) const
{
	double xi  = local_coordinates[0];
	double eta = local_coordinates[1];
	shape_functions[0] = xi ;
	shape_functions[1] = eta ;
	shape_functions[2] = (1.0 - xi - eta) ;
}

void
TriaL3FaceT::ComputeShapeFunctions
(const double* local_coordinates, dMatrixT& shape_functions) const
{
	shape_functions = 0.0;
	dArrayT shape_f(4);
	ComputeShapeFunctions(local_coordinates, shape_f);
	shape_functions(0,0) = shape_f[0];
	shape_functions(1,1) = shape_f[0];
	shape_functions(2,2) = shape_f[0];
	shape_functions(3,0) = shape_f[1];
	shape_functions(4,1) = shape_f[1];
	shape_functions(5,2) = shape_f[1];
	shape_functions(6,0) = shape_f[2];
	shape_functions(7,1) = shape_f[2];
	shape_functions(8,2) = shape_f[2];
}

void
TriaL3FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, dArrayT& shape_functions) const
{
}

void
TriaL3FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, dMatrixT& shape_functions) const
{
}


double
TriaL3FaceT::ComputeJacobian (const double* local_coordinates) const
{
	//HACK
	return 1.0;
}

bool
TriaL3FaceT::Projection 
(ContactNodeT* node, const dArrayT& parameters)  const
{
	double tol_g  = parameters[ContactElementT::kGapTol];
	double tol_xi = parameters[ContactElementT::kXiTol];

	const double* nm = node->Normal();
	/* check normal opposition */
	if ( Dot(nm,fnormal) < 0.0 ) {
	  const double* x0 = node->Position();
	  /* compute local coordinates */
	  double a[3], b[3], c[3];
	  Polynomial(a,b,c);
	  /* components */
	  double a1,b1,c1,a2,b2,c2,x1,x2;
	  const double* t1 = node->Tangent1();
	  const double* t2 = node->Tangent2();
	  x1 = Dot(x0,t1); 
	  x2 = Dot(x0,t2);
	  a1 = Dot(a,t1); b1 = Dot(b,t1); c1 = Dot(c,t1); 
	  a2 = Dot(a,t2); b2 = Dot(b,t2); c2 = Dot(c,t2); 
	  double p0,p1,p2,/*p3,*/m0,m1,m2;//,m3;
	  /*difference*/
	  m0 = a1 - a2 - x1 + x2;
	  m1 = b1 - b2; m2 = c1 - c2; 
	  /*average*/
	  p0 = a1 + a2 - x1 - x2;
	  p1 = b1 + b2; p2 = c1 + c2; 
	  /* reduced equation for xi */
	  double lin = p2*m1 - p1*m2;
	  double con = p2*m0 - p0*m2;
	  double xi[2];
	  xi[0] = con/lin;
	  con = p0*m1 - p1*m0;
	  xi[1] = con/lin;
/* "override" base class */
//	  if( CheckLocalCoordinates(xi,tol_xi) ) { 
  	  if( (xi[0]                > -tol_xi) 
       && (xi[1]                > -tol_xi) 
       && ((1.0 - xi[0] -xi[1]) > -tol_xi) ) { 
	    double a3,b3,c3,/*d3,*/x3;
	    x3 = Dot(x0,nm);
	    a3 = Dot(a,nm); b3 = Dot(b,nm); c3 = Dot(c,nm);
	    /* compute gap */
	    double g =  a3 + b3*xi[0] + c3*xi[1]+ - x3;
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
TriaL3FaceT::LocalBasis  
(double* normal, double* tangent1, double* tangent2) const
{
	double t2[3];
	/* calculate (approx) face tangent */
	Diff(fx[0],fx[1],t2); 	
	/* calculate tangents */
	Cross(normal,t2,tangent1);
	Normalize(tangent1);
	Cross(normal,tangent1,tangent2);
	Normalize(tangent2);
}

void
TriaL3FaceT::Quadrature
(dArray2DT& points, dArrayT& weights) const
{
	points = fIntegrationPoints;//this is dangerous
	for (int i = 0; i < fIntegrationPoints.MajorDim(); i++) {
		weights[i] = ComputeJacobian(points(i));
	}
}

/* suppress CW warning messages */
#pragma warn_unusedarg reset
