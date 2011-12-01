/* $Id: QuadL4FaceT.cpp,v 1.28 2011-12-01 20:38:01 beichuan Exp $ */
#include "QuadL4FaceT.h"

#include "ContactElementT.h"
#include "ContactNodeT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include <cmath>

/* vector functions */
#include "vector3D.h"

using namespace Tahoe;

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
		fIntegrationPoints.Dimension(4,2);
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
	fGeometryType = GeometryT::kQuadrilateral;
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
QuadL4FaceT::ComputeCentroid(double* centroid) const
{
	Ave(fx[0],fx[1],fx[2],fx[3],centroid);
}

double
QuadL4FaceT::ComputeRadius(void) const
{
	double diagonal[3];
	Diff (fx[0],fx[2],diagonal);
	double radius = 0.5*Magnitude(diagonal);
	return radius;
}

void
QuadL4FaceT::NodeNormal(int local_node_number,double* normal) const
{ /* computes (unnormalized) outward normal at vertex node */
	int curr = local_node_number;
	int prev = Prev(local_node_number);	
	int next = Next(local_node_number);	
	double t1[3], t2[3];
	Diff(fx[next],fx[curr],t1);
	Diff(fx[prev],fx[curr],t2);
	Cross(t1,t2,normal); 
}

void
QuadL4FaceT::CalcFaceNormal(void)
{ /* compute face average normal */
#if 0
	double e1[3], e2[3], e3[3], e4[3], v1[3], v2[3];
	Add(fx[0],fx[1],e1);
	Add(fx[1],fx[2],e2);
	Add(fx[2],fx[3],e3);
	Add(fx[3],fx[0],e4);
	Diff(e4,e2,v1);
	Diff(e1,e3,v2);
	Cross(v1,v2,fnormal);
//Normalize(&fnormal);
#endif
	double a[3], b[3], c[3], d[3];
    Polynomial(a,b,c,d);
	Cross(b,c,fnormal);
}

void
QuadL4FaceT::ComputeNormal
(const double* local_coordinates,double* normal) const
{
	double a[3], b[3], c[3], d[3];
	Polynomial(a,b,c,d);
	double v[3];
	double xi  = local_coordinates[0];
	double eta = local_coordinates[1];
	Add(xi,b,-eta,c,v);
	Cross(v,d,normal);
	double fn[3];
	Cross(b,c,fn);
	Add(normal,fn,normal);
	Normalize(normal);
}

void
QuadL4FaceT::ComputeTangent1
(const double* local_coordinates,double* tangent1) const
{
	double a[3], b[3], c[3], d[3];
	Polynomial(a,b,c,d);
	double eta = local_coordinates[1];
	Add(1.0,b,eta,d,tangent1);
}

void
QuadL4FaceT::ComputeTangent2
(const double* local_coordinates,double* tangent2) const
{
	double a[3], b[3], c[3], d[3];
	Polynomial(a,b,c,d);
	double xi  = local_coordinates[0];
	Add(1.0,c,xi,d,tangent2);
}

void
QuadL4FaceT::InterpolatePosition
(const double* local_coordinates, double* x)
const
{
        dArrayT shape_f(4);
        ComputeShapeFunctions (local_coordinates, shape_f);
        x[0] = shape_f[0]*fx[0][0]
             + shape_f[1]*fx[1][0]
             + shape_f[2]*fx[2][0]
             + shape_f[3]*fx[3][0];
        x[1] = shape_f[0]*fx[0][1]
             + shape_f[1]*fx[1][1]
             + shape_f[2]*fx[2][1]
             + shape_f[3]*fx[3][1];
        x[2] = shape_f[0]*fx[0][2]
             + shape_f[1]*fx[1][2]
             + shape_f[2]*fx[2][2]
             + shape_f[3]*fx[3][2];
}

double
QuadL4FaceT::Interpolate
(const double* local_coordinates, dArrayT& nodal_values)
const
{
        dArrayT shape_f(4);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*nodal_values[0]
                     + shape_f[1]*nodal_values[1] 
                     + shape_f[2]*nodal_values[2] 
                     + shape_f[3]*nodal_values[3];
        return value;
}

double
QuadL4FaceT::Interpolate
(const double* local_coordinates, ArrayT<double*>& nodal_values)
const
{
        dArrayT shape_f(4);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*(*nodal_values[0])
                     + shape_f[1]*(*nodal_values[1])
                     + shape_f[2]*(*nodal_values[2])
                     + shape_f[3]*(*nodal_values[3]);
        return value;
}

void
QuadL4FaceT::InterpolateVector
(const double* local_coordinates, dArray2DT& nodal_vectors, double* vector)
const
{
        dArrayT shape_f(4);
        ComputeShapeFunctions (local_coordinates, shape_f);
        vector[0] = shape_f[0]*nodal_vectors(0)[0]
                  + shape_f[1]*nodal_vectors(1)[0] 
                  + shape_f[2]*nodal_vectors(2)[0] 
                  + shape_f[3]*nodal_vectors(3)[0];
        vector[1] = shape_f[0]*nodal_vectors(0)[1]
                  + shape_f[1]*nodal_vectors(1)[1] 
                  + shape_f[2]*nodal_vectors(2)[1] 
                  + shape_f[3]*nodal_vectors(3)[1];
        vector[2] = shape_f[0]*nodal_vectors(0)[2]
                  + shape_f[1]*nodal_vectors(1)[2] 
                  + shape_f[2]*nodal_vectors(2)[2] 
                  + shape_f[3]*nodal_vectors(3)[2];
}

void
QuadL4FaceT::ComputeShapeFunctions 
(const double* local_coordinates, dArrayT& shape_functions) const
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
	shape_functions(9,0) = shape_f[3];
	shape_functions(10,1)= shape_f[3];
	shape_functions(11,2)= shape_f[3];
}

void
QuadL4FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, 
dArrayT& shape_derivatives1, dArrayT& shape_derivatives2 ) const
{
    double xi  = local_coordinates[0];
    double eta = local_coordinates[1];
    shape_derivatives1[0] = -0.25 * (1.0 - eta) ;
    shape_derivatives1[1] =  0.25 * (1.0 - eta) ;
    shape_derivatives1[2] =  0.25 * (1.0 + eta) ;
    shape_derivatives1[3] = -0.25 * (1.0 + eta) ;

    shape_derivatives2[0] = -0.25 * (1.0 - xi ) ;
    shape_derivatives2[1] = -0.25 * (1.0 + xi ) ;
    shape_derivatives2[2] =  0.25 * (1.0 + xi ) ;
    shape_derivatives2[3] =  0.25 * (1.0 - xi ) ;
}

void
QuadL4FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, 
dMatrixT& shape_derivatives1, dMatrixT& shape_derivatives2) const
{
    dArrayT shape_d1(4), shape_d2(4);
    ComputeShapeFunctionDerivatives(local_coordinates, shape_d1, shape_d2);
    shape_derivatives1 = 0.0;
    shape_derivatives1(0,0)  = shape_d1[0];
    shape_derivatives1(1,1)  = shape_d1[0];
    shape_derivatives1(2,2)  = shape_d1[0];
    shape_derivatives1(3,0)  = shape_d1[1];
    shape_derivatives1(4,1)  = shape_d1[1];
    shape_derivatives1(5,2)  = shape_d1[1];
    shape_derivatives1(6,0)  = shape_d1[2];
    shape_derivatives1(7,1)  = shape_d1[2];
    shape_derivatives1(8,2)  = shape_d1[2];
    shape_derivatives1(9,0)  = shape_d1[3];
    shape_derivatives1(10,1) = shape_d1[3];
    shape_derivatives1(11,2) = shape_d1[3];
    shape_derivatives2 = 0.0;
    shape_derivatives2(0,0)  = shape_d2[0];
    shape_derivatives2(1,1)  = shape_d2[0];
    shape_derivatives2(2,2)  = shape_d2[0];
    shape_derivatives2(3,0)  = shape_d2[1];
    shape_derivatives2(4,1)  = shape_d2[1];
    shape_derivatives2(5,2)  = shape_d2[1];
    shape_derivatives2(6,0)  = shape_d2[2];
    shape_derivatives2(7,1)  = shape_d2[2];
    shape_derivatives2(8,2)  = shape_d2[2];
    shape_derivatives2(9,0)  = shape_d2[3];
    shape_derivatives2(10,1) = shape_d2[3];
    shape_derivatives2(11,2) = shape_d2[3];
}


double
QuadL4FaceT::ComputeJacobian (const double* local_coordinates) const
{
    double a[3], b[3], c[3], d[3];
    Polynomial(a,b,c,d);
    double v[3];
    double xi  = local_coordinates[0];
    double eta = local_coordinates[1];
    Add(xi,b,-eta,c,v);
	double normal[3];
    Cross(v,d,normal);
    double fn[3];
    Cross(b,c,fn);
    Add(normal,fn,normal);
    return Magnitude(normal);
}

bool
QuadL4FaceT::Projection 
(ContactNodeT* node, const dArrayT& parameters)  const
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
	  double qua = p3*m1 - p1*m3;                 // "a"
	  double lin = p2*m1 - p1*m2 + p3*m0 - p0*m3; // "b"
	  double con = p2*m0 - p0*m2;                 // "c"
	  double xi[2];
	  if (fabs(qua) < kTol_Quad) { 
        xi[0] = -con/lin; }
	  else {
		double b2a = 0.5*lin/qua;
		double discrim = lin*lin - 4.0*con*qua;
		if (discrim < 0.0) { return 0; }
		else if (b2a > kTol_One ) 
		  { xi[0] = -b2a + sqrt(b2a*b2a - con/qua); }
		else if (b2a <-kTol_One ) 
		  { xi[0] = -b2a - sqrt(b2a*b2a - con/qua); }
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
	    } // CheckGap failed
	  } // CheckLocalCoordinates failed
	} // Normal opposition failed
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
	points = fIntegrationPoints;//this is dangerous

// error checking
	if ( fIntegrationPoints.MajorDim() != 4) throw;

	double jac[4];
	for (int i = 0; i < fIntegrationPoints.MajorDim(); i++) {
		jac[i] = ComputeJacobian(points(i));
	}
	
//.. weights for 8-point serendipity condensed on nodes
	weights[0] = (jac[3]+jac[1]+jac[0])/3.0;
	weights[1] = (jac[0]+jac[2]+jac[1])/3.0;
	weights[2] = (jac[1]+jac[3]+jac[2])/3.0;
	weights[3] = (jac[2]+jac[0]+jac[3])/3.0;
}

