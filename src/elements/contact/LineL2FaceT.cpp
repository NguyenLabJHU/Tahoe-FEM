/* $Id: LineL2FaceT.cpp,v 1.30 2003-11-21 22:54:34 paklein Exp $ */
#include "LineL2FaceT.h"

#include "ContactElementT.h"
#include "ContactNodeT.h"
#include "dArrayT.h"
#include "dMatrixT.h"

/* vector functions */
#include "vector2D.h"

using namespace Tahoe;

dArray2DT LineL2FaceT::fIntegrationPoints;

LineL2FaceT::LineL2FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
	FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	fNumVertexNodes = 2;
	if (!fIntegrationPoints.IsAllocated()){
        	fIntegrationPoints.Dimension(2,1);
        	double* ip;
        	ip = fIntegrationPoints(0);
        	ip[0] = -1.0 ;
        	ip = fIntegrationPoints(1);   
        	ip[0] =  1.0 ; 
	}
	fGeometryType = GeometryT::kLine;

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
LineL2FaceT::ComputeCentroid(double* centroid) const
{
	Ave(fx[0],fx[1],centroid); 
}

double
LineL2FaceT::ComputeRadius(void) const
{
	double diagonal[2];
	Diff (fx[0],fx[1],diagonal);
	double radius = 0.5* Mag(diagonal);
	return radius;
}

void
LineL2FaceT::NodeNormal(int local_node_number, double* normal) const
{ // could just retrieve face normal
#pragma unused(local_node_number)
	/* right to left */
	double t1[2];
	Diff(fx[0],fx[1],t1);
	RCross(t1,normal);
}

void
LineL2FaceT::CalcFaceNormal(void)
{
	/* right to left */
	double t1[2];
	Diff(fx[0],fx[1],t1);
	RCross(t1,fnormal);
}


void
LineL2FaceT::ComputeNormal
(const double* local_coordinates, double* normal) const
{
#pragma unused(local_coordinates)
	double t1[2];
	Diff(fx[0],fx[1],t1);
	/* this assumes a CW parameterization of the boundary */
	RCross(t1,normal);
	Normalize(normal);
}

void
LineL2FaceT::ComputeTangent1
(const double* local_coordinates,double* tangent1) const
{
#pragma unused(local_coordinates)

	Diff(fx[0],fx[1],tangent1);
	Scale(0.5,tangent1);
}


void
LineL2FaceT::ComputeShapeFunctions
(const double* local_coordinates, dArrayT& shape_functions) const
{
	double xi  = local_coordinates[0];
	shape_functions[0] = 0.5 * (1.0 - xi );
	shape_functions[1] = 0.5 * (1.0 + xi );
}

void
LineL2FaceT::ComputeShapeFunctions
(const double* local_coordinates, dMatrixT& shape_functions) const
{
// Move this to ContactElement ??
	shape_functions = 0.0;
	dArrayT shape_f(2);
	ComputeShapeFunctions(local_coordinates, shape_f);
	int i = 0, j, k = 0;
	if (shape_functions.Rows()%shape_functions.Cols() != 0
	 || shape_functions.Rows()/shape_functions.Cols() != 2) 
		throw ExceptionT::kGeneralFail;
	while (i < shape_functions.Rows()) {
		for (j = 0; j < shape_functions.Cols(); j++) {
			shape_functions(i++,j) = shape_f[k];
		}
		k++;
	}
#if 0
        shape_functions(0,0) = shape_f[0];
        shape_functions(1,1) = shape_f[0];
        shape_functions(2,0) = shape_f[1];
        shape_functions(3,1) = shape_f[1];
#endif
}

void
LineL2FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, dArrayT& shape_derivatives) const
{
	double xi  = local_coordinates[0];
	shape_derivatives[0] = -0.5 ;
	shape_derivatives[1] =  0.5 ;

}

void
LineL2FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, dMatrixT& shape_derivatives) const
{
	shape_derivatives = 0.0;
	dArrayT shape_d(2);
	ComputeShapeFunctionDerivatives(local_coordinates, shape_d);
        shape_derivatives(0,0) = shape_d[0];
        shape_derivatives(1,1) = shape_d[0];
        shape_derivatives(2,0) = shape_d[1];
        shape_derivatives(3,1) = shape_d[1];
}

void
LineL2FaceT::InterpolatePosition
(const double* local_coordinates, double* x)
const
{
        dArrayT shape_f(2);
        ComputeShapeFunctions (local_coordinates, shape_f);
        x[0] = shape_f[0]*fx[0][0]
             + shape_f[1]*fx[1][0];
        x[1] = shape_f[0]*fx[0][1]
             + shape_f[1]*fx[1][1];
}



double
LineL2FaceT::Interpolate
(const double* local_coordinates, dArrayT& nodal_values) const
{
	dArrayT shape_f(2);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*nodal_values[0]
                     + shape_f[1]*nodal_values[1];
	return value;
}

double
LineL2FaceT::Interpolate
(const double* local_coordinates, ArrayT<double*>& nodal_values) const
{
    dArrayT shape_f(2);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*(*nodal_values[0])
                     + shape_f[1]*(*nodal_values[1]);
    return value;
}



void
LineL2FaceT::InterpolateVector
(const double* local_coordinates, dArray2DT& nodal_vectors, double* vector) 
const
{
	dArrayT shape_f(2);
	ComputeShapeFunctions (local_coordinates, shape_f);
	vector[0] = shape_f[0]*nodal_vectors(0)[0] 
	          + shape_f[1]*nodal_vectors(1)[0];
	vector[1] = shape_f[0]*nodal_vectors(0)[1] 
	          + shape_f[1]*nodal_vectors(1)[1];
}


double
LineL2FaceT::ComputeJacobian (const double* local_coordinates) const
{
#pragma unused(local_coordinates)

	double t1[2];
	Diff(fx[0],fx[1],t1);
	return 0.5*Mag(t1);
}

bool
LineL2FaceT::Projection 
(ContactNodeT* node, const dArrayT& parameters)  const
{
        double tol_g  = parameters[ContactElementT::kGapTol];
        double tol_xi = parameters[ContactElementT::kXiTol];

        const double* nm = node->Normal();
        /* check normal opposition */
//cout << node->Tag() << " dot " << Dot(nm,fnormal) ;
        if ( Dot(nm,fnormal) < 0.0 ) {
          const double* x0 = node->Position();
          /* compute local coordinates */
          double a[2], b[2];
          Polynomial(a,b);
          /* components */
          double a1,b1,x1;
          const double* t1 = node->Tangent1();
          x1 = Dot(x0,t1);
          a1 = Dot(a,t1); b1 = Dot(b,t1);
	  double xi;
	  xi = (x1 - a1)/b1;
//cout << " xi " << xi ;
          if( CheckLocalCoordinates(xi,tol_xi) ) {
            double x3 = Dot(x0,nm);
            double a3 = Dot(a,nm);
            double b3 = Dot(b,nm);
            /* compute gap */
            double g =  a3 + b3*xi - x3;
//cout << " g " << g ;
            if (CheckGap(g,tol_g) ) {
                /*assign opposite (chooses closest)*/
                bool isbetter = node->AssignOpposing(fSurface,*this,&xi,g);
//cout << " YES \n";
                return isbetter;
            }
          }
        }
//cout << '\n';
        return 0;


}

void
LineL2FaceT::LocalBasis
(double* normal, double* tangent1) const
{
	/* calculate face tangent */
  	Diff(fx[0],fx[1],tangent1);
	Proj(tangent1, normal, tangent1);
	Normalize(tangent1);

}

void
LineL2FaceT::Quadrature
(dArray2DT& points, dArrayT& weights) const
{
        points = fIntegrationPoints; // this is dangerous
        for (int i = 0; i < fIntegrationPoints.MajorDim(); i++) {
                weights[i] = ComputeJacobian(points(i));
        }
}

