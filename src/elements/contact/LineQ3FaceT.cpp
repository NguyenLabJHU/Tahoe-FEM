/* $Id: LineQ3FaceT.cpp,v 1.15 2003-11-21 22:54:34 paklein Exp $ */
#include "LineQ3FaceT.h"

#include "ContactElementT.h"
#include "ContactNodeT.h"
#include "dArrayT.h"
#include "dMatrixT.h"

/* vector functions */
#include "vector2D.h"

/* parameters */

using namespace Tahoe;

dArray2DT LineQ3FaceT::fIntegrationPoints;
static const double kTol_Zero  = 0.00000001;

LineQ3FaceT::LineQ3FaceT
(SurfaceT& surface, dArray2DT& surface_coordinates, 
int number_of_face_nodes, int* connectivity):
	FaceT(surface,surface_coordinates,
	number_of_face_nodes,connectivity)
{
	fNumVertexNodes = 2;
	if (!fIntegrationPoints.IsAllocated()){
        	fIntegrationPoints.Dimension(3,1);
        	double* ip;
        	ip = fIntegrationPoints(0);
        	ip[0] = -1.0 ;
        	ip = fIntegrationPoints(1);   
        	ip[0] =  1.0 ; 
        	ip = fIntegrationPoints(2);   
        	ip[0] =  1.0 ; 
	}
	fGeometryType = GeometryT::kLine;

}

LineQ3FaceT::~LineQ3FaceT (void)
{
}

void
LineQ3FaceT::Initialize(void)
{
        for (int i = 0; i < fConnectivity.Length(); i++) {
                fx[i] = fSurfaceCoordinates(fConnectivity[i]);
        }
}


void
LineQ3FaceT::ComputeCentroid(double* centroid) const
{
	centroid[0] = fx[2][0];
	centroid[1] = fx[2][1];
}

double
LineQ3FaceT::ComputeRadius(void) const
{
	double diagonal[2];
	Diff (fx[0],fx[1],diagonal);
	double radius = 0.5* Mag(diagonal);
	return radius;
}

void
LineQ3FaceT::NodeNormal(int local_node_number, double* normal) const
{
	double t1[2];
	switch(local_node_number)
	{
		case 0: // xi -1
			t1[0] = 0.0; t1[1] = 0.0;
			Add(t1,-1.5,fx[0]);
			Add(t1,-0.5,fx[1]);
			Add(t1, 2.0,fx[2]);
			break;
		case 1: // xi +1
			t1[0] = 0.0; t1[1] = 0.0;
			Add(t1, 0.5,fx[0]);
			Add(t1, 1.5,fx[1]);
			Add(t1,-2.0,fx[2]);
			break;
		case 2: // xi  0
			Diff (fx[0],fx[1],t1);
			break;
	}
	RCross(t1,normal);
	// normalize??
}

void
LineQ3FaceT::CalcFaceNormal(void)
{
	/* this assumes a CW parameterization of the boundary */
        /* right (-1) to left (+1) */
	double t1[2];
        Diff(fx[0],fx[1],t1);
        RCross(t1,fnormal);
}


void
LineQ3FaceT::ComputeNormal
(const double* local_coordinates, double* normal) const
{
	double t1[2];
	ComputeTangent1(local_coordinates,t1);
	/* this assumes a CW parameterization of the boundary */
	RCross(t1,normal);
	Normalize(normal);
}

void
LineQ3FaceT::ComputeTangent1
(const double* local_coordinates,double* tangent1) const
{
	double xi  = local_coordinates[0];
	tangent1[0] = 0.0;	// need a '=' operator
	tangent1[1] = 0.0;	
	Add(tangent1,(xi-0.5),fx[0]);
	Add(tangent1,(xi+0.5),fx[1]);
	Add(tangent1,-2.0*xi, fx[2]);
}

void
LineQ3FaceT::ComputeShapeFunctions
(const double* local_coordinates, dArrayT& shape_functions) const
{
	double xi  = local_coordinates[0];
	shape_functions[0] = 0.5 * xi * ( xi - 1.0 );
	shape_functions[1] = 0.5 * xi * ( xi + 1.0 );
	shape_functions[2] = 1.0 - xi * xi ;
}

void
LineQ3FaceT::ComputeShapeFunctions
(const double* local_coordinates, dMatrixT& shape_functions) const
{
	shape_functions = 0.0;
	dArrayT shape_f(3);
	ComputeShapeFunctions(local_coordinates, shape_f);
	shape_functions(0,0) = shape_f[0];
	shape_functions(1,1) = shape_f[0];
	shape_functions(2,0) = shape_f[1];
	shape_functions(3,1) = shape_f[1];
	shape_functions(4,0) = shape_f[2];
	shape_functions(5,1) = shape_f[2];
}

void
LineQ3FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, dArrayT& shape_derivatives) const
{
	double xi  = local_coordinates[0];
	shape_derivatives[0] =  xi - 0.5 ;
	shape_derivatives[1] =  xi + 0.5 ;
	shape_derivatives[2] = - 2.0 * xi ;
}

void
LineQ3FaceT::ComputeShapeFunctionDerivatives
(const double* local_coordinates, dMatrixT& shape_derivatives) const
{
	shape_derivatives = 0.0;
	dArrayT shape_d(3);
	ComputeShapeFunctions(local_coordinates, shape_d);
	shape_derivatives(0,0) = shape_d[0];
	shape_derivatives(1,1) = shape_d[0];
	shape_derivatives(2,0) = shape_d[1];
	shape_derivatives(3,1) = shape_d[1];
	shape_derivatives(4,0) = shape_d[2];
	shape_derivatives(5,1) = shape_d[2];
}

void
LineQ3FaceT::InterpolatePosition
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
}


double
LineQ3FaceT::Interpolate
(const double* local_coordinates, dArrayT& nodal_values)
const
{
	dArrayT shape_f(3);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*nodal_values[0]
                     + shape_f[1]*nodal_values[1] 
                     + shape_f[2]*nodal_values[2];
	return value;
}

double
LineQ3FaceT::Interpolate
(const double* local_coordinates, ArrayT<double*>& nodal_values)
const
{
	dArrayT shape_f(3);
        ComputeShapeFunctions (local_coordinates, shape_f);
        double value = shape_f[0]*(*nodal_values[0])
                     + shape_f[1]*(*nodal_values[1]) 
                     + shape_f[2]*(*nodal_values[2]);
	return value;
}

void
LineQ3FaceT::InterpolateVector
(const double* local_coordinates, dArray2DT& nodal_vectors, double* vector)
const
{
	dArrayT shape_f(3);
        ComputeShapeFunctions (local_coordinates, shape_f);
	vector[0] = shape_f[0]*nodal_vectors(0)[0] 
	          + shape_f[1]*nodal_vectors(1)[0] 
	          + shape_f[2]*nodal_vectors(2)[0];
	vector[1] = shape_f[0]*nodal_vectors(0)[1] 
	          + shape_f[1]*nodal_vectors(1)[1] 
	          + shape_f[2]*nodal_vectors(2)[1];
}



double
LineQ3FaceT::ComputeJacobian (const double* local_coordinates) const
{
#pragma unused(local_coordinates)

	//HACK
	// mag of tangent
	return 1.0;
}

bool
LineQ3FaceT::Projection 
(ContactNodeT* node, const dArrayT& parameters)  const
{
        double tol_g  = parameters[ContactElementT::kGapTol];
        double tol_xi = parameters[ContactElementT::kXiTol];

        const double* nm = node->Normal();
        /* check normal opposition */
        if ( Dot(nm,fnormal) < 0.0 ) {
          const double* x0 = node->Position();
          /* compute local coordinates */
          double a[2], b[2], c[2];
          Polynomial(a,b,c);
          /* components */
          double a1,b1,c1,x1;
          const double* t1 = node->Tangent1();
          x1 = Dot(x0,t1);
          a1 = Dot(a,t1); b1 = Dot(b,t1); c1 = Dot(c,t1);
	  double xi;
	  if (fabs(c1) < kTol_Zero) {
	    xi = (x1 - a1)/b1;
	  } else {
		double discrim = b1*b1 - 4.0*a1*c1;
		if ( discrim < 0.0 ) { // no intersection
		  xi = (x1 - a1)/b1;
		}
		discrim = 
		    -0.5*(b1 + ((b1>0.0) ? sqrt(discrim) : -sqrt(discrim)));
		double root_1 = discrim/a1;
		double root_2 = c1/discrim;
		// pick xi closer to zero
		xi = (fabs(root_1) < fabs(root_2)) ? root_1 : root_2;
	  }
          if( CheckLocalCoordinates(xi,tol_xi) ) {
            double x3 = Dot(x0,nm);
            double a3 = Dot(a,nm);
            double b3 = Dot(b,nm);
            /* compute gap */
            double g =  a3 + b3*xi - x3;
            if (CheckGap(g,tol_g) ) {
                /*assign opposite (chooses closest)*/
                bool isbetter = node->AssignOpposing(fSurface,*this,&xi,g);
                return isbetter;
            }
          }
        }
        return 0;


}

void
LineQ3FaceT::LocalBasis
(double* normal, double* tangent1) const
{
	/* calculate face tangent */
        Diff(fx[0],fx[1],tangent1);
        Proj(tangent1, normal, tangent1);
        Normalize(tangent1);
	
 	/* NOTE: nothing is done with tangent2 (NULL) */

}

void
LineQ3FaceT::Quadrature
(dArray2DT& points, dArrayT& weights) const
{
        points = fIntegrationPoints;
        for (int i = 0; i < fIntegrationPoints.MajorDim(); i++) {
                weights[i] = ComputeJacobian(points(i));
        }
}

