/* $Id: LineL2FaceT.cpp,v 1.12 2001-04-30 19:30:19 rjones Exp $ */

#include "LineL2FaceT.h"
#include "FaceT.h"

#include "ContactElementT.h"
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
        fIntegrationPoints.Allocate(2,1);
        double* ip;
        ip = fIntegrationPoints(0);
        ip[0] = -1.0 ;
        ip = fIntegrationPoints(1);   
        ip[0] =  1.0 ; 

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
LineL2FaceT::ComputeCentroid(double& centroid) const
{
	Ave(fx[0],fx[1],&centroid); 
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
LineL2FaceT::NodeNormal(int local_node_number, double& normal) const
{
	int c = local_node_number;
	int n = Next(c);
	double t1[2];
	/* get sense of left and right */
	(n > c) ? Diff(fx[c],fx[n],t1) : Diff(fx[n],fx[c],t1);
	RCross(t1,&normal);
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
LineL2FaceT::ComputeNormal(double* local_coordinates, double& normal) const
{
	double t1[2];
	Diff(fx[0],fx[1],t1);
	/* this assumes a CW parameterization of the boundary */
	RCross(t1,&normal);
	Normalize(&normal);
}

void
LineL2FaceT::ComputeShapeFunctions
(double* local_coordinates, dArrayT& shape_functions) const
{
	double xi  = local_coordinates[0];
	shape_functions[0] = 0.5 * (1.0 - xi );
	shape_functions[1] = 0.5 * (1.0 + xi );
}

void
LineL2FaceT::ComputeShapeFunctions
(double* local_coordinates, dMatrixT& shape_functions) const
{
	dArrayT shape_f;
	ComputeShapeFunctions(local_coordinates, shape_f);
// MORE
}

#if 0
void
LineL2FaceT::Interpolate
(double* local_coordinates, dArrayT& nodal_values, double value);
{
        ComputeShapeFunctions
                (dArrayT& local_coordinates, dArrayT& shape_functions);
        value     = shape_function[0]*nodal_value[0]
                  + shape_function[1]*nodal_value[1];
}

void
LineL2FaceT::InterpolateVector
(double* local_coordinates, dArray2DT& nodal_vectors, double& vector);
{
	ComputeShapeFunctions
		(dArrayT& local_coordinates, dArrayT& shape_functions);
	vector[0] = shape_function[0]*nodal_vector[0][0] 
	          + shape_function[1]*nodal_vector[1][0];
	vector[1] = shape_function[0]*nodal_vector[0][1] 
	          + shape_function[1]*nodal_vector[1][1];
}
#endif


double
LineL2FaceT::ComputeJacobian (double* local_coordinates) const
{
	//HACK
	return 1.0;
}

bool
LineL2FaceT::Projection 
(ContactNodeT* node, dArrayT& parameters)  const
{
        double tol_g  = parameters[ContactElementT::kGapTol];
        double tol_xi = parameters[ContactElementT::kXiTol];

        const double* nm = node->Normal();
        /* check normal opposition */
        if ( Dot(nm,fnormal) < 0.0 ) {
          const double* x0 = node->Position();
          /* compute local coordinates */
          double a[3], b[3];
          Polynomial(a,b);
          /* components */
          double a1,b1,x1;
          const double* t1 = node->Tangent1();
          x1 = Dot(x0,t1);
          a1 = Dot(a,t1); b1 = Dot(b,t1);
	  double xi;
	  xi = (x1 - a1)/b1;
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
LineL2FaceT::LocalBasis
(double* normal, double* tangent1, double* tangent2) const
{
	/* calculate face tangent */
        Diff(fx[0],fx[1],tangent1);
        Proj(tangent1, normal, tangent1);
        Normalize(tangent1);
	
 	/* NOTE: nothing is done with tangent2 (NULL) */

}

void
LineL2FaceT::Quadrature
(dArray2DT& points, dArrayT& weights) const
{
        points = fIntegrationPoints;
        for (int i = 0; i < fIntegrationPoints.Length(); i++) {
                weights[i] = ComputeJacobian(points(i));
        }
}

