/* $Id: LineL2FaceT.h,v 1.3 2001-04-11 14:48:58 rjones Exp $ */

#ifndef _LINEL2_FACE_T_H_
#define _LINEL2_FACE_T_H_

/* base class */
#include "FaceT.h"

/* direct members */

/* forward declarations */
class SurfaceT;
class iArrayT;
class dArrayT;
class dMatrixT;

/*  connectivity
 *  1--2
 */

class LineL2FaceT : public FaceT
{
public:

        /* constructor */
        LineL2FaceT
		(SurfaceT& surface,
		dArray2DT& surface_coordinates,
		int num_face_nodes,
		int* connectivity);

        /* destructor */
        ~LineL2FaceT(void);

        /* initialization after construction */
        void Initialize(void);

        /* geometric computation */
        void ComputeCentroid(double& centroid); 
	double ComputeRadius();
        void ComputeNormal(dArrayT& local_coordinates, double& normal); 
#if 0
        void ComputeTangents // ?????????
		(double& local_coordinates, double& tangent1,double& tangent2); 
#endif
	void ComputeShapeFunctions
		(dArrayT& local_coordinates, dArrayT& shape_functions);
	void ComputeShapeFunctions
		(dArrayT& local_coordinates, dMatrixT& shape_functions);
#if 0
	void ComputeShapeFunctionDerivatives
		(ArrayT& local_coordinates, ArrayT& shape_derivatives);
	void ComputeShapeFunctionDerivatives
		(ArrayT& local_coordinates, MatrixT& shape_derivatives);
#endif
	double ComputeJacobian (dArrayT& local_coordinates);
        bool Projection
		(double& point, double& normal, 
		dArrayT& local_coordinates, double gap); 
#if 0
        bool Projection
		(double& point, 
		double& local_coordinates, double gap); 
#endif
protected:

private:
	/* nodal coordinates */
	double*  fx[2];

};

#endif /* _LINEL2_FACE_T_H_ */

