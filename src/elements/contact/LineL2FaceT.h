/* $Id: LineL2FaceT.h,v 1.6 2001-04-23 17:50:26 rjones Exp $ */

#ifndef _LINEL2_FACE_T_H_
#define _LINEL2_FACE_T_H_

/* base class */
#include "FaceT.h"

/* direct members */

/* forward declarations */

/*  connectivity
 * L  1--2  R  (outward normal up)
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
        void NodeNormal(int local_node_number, double& normal); 
	void FaceNormal(void);
	void LocalBasis(double* normal, double* tangent1, double* tangent2);
	void ComputeShapeFunctions
		(dArrayT& local_coordinates, dArrayT& shape_functions);
	void ComputeShapeFunctions
		(dArrayT& local_coordinates, dMatrixT& shape_functions);
	double ComputeJacobian (dArrayT& local_coordinates);
        bool Projection (ContactNodeT* node, dArrayT& parameters);
protected:

private:
	/* nodal coordinates */
	double*  fx[2];

};

#endif /* _LINEL2_FACE_T_H_ */

