/* $Id: QuadL4FaceT.h,v 1.5 2001-04-19 23:47:01 rjones Exp $ */

#ifndef _QUADL4_FACE_T_H_
#define _QUADL4_FACE_T_H_

/* base class */
#include "FaceT.h"

/* direct members */

/* forward declarations */

/*  connectivity
 *  4--3
 *  |  |    (outward normal out-of-plane)
 *  1--2
 */

class QuadL4FaceT : public FaceT
{
public:

        /* constructor */
        QuadL4FaceT
		(SurfaceT& surface,
		dArray2DT& surface_coordinates,
		int num_face_nodes,	
		int* connectivity);

        /* destructor */
        ~QuadL4FaceT(void);

        /* initialization after construction */
        void Initialize(void);

        /* geometric computation */
        void ComputeCentroid(double& centroid); 
	double ComputeRadius();
        void ComputeNormal(dArrayT& local_coordinates, double& normal);
        void NodeNormal(int local_node_number, double& normal);
        void FaceNormal(void);
        void ComputeShapeFunctions
                (dArrayT& local_coordinates, dArrayT& shape_functions);
        void ComputeShapeFunctions
                (dArrayT& local_coordinates, dMatrixT& shape_functions);
        double ComputeJacobian (dArrayT& local_coordinates);
        bool Projection (ContactNodeT* node, dArrayT& parameters);

protected:

private:
        /* nodal coordinates */
        double* fx[4];


};

#endif /* _QUADL4_FACE_T_H_ */

