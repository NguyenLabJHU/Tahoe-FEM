/* $Id: QuadL4FaceT.h,v 1.6 2001-04-23 17:50:27 rjones Exp $ */

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
	void LocalBasis 
		(double* normal, double* tangent1, double* tangent2);

        void ComputeShapeFunctions
                (dArrayT& local_coordinates, dArrayT& shape_functions);
        void ComputeShapeFunctions
                (dArrayT& local_coordinates, dMatrixT& shape_functions);
        double ComputeJacobian (dArrayT& local_coordinates);
        bool Projection (ContactNodeT* node, dArrayT& parameters);

	inline void Polynomial(double* a, double* b, double* c, double* d);

protected:

private:
        /* nodal coordinates */
        double* fx[4];
	/* workspace */
	double a[3], b[3], c[3], d[3];
	double xi[2];


};

inline void
QuadL4FaceT::Polynomial
(double* a, double* b, double* c, double* d)
{	/* const term */
        a[0] = 0.25*( fx[0][0]+fx[1][0]+fx[2][0]+fx[3][0]);
        a[1] = 0.25*( fx[0][1]+fx[1][1]+fx[2][1]+fx[3][1]);
        a[2] = 0.25*( fx[0][2]+fx[1][2]+fx[2][2]+fx[3][2]);
	/* xi term */ 
        b[0] = 0.25*(-fx[0][0]+fx[1][0]+fx[2][0]-fx[3][0]);
        b[1] = 0.25*(-fx[0][1]+fx[1][1]+fx[2][1]-fx[3][1]);
        b[2] = 0.25*(-fx[0][2]+fx[1][2]+fx[2][2]-fx[3][2]);
	/* eta term */
        c[0] = 0.25*(-fx[0][0]-fx[1][0]+fx[2][0]+fx[3][0]);
        c[1] = 0.25*(-fx[0][1]-fx[1][1]+fx[2][1]+fx[3][1]);
        c[2] = 0.25*(-fx[0][2]-fx[1][2]+fx[2][2]+fx[3][2]);
	/* xi eta term */ 
        d[0] = 0.25*( fx[0][0]-fx[1][0]+fx[2][0]-fx[3][0]);
        d[1] = 0.25*( fx[0][1]-fx[1][1]+fx[2][1]-fx[3][1]);
        d[2] = 0.25*( fx[0][2]-fx[1][2]+fx[2][2]-fx[3][2]);
}

#endif /* _QUADL4_FACE_T_H_ */

