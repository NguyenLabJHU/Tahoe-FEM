/* $Id: LineL2FaceT.h,v 1.10 2001-04-25 17:26:44 rjones Exp $ */

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
        void ComputeCentroid(double& centroid) const; 
	double ComputeRadius(void) const;
        void ComputeNormal(dArrayT& local_coordinates, double& normal) const; 
        void NodeNormal(int local_node_number, double& normal) const; 
	void CalcFaceNormal(void);
	void LocalBasis
		(double* normal, double* tangent1, double* tangent2) const;
	void ComputeShapeFunctions
		(dArrayT& local_coordinates, dArrayT& shape_functions) const;
	void ComputeShapeFunctions
		(dArrayT& local_coordinates, dMatrixT& shape_functions) const;
	double ComputeJacobian (dArrayT& local_coordinates) const;
        bool Projection (ContactNodeT* node, dArrayT& parameters) const ;
        inline void Polynomial
                (double* a, double* b) const ;

protected:

private:
	/* nodal coordinates */
	double*  fx[2];

};

inline void
LineL2FaceT::Polynomial
(double* a, double* b) const
{       /* const term */
        a[0] = 0.5*( fx[0][0]+fx[1][0]);
        a[1] = 0.5*( fx[0][1]+fx[1][1]);
        /* xi term */
        b[0] = 0.5*(-fx[0][0]+fx[1][0]);
        b[1] = 0.5*(-fx[0][1]+fx[1][1]);
}


#endif /* _LINEL2_FACE_T_H_ */

