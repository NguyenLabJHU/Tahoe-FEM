/* $Id: QuadL4FaceT.h,v 1.21 2003-11-21 22:54:35 paklein Exp $ */
#ifndef _QUADL4_FACE_T_H_
#define _QUADL4_FACE_T_H_

/* base class */
#include "FaceT.h"

namespace Tahoe {

/**  connectivity
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
        void ComputeCentroid(double* centroid) const; 
		double ComputeRadius() const;
        void ComputeNormal
		(const double* local_coordinates, double* normal) const;
        void ComputeTangent1 
		(const double* local_coordinates, double* tangent1) const;
        void ComputeTangent2 
		(const double* local_coordinates, double* tangent2) const;
        void NodeNormal(int local_node_number, double* normal) const;
        void CalcFaceNormal(void);
		void LocalBasis 
			(double* normal, double* tangent1, double* tangent2) const;
        void ComputeShapeFunctions
            	(const double* local_coordinates, dArrayT& shape_functions) 
				const;
        void ComputeShapeFunctions
          	(const double* local_coordinates, dMatrixT& shape_functions) 
			const;
        void ComputeShapeFunctionDerivatives
                (const double* local_coordinates, 
				dArrayT& shape_derivatives1,  dArrayT& shape_derivatives2) 
				const;
        void ComputeShapeFunctionDerivatives
                (const double* local_coordinates, 
				dMatrixT& shape_derivatives1, dMatrixT& shape_derivatives2) 
				const;
        void InterpolatePosition
                (const double* local_coordinates, double* x) const;
        double Interpolate 
		(const double* local_coordinates, dArrayT& nodal_values) const;
        double Interpolate
        (const double* local_coordinates, ArrayT<double*>& nodal_values) const ;
        void InterpolateVector 
		(const double* local_coordinates, dArray2DT& nodal_values, 
		double* vector) const;
        double ComputeJacobian (const double* local_coordinates) const;
        bool Projection (ContactNodeT* node, const dArrayT& parameters)  const;
	inline void Polynomial 
		(double* a, double* b, double* c, double* d) const ;
        void Quadrature
                (dArray2DT& points, dArrayT& weights) const;


protected:

private:
        /* nodal coordinates */
        const double* fx[4];
	
	/* integration points */
	static dArray2DT fIntegrationPoints;
};

inline void
QuadL4FaceT::Polynomial
(double* a, double* b, double* c, double* d) const
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
	/* xi.eta term */ 
        d[0] = 0.25*( fx[0][0]-fx[1][0]+fx[2][0]-fx[3][0]);
        d[1] = 0.25*( fx[0][1]-fx[1][1]+fx[2][1]-fx[3][1]);
        d[2] = 0.25*( fx[0][2]-fx[1][2]+fx[2][2]-fx[3][2]);
}

} // namespace Tahoe 
#endif /* _QUADL4_FACE_T_H_ */

