/* $Id: LineL2FaceT.h,v 1.24 2003-11-21 22:54:34 paklein Exp $ */
#ifndef _LINEL2_FACE_T_H_
#define _LINEL2_FACE_T_H_

/* base class */
#include "FaceT.h"

namespace Tahoe {

/**  connectivity
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
	void ComputeCentroid(double* centroid) const; 
	double ComputeRadius(void) const;
	void ComputeNormal
		(const double* local_coordinates, double* normal) const; 
	void ComputeTangent1 
		(const double* local_coordinates, double* tangent1) const;
	void NodeNormal(int local_node_number, double* normal) const; 
	void CalcFaceNormal(void);
	void LocalBasis
		(double* normal, double* tangent1) const;
	void ComputeShapeFunctions
	  	(const double* local_coordinates, dArrayT& shape_functions) 
		const;
	void ComputeShapeFunctions
	  	(const double* local_coordinates, dMatrixT& shape_functions) 
		const;
	void ComputeShapeFunctionDerivatives
                (const double* local_coordinates, dArrayT& shape_derivatives) 
		const;
	void ComputeShapeFunctionDerivatives
                (const double* local_coordinates, dMatrixT& shape_derivatives) 
		const;
	void InterpolatePosition 
		(const double* local_coordinates, double* x) const;
	double Interpolate 
		(const double* local_coordinates, dArrayT& nodal_values) const;
    double Interpolate
        (const double* local_coordinates, ArrayT<double*>& nodal_values) const;
	void InterpolateVector 
		(const double* local_coordinates, dArray2DT& nodal_vectors, 
		double* vector) const;
	double ComputeJacobian (const double* local_coordinates) const;
	bool Projection (ContactNodeT* node, const dArrayT& parameters) const ;
	void Quadrature (dArray2DT& points, dArrayT& weights) const;
		// points should be const

	inline void Polynomial
                (double* a, double* b) const ;
protected:

private:
	/* nodal coordinates */
	const double*  fx[2];
	
	/* integration points */  
	static dArray2DT fIntegrationPoints;

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


} // namespace Tahoe 
#endif /* _LINEL2_FACE_T_H_ */

