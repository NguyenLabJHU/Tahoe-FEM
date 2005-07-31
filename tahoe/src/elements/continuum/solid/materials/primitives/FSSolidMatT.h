/* $Id: FSSolidMatT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Defines the interface large strain materials which account             */
/* for thermal strains with the multiplicative split:                     */
/* F_tot = F_(the rest)*F_thermal                                         */
/* Note: At this point, no optimizations are added to compute             */
/* the inverse of F_thermal only once per time step.                      */

#ifndef _FD_STRUCT_MAT_T_H_
#define _FD_STRUCT_MAT_T_H_

/* base class */
#include "FDContinuumT.h"
#include "StructuralMaterialT.h"

/* forward declarations */
class ShapeFunctionT;

class FSSolidMatT: protected FDContinuumT, public StructuralMaterialT
{
public:

	/* constructor */
	FSSolidMatT(ifstreamT& in, const ElasticT& element);

	/* I/O functions */
	virtual void PrintName(ostream& out) const;

	/* required parameter flags */
	virtual bool NeedDisp(void) const;

	/* the shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;

	/* initialization */
	virtual void Initialize(void);

	/* strains/deformation measures */
	const dMatrixT& F(void); // deformation gradient
	const dMatrixT& F(const LocalArrayT& disp); 	
	const dSymMatrixT& C(void); // right stretch
	const dSymMatrixT& b(void); // left stretch
	const dSymMatrixT& E(void); // Green-Lagrange strain
	
	/* general spatial gradients */

	/* Test for localization using "current" values for Cauchy
	 * stress and the spatial tangent moduli. Returns 1 if the
	 * determinant of the acoustic tensor is negative and returns
	 * the normal for which the determinant is minimum. Returns 0
	 * of the determinant is positive. */
	virtual int IsLocalized(dArrayT& normal);

	/* apply pre-conditions at the current time step: compute
	 * thermal dilatation correction */
	virtual void InitStep(void);
	
protected:

	/* return the acoustical tensor and wave speeds */
	virtual const dSymMatrixT& AcousticalTensor(const dArrayT& normal);

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);

	/* acoustical tensor routines */
	void ComputeQ_2D(const dMatrixT& CIJKL, const dSymMatrixT& SIJ,
		const dMatrixT& FkK, const dArrayT& N, dSymMatrixT& Q) const;
	void ComputeQ_3D(const dMatrixT& CIJKL, const dSymMatrixT& SIJ,
		const dMatrixT& FkK, const dArrayT& N, dSymMatrixT& Q) const;

private:

	/* shape functions */
	const ShapeFunctionT& fShapes;
	
	/* nodal displacements */
	const LocalArrayT& fLocDisp;

	/* work space */
	dSymMatrixT fQ;  // return value
	dMatrixT fGradU; // displacement gradient matrix

	/* multiplicative thermal dilatation F */
	dMatrixT fFtherminverse;		
};

/* inlines */
inline const ShapeFunctionT& FSSolidMatT::ShapeFunction(void) const { return fShapes; }

#endif /* _FD_STRUCT_MAT_T_H_ */
