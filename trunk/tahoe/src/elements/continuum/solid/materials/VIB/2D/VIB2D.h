/* $Id: VIB2D.h,v 1.6 2003-01-29 07:34:52 paklein Exp $ */
/* created: paklein (04/09/1997) */
#ifndef _VIB_2D_H_
#define _VIB_2D_H_

/* base classes */
#include "NL_E_Mat2DT.h"
#include "VIB_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;
class CirclePointsT;

/** 2D VIB solver */
class VIB2D: public NL_E_Mat2DT, public VIB_E_MatT
{
public:

/* point generator codes */
	enum GeneratorCodeT {kEvenSpace = 0,
	                     KGaussRule = 1};

	/* constructor */
	VIB2D(ifstreamT& in, const FSMatSupportT& support);

	/* destructor */
	virtual ~VIB2D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;	

	/* set angle offset - for testing onset of amorphous behavior */
	void SetAngle(double angleoffset);

//TEMP - microscopic test of stability
	/* compute F.d for perturbation of magnidute epsilon along
	 * evenly spaced intervals of number AtoF.Length() for the
	 * current state of deformation */
	void Perturb(dArrayT& dU, double eps);

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
	
private:

	/* integration point generator */
	CirclePointsT*	fCircle;	

};

} // namespace Tahoe 
#endif /* _VIB_2D_H_ */
