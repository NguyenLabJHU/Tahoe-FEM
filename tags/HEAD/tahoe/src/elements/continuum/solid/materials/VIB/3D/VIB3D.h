/* $Id: VIB3D.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (04/20/1997)                                          */
/* 3D isotropic VIB solver.                                               */

#ifndef _VIB_3D_H_
#define _VIB_3D_H_

/* base class */
#include "NL_E_MatT.h"
#include "VIB_E_MatT.h"
#include "SpherePointsT.h"

/* forward declarations */
class dMatrixT;

class VIB3D: public NL_E_MatT, public VIB_E_MatT
{
public:

	/* constructor */
	VIB3D(ifstreamT& in, const ElasticT& element);

	/* destructor */
	~VIB3D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;	

	/* set angle offset - for testing onset of amorphous behavior
	 * Angles given in degrees */
	void SetAngles(double phi, double theta);

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
		
private:

	/* integration point generator */
	SpherePointsT*	fSphere;	
		
};

#endif /* _VIB_3D_H_ */
