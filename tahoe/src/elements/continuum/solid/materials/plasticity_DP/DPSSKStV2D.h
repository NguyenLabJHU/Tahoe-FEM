/* $Id: DPSSKStV2D.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: myip (06/01/1999)                                             */

#ifndef _DP_SS_KSTV_2D_H_
#define _DP_SS_KSTV_2D_H_

/* base class */
#include "Material2DT.h"
#include "DPSSKStV.h"

class DPSSKStV2D: public DPSSKStV, public Material2DT
{
public:

	/* constructor */
	DPSSKStV2D(ifstreamT& in, const ElasticT& element);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

private:

	/* return values */
	dSymMatrixT	fStress2D;
	dMatrixT	fModulus2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

#endif /* _DP_SS_KSTV_2D_H_ */
