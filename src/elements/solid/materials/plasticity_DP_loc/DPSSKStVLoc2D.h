/* $Id: DPSSKStVLoc2D.h,v 1.1 2004-03-20 23:35:32 raregue Exp $ */
/* created: myip (06/01/1999) */
#ifndef _DP_SS_KSTV_LOC_2D_H_
#define _DP_SS_KSTV_LOC_2D_H_

/* base class */
#include "Material2DT.h"
#include "DPSSKStVLoc.h"

namespace Tahoe {

class DPSSKStVLoc2D: public DPSSKStVLoc, public Material2DT
{
  public:

	/* constructor */
	DPSSKStVLoc2D(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
				const dSymMatrixT& totalstrain, 
				const ElementCardT& element, int ip);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& cdisc_ijkl(void);
  	
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

} // namespace Tahoe 
#endif /* _DP_SS_KSTV_LOC_2D_H_ */
