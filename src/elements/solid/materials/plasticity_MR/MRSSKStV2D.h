/* created: Majid T. Manzari (04/16/2003) */
#ifndef _MR_SS_KSTV_2D_H_
#define _MR_SS_KSTV_2D_H_

/* base class */
//#include "Material2DT.h"
#include "MRSSKStV.h"

namespace Tahoe {

class MRSSKStV2D: public MRSSKStV//, public Material2DT
{
  public:

	/* constructor */
	MRSSKStV2D(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* returns 3D strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
	
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
#endif /* _MR_SS_KSTV_2D_H_ */