/* created: Majid T. Manzari (04/16/2003) */
#ifndef _GRAD_MR_SS_KSTV_2D_H_
#define _GRAD_MR_SS_KSTV_2D_H_

/* base class */
#include "Material2DT.h"
#include "GRAD_MRSSKStV.h"

namespace Tahoe 
{

class GRAD_MRSSKStV2D: public GRAD_MRSSKStV, public Material2DT
{
  public:

	/* constructor */
	GRAD_MRSSKStV2D(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returns 3D strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
		
	/* returns 3D gradient of strain (3D) */
	virtual const dSymMatrixT& GradElasticStrain(
                const dSymMatrixT& del2_totalstrain, 
		const ElementCardT& element, int ip);
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& cdisc_ijkl(void);
  	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);
	
	/* yield function */
	virtual const double& Yield_Function(void);

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

  private:
  
  	/* return values */
  	dSymMatrixT	fStress2D;
  	dMatrixT	fModulus2D;
  	double      fYieldFunction2D; //yield function

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _GRAD_MR_SS_KSTV_2D_H_ */