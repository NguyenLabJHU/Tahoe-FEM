#ifndef _GEOMODEL_SS_2D_T_H_
#define _GEOMODEL_SS_2D_T_H_

/* base class */
#include "GeoModelSST.h"

#include "SSSolidMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

class GeoModelSS2DT: public GeoModelSST
{
public:

	/* constructor */
	GeoModelSS2DT(void);

	/* initialization */
	//virtual void Initialize(void);
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
			const dSymMatrixT& totalstrain, 
			const ElementCardT& element, int ip);

	/* moduli */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& ce_ijkl(void);
	virtual const dMatrixT& con_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
	virtual const dMatrixT& con_perfplas_ijkl(void);
        
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

private:

	// pointer to material support
	const SSEnhLocMatSupportT* fSSEnhLocMatSupport;
	  
	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT fModulus2D;
	dMatrixT fModulusElas2D;
	dMatrixT fModulusPerfPlas2D;
	dMatrixT fModulusContinuum2D;
	dMatrixT fModulusContinuumPerfPlas2D;        

	/* work space */
	dSymMatrixT fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _GEOMODEL_SS_2D_T_H_ */
