/* $Id: FossumSSIso2DT.h,v 1.12 2004-08-04 02:28:46 cfoster Exp $ */
#ifndef _FOSSUM_SS_ISO_2D_T_H_
#define _FOSSUM_SS_ISO_2D_T_H_

/* base class */
//#include "Material2DT.h"
#include "FossumSSIsoT.h"

#include "SSSolidMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

class FossumSSIso2DT: public FossumSSIsoT//, public Material2DT
{
public:

	/* constructor */
	FossumSSIso2DT(void);

	/* initialization */
	//virtual void Initialize(void);
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	//virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
			const dSymMatrixT& totalstrain, 
			const ElementCardT& element, int ip);

	/* moduli */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& con_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
	virtual const dMatrixT& con_perfplas_ijkl(void);
        
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

private:
  
	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT fModulus2D;
	dMatrixT fModulusPerfPlas2D;
	dMatrixT fModulusContinuum2D;
	dMatrixT fModulusContinuumPerfPlas2D;        

	/* work space */
	dSymMatrixT fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _FOSSUM_SS_ISO_2D_T_H_ */
