/* $Id: GradSSSolidMatT.h,v 1.14 2004-08-05 17:26:11 rdorgan Exp $ */
#ifndef _GRAD_SS_SOLID_MAT_T_H_
#define _GRAD_SS_SOLID_MAT_T_H_

/* base class */
#include "SSSolidMatT.h"

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class GradSSMatSupportT;

/** defines the interface for gradient dependent small strain continuum materials */
class GradSSSolidMatT: public SSSolidMatT
{
public:

	/** constructor */
	GradSSSolidMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetGradSSMatSupport(const GradSSMatSupportT* support);

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);
	
	/** \name field */
	/*@{*/
	const double& Lambda(void) const;
	const double& Lambda(int ip) const;
	const double& Lambda_last(void) const;
	const double& Lambda_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const double& GradLambda(void) const;
	const double& GradLambda(int ip) const;
	const double& GradLambda_last(void) const;
	const double& GradLambda_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const double& LapLambda(void) const;
	const double& LapLambda(int ip) const;
	const double& LapLambda_last(void) const;
	const double& LapLambda_last(int ip) const;
	/*@}*/

	/** \name spatial description */
	/*@{*/
	virtual const dMatrixT& odm_bh_ij(void) = 0;
	virtual const dMatrixT& odm_hb_ij(void) = 0;
	virtual const dMatrixT& gm_hh(void) = 0;
	virtual const dMatrixT& gm_hp(void) = 0;
	virtual const dMatrixT& gm_hq(void) = 0;
	virtual double yc(void) = 0;
	virtual int weakened(void) = 0;
	/*@}*/
	
	/** return the plastic multiplier in the material at the current integration point. 
	 * Returns the plastic multiplier. */
	virtual void PMultiplier(dSymMatrixT& pmultiplier) { pmultiplier = Lambda(); };
	virtual void GradPMultiplier(dSymMatrixT& gradpmultiplier) { gradpmultiplier = GradLambda(); };
	virtual void LapPMultiplier(dSymMatrixT& lappmultiplier) { lappmultiplier = LapLambda(); };

protected:

	/** gradient small strain material support */
	const GradSSMatSupportT* fGradSSMatSupport;
};

} // namespace Tahoe 
#endif /* _GRAD_SS_SOLID_MAT_T_H_ */
