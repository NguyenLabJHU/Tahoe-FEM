/* $Id: GradSSSolidMatT.h,v 1.15 2004-09-02 18:25:04 rdorgan Exp $ */
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
class ElementCardT;

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
	const dMatrixT& Lambda(void) const;
	const dMatrixT& Lambda(int ip) const;
	const dMatrixT& Lambda_last(void) const;
	const dMatrixT& Lambda_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const dMatrixT& GradLambda(void) const;
	const dMatrixT& GradLambda(int ip) const;
	const dMatrixT& GradLambda_last(void) const;
	const dMatrixT& GradLambda_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const dMatrixT& LapLambda(void) const;
	const dMatrixT& LapLambda(int ip) const;
	const dMatrixT& LapLambda_last(void) const;
	const dMatrixT& LapLambda_last(int ip) const;
	/*@}*/

	/** \name spatial description */
	/*@{*/
	virtual const dMatrixT& odm_bh_ij(void) = 0;
	virtual const dMatrixT& odm_hb_ij(void) = 0;
	virtual const dMatrixT& gm_hh(void) = 0;
	virtual const dMatrixT& gm_hp(void) = 0;
	virtual const dMatrixT& gm_hq(void) = 0;
	virtual const dSymMatrixT& n_ij(void) = 0;
	virtual double yc(void) = 0;
	virtual double ys(void) = 0;
	virtual int weakened(void) = 0;
	virtual void UpdateWeakened(const ElementCardT& element, int ip) = 0;
	virtual void ResetWeakened(const ElementCardT& element, int ip) = 0;
	/*@}*/
	
	/** return the plastic multiplier in the material at the current integration point. */
	virtual void PMultiplier(dMatrixT& pmultiplier) { pmultiplier = Lambda(); };
	virtual void GradPMultiplier(dMatrixT& gradpmultiplier) { gradpmultiplier = GradLambda(); };
	virtual void LapPMultiplier(dMatrixT& lappmultiplier) { lappmultiplier = LapLambda(); };

protected:

	/** gradient small strain material support */
	const GradSSMatSupportT* fGradSSMatSupport;
};

} // namespace Tahoe 
#endif /* _GRAD_SS_SOLID_MAT_T_H_ */
