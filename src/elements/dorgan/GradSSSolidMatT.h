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
class DorganVoyiadjisMarinT;

/** defines the interface for gradient dependent small strain continuum materials */
class GradSSSolidMatT: public SSSolidMatT
{
public:

	/** constructor */
	GradSSSolidMatT(ifstreamT& in, const GradSSMatSupportT& support);

	/** destructor */
	~GradSSSolidMatT(void);

	virtual void Initialize(void);

	/* I/O functions */
	virtual void PrintName(ostream& out) const;

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);

	/** \name isotropic hardening */
	/*@{*/
	double& R(void);
	double& R(int ip);
	/*@}*/

	/** \name isotropic hardening from the end of the previous time step */
	/*@{*/
	double& R_last(void);
	double& R_last(int ip);
	/*@}*/

	/** \name Laplacian isotropic hardening */
	/*@{*/
	double& LaplacianR(void);
	double& LaplacianR(int ip);
	/*@}*/

	/** \name Laplacian isotropic hardening from the end of the previous time step */
	/*@{*/
	double& LaplacianR_last(void);
	double& LaplacianR_last(int ip);
	/*@}*/

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void) = 0;

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void) = 0;

	/** off diagonal moduli */
	virtual const dMatrixT& odm_ij(void) = 0;

	/** gradient dependent moduli */
	virtual const dMatrixT& gm(void) = 0;

	/** yield criteria moduli */
	virtual double yc(void) = 0;

	/** incremental change in R_bar */
	virtual double del_RBar(void) = 0;
	/*@}*/

protected:

	/** small strain material support */
	const GradSSMatSupportT& fGradSSMatSupport;
};

} // namespace Tahoe 
#endif /* _GRAD_SS_SOLID_MAT_T_H_ */
