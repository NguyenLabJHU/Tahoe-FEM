/* $Id: GradSSSolidMatT.h,v 1.3 2003-09-29 19:58:57 rdorgan Exp $ */
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

	/** number of degrees of freedom for iso_hard field (per node) in the host
	 * element group. */
	int NumDOF_R(void) const;

	/** number of total degrees of freedom (per node) in the host
	 * element group. */
	int NumDOF_Total(void) const;

	/** the total number of integration points per element in the
	 * host element group for iso_hard field. */
	int NumIP_R(void) const;

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

	/** incremental change in R */
	virtual double del_R(void) = 0;

	/** incremental change in LaplacianR */
	virtual double del_LaplacianR(void) = 0;
	/*@}*/

protected:

	/** small strain material support */
	const GradSSMatSupportT& fGradSSMatSupport;

	/** number of degrees of freedom for iso_hard field*/
	int fNumDOF_R;

	/** total number of degrees of freedom */
	int fNumDOF_Total;
	
	/** number of integration points for iso_hard field */
	int fNumIP_R;
};

/* inlines */
inline int GradSSSolidMatT::NumDOF_R(void) const { return fNumDOF_R; }
inline int GradSSSolidMatT::NumDOF_Total(void) const { return fNumDOF_Total; }
inline int GradSSSolidMatT::NumIP_R(void) const { return fNumIP_R; }

} // namespace Tahoe 
#endif /* _GRAD_SS_SOLID_MAT_T_H_ */
