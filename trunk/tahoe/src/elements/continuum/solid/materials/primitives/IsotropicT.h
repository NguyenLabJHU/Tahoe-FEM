/* $Id: IsotropicT.h,v 1.7 2003-08-16 01:31:48 rdorgan Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _ISOTROPIC_T_H_
#define _ISOTROPIC_T_H_

#include "Environment.h"

/* direct members */
#include "Material2DT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;
class ifstreamT;

class IsotropicT
{
public:

	/* constructor */
	IsotropicT(ifstreamT& in);
	IsotropicT(void);

	/* set moduli */
	void Set_E_nu(double E, double nu);
	void Set_mu_kappa(double mu, double kappa);
	void Set_PurePlaneStress_mu_lambda(double mu, double lambda);
	
	/* accessors */
	double Young(void) const;
	double Poisson(void) const;
	double Mu(void) const;
	double Kappa(void) const;
	double Lambda(void) const;
	
	/* print parameters */
	void Print(ostream& out) const;

protected:

	/* compute isotropic moduli tensor */
	void ComputeModuli(dMatrixT& moduli) const;
	void ComputeModuli2D(dMatrixT& moduli, Material2DT::ConstraintOptionT constraint) const;
	void ComputeModuli1D(dMatrixT& moduli) const;

	/* scale factor for constrained dilatation */
	double DilatationFactor2D(Material2DT::ConstraintOptionT constraint) const;   	

private:

	/* moduli */
	double fYoung;
	double fPoisson;
	double fMu;
	double fKappa;
	double fLambda;
};

/* inline functions */
inline double IsotropicT::Young(void) const { return fYoung; }
inline double IsotropicT::Poisson(void) const { return fPoisson; }
inline double IsotropicT::Mu(void) const { return fMu; }
inline double IsotropicT::Kappa(void) const { return fKappa; }
inline double IsotropicT::Lambda(void) const { return fLambda; }
} // namespace Tahoe 
#endif /* _ISOTROPIC_T_H_ */
