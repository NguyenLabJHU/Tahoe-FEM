/* $Id: IsotropicT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _ISOTROPIC_T_H_
#define _ISOTROPIC_T_H_

#include "Environment.h"

/* forward declarations */
class dMatrixT;
#include "ios_fwd_decl.h"
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
	
	/* accessors */
	double Young(void) const;
	double Poisson(void) const;
	
	/* returns the Lame constants (calculated from E, nu) */
	void Lame(double& mu, double& lambda) const;

	/* shear and bulk moduli */
	void MuAndKappa(double& mu, double& kappa) const;
	double Mu(void) const;
	double Lambda(void) const;

	/* print parameters */
	void Print(ostream& out) const;

protected:

	/* compute isotropic moduli tensor */
	void ComputeModuli(dMatrixT& moduli, double mu, double lambda) const;

private:

	double	fYoung;
	double	fPoisson;
		
};

/* inline functions */
inline double IsotropicT::Mu(void) const
{
	return 0.5*fYoung/(1.0 + fPoisson);
}

inline double IsotropicT::Lambda(void) const
{
	return fYoung*fPoisson/((1.0 + fPoisson)*(1.0 - 2.0*fPoisson));
}

#endif /* _ISOTROPIC_T_H_ */
