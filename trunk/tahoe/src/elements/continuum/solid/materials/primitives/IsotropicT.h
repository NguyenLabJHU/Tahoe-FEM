/* $Id: IsotropicT.h,v 1.8 2004-07-15 08:29:19 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _ISOTROPIC_T_H_
#define _ISOTROPIC_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "SolidMaterialT.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;
class ifstreamT;

class IsotropicT: virtual public ParameterInterfaceT
{
public:

	/* constructor */
	IsotropicT(ifstreamT& in);
	IsotropicT(void);

	/** \name set moduli */
	/*@{*/
	void Set_E_nu(double E, double nu);
	void Set_mu_kappa(double mu, double kappa);
	void Set_PurePlaneStress_mu_lambda(double mu, double lambda);
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	double Young(void) const;
	double Poisson(void) const;
	double Mu(void) const;
	double Kappa(void) const;
	double Lambda(void) const;
	/*@}*/
	
	/* print parameters */
	void Print(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* compute isotropic moduli tensor */
	void ComputeModuli(dMatrixT& moduli) const;
	void ComputeModuli2D(dMatrixT& moduli, SolidMaterialT::ConstraintT constraint) const;
	void ComputeModuli1D(dMatrixT& moduli) const;

	/* scale factor for constrained dilatation */
	double DilatationFactor2D(SolidMaterialT::ConstraintT constraint) const;   	

private:

	/** \name moduli */
	/*@{*/
	double fYoung;
	double fPoisson;
	double fMu;
	double fKappa;
	double fLambda;
	/*@}*/
};

/* inline functions */
inline double IsotropicT::Young(void) const { return fYoung; }
inline double IsotropicT::Poisson(void) const { return fPoisson; }
inline double IsotropicT::Mu(void) const { return fMu; }
inline double IsotropicT::Kappa(void) const { return fKappa; }
inline double IsotropicT::Lambda(void) const { return fLambda; }
} // namespace Tahoe 
#endif /* _ISOTROPIC_T_H_ */
