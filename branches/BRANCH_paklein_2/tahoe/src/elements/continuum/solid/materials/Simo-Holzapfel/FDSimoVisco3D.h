/* $Id: FDSimoVisco3D.h,v 1.2.6.1 2002-10-28 06:49:07 paklein Exp $ */
/* created:   TDN (5/31/2001) */
#ifndef _FD_SIMO_VISCO3D_H_
#define _FD_SIMO_VISCO3D_H_
 
#include "FDSimoViscoBaseT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

class FDSimoVisco3D: public FDSimoViscoBaseT
{
	public:

	/*constructor*/
	FDSimoVisco3D(ifstreamT& in, const FDMatSupportT& support);

	/*print parameters*/
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;
		
	/* spatial description */
	const dMatrixT& c_ijkl(void); // spatial tangent moduli
	const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	const dMatrixT& C_IJKL(void); // material tangent moduli
	const dSymMatrixT& S_IJ(void); // PK2 stress

	protected:

	enum Spring {kEquilibrium = 0, kNonEquilibrium = 1};

	virtual double Phi(const dMatrixT& Fbar, const double& J,
			   const int SpringType) = 0;
		
	virtual void Sig_dev(const dMatrixT& Fbar, const double& J,
		      dSymMatrixT& stress, const int SpringType) = 0;
	virtual void Devmod(const dMatrixT& Fbar, const double& J, 
			    dMatrixT& mod, const int SpringType) = 0;
	virtual double dUdJ(const double& J, const int SpringType) = 0;
	virtual double ddUddJ(const double& J, const int SpringType) = 0;

	protected:	 

	/*volumetric/deviatoric deformation measures*/
	double fJ;	
	dMatrixT fF;
	double fJ_E;
	double fJ_I;
	dMatrixT fFbar_E;
	dMatrixT fFbar_I;

	/*stress/modulus*/
	dMatrixT fModulus;
	dSymMatrixT fStress;

 	 /*relaxation times*/
	double ftauS;
	double ftauB;

	/* exp(-a* dt/tau)*/
	double falphaS;
	double falphaB;
	double fbetaS;
	double fbetaB;

	const double fthird;
};
}
#endif /*_FD_SIMO_VISCO3D_H_*/
