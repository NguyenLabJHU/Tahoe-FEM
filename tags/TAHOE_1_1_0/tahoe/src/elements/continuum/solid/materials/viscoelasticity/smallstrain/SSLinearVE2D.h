/* $Id: SSLinearVE2D.h,v 1.1 2003-04-05 20:05:37 thao Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_LINEAR_VE_2D_H_
#define _SS_LINEAR_VE_2D_H_

#include "SSViscoelasticityT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class IsotropicT;

/** base class for standard solid Kirchhoff St. Venant constitutive models 
 * constitutive law */
class SSLinearVE2D: public SSViscoelasticityT
{
	public:
	
	/*constructor*/
	SSLinearVE2D(ifstreamT& in, const SSMatSupportT& support);
		
	/*print parameters*/
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;

	virtual double StrainEnergyDensity(void);

	/* spatial description */ 
	const dMatrixT& c_ijkl(void); // spatial tangent moduli 
	const dSymMatrixT& s_ij(void); // Cauchy stress 
	
	/* material description */ 
	const dMatrixT& C_IJKL(void); // material tangent moduli 
	const dSymMatrixT& S_IJ(void); // PK2 stress 

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	 
    protected: 
	
	/*1/3*/
	const double fthird;
        /*strain energy potentials*/ 
	dArrayT fMu;
	dArrayT fKappa;

	/*stress/modulus*/ 
	dMatrixT fModulus; 
	dSymMatrixT fStress;
    
	/*work spaces*/
	dMatrixT fModMat; 
	dSymMatrixT fStress3D;
	dSymMatrixT fStrain3D;
 
	/*relaxation times*/ 
	double ftauS; 
	double ftauB; 

	/* exp(-a* dt/tau)*/ 
	double falphaS; 
	double falphaB; 
	double fbetaS; 
	double fbetaB; 
};
}
#endif  /* _SS_LINEAR_VE_2D_H_ */
