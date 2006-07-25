/* $Id: SSLinearVE2D.h,v 1.3 2004-07-15 08:29:34 paklein Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_LINEAR_VE_2D_H_
#define _SS_LINEAR_VE_2D_H_

/* base class */
#include "SSViscoelasticityT.h"

namespace Tahoe {

/** base class for standard solid Kirchhoff St. Venant constitutive models 
 * constitutive law */
class SSLinearVE2D: public SSViscoelasticityT
{
	public:
	
	/** constructor */
	SSLinearVE2D(void);
		
	virtual double StrainEnergyDensity(void);

	/* spatial description */ 
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli 
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress 
	
	/* material description */ 
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli 
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress 

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	 
    protected: 
	
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
