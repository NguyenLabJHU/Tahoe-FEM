/* $Id: SSSV_KStV2D.h,v 1.6 2004-07-15 08:28:50 paklein Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_SV_KStV_2D_H_
#define _SS_SV_KStV_2D_H_

#include "SSSimoViscoT.h"
#include "IsotropicT.h"
//#include "Material2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** base class for standard solid Kirchhoff St. Venant constitutive models 
 * constitutive law */
class SSSV_KStV2D: public SSSimoViscoT, public IsotropicT//, public Material2DT
{
	public:
	
	/*constructor*/
	SSSV_KStV2D(ifstreamT& in, const SSMatSupportT& support);

	virtual double StrainEnergyDensity(void);

	/* spatial description */ 
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli 
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress 
 
	/* material description */ 
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli 
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress 

	/*return internal dissipation variables*/
	virtual const dArrayT& InternalStressVars(void);
	virtual const dArrayT& InternalStrainVars(void);
			
	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

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
#endif  /* _SS_SV_KStV_2D_H_ */
