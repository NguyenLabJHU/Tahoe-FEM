/* $Id: SSSV_KStV3D.h,v 1.2 2003-04-05 20:38:08 thao Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_SV_KStV_3D_H_
#define _SS_SV_KStV_3D_H_

#include "SSSimoViscoT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class IsotropicT;

/** base class for standard solid Kirchhoff St. Venant constitutive models 
 * constitutive law */
class SSSV_KStV3D: public SSSimoViscoT
{
	public:
	
	/*constructor*/
	SSSV_KStV3D(ifstreamT& in, const SSMatSupportT& support);
		
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

	/*strain*/
	dSymMatrixT fe;

        /*stress/modulus*/ 
        dMatrixT fModulus; 
        dSymMatrixT fStress; 

	    dMatrixT fModMat;

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
#endif  /* _SS_SV_KStV_3D_H_ */
