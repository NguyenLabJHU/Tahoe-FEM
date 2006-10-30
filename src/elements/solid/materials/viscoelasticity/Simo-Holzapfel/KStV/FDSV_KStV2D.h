/* $Id: FDSV_KStV2D.h,v 1.1 2006-10-30 23:32:06 thao Exp $ */
/* created:   TDN (5/31/2001) */
#ifndef _FD_SV_KStV_2D_H_
#define _FD_SV_KStV_2D_H_

#include "FDSimoViscoBaseT.h"
//#include "Material2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class IsotropicT;

/** base class for standard solid Kirchhoff St. Venant constitutive models 
 * constitutive law */
class FDSV_KStV2D: public FDSimoViscoBaseT//, public Material2DT
{
	public:
	
	/*constructor*/
	FDSV_KStV2D(ifstreamT& in, const FSMatSupportT& support);

        /* spatial description */ 
        virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli 
        virtual const dSymMatrixT& s_ij(void); // Cauchy stress 
 
        /* material description */ 
        virtual const dMatrixT& C_IJKL(void); // material tangent moduli 
        virtual const dSymMatrixT& S_IJ(void); // PK2 stress 
 
        protected: 
	/*1/3*/
	const double fthird;

        /*strain energy potentials*/ 
	dArrayT fMu;
	dArrayT fKappa;

	/*strain*/
	dSymMatrixT fE;

        /*stress/modulus*/ 
        dMatrixT fModulus; 
        dMatrixT fModMat;
        dSymMatrixT fStress; 
 
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
#endif  /* _FD_SV_KStV_2D_H_ */
