/* $Id: GRAD_MRSSKStV2D.h,v 1.7 2005-08-05 22:26:12 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   Gradient Enhanced MR Model
*/
#ifndef _GRAD_MR_SS_KSTV_2D_H_
#define _GRAD_MR_SS_KSTV_2D_H_

/* base class */
#include "GRAD_MRSSKStV.h"

namespace Tahoe 
{

class GRAD_MRSSKStV2D: public GRAD_MRSSKStV
{
  public:

	/* constructor */
	GRAD_MRSSKStV2D(void);
	
	/* returns 3D strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
		
	/* returns 3D laplacian of strain (3D) */
	virtual const dSymMatrixT& LapElasticStrain(
                const dSymMatrixT& laptotalstrain, 
		const ElementCardT& element, int ip);
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
	
  	/*@{*/
	virtual const dMatrixT& c_UU1_ijkl(void);
	virtual const dMatrixT& c_UU2_ijkl(void);
	virtual const dMatrixT& c_ULam1_ij(void);
	virtual const dMatrixT& c_ULam2_ij(void);
	virtual const dMatrixT& c_LamU1_ij(void);
	virtual const dMatrixT& c_LamU2_ij(void);
	virtual const dMatrixT& c_LamLam1(void);
	virtual const dMatrixT& c_LamLam2(void);
	/*@{*/
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);
	
	/* yield function */
	virtual const double& YieldF(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

  private:

	/** \name 2D/3D dimension transformation for coupled moduli 
		C_ULambda1, C_ULambda2, C_LambdaU1, C_LambdaU2 */
	/*@{*/
	/** fill undefined elements with zeroes */
	void ReduceOffDiagonalModulus(const dMatrixT& mat3D, dMatrixT& mat2D);
	/*@}*/
	
  private:
  
  	/* return values */
  	dSymMatrixT	fStress2D;
  	dMatrixT	fModulus2D, fModulusPerfPlas2D;
  	dMatrixT    fModulusUU1_2D, fModulusUU2_2D;
    dMatrixT    fModulusULam1_2D, fModulusULam2_2D;
    dMatrixT    fModulusLamU1_2D, fModulusLamU2_2D;
    dMatrixT    fModulusLamLam1_2D, fModulusLamLam2_2D;
    dMatrixT    fTemp2DA, fTemp2DB;
  	double      fYieldFunction2D; //yield function

	/* work space */
	dSymMatrixT	fTotalStrain3D, fLapTotalStrain3D;
};

inline void GRAD_MRSSKStV2D::ReduceOffDiagonalModulus(const dMatrixT& mat3D, dMatrixT& mat2D)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (mat2D.Rows() != 3 || mat2D.Cols() != 1) throw ExceptionT::kGeneralFail;
	if (mat3D.Rows() != 6 || mat3D.Cols() != 1) throw ExceptionT::kGeneralFail;
#endif

	double* p = mat2D.Pointer();
	*p++ = mat3D[0]; // 1,1
	*p++ = mat3D[1]; // 2,1
	*p   = mat3D[5]; // 3,1
}

} // namespace Tahoe 
#endif /* _GRAD_MR_SS_KSTV_2D_H_ */