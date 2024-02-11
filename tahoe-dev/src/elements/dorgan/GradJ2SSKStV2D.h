/* $Id: GradJ2SSKStV2D.h,v 1.1 2004-09-02 18:25:04 rdorgan Exp $ */
#ifndef _GRAD_J2_SS_KSTV_2D_H_
#define _GRAD_J2_SS_KSTV_2D_H_

/* base classes */
#include "GradJ2SSKStV.h"

namespace Tahoe {

class GradJ2SSKStV2D: public GradJ2SSKStV
{
public:

	/** constructor */
	GradJ2SSKStV2D(void);

	/** returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);
	
	/** modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/** off diagonal moduli for Kar */
	virtual const dMatrixT& odm_bh_ij(void);
	
	/** off diagonal moduli for Kra */
	virtual const dMatrixT& odm_hb_ij(void);
	
	/** modulus for second term in Krr */
	virtual const dMatrixT& gm_hp(void);
	
	/** stress */
	virtual const dSymMatrixT& s_ij(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name 2D/3D dimension transformation for GradientModuli_hp */
	/*@{*/
	/** fill undefined elements with zeroes */
	void ReduceGradientModuli_hp(const dMatrixT& mat3D, dMatrixT& mat2D);
	void ReduceOffDiagonalModulus_bh(const dMatrixT& mat3D, dMatrixT& mat2D);
	/*@}*/

private:

	/* return values */
	dSymMatrixT	fStress2D;
	dMatrixT	fModulus2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
	dMatrixT fTensorTemp1;
	dMatrixT fOffDiagonalModulus_bh_2D, fOffDiagonalModulus_hb_2D;
	dMatrixT fGradientModulus_hp_2D;
};

inline void GradJ2SSKStV2D::ReduceGradientModuli_hp(const dMatrixT& mat3D, dMatrixT& mat2D)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (mat2D.Rows() != 1 || mat2D.Cols() != 2) throw ExceptionT::kGeneralFail;
	if (mat3D.Rows() != 1 || mat3D.Cols() != 3) throw ExceptionT::kGeneralFail;
#endif

	double* p = mat2D.Pointer();
	*p++ = mat3D[0];
	*p = mat3D[1];
}
 
inline void GradJ2SSKStV2D::ReduceOffDiagonalModulus_bh(const dMatrixT& mat3D, dMatrixT& mat2D)
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
#endif /* _GRAD_J2_SS_KSTV_2D_H_ */
