/* $Id: FSSolidMatT.h,v 1.17 2004-08-01 20:40:42 paklein Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _FD_STRUCT_MAT_T_H_
#define _FD_STRUCT_MAT_T_H_

/* base class */
#include "SolidMaterialT.h"
#include "TensorTransformT.h"

/* direct members */
#include "FSMatSupportT.h"

namespace Tahoe {

/* forward declarations */
class FiniteStrainT;

/** base class for finite deformation constitutive models. The interface
 * provides access to the element-computed deformation as well as
 * functions to compute the Green Lagragian strain tensor \b E
 * and the stretch tensors \b C and \b b. The class provides support
 * for a multiplicative thermal strain:\n
 * F_total = F_mechanical F_thermal\n
 * where the total deformation gradient is available through
 * FSSolidMatT::F_total, the "mechanical" part of the deformation
 * gradient is available through FSSolidMatT::F_total, and the
 * \a inverse of the thermal deformation gradient is available
 * through FSSolidMatT::F_thermal_inverse. */
class FSSolidMatT: public SolidMaterialT, protected TensorTransformT
{
public:

	/** constructor */
	FSSolidMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetFSMatSupport(const FSMatSupportT* support);

	/** finite strain materials support */
	const FSMatSupportT& FSMatSupport(void) const;

	/** \name tangent moduli 
	 * FSSolidMatT::c_ijkl computes the tangent moduli in the spatial representation
	 * using the finite difference approximation developed by Miehe, CMAME \b 134, 1996.
	 * FSSolidMatT::C_IJKL calls FSSolidMatT::c_ijkl and pulls the result back to the
	 * material representation. */
	/*@{*/
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& C_IJKL(void);
	/*@}*/

	/** compute the 2nd Piola-Kirchhoff stress by pulling back the result computed with
	 * SolidMaterialT::s_ij */
	const dSymMatrixT& S_IJ(void);
	
	/** test for localization. check for bifurvation using current
	 * Cauchy stress and the spatial tangent moduli.
	 * \param normal orientation of the localization if localized
	 * \return 1 if the determinant of the acoustical tensor is negative
	 * or 0 if the determinant is positive. */
	virtual int IsLocalized(dArrayT& normal);

	/** initialize current step. compute thermal dilatation */
	virtual void InitStep(void);

	/** close current step. store thermal dilatation */
	virtual void CloseStep(void);

	/** required parameter flag. Indicates whether the constitutive model
	 * requires calculation of the deformation gradient.
	 * \return true by default. */
	virtual bool Need_F(void) const { return true; };

	/** required parameter flag. Indicates whether the constitutive model
	 * requires the deformation gradient from the previous time increment.
	 * \return false by default. */
	virtual bool Need_F_last(void) const { return false; };

	/** total deformation gradient. \note This function is on its
	 * way out. Use FSSolidMatT::F_total */
	const dMatrixT& F(void) const; 

	/** total deformation gradient at the given integration point. \note This 
	 * function is on its way out. Use FSSolidMatT::F_total */
	const dMatrixT& F(int ip) const; 

	/** total deformation gradient */
	const dMatrixT& F_total(void) const; 

	/** total deformation gradient at the given integration point */
	const dMatrixT& F_total(int ip) const;

	/** mechanical part of the deformation gradient. The part of the
	 * deformation gradient not associated with an imposed thermal
	 * strain. */
	const dMatrixT& F_mechanical(void); 

	/** mechanical part of the deformation gradient at the given integration
	 * point. The part of the deformation gradient not associated with an 
	 * imposed thermal strain. */
	const dMatrixT& F_mechanical(int ip);

	/** total deformation gradient from end of previous step */
	const dMatrixT& F_total_last(void) const; 

	/** total deformation gradient at the given integration point 
	 * from end of previous step */
	const dMatrixT& F_total_last(int ip) const;
	
	/** inverse of the deformation gradient associated with the
	 * imposed thermal strain */
	const dMatrixT& F_thermal_inverse(void) const { return fF_therm_inv; };

	/** inverse of the deformation gradient associated with the
	 * imposed thermal strain from the previous time increment. */
	const dMatrixT& F_thermal_inverse_last(void) const { return fF_therm_inv_last; };

	/** mechanical part of the deformation gradient from the previous
	 * time increment. The part of the deformation gradient not associated 
	 * with an imposed thermal strain. */
	const dMatrixT& F_mechanical_last(void); 

	/** mechanical part of the deformation gradient from the previous
	 * time increment at the given integration point. The part of the 
	 * deformation gradient not associated 
	 * with an imposed thermal strain. */
	const dMatrixT& F_mechanical_last(int ip);

	/** return the strain in the material at the current integration point. 
	 * Returns the Green-Lagrangian strain. */
	virtual void Strain(dSymMatrixT& strain) { Compute_E(F_mechanical(), strain);}
	virtual void Stretch(dSymMatrixT& stretch) {Compute_C(F_mechanical(), stretch);}

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** enum for use by derived classes to to track last stress called */
	enum FrameT {kNone, kMaterial, kSpatial};

	/** left stretch tensor.
	 * \param F the deformation gradient
	 * \param b return value */
	void Compute_b(const dMatrixT& F, dSymMatrixT& b) const;

	/** right stretch tensor. 
	 * \param F the deformation gradient
	 * \param C return value */
	void Compute_C(const dMatrixT& F, dSymMatrixT& C) const;

	/** Green-Lagrangian strain.
	 * \param F the deformation gradient
	 * \param E return value */
	void Compute_E(const dMatrixT& F, dSymMatrixT& E) const;

	/** left stretch tensor.
	 * \note this version is being removed. Use the version
	 * which requires the deformation to be passed in
	 * \param b return value */
	void Compute_b(dSymMatrixT& b) const;

	/** right stretch tensor. 
	 * \note this version of is being removed. Use the version
	 * which requires the deformation to be passed in
	 * \param C return value */
	void Compute_C(dSymMatrixT& C) const;

	/** Green-Lagrangian strain. 
	 * \note this version is being removed. Use the version
	 * which requires the deformation to be passed in
	 * \param E return value */
	void Compute_E(dSymMatrixT& E) const;

	/** acoustical tensor.
	 * \param normal wave propagation direction
	 * \return acoustical tensor */
	virtual const dSymMatrixT& AcousticalTensor(const dArrayT& normal);

	/** finite strain element group.
	 * allows access to all const functions of the finite strain element
	 * class that are not currently supported with wrappers. \note this
	 * method is not guaranteed to be supported. If no FiniteStrainT is
	 * available, this function will return NULL.
	 * \return a const pointer to the supporting element group */
	const FiniteStrainT* FiniteStrain(void) const { return fFSMatSupport->FiniteStrain(); };

private:

	/* set inverse of thermal transformation - return true if active */
	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);

	/** compute acoustical tensor in 2D.
	 * \param CIJKL material tangent modulus
	 * \param SIJ 2nd Piola-Kirchhoff stress 
	 * \param FkK deformation gradient
	 * \param N wave propogation direction
	 * \param Q resulting acoustical tensor */
	void ComputeQ_2D(const dMatrixT& CIJKL, const dSymMatrixT& SIJ,
		const dMatrixT& FkK, const dArrayT& N, dSymMatrixT& Q) const;

	/** compute acoustical tensor in 3D.
	 * \param CIJKL material tangent modulus
	 * \param SIJ 2nd Piola-Kirchhoff stress 
	 * \param FkK deformation gradient
	 * \param N wave propogation direction
	 * \param Q resulting acoustical tensor */
	void ComputeQ_3D(const dMatrixT& CIJKL, const dSymMatrixT& SIJ,
		const dMatrixT& FkK, const dArrayT& N, dSymMatrixT& Q) const;

protected:

	/** support for finite strain materials */
	const FSMatSupportT* fFSMatSupport;

	/** \name return values */
	/*@{*/
	dSymMatrixT fStress;
	dMatrixT fModulus;
	/*@}*/

private:

	/** return value for FSSolidMatT::AcousticalTensor */
	dSymMatrixT fQ;  

	/** inverse of the multiplicative thermal deformation gradient */
	dMatrixT fF_therm_inv;

	/** inverse of the multiplicative thermal deformation gradient
	 * from the previous time increment. */
	dMatrixT fF_therm_inv_last;
	
	/** return value. Used as the return value of the mechanical part
	 * of the deformation gradient, if there are thermal strain. Otherwise,
	 * is unused. */
	dMatrixT fF_mechanical;
	
	/** true if temperature field found during FSSolidMatT::Initialize */
	bool fTemperatureField;
	dArrayT fTemperature;
	
	/** \name FSSolidMatT::c_ijkl work space */
	/*@{*/
	dMatrixT F_0_;
	dArrayT vec_;
	dSymMatrixT stress_;
	/*@}*/
};

/* finite strain materials support */
inline const FSMatSupportT& FSSolidMatT::FSMatSupport(void) const
{ 
#if __option(extended_errorcheck)
	if (!fFSMatSupport) 
		ExceptionT::GeneralFail("FSSolidMatT::FSMatSupport", "pointer not set");
#endif

	return *fFSMatSupport; 
}

} /* namespace Tahoe */

#endif /* _FD_STRUCT_MAT_T_H_ */
