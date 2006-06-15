/* $Id: SSEnhLocMatSupportT.h,v 1.2 2006-06-15 18:07:17 regueiro Exp $ */
#ifndef _SS_ENH_LOC_MAT_SUPPORT_T_H_
#define _SS_ENH_LOC_MAT_SUPPORT_T_H_

/* base class */
#include "SSMatSupportT.h"

#include "Array2DT.h"

#include "DevelopmentElementsConfig.h"

namespace Tahoe {

/* forward declarations */
#ifdef ENHANCED_STRAIN_LOC_DEV
class SmallStrainEnhLocT;
#endif

/** support for the small strain embedded discontinuity Tahoe materials classes */
class SSEnhLocMatSupportT: public SSMatSupportT
{
public:

	/** constructor */
	SSEnhLocMatSupportT(int ndof, int nip);

	/** destructor */
	~SSEnhLocMatSupportT(void);
	
	/** set source for the stress */
	void SetElementStress(const ArrayT<dSymMatrixT>* stress_List);
	void SetElementStress(const Array2DT<dSymMatrixT>* elementstress_List);
	
	/** set source for the loc_flag */
	void SetElementLocFlag(const int* loc_flag);
	void SetElementLocFlag(const iArrayT* elementloc_flag);

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	#ifdef ENHANCED_STRAIN_LOC_DEV	
	const SmallStrainEnhLocT* SmallStrainEnhLoc(void) const { return fSmallStrainEnhLoc; };
	#endif

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/
	
	// functions to access data in SmallStrainEnhLocT
	const dSymMatrixT& ElementStress(void) const;
	const dSymMatrixT& ElementStress(int ip) const;
	const dSymMatrixT& ElementStress(int elem, int ip) const;
	const int& ElementLocflag(void) const;
	const int& ElementLocflag(int elem) const;

  private:

  	/** \name return values */
	/*@{*/
  	const ArrayT<dSymMatrixT>* fStress_List;
  	
  	const Array2DT<dSymMatrixT>* fElementStress_List;
  	
  	const int* fLocFlag;
  	
  	const iArrayT* fElemLocFlag;
	/*@}*/
	
	/** pointer to the small strain embedded strong discontinuity element */
	#ifdef ENHANCED_STRAIN_LOC_DEV
	const SmallStrainEnhLocT* fSmallStrainEnhLoc;
	#endif

};

/* inlines */
inline const dSymMatrixT& SSEnhLocMatSupportT::ElementStress(void) const
{
	if (!fStress_List) throw ExceptionT::kGeneralFail;
	return (*fStress_List)[CurrIP()];
}

inline const dSymMatrixT& SSEnhLocMatSupportT::ElementStress(int ip) const
{
	if (!fStress_List) throw ExceptionT::kGeneralFail;
	return (*fStress_List)[ip];
}

inline const dSymMatrixT& SSEnhLocMatSupportT::ElementStress(int elem, int ip) const
{
	if (!fElementStress_List) throw ExceptionT::kGeneralFail;
	return (*fElementStress_List)[elem,ip];
}

inline const int& SSEnhLocMatSupportT::ElementLocflag(void) const
{
	if (!fLocFlag) throw ExceptionT::kGeneralFail;
	return *fLocFlag;
}

inline const int& SSEnhLocMatSupportT::ElementLocflag(int elem) const
{
	if (!fElemLocFlag) throw ExceptionT::kGeneralFail;
	return (*fElemLocFlag)[elem];
}



} /* namespace Tahoe */
#endif /* _SS_ENH_LOC_MAT_SUPPORT_T_H_ */
