/* $Id: SSMatSupportT.h,v 1.4 2004-07-15 08:28:22 paklein Exp $ */
#ifndef _SS_MAT_SUPPORT_T_H_
#define _SS_MAT_SUPPORT_T_H_

/* base class */
#include "SolidMatSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class SmallStrainT;

/** support for the small strain Tahoe materials classes */
class SSMatSupportT: public SolidMatSupportT
{
public:

	/** constructor */
	SSMatSupportT(int ndof, int nip);

	/** destructor */
	~SSMatSupportT(void);

	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/

	/** total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;

	/** set source for the strain */
	void SetLinearStrain(const ArrayT<dSymMatrixT>* strain_List);

	/** set source for the strain from the end of the previous time step */
	void SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List);
	/*@}*/

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const SmallStrainT* SmallStrain(void) const { return fSmallStrain; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/

  private:

  	/** \name return values */
	/*@{*/
  	const ArrayT<dSymMatrixT>* fStrain_List;
  	const ArrayT<dSymMatrixT>* fStrain_last_List;
	/*@}*/

  	/** pointer to the small strain element */
	const SmallStrainT* fSmallStrain;	
};

/* inlines */
inline const dSymMatrixT& SSMatSupportT::LinearStrain(void) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_List)[CurrIP()]; 
}

inline const dSymMatrixT& SSMatSupportT::LinearStrain(int ip) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_List)[ip]; 
}

inline const dSymMatrixT& SSMatSupportT::LinearStrain_last(void) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_last_List)[CurrIP()]; 
}

inline const dSymMatrixT& SSMatSupportT::LinearStrain_last(int ip) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_last_List)[ip]; 
}

} /* namespace Tahoe */
#endif /* _SS_MAT_SUPPORT_T_H_ */
