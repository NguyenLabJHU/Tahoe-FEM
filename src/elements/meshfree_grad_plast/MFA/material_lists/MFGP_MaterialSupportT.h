/* $Id: MFGP_MaterialSupportT.h,v 1.2 2005-01-14 00:12:25 raregue Exp $ */
#ifndef _MFGP_MATERIAL_SUPPORT_T_H_
#define _MFGP_MATERIAL_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "LocalArrayT.h"
#include "AutoArrayT.h"
#include "ElementCardT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MFGP_AssemblyT;
class ElementCardT;

/** support for the Tahoe materials classes. */
class MFGP_MaterialSupportT: public MaterialSupportT
{
public:

	/** constructor */
	MFGP_MaterialSupportT(int ndof, int nip);

	/** destructor */
	virtual ~MFGP_MaterialSupportT(void);

	/** the low-level communicator only including processes with non-zero numbers
	 * of elements, or NULL if it doesn't exist */
	const CommunicatorT* GroupCommunicator(void) const { return fGroupCommunicator; };
	
	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const MFGP_AssemblyT* MFGP_Assembly(void) const;

	/** solver iteration number for the group set with MaterialSupportT::SetGroup */
	const int& GroupIterationNumber(void) const;

	/** return the current element.  If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return -1 */
	int CurrElementNumber(void) const;

	/** return the specified card.  If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return NULL */
	//ElementCardT* ElementCard(int card);
	ElementCardT* ElementCard(int card) const;

	/** return the current card.  If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return NULL */
	ElementCardT* CurrentElement(void) const;

	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;
	/*@}*/

	/** \name set host code information */
	/*@{*/
	/** set the element group pointer */
	virtual void SetMFGPAssembly(const MFGP_AssemblyT* p);

	/** set pointer local array */
	virtual void SetLocalArray(const LocalArrayT& array);

	/** set the source for element cards */
	void SetElementCards(AutoArrayT<ElementCardT>* element_cards);

	/** set solver group */
	void SetGroup(int group) { fGroup = group; };
	/*@}*/
	//********************************************
	// SSMatSupportT
	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/

	/** total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;
	
	/** \name laplacian of total strain */
	/*@{*/
	const dSymMatrixT& LapLinearStrain(void) const;
	const dSymMatrixT& LapLinearStrain(int ip) const;
	/*@}*/

	/** laplacian of total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LapLinearStrain_last(void) const;
	const dSymMatrixT& LapLinearStrain_last(int ip) const;

	/** set source for the strain */
	void SetLinearStrain(const ArrayT<dSymMatrixT>* strain_List);

	/** set source for the strain from the end of the previous time step */
	void SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List);
	
	/** set source for the laplacian of strain */
	void SetLapLinearStrain(const ArrayT<dSymMatrixT>* strain_List);

	/** set source for the laplacian of strain from the end of the previous time step */
	void SetLapLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List);
	/*@}*/
	//********************************************
  private:
  
	/** communicator including only processes with non-zero numbers of elements */
	const CommunicatorT* fGroupCommunicator;

	/** pointer to element card information */
	AutoArrayT<ElementCardT>* fElementCards;	
  
  	/** pointer to the continuum element */
  	const MFGP_AssemblyT* fMFGPAssembly;

	/** solver group for MaterialSupportT::fContinuumElement */
	int fGroup;

	/** \name pointers to local arrays */
	/*@{*/
	//const LocalArrayT* fInitCoords;
	//const LocalArrayT* fDisp;
	// SolidMatSupportT
	const LocalArrayT* fLastDisp;
	const LocalArrayT* fVel;
	const LocalArrayT* fAcc;
	// SSMatSupportT
	const ArrayT<dSymMatrixT>* fStrain_List;
  	const ArrayT<dSymMatrixT>* fStrain_last_List;
  	const ArrayT<dSymMatrixT>* fLapStrain_List;
  	const ArrayT<dSymMatrixT>* fLapStrain_last_List;
	/*@}*/
};

/* inlines functions */
inline const MFGP_AssemblyT* MFGP_MaterialSupportT::MFGP_Assembly(void) const
{
	return fMFGPAssembly;
}

/* set the source for element cards */
inline void MFGP_MaterialSupportT::SetElementCards(AutoArrayT<ElementCardT>* element_cards)
{
	fElementCards = element_cards;
}
//*********************************
// SSMatSupportT
inline const dSymMatrixT& MFGP_MaterialSupportT::LinearStrain(void) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LinearStrain(int ip) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_List)[ip]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LinearStrain_last(void) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_last_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LinearStrain_last(int ip) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_last_List)[ip]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LapLinearStrain(void) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LapLinearStrain(int ip) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_List)[ip]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LapLinearStrain_last(void) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_last_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGP_MaterialSupportT::LapLinearStrain_last(int ip) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_last_List)[ip]; 
}
//*********************************
} /* namespace Tahoe */

#endif /* _MFGP_MATERIAL_SUPPORT_T_H_ */
