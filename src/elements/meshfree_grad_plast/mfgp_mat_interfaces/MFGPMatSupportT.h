/* $Id: MFGPMatSupportT.h */
#ifndef _MFGP_MAT_SUPPORT_T_H_
#define _MFGP_MAT_SUPPORT_T_H_

/* base class */
#include "BasicSupportT.h"

/* direct members */
#include "LocalArrayT.h"
#include "AutoArrayT.h"
#include "ElementCardT.h"

#include "dArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MFGPAssemblyT;
class ElementCardT;

/** support for the MFGP materials classes. */
class MFGPMatSupportT: public BasicSupportT
{
public:

	/** constructor */
	MFGPMatSupportT(int ndof, int nip);
	MFGPMatSupportT(int ndof_displ, int ndof_plast, int nip_displ, int nip_plast);

	/** destructor */
	virtual ~MFGPMatSupportT(void);

	/** \name dimensions */
	/*@{*/
	/** number of degrees of freedom (per node) */
	int NumDOF(void) const { return fNumDOF_displ; };
	int NumDOF_PM(void) const { return fNumDOF_plast; };

	/** stress evaluation points per element */
	int NumIP(void) const { return fNumIP_displ; };
	/*@}*/

	/** the low-level communicator only including processes with non-zero numbers
	 * of elements, or NULL if it doesn't exist */
	const CommunicatorT* GroupCommunicator(void) const { return fGroupCommunicator; };

	/** \name run time status */
	/*@{*/

	/** current stress evaluation point within the element. If
	 * no source for the current point is set using 
	 * MaterialSupportT::SetCurrIP, will return 0. */
	int CurrIP(void) const;

	/** set the source for the current stress evaluation point */
	void SetCurrIP(const int& curr_ip);
	/*@}*/
	
	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const MFGPAssemblyT* MFGPAssembly(void) const;

	/** solver iteration number for the group set with MaterialSupportT::SetGroup */
	const int& GroupIterationNumber(void) const;

	/** return the number of elements. If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return 0 */
	int NumElements(void) const;

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
	virtual void SetMFGPAssembly(const MFGPAssemblyT* p);

	/** set pointer local array */
	virtual void SetLocalArray(const LocalArrayT& array);

	/** set the source for element cards */
	void SetElementCards(AutoArrayT<ElementCardT>* element_cards);

	/** set solver group */
	void SetGroup(int group) { fGroup = group; };
	/*@}*/
	
	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/
	
	/** \name laplacian of total strain */
	/*@{*/
	const dSymMatrixT& LapLinearStrain(void) const;
	const dSymMatrixT& LapLinearStrain(int ip) const;
	/*@}*/

	/** total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;
	/*@}*/
	
	/** laplacian of total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LapLinearStrain_last(void) const;
	const dSymMatrixT& LapLinearStrain_last(int ip) const;
	/*@}*/
	
	/*@}*/
	/** set source for the strain */
	void SetLinearStrain(const ArrayT<dSymMatrixT>* strain_List);

	/** set source for the strain from the end of the previous time step */
	void SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List);
	/*@}*/
	
	/*@}*/
	/** set source for the laplacian of strain */
	void SetLapLinearStrain(const ArrayT<dSymMatrixT>* lapstrain_List);

	/** set source for the laplacian of strain from the end of the previous time step */
	void SetLapLinearStrain_last(const ArrayT<dSymMatrixT>* lapstrain_last_List);
	/*@}*/
	
	/** \name plastic multiplier */
	/*@{*/
	const dArrayT& LambdaPM(void) const;
	const dArrayT& LambdaPM(int ip) const;
	/*@}*/
	
	/** \name laplacian of plastic multiplier */
	/*@{*/
	const dArrayT& LapLambdaPM(void) const;
	const dArrayT& LapLambdaPM(int ip) const;
	/*@}*/

	/** plastic multiplier from the end of the previous time step */
	/*@{*/
	const dArrayT& LambdaPM_last(void) const;
	const dArrayT& LambdaPM_last(int ip) const;
	/*@}*/
	
	/** laplacian of plastic multiplier from the end of the previous time step */
	/*@{*/
	const dArrayT& LapLambdaPM_last(void) const;
	const dArrayT& LapLambdaPM_last(int ip) const;
	/*@}*/
	
	/*@}*/
	/** set source for the plastic multiplier */
	void SetLambdaPM(const ArrayT<dArrayT>* lambda_List);

	/** set source for the plastic multiplier from the end of the previous time step */
	void SetLambdaPM_last(const ArrayT<dArrayT>* lambda_last_List);
	/*@}*/
	
	/*@}*/
	/** set source for the laplacian of plastic multiplier */
	void SetLapLambdaPM(const ArrayT<dArrayT>* laplambda_List);

	/** set source for the laplacian of plastic multiplier from the end of the previous time step */
	void SetLapLambdaPM_last(const ArrayT<dArrayT>* laplambda_last_List);
	/*@}*/

  private:
  
  	/** \name dimensions */
  	/*@{*/
	/** number of degrees of freedom */
	int fNumDOF_displ;
	int fNumDOF_plast;
	
	/** number of integration points */
	int fNumIP_displ;
	int fNumIP_plast;
  	/*@}*/
  	
  	/** source for the current integration point */
	const int* fCurrIP;

	/** communicator including only processes with non-zero numbers of elements */
	const CommunicatorT* fGroupCommunicator;

	/** pointer to element card information */
	AutoArrayT<ElementCardT>* fElementCards;	
  
  	/** pointer to the continuum element */
  	const MFGPAssemblyT* fMFGPAssembly;

	/** solver group for MaterialSupportT::fContinuumElement */
	int fGroup;

	/** \name pointers to local arrays */
	/*@{*/
	const LocalArrayT* fInitCoords;
	const LocalArrayT* fDisp;
	const LocalArrayT* fLastDisp;
	const LocalArrayT* fVel;
	const LocalArrayT* fAcc;
	/*@}*/
	
	/** \name return values */
	/*@{*/
	/* strains */
  	const ArrayT<dSymMatrixT>* fStrain_List;
  	const ArrayT<dSymMatrixT>* fStrain_last_List;
  	const ArrayT<dSymMatrixT>* fLapStrain_List;
  	const ArrayT<dSymMatrixT>* fLapStrain_last_List;
  	
  	/* plastic multiplier */
  	const ArrayT<dArrayT>* fLambda_List;
  	const ArrayT<dArrayT>* fLambda_last_List;
  	const ArrayT<dArrayT>* fLapLambda_List;
  	const ArrayT<dArrayT>* fLapLambda_last_List;
	/*@}*/
};

/* inlines functions */
inline const MFGPAssemblyT* MFGPMatSupportT::MFGPAssembly(void) const
{
	return fMFGPAssembly;
}

/* solver iteration number for the group set with MaterialSupportT::SetGroup */
inline const int& MFGPMatSupportT::GroupIterationNumber(void) const {
	if (fGroup == -1) ExceptionT::GeneralFail("MFGPMatSupportT::GroupIterationNumber", "solver group not set");
	return IterationNumber(fGroup); /* inherited */
}

/* set the source for element cards */
inline void MFGPMatSupportT::SetElementCards(AutoArrayT<ElementCardT>* element_cards)
{
	fElementCards = element_cards;
}

/* return the number of elements */
inline int MFGPMatSupportT::NumElements(void) const
{
	if (fElementCards) 
		return fElementCards->Length();
	else
		return 0;
}

/* return the current element */
inline int MFGPMatSupportT::CurrElementNumber(void) const
{
	if (fElementCards) 
		return fElementCards->Position();
	else
		return -1;
}

/* return the specified card */
inline ElementCardT* MFGPMatSupportT::ElementCard(int card) const
{
	if (fElementCards) 
		return fElementCards->Pointer(card);
	else
		return NULL;
}

/* return the current */
inline ElementCardT* MFGPMatSupportT::CurrentElement(void) const
{
	if (fElementCards && fElementCards->InRange()) 
		return &(fElementCards->Current());
	else
		return NULL;
}

inline int MFGPMatSupportT::CurrIP(void) const {
	if (fCurrIP) return *fCurrIP;
	else return 0;
}

inline void MFGPMatSupportT::SetCurrIP(const int& curr_ip) { fCurrIP = &curr_ip; }

/* inlines */
/* return strains */
inline const dSymMatrixT& MFGPMatSupportT::LinearStrain(void) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LinearStrain(int ip) const
{
	if (!fStrain_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_List)[ip]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LinearStrain_last(void) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_last_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LinearStrain_last(int ip) const
{
	if (!fStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fStrain_last_List)[ip]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LapLinearStrain(void) const
{
	if (!fLapStrain_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LapLinearStrain(int ip) const
{
	if (!fLapStrain_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_List)[ip]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LapLinearStrain_last(void) const
{
	if (!fLapStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_last_List)[CurrIP()]; 
}

inline const dSymMatrixT& MFGPMatSupportT::LapLinearStrain_last(int ip) const
{
	if (!fLapStrain_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapStrain_last_List)[ip]; 
}

/* return plastic multiplier */
inline const dArrayT& MFGPMatSupportT::LambdaPM(void) const
{
	if (!fLambda_List) throw ExceptionT::kGeneralFail;
	return (*fLambda_List)[CurrIP()]; 
}

inline const dArrayT& MFGPMatSupportT::LambdaPM(int ip) const
{
	if (!fLambda_List) throw ExceptionT::kGeneralFail;
	return (*fLambda_List)[ip]; 
}

inline const dArrayT& MFGPMatSupportT::LambdaPM_last(void) const
{
	if (!fLambda_last_List) throw ExceptionT::kGeneralFail;
	return (*fLambda_last_List)[CurrIP()]; 
}

inline const dArrayT& MFGPMatSupportT::LambdaPM_last(int ip) const
{
	if (!fLambda_last_List) throw ExceptionT::kGeneralFail;
	return (*fLambda_last_List)[ip]; 
}

inline const dArrayT& MFGPMatSupportT::LapLambdaPM(void) const
{
	if (!fLapLambda_List) throw ExceptionT::kGeneralFail;
	return (*fLapLambda_List)[CurrIP()]; 
}

inline const dArrayT& MFGPMatSupportT::LapLambdaPM(int ip) const
{
	if (!fLapLambda_List) throw ExceptionT::kGeneralFail;
	return (*fLapLambda_List)[ip]; 
}

inline const dArrayT& MFGPMatSupportT::LapLambdaPM_last(void) const
{
	if (!fLapLambda_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapLambda_last_List)[CurrIP()]; 
}

inline const dArrayT& MFGPMatSupportT::LapLambdaPM_last(int ip) const
{
	if (!fLapLambda_last_List) throw ExceptionT::kGeneralFail;
	return (*fLapLambda_last_List)[ip]; 
}

} /* namespace Tahoe */

#endif /* _MFGP_MAT_SUPPORT_T_H_ */
