/* $Id: MFGP_MaterialSupportT.cpp */

#include "MFGP_MaterialSupportT.h"
#include "ElementsConfig.h"

#include "MFGP_AssemblyT.h"
#include "ElementSupportT.h"

using namespace Tahoe;

/* constructor */
MFGP_MaterialSupportT::MFGP_MaterialSupportT(int ndof, int nip):
	MaterialSupportT(ndof, nip),
	/* multiprocessor information */
	fGroupCommunicator(NULL),
	fElementCards(NULL),
	fMFGPAssembly(NULL),
	fGroup(-1),
	fLastDisp(NULL),
	fVel(NULL),
	fAcc(NULL),
	fStrain_List(NULL),
	fStrain_last_List(NULL),
	fLapStrain_List(NULL),
	fLapStrain_last_List(NULL)
{ 

}
 
/* destructor */
MFGP_MaterialSupportT::~MFGP_MaterialSupportT(void) { }

void MFGP_MaterialSupportT::SetMFGPAssembly(const MFGP_AssemblyT* p)
{
	fMFGPAssembly = p;
#ifdef MESHFREE_GRAD_PLAST_DEV
	if (fMFGPAssembly)
	{
		const ElementSupportT& element_support = fMFGPAssembly->ElementSupport();
		
		/* set FEManagerT */
		const FEManagerT& fe_man = element_support.FEManager();
		SetFEManager(&fe_man);
		
		fGroupCommunicator = &(fMFGPAssembly->GroupCommunicator());
	}
	else 
	{
		SetFEManager(NULL);
		fGroupCommunicator = NULL;
		fGroup = -1;
	}
#endif
}

/* return a pointer the specified local array */
const LocalArrayT* MFGP_MaterialSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	switch (t)
	{
		case LocalArrayT::kLastDisp:
			return fLastDisp;
	
		case LocalArrayT::kVel:
			return fVel;

		case LocalArrayT::kAcc:
			return fAcc;

		default:
			/* inherited */
			return MaterialSupportT::LocalArray(t);
	}
}

/* set pointer */
void MFGP_MaterialSupportT::SetLocalArray(const LocalArrayT& array)
{
	switch (array.Type())
	{
		case LocalArrayT::kLastDisp:
			fLastDisp = &array;
			break;
	
		case LocalArrayT::kVel:
			fVel = &array;
			break;

		case LocalArrayT::kAcc:
			fAcc = &array;
			break;

		default:
			/* inherited */
			MaterialSupportT::SetLocalArray(array);
	}
}
//***************************************
//SSMatSupportT
/* set source for the strain */
void MFGP_MaterialSupportT::SetLinearStrain(const ArrayT<dSymMatrixT>* strain_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!strain_List) throw ExceptionT::kGeneralFail;
	if (strain_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fStrain_List = strain_List;
}

/** set source for the strain from the end of the previous time step */
void MFGP_MaterialSupportT::SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!strain_last_List) throw ExceptionT::kGeneralFail;
	if (strain_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fStrain_last_List = strain_last_List;
}
void MFGP_MaterialSupportT::SetLapLinearStrain(const ArrayT<dSymMatrixT>* strain_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!strain_List) throw ExceptionT::kGeneralFail;
	if (strain_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fStrain_List = strain_List;
}

/** set source for the strain from the end of the previous time step */
void MFGP_MaterialSupportT::SetLapLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!strain_last_List) throw ExceptionT::kGeneralFail;
	if (strain_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fStrain_last_List = strain_last_List;
}
//***************************************