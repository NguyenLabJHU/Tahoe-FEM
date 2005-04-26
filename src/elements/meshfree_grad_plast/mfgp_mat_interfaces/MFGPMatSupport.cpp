/* $Id: MFGPMatSupportT.cpp */
#include "MFGPMatSupportT.h"
#include "ElementsConfig.h"
#include "MFGP_AssemblyT.h"
#include "ElementSupportT.h"

using namespace Tahoe;

/* constructor */
MFGPMatSupportT::MFGPMatSupportT(int ndof, int nip):
	fNumDOF(ndof),
	fNumIP(nip),
	fCurrIP(NULL),

	/* multiprocessor information */
	fGroupCommunicator(NULL),
	fElementCards(NULL),
	fMFGPAssembly(NULL),
	fGroup(-1),
	fInitCoords(NULL),
	fDisp(NULL)
{ 

}
 
/* destructor */
MFGPMatSupportT::~MFGPMatSupportT(void) { }

void MFGPMatSupportT::SetMFGPAssembly(const MFGP_AssemblyT* p)
{
	fMFGPAssembly = p;
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
}

/* return a pointer the specified local array */
const LocalArrayT* MFGPMatSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	switch (t)
	{
		case LocalArrayT::kInitCoords:
			return fInitCoords;
	
		case LocalArrayT::kDisp:
			return fDisp;

		default:
			return NULL;
	}
}

/* set pointer */
void MFGPMatSupportT::SetLocalArray(const LocalArrayT& array)
{
	switch (array.Type())
	{
		case LocalArrayT::kInitCoords:
			fInitCoords = &array;
			break;
		case LocalArrayT::kDisp:
			fDisp = &array;
			break;
		default:
			ExceptionT::GeneralFail("MaterialSupportT::LocalArray",
				"unrecognized array type: %d", array.Type());
	}
}

/* set source for the strain */
void MFGPMatSupportT::SetLinearStrain(const ArrayT<dSymMatrixT>* strain_List)
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
void MFGPMatSupportT::SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List)
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
	
/* set source for the laplacian of strain */
void MFGPMatSupportT::SetLapLinearStrain(const ArrayT<dSymMatrixT>* lapstrain_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!lapstrain_List) throw ExceptionT::kGeneralFail;
	if (lapstrain_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*lapstrain_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fLapStrain_List = lapstrain_List;
}

/** set source for the strain from the end of the previous time step */
void MFGPMatSupportT::SetLapLinearStrain_last(const ArrayT<dSymMatrixT>* lapstrain_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!lapstrain_last_List) throw ExceptionT::kGeneralFail;
	if (lapstrain_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*lapstrain_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fLapStrain_last_List = lapstrain_last_List;
}

/* set source for the plastic multiplier */
void MFGPMatSupportT::SetLambda(const ArrayT<dArrayT>* lambda_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!lambda_List) throw ExceptionT::kGeneralFail;
	if (lambda_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*lambda_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fLambda_List = lambda_List;
}

/** set source for the plastic multiplier from the end of the previous time step */
void MFGPMatSupportT::SetLambda_last(const ArrayT<dArrayT>* lambda_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!lambda_last_List) throw ExceptionT::kGeneralFail;
	if (lambda_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*lambda_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fLambda_last_List = lambda_last_List;
}
	
/* set source for the laplacian of plastic multiplier */
void MFGPMatSupportT::SetLapLambda(const ArrayT<dArrayT>* laplambda_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!laplambda_List) throw ExceptionT::kGeneralFail;
	if (laplambda_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*laplambda_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fLapLambda_List = laplambda_List;
}

/** set source for the plastic multiplier from the end of the previous time step */
void MFGPMatSupportT::SetLapLambda_last(const ArrayT<dArrayT>* laplambda_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!laplambda_last_List) throw ExceptionT::kGeneralFail;
	if (laplambda_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*laplambda_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fLapLambda_last_List = laplambda_last_List;
}
