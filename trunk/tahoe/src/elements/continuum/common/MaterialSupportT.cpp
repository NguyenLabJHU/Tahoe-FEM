/* $Id: MaterialSupportT.cpp,v 1.10 2004-07-15 08:26:14 paklein Exp $ */
#include "MaterialSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "ContinuumElementT.h"
#include "ElementSupportT.h"
#endif

using namespace Tahoe;

/* constructor */
MaterialSupportT::MaterialSupportT(int ndof, int nip):
//	fNumSD(nsd),
	fNumDOF(ndof),
	fNumIP(nip),
	fCurrIP(NULL),

	/* multiprocessor information */
	fGroupCommunicator(NULL),
	fElementCards(NULL),
	fContinuumElement(NULL),
	fGroup(-1),
	fInitCoords(NULL),
	fDisp(NULL)
{ 

}
 
/* destructor */
MaterialSupportT::~MaterialSupportT(void) { }

void MaterialSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	fContinuumElement = p;
#ifdef CONTINUUM_ELEMENT
	if (fContinuumElement)
	{
		const ElementSupportT& element_support = fContinuumElement->ElementSupport();
		
		/* set FEManagerT */
		const FEManagerT& fe_man = element_support.FEManager();
		SetFEManager(&fe_man);
		
		fGroupCommunicator = &(fContinuumElement->GroupCommunicator());
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
const LocalArrayT* MaterialSupportT::LocalArray(LocalArrayT::TypeT t) const
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
void MaterialSupportT::SetLocalArray(const LocalArrayT& array)
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

/* interpolate the given field to the current integration point */
bool MaterialSupportT::Interpolate(const LocalArrayT& u, dArrayT& u_ip) const
{
#ifdef CONTINUUM_ELEMENT
	if (!fContinuumElement) 
	{
		u_ip = 0.0;
		return false;
	}
	else
	{
		fContinuumElement->IP_Interpolate(u, u_ip);
		return true;
	}
#else
#pragma unused(u)
	u_ip = 0.0;
	return false;
#endif
}

/* interpolate the given field to the given integration point */
bool MaterialSupportT::Interpolate(const LocalArrayT& u, dArrayT& u_ip, int ip) const
{
#ifdef CONTINUUM_ELEMENT
	if (!fContinuumElement) 
	{
		u_ip = 0.0;
		return false;
	}
	else
	{
		fContinuumElement->IP_Interpolate(u, u_ip, ip);
		return true;
	}
#else
#pragma unused(u)
#pragma unused(ip)
	u_ip = 0.0;
	return false;
#endif
}
