/* $Id: MaterialSupportT.cpp,v 1.5 2003-01-27 07:00:28 paklein Exp $ */
#include "MaterialSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "ContinuumElementT.h"
#endif

using namespace Tahoe;

/* constructor */
MaterialSupportT::MaterialSupportT(int nsd, int ndof, int nip):
	fNumSD(nsd),
	fNumDOF(ndof),
	fNumIP(nip),
	fRunState(NULL),

	/* sources for run time information */
	fCurrIP(NULL),
	fIterationNumber(NULL),
	fTime(NULL),
	fTimeStep(NULL),
	fStepNumber(NULL),
	fNumberOfSteps(NULL),

	fElementCards(NULL),
	fContinuumElement(NULL),
	fInitCoords(NULL),
	fDisp(NULL)
{ 

}
 
/* destructor */
MaterialSupportT::~MaterialSupportT(void) { }

/* return a pointer to the specified LoadTime function */
const ScheduleT* MaterialSupportT::Schedule(int num) const
{
#ifdef CONTINUUM_ELEMENT
	if (fContinuumElement) 
		return fContinuumElement->Schedule(num);
	else
		return NULL;
#else
#pragma unused(num)
	return NULL;
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
