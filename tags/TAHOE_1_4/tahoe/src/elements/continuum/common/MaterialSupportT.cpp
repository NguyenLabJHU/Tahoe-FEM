/* $Id: MaterialSupportT.cpp,v 1.9 2004-04-08 23:37:16 paklein Exp $ */
#include "MaterialSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "ContinuumElementT.h"
#include "ElementSupportT.h"
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

	/* multiprocessor information */
	fSize(1),
	fRank(0),
	fCommunicator(NULL),
	fGroupCommunicator(NULL),

	fElementCards(NULL),
	fContinuumElement(NULL),
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
		fSize = element_support.Size();
		fRank = element_support.Rank();
		fCommunicator = &(element_support.Communicator());
		fGroupCommunicator = &(fContinuumElement->GroupCommunicator());
	}
	else {
		fCommunicator = NULL;
		fGroupCommunicator = NULL;
	}
#endif
}

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

/* the parameters stream */
ifstreamT& MaterialSupportT::Input(void) const
{
	const char caller[] = "MaterialSupportT::Input";

#ifdef CONTINUUM_ELEMENT
	if (!fContinuumElement) ExceptionT::GeneralFail(caller, "continuum element not defined");
	return fContinuumElement->ElementSupport().Input();
#else
	ExceptionT::GeneralFail(caller, "requires option CONTINUUM_ELEMENT");
	ifstreamT* dummy;
	return &dummy;
#endif
}

/* the echo file */
ofstreamT& MaterialSupportT::Output(void) const
{
	const char caller[] = "MaterialSupportT::Output";

#ifdef CONTINUUM_ELEMENT
	if (!fContinuumElement) ExceptionT::GeneralFail(caller, "continuum element not defined");
	return fContinuumElement->ElementSupport().Output();
#else
	ExceptionT::GeneralFail(caller, "requires option CONTINUUM_ELEMENT");
	ofstreamT* dummy;
	return &dummy;
#endif
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
