/* $Id: MaterialSupportT.cpp,v 1.4 2002-11-15 02:46:33 paklein Exp $ */
#include "MaterialSupportT.h"
#include "ContinuumElementT.h"

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
	fContinuumElement(NULL)
{ 

}
 
/* destructor */
MaterialSupportT::~MaterialSupportT(void)
{

}

/* return a pointer to the specified LoadTime function */
const ScheduleT* MaterialSupportT::Schedule(int num) const
{
	if (fContinuumElement) 
		return fContinuumElement->Schedule(num);
	else
		return NULL;
}

/* return a pointer the specified local array */
const LocalArrayT* MaterialSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	/* no source */
	if (!fContinuumElement) return NULL;

	switch (t)
	{
		case LocalArrayT::kInitCoords:
			return &(fContinuumElement->InitialCoordinates());
	
		case LocalArrayT::kDisp:
			return &(fContinuumElement->Displacements());

		default:
			return NULL;
	}
}

/* interpolate the given field to the current integration point */
bool MaterialSupportT::Interpolate(const LocalArrayT& u, dArrayT& u_ip) const
{
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
}

/* interpolate the given field to the given integration point */
bool MaterialSupportT::Interpolate(const LocalArrayT& u, dArrayT& u_ip, int ip) const
{
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
}
