/* $Id: MaterialSupportT.cpp,v 1.2.8.2 2002-10-28 06:49:15 paklein Exp $ */
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
