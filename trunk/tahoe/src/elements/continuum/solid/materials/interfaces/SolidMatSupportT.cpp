/* $Id: SolidMatSupportT.cpp,v 1.2 2002-11-14 17:06:21 paklein Exp $ */
#include "SolidMatSupportT.h"
#include "ElasticT.h"

using namespace Tahoe;

/* constructor */
SolidMatSupportT::SolidMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fElastic(NULL)
{

}

/* set the element group pointer */
void SolidMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

	/* cast to elastic element pointer */
	fElastic = dynamic_cast<const ElasticT*>(p);
}

/* return a pointer the specified local array */
const LocalArrayT* SolidMatSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	/* quick exit to inherited */
	if (!fElastic) return MaterialSupportT::LocalArray(t);

	switch (t)
	{
		case LocalArrayT::kLastDisp:
			return &(fElastic->LastDisplacements());
	
		case LocalArrayT::kVel:
			return &(fElastic->Velocities());

		case LocalArrayT::kAcc:
			return &(fElastic->Accelerations());

		default:
			/* inherited */
			return MaterialSupportT::LocalArray(t);
	}
}

/* nodal temperatures */
const LocalArrayT* SolidMatSupportT::Temperatures(void) const
{
	if (!fElastic)
		return NULL;
	else
		return fElastic->Temperatures();
}

/* nodal temperatures from the last time step */
const LocalArrayT* SolidMatSupportT::LastTemperatures(void) const
{
	if (!fElastic)
		return NULL;
	else
		return fElastic->LastTemperatures();
}
