/* $Id: SSMatSupportT.cpp,v 1.1.2.1 2002-10-28 06:49:15 paklein Exp $ */
#include "SSMatSupportT.h"
#include "SmallStrainT.h"

using namespace Tahoe;

/* constructor */
SSMatSupportT::SSMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fSmallStrain(NULL)
{


}
 
/* destructor */
SSMatSupportT::~SSMatSupportT(void)
{

}

/* set the element group pointer */
void SSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

	/* cast to small strain pointer */
	fSmallStrain = dynamic_cast<const SmallStrainT*>(p);
}

/* return a pointer the specified local array */
const LocalArrayT* SSMatSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	/* quick exit to inherited */
	if (!fSmallStrain) return MaterialSupportT::LocalArray(t);

	switch (t)
	{
		case LocalArrayT::kLastDisp:
			return &(fSmallStrain->LastDisplacements());
	
		case LocalArrayT::kVel:
			return &(fSmallStrain->Velocities());

		case LocalArrayT::kAcc:
			return &(fSmallStrain->Accelerations());

		default:
			/* inherited */
			return MaterialSupportT::LocalArray(t);
	}
}
