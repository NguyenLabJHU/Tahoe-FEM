/* $Id: FSMatSupportT.cpp,v 1.1.2.1 2002-10-28 06:49:15 paklein Exp $ */
#include "FSMatSupportT.h"
#include "FiniteStrainT.h"

using namespace Tahoe;

/* constructor */
FSMatSupportT::FSMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fFiniteStrain(NULL)
{


}
 
/* destructor */
FSMatSupportT::~FSMatSupportT(void)
{

}

/* set the element group pointer */
void FSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

	/* cast to finite strain pointer */
	fFiniteStrain = dynamic_cast<const FiniteStrainT*>(p);
}

/* return a pointer the specified local array */
const LocalArrayT* FSMatSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	/* quick exit to inherited */
	if (!fFiniteStrain) return MaterialSupportT::LocalArray(t);

	switch (t)
	{
		case LocalArrayT::kLastDisp:
			return &(fFiniteStrain->LastDisplacements());
	
		case LocalArrayT::kVel:
			return &(fFiniteStrain->Velocities());

		case LocalArrayT::kAcc:
			return &(fFiniteStrain->Accelerations());

		default:
			/* inherited */
			return MaterialSupportT::LocalArray(t);
	}
}
