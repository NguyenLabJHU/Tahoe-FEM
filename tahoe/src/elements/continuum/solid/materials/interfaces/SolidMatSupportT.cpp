/* $Id: SolidMatSupportT.cpp,v 1.4 2003-01-29 07:34:57 paklein Exp $ */
#include "SolidMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "SolidElementT.h"
#endif

using namespace Tahoe;

/* constructor */
SolidMatSupportT::SolidMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fSolidElement(NULL),
	fLastDisp(NULL),
	fVel(NULL),
	fAcc(NULL),
	fTemperatures(NULL),
	fLastTemperatures(NULL)
{

}

/* set the element group pointer */
void SolidMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to elastic element pointer */
	fSolidElement = dynamic_cast<const SolidElementT*>(p);
#endif
}

/* return a pointer the specified local array */
const LocalArrayT* SolidMatSupportT::LocalArray(LocalArrayT::TypeT t) const
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

void SolidMatSupportT::SetLocalArray(const LocalArrayT& array)
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
