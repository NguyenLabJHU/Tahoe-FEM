/* $Id: SolidMatSupportT.cpp,v 1.3 2003-01-27 07:00:28 paklein Exp $ */
#include "SolidMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "ElasticT.h"
#endif

using namespace Tahoe;

/* constructor */
SolidMatSupportT::SolidMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fElastic(NULL),
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
	fElastic = dynamic_cast<const ElasticT*>(p);
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
