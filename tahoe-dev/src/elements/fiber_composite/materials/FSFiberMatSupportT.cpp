/* $Id: FSFiberMatSupportT.cpp,v 1.3 2011-09-16 21:00:28 thao Exp $ */
#include "FSFiberMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "UpLagFiberCompT.h"
#endif

using namespace Tahoe;

/* constructor */
FSFiberMatSupportT::FSFiberMatSupportT(int ndof, int nip):
	FSMatSupportT(ndof, nip),
	fFiber_list(NULL),
	fFiberDispersion_list(NULL)
{

}

/* set source for fiber vector list */
void FSFiberMatSupportT::SetFibers(const ArrayT<dArray2DT>* Fiber_list)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
	/* keep pointer */
	fFiber_list = Fiber_list;
}

/* set source for fiber vector list */
void FSFiberMatSupportT::SetFibersDispersion(const ArrayT<dArrayT>* FiberDispersion_list)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
	/* keep pointer */
	fFiberDispersion_list = FiberDispersion_list;
}

/* set the element group pointer */
void FSFiberMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	FSMatSupportT::SetContinuumElement(p);

/*
#ifdef CONTINUUM_ELEMENT
	fFiberElement = TB_DYNAMIC_CAST(const UpLagFiberCompT*, p);
#endif
*/
}
