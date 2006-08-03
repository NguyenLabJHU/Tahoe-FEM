/* $Id: FSFiberMatSupportT.cpp,v 1.1 2006-08-03 01:10:41 thao Exp $ */
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
	fFiberElement(NULL)
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

/* set the element group pointer */
void FSFiberMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	SolidMatSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to finite strain pointer */
	fFiberElement = TB_DYNAMIC_CAST(const UpLagFiberCompT*, p);
#endif
}
