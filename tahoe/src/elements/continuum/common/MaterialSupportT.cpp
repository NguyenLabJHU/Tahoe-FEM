/* $Id: MaterialSupportT.cpp,v 1.2.8.1 2002-10-26 16:24:21 paklein Exp $ */
#include "MaterialSupportT.h"

using namespace Tahoe;

/* constructor */
MaterialSupportT::MaterialSupportT(int nsd, int ndof, int nip):
	fNumSD(nsd),
	fNumDOF(ndof),
	fNumIP(nip),
	fCurrIP(NULL),
	fIterationNumber(NULL),
	fContinuumElement(NULL) 
{ 

}
 
/* destructor */
MaterialSupportT::~MaterialSupportT(void)
{

}
