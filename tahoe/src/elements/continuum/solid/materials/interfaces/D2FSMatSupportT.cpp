/* $Id: D2FSMatSupportT.cpp,v 1.1.2.1 2002-10-28 06:49:15 paklein Exp $ */
#include "D2FSMatSupportT.h"

using namespace Tahoe;

/* constructor */
D2FSMatSupportT::D2FSMatSupportT(int nsd, int ndof, int nip):
	FSMatSupportT(nsd, ndof, nip),
	fD2MeshFreeFDElastic(NULL)
{


}
 
/* destructor */
D2FSMatSupportT::~D2FSMatSupportT(void)
{

}
