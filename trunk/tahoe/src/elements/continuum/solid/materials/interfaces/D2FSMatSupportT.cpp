/* $Id: D2FSMatSupportT.cpp,v 1.2 2002-11-14 17:06:21 paklein Exp $ */
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
