/* $Id: FieldSupportT.cpp,v 1.3 2002-07-02 19:56:26 cjkimme Exp $ */
#include "FieldSupportT.h"
#include "FEManagerT.h"


using namespace Tahoe;

void FieldSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const
{
	/* wrapper */
	fFEManager.AssembleLHS(group, elMat, eqnos);
}

void FieldSupportT::AssembleRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const
{
	/* wrapper */
	fFEManager.AssembleRHS(group, elRes, eqnos);
}

ifstreamT& FieldSupportT::Input(void) const
{
	/* wrapper */
	return fFEManager.Input();
}

ofstreamT& FieldSupportT::Output(void) const
{
	/* wrapper */
	return fFEManager.Output();
}
