/* $Id: FieldSupportT.cpp,v 1.1.2.1 2002-04-30 00:07:12 paklein Exp $ */
#include "FieldSupportT.h"
#include "FEManagerT.h"

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
