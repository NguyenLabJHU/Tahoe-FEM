/* $Id: FieldSupportT.cpp,v 1.4 2003-08-18 03:46:37 paklein Exp $ */
#include "FieldSupportT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"

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

/* construct new KBC controller */
KBC_ControllerT* FieldSupportT::NewKBC_Controller(FieldT& field, int code) const
{
	/* node manager */
	return fNodeManager.NewKBC_Controller(field, code);
}

FBC_ControllerT* FieldSupportT::NewFBC_Controller(FieldT& field, int code) const
{
	/* node manager */
	return fNodeManager.NewFBC_Controller(field, code);
}
