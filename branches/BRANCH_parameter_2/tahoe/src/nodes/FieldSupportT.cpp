/* $Id: FieldSupportT.cpp,v 1.4.16.2 2004-03-22 18:39:39 paklein Exp $ */
#include "FieldSupportT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"

using namespace Tahoe;

/* the model */
ModelManagerT& FieldSupportT::ModelManager(void) const
{
	/* check pointer */
	ModelManagerT* model = fFEManager.ModelManager();
	if (!model)
		ExceptionT::GeneralFail("FieldSupportT::ModelManager", "not defined");
	return *model;
}

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

/* return a pointer to the specified schedule */
const ScheduleT* FieldSupportT::Schedule(int num) const
{
	return fNodeManager.Schedule(num);
}
