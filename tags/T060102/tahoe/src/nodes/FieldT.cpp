/* $Id: FieldT.cpp,v 1.1 2002-04-21 07:13:33 paklein Exp $ */
#include "FieldT.h"
#include "nControllerT.h"
#include "KBC_ControllerT.h"
#include "FBC_ControllerT.h"

/* constructor */
FieldT::FieldT(const StringT& name, nControllerT& controller, int ndof):
	fName(name),
	fnController(controller),
	fField(fnController.Order()),
	fField_last(fnController.Order()),
	fEqnos(0, ndof)
{
	/* set default labels */
	fLabels.Dimension(ndof);
	for (int i = 0; i < fLabels.Length(); i++)
		fLabels[i].Append("D_", i+1);
		
	/* set integrator */
	for (int i = 0; i < fField.Length(); i++)
		fnController.RegisterField(fField[i], i);
}

/* set field labels */
void FieldT::SetLabels(const ArrayT<StringT>& labels)
{
	fLabels = labels;
}

/* set number of nodes */
void FieldT::Dimension(int nnd)
{
	/* number of degrees of freedom */
	int ndof = fEqnos.MinorDim();

	/* dimension field */
	for (int i = 0; i < fField.Length(); i++)
	{
		fField[i].Dimension(nnd, ndof);
		fField_last[i].Dimension(fField[i]);
	}

	/* dimension equations array */
	fEqnos.Dimension(nnd, ndof);
}

/* beginning of time series */
void FieldT::InitialCondition(void)
{
	fField_last = fField;
}

/* apply predictor to all degrees of freedom */
void FieldT::InitStep(void)
{
	/* predictor */
	fnController.Predictor();
}

/* update history */
void FieldT::CloseStep(void)
{

}

/* update the active degrees of freedom */
void FieldT::Update(const dArrayT& update, int eq_start, int eq_stop)
{
	/* corrector */
	fnController.Corrector(fEqnos, update, eq_start, eq_stop);
}

/* reset displacements (and configuration to the last known solution) */
void FieldT::ResetStep(void)
{
	fField = fField_last;
}
