/* $Id: ContinuumMaterialT.cpp,v 1.2 2001-02-20 00:23:20 paklein Exp $ */
/* created: paklein (11/20/1996)                                          */

#include "ContinuumMaterialT.h"
#include "ContinuumElementT.h"
#include "ShapeFunctionT.h"
#include "ArrayT.h"
#include "StringT.h"

/* constructor */
ContinuumMaterialT::ContinuumMaterialT(const ContinuumElementT& element):
	fContinuumElement(element),
	fNumDOF(element.NumDOF()),
	fNumIP(element.NumIP()),
	fCurrIP(element.CurrIP())
{

}

/* destructor */
ContinuumMaterialT::~ContinuumMaterialT(void) { }

/* element card data */
int ContinuumMaterialT::NumElements(void) const
{
	return fContinuumElement.NumElements();
}

int ContinuumMaterialT::CurrElementNumber(void) const
{
	return fContinuumElement.CurrElementNumber();
}

ElementCardT& ContinuumMaterialT::ElementCard(int card) const
{
	return fContinuumElement.ElementCard(card);
}

ElementCardT& ContinuumMaterialT::CurrentElement(void) const
{
	return fContinuumElement.CurrentElement();
}

/* initialization */
void ContinuumMaterialT::Initialize(void)
{
/* do nothing */
}

/* storage initialization */
bool ContinuumMaterialT::NeedsPointInitialization(void) const { return false; }
void ContinuumMaterialT::PointInitialize(void) { /* nothing to do */ }

/* form of tangent matrix */
GlobalT::SystemTypeT ContinuumMaterialT::TangentType(void) const
{
	/* symmetric by default */
	return GlobalT::kSymmetric;
}

/* apply pre-conditions at the current time step */
void ContinuumMaterialT::InitStep(void) { }

/* update/reset internal variables */
void ContinuumMaterialT::UpdateHistory(void) { }
void ContinuumMaterialT::ResetHistory(void) { }

/* print parameters */
void ContinuumMaterialT::Print(ostream& out) const
{
	PrintName(out);
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ContinuumMaterialT::NumOutputVariables(void) const { return 0; }
void ContinuumMaterialT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Free();	
}
void ContinuumMaterialT::ComputeOutput(dArrayT& output)
{
#pragma unused(output)
}

/* returns true if two materials have compatible output variables */
bool ContinuumMaterialT::CompatibleOutput(const ContinuumMaterialT& m1, 
	const ContinuumMaterialT& m2)
{
	/* number of variables */
	if (m1.NumOutputVariables() != m2.NumOutputVariables())
		return false;
	/* labels */
	else
	{
		ArrayT<StringT> labels1, labels2;
		m1.OutputLabels(labels1);
		m2.OutputLabels(labels2);
		for (int i = 0; i < labels1.Length(); i++)
			if (labels1[i] != labels2[i])
				return false;

		/* compatible if execution false through */
		return true;
	}
}	

/***********************************************************************
* Protected
***********************************************************************/

void ContinuumMaterialT::PrintName(ostream& out) const
{
	out << " Material name:\n";
}
