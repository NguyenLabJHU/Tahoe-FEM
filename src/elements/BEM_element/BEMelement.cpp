/* $Id: BEMelement.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: AFLP (02/28/1998)                                             */

#include "BEMelement.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "StringT.h"

BEMelement::BEMelement(FEManagerT& fe_manager, const StringT& infile):
	ElementBaseT(fe_manager),
	fInfile(infile)
{
	/* reset format */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);

	/* read in all parameters */

	//need to define
	// fNumElements -> just 1
	// fNumDOF
	// fNumSD
	// fNumElemNodes -> (not 3!!!!! )
}

/* destructor */
BEMelement::~BEMelement(void) {	}

/* form of tangent matrix */
GlobalT::SystemTypeT BEMelement::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* solution calls */
void BEMelement::AddNodalForce(int node, dArrayT& force)
{
#pragma unused(node)
#pragma unused(force)
}

/* returns the energy as defined by the derived class types */
double BEMelement::InternalEnergy(void)
{
	return 0.0;
}


/* output */
void BEMelement::RegisterOutput(void)
{

}

void BEMelement::WriteOutput(IOBaseT::OutputModeT mode)
{
#pragma unused(mode)
}

void BEMelement::SendOutput(int kincode)
{
#pragma unused(kincode)
}

/***********************************************************************
* Protected
***********************************************************************/

/* called by FormRHS and FormLHS */
void BEMelement::LHSDriver(void)
{

}

void BEMelement::RHSDriver(void)
{

}

/* print element group data */
void BEMelement::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	// output your control data
}
