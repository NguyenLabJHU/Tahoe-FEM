/* $Id: BEMelement.cpp,v 1.3 2002-07-02 19:55:13 cjkimme Exp $ */
/* created: AFLP (02/28/1998) */

#include "BEMelement.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"
#include "StringT.h"


using namespace Tahoe;

BEMelement::BEMelement(const ElementSupportT& support, const FieldT& field, 
	const StringT& infile):
	ElementBaseT(support, field),
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
void BEMelement::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
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
