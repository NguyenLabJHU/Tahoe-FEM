/* $Id: GradSSSolidMatT.cpp,v 1.10 2004-07-15 08:28:12 paklein Exp $ */ 
#include "GradSSSolidMatT.h"
#include <iostream.h>
#include "GradSSMatSupportT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatT::GradSSSolidMatT(ifstreamT& in, const GradSSMatSupportT& support):
	ParameterInterfaceT("gradient_small_strain_solid_material"),
//	SSSolidMatT(support),
	fGradSSMatSupport(&support),

	fNumDOF_Field(support.NumDOF_Field()),
	fNumDOF_Total(support.NumDOF() + fNumDOF_Field),

	fNumIP_Field(support.NumIP_Field())
{

}

GradSSSolidMatT::GradSSSolidMatT(void):
	ParameterInterfaceT("gradient_small_strain_solid_material"),
	fGradSSMatSupport(NULL)
{

}

/* destructor */
GradSSSolidMatT::~GradSSSolidMatT(void) { }

/* field */
const double& GradSSSolidMatT::Field(void) const
{
	return fGradSSMatSupport->LinearField(); 
}

const double& GradSSSolidMatT::Field(int ip) const
{
	return fGradSSMatSupport->LinearField(ip); 
}

/* field from the end of the previous time step */
const double& GradSSSolidMatT::Field_last(void) const
{
	return fGradSSMatSupport->LinearField_last(); 
}

const double& GradSSSolidMatT::Field_last(int ip) const
{
	return fGradSSMatSupport->LinearField_last(ip); 
}

/* gradient field */
const double& GradSSSolidMatT::GradField(void) const
{
	return fGradSSMatSupport->LinearGradField(); 
}

const double& GradSSSolidMatT::GradField(int ip) const
{
	return fGradSSMatSupport->LinearGradField(ip); 
}

/* gradient field from the end of the previous time step */
const double& GradSSSolidMatT::GradField_last(void) const
{
	return fGradSSMatSupport->LinearGradField_last(); 
}

const double& GradSSSolidMatT::GradField_last(int ip) const
{
	return fGradSSMatSupport->LinearGradField_last(ip); 
}

/* Laplacian field */
const double& GradSSSolidMatT::LapField(void) const
{
	return fGradSSMatSupport->LinearLapField(); 
}

const double& GradSSSolidMatT::LapField(int ip) const
{
	return fGradSSMatSupport->LinearLapField(ip); 
}

/* Laplacian field from the end of the previous time step */
const double& GradSSSolidMatT::LapField_last(void) const
{
	return fGradSSMatSupport->LinearLapField_last(); 
}

const double& GradSSSolidMatT::LapField_last(int ip) const
{
	return fGradSSMatSupport->LinearLapField_last(ip); 
}

/* apply pre-conditions at the current time step */
void GradSSSolidMatT::InitStep(void)
{
	/* inherited */
	SSSolidMatT::InitStep();
}
