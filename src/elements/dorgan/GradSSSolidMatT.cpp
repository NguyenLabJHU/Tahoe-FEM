/* $Id: GradSSSolidMatT.cpp,v 1.8 2004-04-23 18:44:36 rdorgan Exp $ */ 
#include "GradSSSolidMatT.h"
#include <iostream.h>
#include "GradSSMatSupportT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatT::GradSSSolidMatT(ifstreamT& in, const GradSSMatSupportT& support):
	SSSolidMatT(in, support),
	fGradSSMatSupport(&support),

	fNumDOF_Field(support.NumDOF_Field()),
	fNumDOF_Total(support.NumDOF() + fNumDOF_Field),

	fNumIP_Field(support.NumIP_Field())
{
	SetName("gradient_small_strain_solid_material");
}

GradSSSolidMatT::GradSSSolidMatT(void):
	fGradSSMatSupport(NULL)
{
	SetName("gradient_small_strain_solid_material");
}

/* destructor */
GradSSSolidMatT::~GradSSSolidMatT(void) { }

/* initialization */
void GradSSSolidMatT::Initialize(void)
{
	/* inherited */
	SSSolidMatT::Initialize();
}

/* I/O */
void GradSSSolidMatT::PrintName(ostream& out) const
{
	/* inherited */
	SSSolidMatT::PrintName(out);
	
	out << "    Multi Field solution\n";
}

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

/* Laplacian field */
const double& GradSSSolidMatT::LaplacianField(void) const
{
	return fGradSSMatSupport->LinearLaplacianField(); 
}

const double& GradSSSolidMatT::LaplacianField(int ip) const
{
	return fGradSSMatSupport->LinearLaplacianField(ip); 
}

/* Laplacian field from the end of the previous time step */
const double& GradSSSolidMatT::LaplacianField_last(void) const
{
	return fGradSSMatSupport->LinearLaplacianField_last(); 
}

const double& GradSSSolidMatT::LaplacianField_last(int ip) const
{
	return fGradSSMatSupport->LinearLaplacianField_last(ip); 
}

/* apply pre-conditions at the current time step */
void GradSSSolidMatT::InitStep(void)
{
	/* inherited */
	SSSolidMatT::InitStep();
}
