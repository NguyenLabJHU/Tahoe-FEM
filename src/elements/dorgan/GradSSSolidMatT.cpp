/* $Id: GradSSSolidMatT.cpp,v 1.11 2004-07-20 23:16:50 rdorgan Exp $ */ 
#include "GradSSSolidMatT.h"
#include "GradSSMatSupportT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<GradSSSolidMatT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<GradSSSolidMatT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
GradSSSolidMatT::GradSSSolidMatT(void):
	ParameterInterfaceT("grad_small_strain_solid_material"),
	fGradSSMatSupport(NULL)
{

}

/* set the material support or pass NULL to clear */
void GradSSSolidMatT::SetGradSSMatSupport(const GradSSMatSupportT* support)
{
	/* set inherited material support */
	SetMaterialSupport(support);

	fGradSSMatSupport = support;
}

/* field */
const double& GradSSSolidMatT::Lambda(void) const
{
	return fGradSSMatSupport->LinearPMultiplier(); 
}

const double& GradSSSolidMatT::Lambda(int ip) const
{
	return fGradSSMatSupport->LinearPMultiplier(ip); 
}

/* field from the end of the previous time step */
const double& GradSSSolidMatT::Lambda_last(void) const
{
	return fGradSSMatSupport->LinearPMultiplier_last(); 
}

const double& GradSSSolidMatT::Lambda_last(int ip) const
{
	return fGradSSMatSupport->LinearPMultiplier_last(ip); 
}

/* gradient field */
const double& GradSSSolidMatT::GradLambda(void) const
{
	return fGradSSMatSupport->LinearGradPMultiplier(); 
}

const double& GradSSSolidMatT::GradLambda(int ip) const
{
	return fGradSSMatSupport->LinearGradPMultiplier(ip); 
}

/* gradient field from the end of the previous time step */
const double& GradSSSolidMatT::GradLambda_last(void) const
{
	return fGradSSMatSupport->LinearGradPMultiplier_last(); 
}

const double& GradSSSolidMatT::GradLambda_last(int ip) const
{
	return fGradSSMatSupport->LinearGradPMultiplier_last(ip); 
}

/* Laplacian field */
const double& GradSSSolidMatT::LapLambda(void) const
{
	return fGradSSMatSupport->LinearLapPMultiplier(); 
}

const double& GradSSSolidMatT::LapLambda(int ip) const
{
	return fGradSSMatSupport->LinearLapPMultiplier(ip); 
}

/* Laplacian field from the end of the previous time step */
const double& GradSSSolidMatT::LapLambda_last(void) const
{
	return fGradSSMatSupport->LinearLapPMultiplier_last(); 
}

const double& GradSSSolidMatT::LapLambda_last(int ip) const
{
	return fGradSSMatSupport->LinearLapPMultiplier_last(ip); 
}

/* apply pre-conditions at the current time step */
void GradSSSolidMatT::InitStep(void)
{
	/* inherited */
	SSSolidMatT::InitStep();
}
