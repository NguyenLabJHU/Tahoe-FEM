/* $Id: NLDiffusionMaterialT.cpp,v 1.3.2.1 2004-01-21 19:09:57 paklein Exp $ */
#include "NLDiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* variation functions */
#include "LinearT.h"

/* constructor */
NLDiffusionMaterialT::NLDiffusionMaterialT(ifstreamT& in, const DiffusionMatSupportT& support):
	ParameterInterfaceT("nonlinear_diffusion"),
	DiffusionMaterialT(in, support),
	fConductivityScaleFunction(NULL),
	fCpScaleFunction(NULL),
	fScaledConductivity(NumSD())
{
	/* parameters in temperature variation in conductivity */
	double A, B;
	in >> A >> B;
	fConductivityScaleFunction = new LinearT(A, B);

	/* parameters in temperature variation in specific heat */
	in >> A >> B;
	fCpScaleFunction = new LinearT(A, B);
}

NLDiffusionMaterialT::NLDiffusionMaterialT(void):
	ParameterInterfaceT("nonlinear_diffusion"),
	fConductivityScaleFunction(NULL),
	fCpScaleFunction(NULL)
{

}

/* destructor */
NLDiffusionMaterialT::~NLDiffusionMaterialT(void)
{
	delete fConductivityScaleFunction;
	delete fCpScaleFunction;
}

/* I/O functions */
void NLDiffusionMaterialT::Print(ostream& out) const
{
	/* inherited */
	DiffusionMaterialT::Print(out);
	
	/* temperature variation functions */
	out << " Temperature variation in conductivity:\n";
	fConductivityScaleFunction->Print(out);
	fConductivityScaleFunction->PrintName(out);
}

void NLDiffusionMaterialT::PrintName(ostream& out) const
{
	/* inherited */
	DiffusionMaterialT::PrintName(out);
	
	out << "    Nonlinear diffusion material\n";
}

/* conductivity */
const dMatrixT& NLDiffusionMaterialT::k_ij(void)
{
	double field = fDiffusionMatSupport->Field();
	fScaledConductivity.SetToScaled(fConductivityScaleFunction->Function(field), fConductivity);
	return fScaledConductivity;
}

/* heat flux */
const dArrayT& NLDiffusionMaterialT::q_i(void)
{
	double scale = -fConductivityScaleFunction->Function(fDiffusionMatSupport->Field());
	fConductivity.Multx(fDiffusionMatSupport->Gradient(), fq_i, scale);
	return fq_i;
}

/* change in heat flux with temperature */
const dArrayT& NLDiffusionMaterialT::dq_i_dT(void)
{
	double scale = -fConductivityScaleFunction->DFunction(fDiffusionMatSupport->Field());
	fConductivity.Multx(fDiffusionMatSupport->Gradient(), fdq_i, scale);
	return fdq_i;
}

/* specific heat */
double NLDiffusionMaterialT::SpecificHeat(void) const
{
	double cp = DiffusionMaterialT::SpecificHeat();
	double scale = fCpScaleFunction->Function(fDiffusionMatSupport->Field());
	return cp*scale;
}

/* change in specific heat with temperature */
double NLDiffusionMaterialT::dCapacity_dT(void) const
{
	double d_cp = fCpScaleFunction->DFunction(fDiffusionMatSupport->Field());
	return fDensity*d_cp;
}
