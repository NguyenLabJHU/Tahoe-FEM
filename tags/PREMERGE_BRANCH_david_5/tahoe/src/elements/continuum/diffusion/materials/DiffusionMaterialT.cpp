/* $Id: DiffusionMaterialT.cpp,v 1.10 2005-01-07 02:16:03 paklein Exp $ */
/* created: paklein (10/02/1999) */
#include "DiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"

#include "StringT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<DiffusionMaterialT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<DiffusionMaterialT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
DiffusionMaterialT::DiffusionMaterialT(void):
	ParameterInterfaceT("linear_diffusion_material"),
	fDiffusionMatSupport(NULL),
	fDensity(0.0),
	fSpecificHeat(0.0)
{

}

/* set support */
void DiffusionMaterialT::SetDiffusionMatSupport(const DiffusionMatSupportT* support)
{
	/* inherited */
	SetMaterialSupport(support);
	fDiffusionMatSupport = support;

	/* dimension */
	int nsd = NumSD();
	fConductivity.Dimension(nsd);
	fq_i.Dimension(nsd);
	fdq_i.Dimension(nsd);
	fdk_ij.Dimension(nsd);

	/* initialize */
	fConductivity = 0.0;
	fq_i = 0.0;
	fdq_i = 0.0;
	fdk_ij = 0.0;
}

/* heat flux */
const dArrayT& DiffusionMaterialT::q_i(void)
{
	/* should be 1 row */
	fConductivity.Multx(fDiffusionMatSupport->Gradient(), fq_i, -1.0);
	return fq_i;
}

/* describe the parameters needed by the interface */
void DiffusionMaterialT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ContinuumMaterialT::DefineParameters(list);

	/* define parameters */
	list.AddParameter(fDensity, "density");
	list.AddParameter(fSpecificHeat, "specific_heat");
	list.AddParameter(ParameterT::Double, "conductivity");
}

/* accept parameter list */
void DiffusionMaterialT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContinuumMaterialT::TakeParameterList(list);

	/* get parameters */
	fDensity = list.GetParameter("density");
	fSpecificHeat = list.GetParameter("specific_heat");
	double k = list.GetParameter("conductivity");
	fConductivity.Identity(k);
}
