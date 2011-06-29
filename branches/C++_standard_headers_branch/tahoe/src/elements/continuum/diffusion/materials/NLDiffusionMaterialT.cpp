/* $Id: NLDiffusionMaterialT.cpp,v 1.5 2005-01-07 02:16:03 paklein Exp $ */
#include "NLDiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* variation functions */
#include "LinearT.h"

/* constructor */
NLDiffusionMaterialT::NLDiffusionMaterialT(void):
	ParameterInterfaceT("nonlinear_diffusion_material"),
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


/* conductivity */
const dMatrixT& NLDiffusionMaterialT::k_ij(void)
{
	double field = fDiffusionMatSupport->Field();
	fScaledConductivity.SetToScaled(fConductivityScaleFunction->Function(field), fConductivity);
	return fScaledConductivity;
}

/* change in conductivity with temperature */
const dMatrixT& NLDiffusionMaterialT::dk_ij(void)
{
	double field = fDiffusionMatSupport->Field();
	fdk_ij.SetToScaled(fConductivityScaleFunction->DFunction(field), fConductivity);
	return fdk_ij;
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

/*information about subordinate parameter lists */
void NLDiffusionMaterialT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	DiffusionMaterialT::DefineSubs(sub_list);
	
	sub_list.AddSub("conductivity_function");
	sub_list.AddSub("specificheat_function");
}

/* return the description of the given inline subordinate parameter list */
void NLDiffusionMaterialT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "NL_diff_mat_function_choice") {
		order = ParameterListT::Choice;

		sub_lists.AddSub("linear_function");
		sub_lists.AddSub("power_law");
		sub_lists.AddSub("cubic_spline");
	}
	else /* inherited */
		DiffusionMaterialT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NLDiffusionMaterialT::NewSub(const StringT& name) const
{
	C1FunctionT* function = C1FunctionT::New(name);
	if (function)
		return function;
	else if (name == "conductivity_function" || name == "specificheat_function") {
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetSubSource(this);
		choice->AddSub("NL_diff_mat_function_choice", ParameterListT::Once, true);
		return choice;
	}	
	else /* inherited */
		return DiffusionMaterialT::NewSub(name);
}

/* accept parameter list */
void NLDiffusionMaterialT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "NLDiffusionMaterialT::TakeParameterList";

	/* inherited */
	DiffusionMaterialT::TakeParameterList(list);

	/* construct temperature dependence functions */
	const ParameterListT& cond_function_choice = list.GetList("conductivity_function");
	const ParameterListT& cond_function = cond_function_choice.GetListChoice(*this, "NL_diff_mat_function_choice");
	fConductivityScaleFunction = C1FunctionT::New(cond_function.Name());
	if (!fConductivityScaleFunction) ExceptionT::GeneralFail(caller, "could not construct %s", cond_function_choice.Name().Pointer());
	fConductivityScaleFunction->TakeParameterList(cond_function);

	const ParameterListT& cp_function_choice = list.GetList("specificheat_function");
	const ParameterListT& cp_function = cp_function_choice.GetListChoice(*this, "NL_diff_mat_function_choice");
	fCpScaleFunction = C1FunctionT::New(cp_function.Name());
	if (!fCpScaleFunction) ExceptionT::GeneralFail(caller, "could not construct %s", cp_function_choice.Name().Pointer());
	fCpScaleFunction->TakeParameterList(cp_function);

	/* dimension work space */
	fScaledConductivity.Dimension(NumSD());
}
