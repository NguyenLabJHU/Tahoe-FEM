/* $Id: NLDiffusionMaterialT.cpp,v 1.3.18.2 2004-06-09 23:17:26 paklein Exp $ */
#include "NLDiffusionMaterialT.h"
#include "DiffusionMatSupportT.h"
#include "ifstreamT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* variation functions */
#include "LinearT.h"

/* constructor */
NLDiffusionMaterialT::NLDiffusionMaterialT(ifstreamT& in, const DiffusionMatSupportT& support):
	ParameterInterfaceT("nonlinear_diffusion_material"),
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
void NLDiffusionMaterialT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "NL_diff_mat_function_choice") {
		order = ParameterListT::Choice;

		sub_sub_list.AddSub("linear_function");
		sub_sub_list.AddSub("power_law");
		sub_sub_list.AddSub("cubic_spline");
	}
	else /* inherited */
		DiffusionMaterialT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NLDiffusionMaterialT::NewSub(const StringT& list_name) const
{
	C1FunctionT* function = C1FunctionT::New(list_name);
	if (function)
		return function;
	else if (list_name == "conductivity_function" || list_name == "specificheat_function") {
		ParameterContainerT* choice = new ParameterContainerT(list_name);
		choice->SetSubSource(this);
		choice->AddSub("NL_diff_mat_function_choice", ParameterListT::Once, true);
		return choice;
	}	
	else /* inherited */
		return DiffusionMaterialT::NewSub(list_name);
}

/* accept parameter list */
void NLDiffusionMaterialT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "NLDiffusionMaterialT::TakeParameterList";

	/* inherited */
	DiffusionMaterialT::TakeParameterList(list);

	/* construct temperature dependence functions */
	const ParameterListT* cond_function = NULL;
	const ParameterListT& cond_function_choice = list.GetList("conductivity_function");
	cond_function = cond_function_choice.ResolveListChoice(*this, "NL_diff_mat_function_choice");
	fConductivityScaleFunction = C1FunctionT::New(cond_function->Name());
	if (!fConductivityScaleFunction)
		ExceptionT::GeneralFail(caller, "could not construct %s", cond_function_choice.Name().Pointer());
	fConductivityScaleFunction->TakeParameterList(*cond_function);

	const ParameterListT& cp_function_choice = list.GetList("specificheat_function");
	cond_function = cp_function_choice.ResolveListChoice(*this, "NL_diff_mat_function_choice");
	fCpScaleFunction = C1FunctionT::New(cond_function->Name());
	if (!fCpScaleFunction)
		ExceptionT::GeneralFail(caller, "could not construct %s", cp_function_choice.Name().Pointer());
	fCpScaleFunction->TakeParameterList(*cond_function);

	/* dimension work space */
	fScaledConductivity.Dimension(NumSD());
}
