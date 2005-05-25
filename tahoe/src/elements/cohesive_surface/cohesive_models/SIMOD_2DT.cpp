/* $Id: SIMOD_2DT.cpp,v 1.2 2005-05-25 17:25:04 paklein Exp $ */
#include "SIMOD_2DT.h"

/* enabled */
#ifdef __SIMOD__

#include "ParameterContainerT.h"
#include <string.h>

/* SIMOD headers */
#include "simod_model_lists.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;

/* list of supported SIMOD models */
const simod_util_spc::model_IDnums SIMODSupported[] = {
	simod_util_spc::xuneedleman
};
const int kNumSIMODSupported = sizeof(SIMODSupported)/sizeof(*SIMODSupported);

/* constructor */
SIMOD_2DT::SIMOD_2DT(void): 
	SurfacePotentialT(knumDOF),
	fSIMOD(NULL),
	finternalVar(NULL)
{
	SetName("SIMOD_2D");
}

/* destructor */
SIMOD_2DT::~SIMOD_2DT(void)
{
	delete fSIMOD;
	delete finternalVar;
}

/* return the number of state variables needed by the model */
int SIMOD_2DT::NumStateVariables(void) const {
	return fSIMOD->num_internal_var();
}

/* surface potential */
double SIMOD_2DT::FractureEnergy(const ArrayT<double>& state)
{
#pragma unused(state)

	return 0.0; 
}

double SIMOD_2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

#pragma unused(state)
#pragma unused(jump_u)

	return 0.0;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& SIMOD_2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

#pragma unused(state)
#pragma unused(sigma)
#pragma unused(qIntegrate)

	/* set size for model with maximum no. of parameters */
//	finternalVar->toFloatAry(state.Pointer());

	/* evaluate the traction */
	double time = 0.0;
	bool active = fSIMOD->calcTractions2d(finternalVar, jump_u[0], jump_u[1], time,
		fTraction[0], fTraction[1]);	

	return fTraction;
}

/* potential stiffness */
const dMatrixT& SIMOD_2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

#pragma unused(jump_u)
#pragma unused(state)
#pragma unused(sigma)

	/* evaluate the stiffness */
	double time = 0.0;
	double tau_t = 0.0, tau_n = 0.0;
	SharedInterfaceModel::a2by2 Cchk = {{0.0,0.0},{0.0,0.0}};
	bool active= fSIMOD->calcC2d(finternalVar, jump_u[0], jump_u[1], time,
		tau_t, tau_n, Cchk);

	/* translate */
	fStiffness(0,0) = Cchk[0][0];
	fStiffness(1,0) = Cchk[1][0];
	fStiffness(0,1) = Cchk[0][1];
	fStiffness(1,1) = Cchk[1][1];
	
	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT SIMOD_2DT::Status(const dArrayT& jump_u, const ArrayT<double>& state)
{
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

#pragma unused(jump_u)

	bool active = fSIMOD->active(finternalVar);
	return (active) ? Precritical : Failed;
}

/* information about subordinate parameter lists */
void SIMOD_2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);

	/* 2D SIMOD model choice */
	sub_list.AddSub("simod_model_choice_2D", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SIMOD_2DT::NewSub(const StringT& name) const
{
	const char caller[] = "SIMOD_2DT::NewSub";
	if (name == "simod_model_choice_2D")
	{
		ParameterContainerT* simod_2D = new ParameterContainerT(name);
		simod_2D->SetListOrder(ParameterListT::Choice);
		simod_2D->SetSubSource(this);
	
		/* supported SIMOD models */
		for (int i = 0; i < kNumSIMODSupported; i++) {
			StringT model_name = "SIMOD_";
			model_name.Append(simod_util_spc::model_names[SIMODSupported[i]].c_str());
			simod_2D->AddSub(model_name);
		}
	
		return simod_2D;
	}
	else if (strncmp("SIMOD_", name, 6) == 0) /* catch all */
	{
		ParameterContainerT* model_params = new ParameterContainerT(name);

		/* create instance */
		string model_name = name.Pointer(6);
		SharedInterfaceModel_spc::SharedInterfaceModel* simod = 
			SharedInterfaceModel_spc::SharedInterfaceModel::factory(model_name);
		if (!simod) ExceptionT::GeneralFail(caller, "error constructing \"%s\"", name.Pointer());
			
		/* define vector of parameters */
		int nv = simod->num_input_parameters();
		for (int i = 0; i < nv; i++) {
			StringT label = "p_";
			label.Append(i+1);
			model_params->AddParameter(ParameterT::Double, label);
		}

		return model_params;
	}
	else /* inherited */
		return SurfacePotentialT::NewSub(name);
}

/* accept parameter list */
void SIMOD_2DT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SIMOD_2DT::TakeParameterList";

	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	/* resolve model choice */
	const ParameterListT* simod_params = list.ListChoice(*this, "simod_model_choice_2D");
	if (!simod_params) ExceptionT::GeneralFail(caller, "could not resolve choice \"simod_model_choice_2D\"");

	/* construct model */
	string model_name = simod_params->Name().Pointer(6);
	fSIMOD = SharedInterfaceModel_spc::SharedInterfaceModel::factory(model_name);
	if (!fSIMOD)
		ExceptionT::GeneralFail(caller, "error constructing \"%s\"", 
			simod_params->Name().Pointer());

	/* extract parameters */
	const ArrayT<ParameterT>& parameters = simod_params->Parameters();
	if (parameters.Length() != fSIMOD->num_input_parameters())
		ExceptionT::GeneralFail(caller, "expecting %d input parameters not %d",
			fSIMOD->num_input_parameters(), parameters.Length());
	dArrayT param_array(parameters.Length());
	for (int i = 0; i < param_array.Length(); i++)
		param_array[i] = parameters[i];

	/* construct parameter input object */
	if (model_name == simod_util_spc::model_names[simod_util_spc::xuneedleman])
	{
		/* define parameter object */
		XuNeedleman::parameters* params = new XuNeedleman::parameters(param_array.Pointer());
		
		/* define model parameters */
		fSIMOD->setParameters(params);
		delete params;
		
		/* construct internal variable array */
		finternalVar = new XuNeedleman::internalVar;
	}	
	else
		ExceptionT::GeneralFail(caller, "unrecognized model name \"%s\"",
			model_name.c_str());
}

#endif /* __SIMOD__ */
