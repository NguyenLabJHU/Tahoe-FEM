/* $Id: SSKStV_Optimize.cpp,v 1.1 2009-04-23 03:03:51 thao Exp $ */
/* created: paklein (06/10/1997) */
#include "SSKStV_Optimize.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;
static const int num_params = 2;

/* constructor */
SSKStV_Optimize::SSKStV_Optimize(void):
	ParameterInterfaceT("ss_optimize_StVenant")
{
}


void SSKStV_Optimize::DefineSubs(SubListT& sub_list) const
{ 
	/*inherited*/
	SSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
}

ParameterInterfaceT* SSKStV_Optimize::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = SSSolidMatT::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "modulus_definition_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetSubSource(this);

		/* set the choices */		
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT E_and_nu("E_and_nu");
		ParameterT E(ParameterT::Double, "Young_modulus");
		E.AddLimit(0.0, LimitT::Lower);
		E_and_nu.AddParameter(E);
		ParameterT Poisson(ParameterT::Double, "Poisson_ratio");
		Poisson.AddLimit(-1.0, LimitT::Lower);
		Poisson.AddLimit( 0.5, LimitT::Upper);
		E_and_nu.AddParameter(Poisson);
		choice->AddSub(E_and_nu);

		return(choice);
	}
	else /* inherited */
		return SSSolidMatT::NewSub(name);
}

void SSKStV_Optimize::DefineParameters(ParameterListT& list) const
{
	/*inherited*/
	SSSolidMatT::DefineParameters(list);
	list.SetDescription("parameter order: [E, Nu]");
}


/* accept parameter list */
void SSKStV_Optimize::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SSKStV_Optimize::TakeParameterList";
	/*inherited*/
	SSKStV::TakeParameterList(list);
	
	flabels.Dimension(num_params);
	fparams.Dimension(num_params);
	
	flabels[0] = "E";
	flabels[1] = "nu";
	
	fparams[0] = fE = IsotropicT::Young();
	fparams[1] =  fnu = IsotropicT::Poisson();

	fParamGrads.Dimension(dSymMatrixT::NumValues(NumSD()), num_params);
	fconstraint_grad.Dimension(num_params);
	fconstraint_grad = 0.0;
}


const dArray2DT& SSKStV_Optimize::ds_ij_dlambda_q(void)
{
	/*calculate delta strain*/
	const dSymMatrixT& strain = e();
	double trace = strain.Trace();
	
	double coeff1 = 1.0/(1.0+fnu);
	double coeff2 = fnu/(1.0 - fnu - 2.0*fnu*fnu);

	fParamGrads(0,0) = coeff1*strain[0] + coeff2*trace;
	fParamGrads(1,0) = coeff1*strain[1] + coeff2*trace;
	fParamGrads(2,0) = coeff1*strain[2] + coeff2*trace;
	
	fParamGrads(3,0) = coeff1*strain[3];
	fParamGrads(4,0) = coeff1*strain[4];
	fParamGrads(5,0) = coeff1*strain[5];
	
	/*derivative with respect to Poisson's Ratio*/
	coeff1 = -fE/(1.0+fnu)/(1.0+fnu);
	coeff2 = fE*(1+2.0*fnu*fnu)/((-1.0+fnu+2.0*fnu*fnu)*(-1.0+fnu+2.0*fnu*fnu));

	fParamGrads(0,1) = coeff1*strain[0] + coeff2*trace;
	fParamGrads(1,1) = coeff1*strain[1] + coeff2*trace;
	fParamGrads(2,1) = coeff1*strain[2] + coeff2*trace;
	
	fParamGrads(3,1) = coeff1*strain[3];
	fParamGrads(4,1) = coeff1*strain[4];
	fParamGrads(5,1) = coeff1*strain[5];
	
//	cout << "stress_grad: "<<fParamGrads<<endl;
	return(fParamGrads);
}
