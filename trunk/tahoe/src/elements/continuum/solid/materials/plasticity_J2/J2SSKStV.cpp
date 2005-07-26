/* $Id: J2SSKStV.cpp,v 1.15 2005-07-26 16:31:24 paklein Exp $ */
/* created: paklein (06/18/1997) */
#include "J2SSKStV.h"
#include "SSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);

/* element output data */
const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {
	"alpha",  // equivalent plastic strain
	   "VM",  // Von Mises stress
	"press"}; // pressure

/* constructor */
J2SSKStV::J2SSKStV(void):
	ParameterInterfaceT("small_strain_StVenant_J2"),
	HookeanMatT(3),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	fElasticIterations(0)
{

}

/* update internal variables */
void J2SSKStV::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) 
		Update(element, NumIP());
}

/* reset internal variables to last converged solution */
void J2SSKStV::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) 
		Reset(element, NumIP());
}

/* modulus */
const dMatrixT& J2SSKStV::c_ijkl(void)
{
	/* elastoplastic correction */
	fModulus.SumOf(HookeanMatT::Modulus(), ModuliCorrection(CurrentElement(), Mu(), NumIP(), CurrIP()));	
	return fModulus;
}

/* stress */
const dSymMatrixT& J2SSKStV::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, NumIP(), ip);

	/* elastic stress */
	HookeanStress(e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	int iteration = fSSMatSupport->GroupIterationNumber();
	if (iteration > fElasticIterations) /* elastic iteration */
		fStress += StressCorrection(e_els, element, Mu(), NumIP(), ip);
	/*	if (CurrElementNumber() == 0) {
	  cout<<"\nIP: "<<CurrIP();
	  cout<<"\nStress: "<<fStress;
	  if (element.IsAllocated()) {
	    cout<< "\nalpha: "<<fInternal[kalpha];
	    cout << "\n: Kalpha: "<<K(fInternal[kalpha]);
	    cout << "\nYield: "<<YieldCondition(fStress, fInternal[kalpha]);
	  }
	  }*/
	return fStress;	
}

/* returns the strain energy density for the specified strain */
double J2SSKStV::StrainEnergyDensity(void)
{
	return HookeanEnergy(ElasticStrain(e(), CurrentElement(), NumIP(), CurrIP()));		
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int J2SSKStV::NumOutputVariables(void) const  { return kNumOutput; }
void J2SSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

const iArrayT& J2SSKStV::InternalDOF(void) const
{
  return(fInternalDOF);
}

const dArrayT& J2SSKStV::InternalStrainVars(void)
{
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) {
    s_ij();
    double* p = fInternalStrainVars.Pointer();

    /* plastic increment */
    double& dgamma = fInternal[kdgamma];

    *p++ = fInternal[kalpha] + sqrt23*dgamma;

    *p++ = fPlasticStrain[0] + dgamma*fUnitNorm[0];
    *p++ = fPlasticStrain[1] + dgamma*fUnitNorm[1];
    *p++ = fPlasticStrain[2] + dgamma*fUnitNorm[2];
    *p++ = fPlasticStrain[3] + dgamma*fUnitNorm[3];
    *p++ = fPlasticStrain[4] + dgamma*fUnitNorm[4];
    *p++ = fPlasticStrain[5] + dgamma*fUnitNorm[5];

    *p++ = fPlasticStrain[0] + dgamma*fUnitNorm[0];
    *p++ = fPlasticStrain[1] + dgamma*fUnitNorm[1];
    *p++ = fPlasticStrain[2] + dgamma*fUnitNorm[2];
    *p++ = fPlasticStrain[3] + dgamma*fUnitNorm[3];
    *p++ = fPlasticStrain[4] + dgamma*fUnitNorm[4];
    *p++ = fPlasticStrain[5] + dgamma*fUnitNorm[5];

    /*    const iArrayT& flags = element.IntegerData();
    if (flags[CurrIP()] == kIsPlastic && CurrElementNumber() == 0) {
      cout << "\nElement: "<<CurrElementNumber()<< "\t IP: "<<CurrIP();
      cout << "\n:Internal Strains: "<<fInternalStrainVars;
      }*/
  }
  else fInternalStrainVars = 0.0;
  return(fInternalStrainVars);
}

const dArrayT& J2SSKStV::InternalStressVars(void)
{
  ElementCardT& element = CurrentElement();
  if (element.IsAllocated()) {
    s_ij();
    double* p = fInternalStressVars.Pointer();
    double& dgamma = fInternal[kdgamma];
    double alpha = fInternal[kalpha];
    alpha += sqrt23*dgamma;

    *p++ = -(K(alpha)-fYield);

    *p++ = -fBeta[0]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
    *p++ = -fBeta[1]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
    *p++ = -fBeta[2]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
    *p++ = -fBeta[3]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
    *p++ = -fBeta[4]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);
    *p++ = -fBeta[5]; //-(2.0*(1.0 - ftheta)*fH_bar*dgamma/3.0, fUnitNorm);

    *p++ = fStress[0];
    *p++ = fStress[1];
    *p++ = fStress[2];
    *p++ = fStress[3];
    *p++ = fStress[4];
    *p++ = fStress[5];

    /*    const iArrayT& flags = element.IntegerData();
    if (flags[CurrIP()] == kIsPlastic && CurrElementNumber() == 0) {
      cout << "\nElement: "<<CurrElementNumber()<< "\t IP: "<<CurrIP();
      cout << "\n:Internal Stress: "<<fInternalStressVars;
      }*/
  }
  else fInternalStressVars = 0.0;
  return(fInternalStressVars);
}

void J2SSKStV::ComputeOutput(dArrayT& output)
{
	/* stress tensor (loads element data and sets fStress) */
	s_ij();

	/* pressure */
	output[2] = fStress.Trace()/3.0;
	
	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		/* plastic strain */
		output[0] = fInternal[kalpha];
		
		/* status flags */
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
			output[0] += sqrt23*fInternal[kdgamma];
	}
	else
		output[0] = 0.0;
}

/* describe the parameters needed by the interface */
void J2SSKStV::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSSolidMatT::DefineParameters(list);
	IsotropicT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
	J2SSC0HardeningT::DefineParameters(list);

	/* number of elastic iterations */
	ParameterT elastic_its(fElasticIterations, "elastic_iterations");
	elastic_its.SetDefault(fElasticIterations);
	list.AddParameter(elastic_its);
}

/* information about subordinate parameter lists */
void J2SSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	J2SSC0HardeningT::DefineSubs(sub_list);
}

/* return the description of the given inline subordinate parameter list */
void J2SSKStV::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	SSSolidMatT::DefineInlineSub(name, order, sub_lists);
	IsotropicT::DefineInlineSub(name, order, sub_lists);
	HookeanMatT::DefineInlineSub(name, order, sub_lists);
	J2SSC0HardeningT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* J2SSKStV::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = NULL;

	/* try each base class */
	sub = SSSolidMatT::NewSub(name);
	if (sub) return sub;

	sub = IsotropicT::NewSub(name);
	if (sub) return sub;

	sub = HookeanMatT::NewSub(name);
	if (sub) return sub;
	
	return J2SSC0HardeningT::NewSub(name);
}

/* accept parameter list */
void J2SSKStV::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
	J2SSC0HardeningT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);

	/* number elastic iterations */
	fElasticIterations = list.GetParameter("elastic_iterations");
	fElasticIterations -= 2;
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set modulus */
void J2SSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}
