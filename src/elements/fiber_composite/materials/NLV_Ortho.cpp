/* $Id: NLV_Ortho.cpp,v 1.1 2008-09-25 13:32:32 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "NLV_Ortho.h"
#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"

#include "ParabolaT.h"
#include "FungType.h"
#include "FungType2.h"
#include "ScaledCsch.h"
#include "LinearExponentialT.h"

const double Pi = acos(-1.0);
static const double third = 1.0/3.0;
const int kNumOutputVar = 10;
static const char* Labels[kNumOutputVar] = {"p1_X", "p1_Y", "p1_Z","p2_X", "p2_Y", "p2_Z","Iv1_m", "Iv4_m", "Iv6_m","J"};

using namespace Tahoe;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
NLV_Ortho::NLV_Ortho(void):
  FSFiberMatViscT(),
  ParameterInterfaceT("nonlinear_viscoelasticity_ortho")
{
	/*reset default*/
	fNumFibProcess = 0;
	fNumMatProcess = 0;
}

/* destructor */
NLV_Ortho::~NLV_Ortho(void) 
{ 
	/*allocated?*/
	double l1 = fPot1_f4.Length();
	double l2 = fPot1_f5.Length();
	double l3 = fVisc1_f.Length();
	for (int i = 0; i < l1; i++){
			delete fPot1_f4[i];
			delete fPot2_f4[i];
	}
	for (int i = 0; i < l2; i++){
			delete fPot1_f5[i];
			delete fPot2_f5[i];
	}
	for (int i = 0; i < l3; i++){
			delete fVisc1_f[i];
			delete fVisc2_f[i];
	}
}


/* free energy density */
double NLV_Ortho::StrainEnergyDensity(void)
{
	/*returns free-energy density*/
	/*matrix contribution*/
	/* stretch */
	const dMatrixT& F = F_mechanical();
	fb.MultAAT(F);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	
	/*equilibrium contribution*/
	double J = sqrt(fC.Det());
	double J23 = pow(J, -2.0*third);
	fEigs *= J23;
	
	double energy =fPot_m[0]->Energy(fEigs, J);
	
	/*nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	for (int i = 0; i < fNumMatProcess; i++)
	{		
		const dSymMatrixT& Cv = fC_v[i];
		fInverse.Inverse(Cv);
		/*be = F.Cv^-1.F^T*/
		fb.MultQBQT(F, fInverse);
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
		fEigs = fSpectralDecompSpat.Eigenvalues();
		double Je = sqrt(fEigs.Product());

		/*deviatoric part*/
		double Je23 = pow(Je,-2.0*third);
		fEigs *= Je23;
		
		energy += fPot_m[i+1]->Energy(fEigs,J);
	}
		
	/*fiber contribution*/
	/*equilibrium contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	double I4 = fFiberStretch[0]; /*Cf11*/
	double I6 = fFiberStretch[1]; /*Cf22*/
	double I5 = fFiberStretch[0]*fFiberStretch[0] + fFiberStretch[4]*fFiberStretch[4] 
				+ fFiberStretch[5]*fFiberStretch[5]; /*I5=Cf11^2+Cf13^2+Cf12^2*/
	double I7 = fFiberStretch[1]*fFiberStretch[1] + fFiberStretch[3]*fFiberStretch[3] 
				+ fFiberStretch[5]*fFiberStretch[5]; 	/*I7=Cf22^2+Cf23^2+Cf12^2*/
/*	if (I4 < 1.0) I4 = 1.0;	
	if (I6 < 1.0) I6 = 1.0;
	if (I5 < 1.0) I5 = 1.0;
	if (I7 < 1.0) I7 = 1.0;
*/	
	energy += fPot1_f4[0]->Function(I4);
	energy += fPot2_f4[0]->Function(I6);
	energy += fPot1_f5[0]->Function(I5);
	energy += fPot2_f5[0]->Function(I7);

	/*non-equilibrium*/
	int j = fNumMatProcess;
	for (int i = 0; i < fNumFibProcess; i++)
	{
		ComputeFiberStretch(fC_v[j], fFiberStretch_v);

		double I4 = fFiberStretch[0]/fFiberStretch_v[0]; /*Cf:M1/Cfv:M1*/
		double I6 = fFiberStretch[1]/fFiberStretch_v[1]; /*Cf:M2/Cfv:M2*/
		
		fInverse.Inverse(fFiberStretch_v);
		double I5 = fFiberStretch[0]*(fFiberStretch[0]*fInverse[0] + fFiberStretch[4]*fInverse[4] + fFiberStretch[5]*fInverse[5])
				+ fFiberStretch[4]*(fFiberStretch[0]*fInverse[4] + fFiberStretch[4]*fInverse[2] + fFiberStretch[5]*fInverse[3]) 
				+ fFiberStretch[5]*(fFiberStretch[0]*fInverse[5] + fFiberStretch[4]*fInverse[3] + fFiberStretch[5]*fInverse[1]);
		I5 /= fFiberStretch_v[0]; /* (C.Cv^-1.C):M1/(Cv:M1) */

		double I7 = fFiberStretch[1]*(fFiberStretch[1]*fInverse[1] + fFiberStretch[3]*fInverse[3] + fFiberStretch[5]*fInverse[5])
				+ fFiberStretch[3]*(fFiberStretch[1]*fInverse[3] + fFiberStretch[3]*fInverse[2] + fFiberStretch[5]*fInverse[4]) 
				+ fFiberStretch[5]*(fFiberStretch[1]*fInverse[5] + fFiberStretch[3]*fInverse[4] + fFiberStretch[5]*fInverse[0]);
		I7 /= fFiberStretch_v[1]; /* (C.Cv^-1.C):M2/(Cv:M2)*/

		energy += fPot1_f4[i+1]->Function(I4);
		energy += fPot2_f4[i+1]->Function(I6);
		energy += fPot1_f5[i+1]->Function(I5);
		energy += fPot2_f5[i+1]->Function(I7);

		j++;
	}
	return energy;
}

int NLV_Ortho::NumOutputVariables() const {
	return kNumOutputVar;
}

void NLV_Ortho::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void NLV_Ortho::ComputeOutput(dArrayT& output)
{
	/*calculates deformed fiber vectors*/
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	
	const double* p_nt = Fibers(0);
	const double* p_is = Fibers(1);
	
	const dMatrixT& F = F_mechanical();
	double* pb = output.Pointer();
	
	/*deformed P1 fiber orientation*/
//	cout << "\nF: "<<F;
//	cout << "\nP1: "<<*p_nt<<"\t"<<*(p_nt+1);
	F.Multx(p_nt, pb);
//	cout << "\np1: "<<*pb<<"\t"<<*(pb+1);
//	cout << "\nJ: "<<F.Det();
	pb += NumSD();
	
	/*deformed P2 fiber orientation*/
	F.Multx(p_is, pb);
	pb += NumSD();
	
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	/*calculate neq. matrix contribution*/
	if (fNumMatProcess > 0)
	{
		const dSymMatrixT& Cv_m = fC_v[0];
		fSpectralDecompSpat.SpectralDecomp_Jacobi(Cv_m, false);	
		const dArrayT& eigenvals = fSpectralDecompSpat.Eigenvalues();
		*pb++ = eigenvals.Max();
	}
	else *pb++ = 1.0;
	if (fNumFibProcess > 0)
	{
		ComputeFiberStretch(fC_v[fNumMatProcess], fFiberStretch_v);
		*pb++ = fFiberStretch_v[0];
		*pb++ = fFiberStretch_v[1];
	}
	else {
		*pb++ = 1.0;
		*pb++ = 1.0;
	}
	output[9] = F.Det();
}

/* information about subordinate parameter lists */
void NLV_Ortho::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSFiberMatViscT::DefineSubs(sub_list);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("eq_fiber1_pot4", ParameterListT::Once);
	sub_list.AddSub("neq_fiber1_pot4", ParameterListT::Any);

	sub_list.AddSub("eq_fiber2_pot4", ParameterListT::Any);
	sub_list.AddSub("neq_fiber2_pot4", ParameterListT::Any);

	sub_list.AddSub("eq_fiber1_pot5", ParameterListT::Once);
	sub_list.AddSub("neq_fiber1_pot5", ParameterListT::Any);

	sub_list.AddSub("eq_fiber2_pot5", ParameterListT::Any);
	sub_list.AddSub("neq_fiber2_pot5", ParameterListT::Any);

	/* choice of viscosity */
	sub_list.AddSub("fiber1_visc_potential", ParameterListT::Any);
	sub_list.AddSub("fiber2_visc_potential", ParameterListT::Any);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NLV_Ortho::NewSub(const StringT& name) const
{

	C1FunctionT* func = NULL;
	if (name == "scaled-csch") 
		func = new ScaledCsch;
	else if (name == "linear_exponential") 
		func = new LinearExponentialT;
	else if (name == "fung-type")
		func = new FungType;
	else if (name == "parabola")
		func = new ParabolaT;
	else if (name == "fung_type_square")
		func = new FungType2;

	if (func)
		return func;

	if (name == "eq_fiber1_pot4" || name == "neq_fiber1_pot4" || name == "eq_fiber1_pot5" || name == "neq_fiber1_pot5" ||
		name == "eq_fiber2_pot4" || name == "neq_fiber2_pot4" || name == "eq_fiber2_pot5" || name == "neq_fiber2_pot5")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		choice->AddSub("parabola");
		choice->AddSub("fung_type");
		choice->AddSub("fung_type_square");

		return(choice);
	}
	else if (name == "fiber1_visc_potential"||name == "fiber2_visc_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		choice->AddSub("scaled-csch");
		choice->AddSub("linear_exponential");

		return(choice);
	}
	else
		return(FSFiberMatViscT::NewSub(name));
}

/* accept parameter list */
void NLV_Ortho::TakeParameterList(const ParameterListT& list)
{
	StringT caller = "TwoFiberViscT::TakeParameterList";
	
	/* inherited */
	FSFiberMatViscT::TakeParameterList(list);
	//	fNumFibStress = dSymMatrixT::NumValues(fNumSD);

	int num_fib1_neq4 =  list.NumLists("neq_fiber1_pot4");
	int num_fib1_neq5 =  list.NumLists("neq_fiber1_pot5");
	int num_fib1_visc =  list.NumLists("fiber1_visc_potential");
	
	int num_fib2_neq4 =  list.NumLists("neq_fiber2_pot4");
	int num_fib2_neq5 =  list.NumLists("neq_fiber2_pot5");
	int num_fib2_visc =  list.NumLists("fiber2_visc_potential");
	
	bool samefibers = false;
	
	if ((num_fib1_neq4 != num_fib1_neq5) &&(num_fib1_neq4 != num_fib1_visc))
		ExceptionT::GeneralFail(caller, 
			"number of viscosity functions does not match number of  nonequilibrium potentials for fiber 1");
	if ((num_fib2_neq4 != num_fib2_neq5) &&(num_fib2_neq4 != num_fib2_visc))
		ExceptionT::GeneralFail(caller, 
			"number of  viscosity functions does not match number of  nonequilibrium potentials for fiber 2");

	int num_fib2_eq4 =  list.NumLists("eq_fiber2_pot4");
	int num_fib2_eq5 =  list.NumLists("eq_fiber2_pot5");
	
	if (num_fib2_eq4 != num_fib2_eq5)
		ExceptionT::GeneralFail(caller, 
			"unequal number of potentials for fiber 2");
	
	if (num_fib2_eq4 ==0)
		samefibers = true;
	else 
		if (num_fib2_neq4 != num_fib2_neq4)
			ExceptionT::GeneralFail(caller, 
				"number of  viscosity functions and  nonequilibrium potentials do not match for fibers 1 and 2");
	
	fNumFibProcess = num_fib1_neq4;
	fVisc1_f.Dimension(fNumFibProcess);
	fPot1_f4.Dimension(fNumFibProcess+1);	
	fPot1_f5.Dimension(fNumFibProcess+1);	
		
	fVisc2_f.Dimension(fNumFibProcess);
	fPot2_f4.Dimension(fNumFibProcess+1);	
	fPot2_f5.Dimension(fNumFibProcess+1);	

	const ParameterListT& fiber1_pot4 = list.GetListChoice(*this, "eq_fiber1_pot4");
	if (fiber1_pot4.Name() == "parabola")
		fPot1_f4[0] = new ParabolaT;
	else if (fiber1_pot4.Name() == "fung-type")
		fPot1_f4[0] = new FungType;
	else if (fiber1_pot4.Name() == "fung_type_Square")
		fPot1_f4[0] = new FungType2;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot1_f4[0]) throw ExceptionT::kOutOfMemory;
	fPot1_f4[0]->TakeParameterList(fiber1_pot4);

	if (samefibers)
	{
		if (fiber1_pot4.Name() == "parabola")
			fPot2_f4[0] = new ParabolaT;
		else if (fiber1_pot4.Name() == "fung-type")
			fPot2_f4[0] = new FungType;
		else if (fiber1_pot4.Name() == "fung_type_Square")
			fPot2_f4[0] = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot2_f4[0]) throw ExceptionT::kOutOfMemory;
		fPot1_f4[0]->TakeParameterList(fiber1_pot4);
	}
	else
	{
		const ParameterListT& fiber2_pot4 = list.GetListChoice(*this, "eq_fiber2_pot4");
		if (fiber2_pot4.Name() == "parabola")
			fPot2_f4[0] = new ParabolaT;
		else if (fiber2_pot4.Name() == "fung-type")
			fPot2_f4[0] = new FungType;
		else if (fiber2_pot4.Name() == "fung_type_Square")
			fPot2_f4[0] = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot2_f4[0]) throw ExceptionT::kOutOfMemory;
		fPot2_f4[0]->TakeParameterList(fiber2_pot4);
	}
	
	const ParameterListT& fiber1_pot5 = list.GetListChoice(*this, "eq_fiber1_pot5");
	if (fiber1_pot5.Name() == "parabola")
		fPot1_f5[0] = new ParabolaT;
	else if (fiber1_pot5.Name() == "fung-type")
		fPot1_f5[0] = new FungType;
	else if (fiber1_pot5.Name() == "fung_type_square")
		fPot1_f5[0] = new FungType2;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot1_f5[0]) throw ExceptionT::kOutOfMemory;
	fPot1_f5[0]->TakeParameterList(fiber1_pot5);
	if (samefibers)
	{
		if (fiber1_pot5.Name() == "parabola")
			fPot2_f5[0] = new ParabolaT;
		else if (fiber1_pot5.Name() == "fung-type")
			fPot2_f5[0] = new FungType;
		else if (fiber1_pot5.Name() == "fung_type_Square")
			fPot2_f5[0] = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot2_f5[0]) throw ExceptionT::kOutOfMemory;
		fPot2_f5[0]->TakeParameterList(fiber1_pot5);
	}
	else
	{
		const ParameterListT& fiber2_pot5 = list.GetListChoice(*this, "eq_fiber2_pot5");
		if (fiber2_pot5.Name() == "parabola")
			fPot2_f5[0] = new ParabolaT;
		else if (fiber2_pot5.Name() == "fung-type")
			fPot2_f5[0] = new FungType;
		else if (fiber2_pot5.Name() == "fung_type_Square")
			fPot2_f5[0] = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot2_f5[0]) throw ExceptionT::kOutOfMemory;
		fPot2_f5[0]->TakeParameterList(fiber2_pot5);
	}

	for (int i = 0; i < fNumFibProcess; i++)
	{

		const ParameterListT& fiber1_potn4 = list.GetListChoice(*this, "neq_fiber1_pot4");
		if (fiber1_potn4.Name() == "parabola")
			fPot1_f4[i+1] = new ParabolaT;
		else if (fiber1_potn4.Name() == "fung-type")
			fPot1_f4[i+1] = new FungType;
		else if	(fiber1_potn4.Name() == "fung_type_square")
			fPot1_f4[i+1] = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot1_f4[i+1]) throw ExceptionT::kOutOfMemory;
		fPot1_f4[i+1]->TakeParameterList(fiber1_potn4);

		if (samefibers)
		{
			if (fiber1_potn4.Name() == "parabola")
				fPot2_f4[i+1] = new ParabolaT;
			else if (fiber1_potn4.Name() == "fung-type")
				fPot2_f4[i+1] = new FungType;
			else if (fiber1_potn4.Name() == "fung_type_Square")
				fPot2_f4[i+1] = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fPot2_f4[i+1]) throw ExceptionT::kOutOfMemory;
			fPot1_f4[i+1]->TakeParameterList(fiber1_potn4);
		}
		else
		{
			const ParameterListT& fiber2_potn4 = list.GetListChoice(*this, "neq_fiber2_pot4");
			if (fiber2_potn4.Name() == "parabola")
				fPot2_f4[i+1] = new ParabolaT;
			else if (fiber2_potn4.Name() == "fung-type")
				fPot2_f4[i+1] = new FungType;
			else if (fiber2_potn4.Name() == "fung_type_Square")
				fPot2_f4[i+1] = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fPot2_f4[i+1]) throw ExceptionT::kOutOfMemory;
			fPot2_f4[i+1]->TakeParameterList(fiber2_potn4);
		}


		const ParameterListT& fiber1_potn5 = list.GetListChoice(*this, "neq_fiber1_pot5");
		if (fiber1_potn5.Name() == "parabola")
			fPot1_f5[i+1] = new ParabolaT;
		else if (fiber1_potn5.Name() == "fung-type")
			fPot1_f5[i+1] = new FungType;
		else if	(fiber1_potn5.Name() == "fung_type_square")
			fPot1_f5[i+1] = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot1_f5[i+1]) throw ExceptionT::kOutOfMemory;
		fPot1_f5[i+1]->TakeParameterList(fiber1_potn5);

		if (samefibers)
		{
			if (fiber1_potn5.Name() == "parabola")
				fPot2_f5[i+1] = new ParabolaT;
			else if (fiber1_potn5.Name() == "fung-type")
				fPot2_f5[i+1] = new FungType;
			else if (fiber1_potn5.Name() == "fung_type_Square")
				fPot2_f5[i+1] = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fPot2_f5[i+1]) throw ExceptionT::kOutOfMemory;
			fPot2_f5[i+1]->TakeParameterList(fiber1_potn5);
		}
		else
		{
			const ParameterListT& fiber2_potn5 = list.GetListChoice(*this, "neq_fiber2_pot5");
			if (fiber2_potn5.Name() == "parabola")
				fPot2_f5[i+1] = new ParabolaT;
			else if (fiber2_potn5.Name() == "fung-type")
				fPot2_f5[i+1] = new FungType;
			else if (fiber2_potn5.Name() == "fung_type_Square")
				fPot2_f5[i+1] = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fPot2_f5[i+1]) throw ExceptionT::kOutOfMemory;
			fPot2_f5[i+1]->TakeParameterList(fiber2_potn5);
		}

		const ParameterListT& fiber1_visc = list.GetListChoice(*this, "fiber1_visc_potential", i);
		if (fiber1_visc.Name() == "linear_exponential")
			fVisc1_f[i] = new LinearExponentialT;
		else if (fiber1_visc.Name() == "scaled-csch")
			fVisc1_f[i] = new ScaledCsch;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fVisc1_f[i]) throw ExceptionT::kOutOfMemory;
		fVisc1_f[i]->TakeParameterList(fiber1_visc);
		
		if (samefibers)
		{
			if (fiber1_visc.Name() == "parabola")
				fVisc2_f[i] = new ParabolaT;
			else if (fiber1_visc.Name() == "fung-type")
				fVisc2_f[i] = new FungType;
			else if (fiber1_visc.Name() == "fung_type_Square")
				fVisc2_f[i] = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fVisc2_f[i]) throw ExceptionT::kOutOfMemory;
			fVisc2_f[i]->TakeParameterList(fiber1_visc);
		}
		else
		{
			const ParameterListT& fiber2_visc = list.GetListChoice(*this, "fiber2_visc_potential", i);
			if (fiber2_visc.Name() == "parabola")
				fVisc2_f[i] = new ParabolaT;
			else if (fiber2_visc.Name() == "fung-type")
				fVisc2_f[i] = new FungType;
			else if (fiber2_visc.Name() == "fung_type_Square")
				fVisc2_f[i] = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fVisc2_f[i]) throw ExceptionT::kOutOfMemory;
			fVisc2_f[i]->TakeParameterList(fiber2_visc);
		}
	}
			
	/*dimension state variable storage arrays*/
	int numprocess = fNumFibProcess + fNumMatProcess;
	SetStateVariables(numprocess);

	/*Dimension work spaces for fiber calculation*/
	fFlowStress.Dimension(fNumSD);
	fResidual.Dimension(fNumFibStress);
	fMod1.Dimension(fNumSD);
	fMod2.Dimension(fNumSD);
	fiK_f.Dimension(fNumFibStress);
	fG.Dimension(fNumFibStress);
	
	int num_inv = 4;
	fI.Dimension(num_inv);
	fdW.Dimension(num_inv);
	fddW.Dimension(num_inv);

}
	
/***********************************************************************
 * Protected
 ***********************************************************************/

/*computes integrated fiber stress in local frame*/
void NLV_Ortho::ComputeFiberStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
			dSymMatrixT& FiberStress, const int pindex)
{
	FiberStress = 0.0;
	if (pindex < 0)  /*eq*/
	{ 
		/*calculate structural invariants*/
		double& I4 = fI[0] = FiberStretch[0]; /*Cf11*/
		double& I6 = fI[2] = FiberStretch[1]; /*Cf22*/
		double& I5 = fI[1] = FiberStretch[0]*FiberStretch[0] + FiberStretch[4]*FiberStretch[4] 
				+ FiberStretch[5]*FiberStretch[5]; /*I5=Cf11^2+Cf13^2+Cf12^2*/
		double& I7 = fI[3] = FiberStretch[1]*FiberStretch[1] + FiberStretch[3]*FiberStretch[3] 
				+ FiberStretch[5]*FiberStretch[5]; 	/*I7=Cf22^2+Cf23^2+Cf12^2*/
		
/*		if (I4 < 1.0) I4 = 1.0;	
		if (I6 < 1.0) I6 = 1.0;
		if (I5 < 1.0) I5 = 1.0;
		if (I7 < 1.0) I7 = 1.0;
*/		
		/*calc 2dW/dI*/
		dWdI(fI, fdW, pindex);
		
		const double& s4 = fdW[0];
		const double& s5 = fdW[1];
		const double& s6 = fdW[2];
		const double& s7 = fdW[3];
				
		FiberStress = 0.0;
		FiberStress[0] = s4 + 2.0*s5*FiberStretch[0]; 
		FiberStress[1] = s6 + 2.0*s7*FiberStretch[1];
		FiberStress[3] = s7*FiberStretch[3];
		FiberStress[4] = s5*FiberStretch[4];
		FiberStress[5] = s5*FiberStretch[5] + s7*FiberStretch[5];
		
	}
	else {
		/*calculate elastic structural invariants*/
		double& I4 = fI[0] = FiberStretch[0]/FiberStretch_v[0]; /*Cf:M1/Cfv:M1*/
		double& I6 = fI[2] = FiberStretch[1]/FiberStretch_v[1]; /*Cf:M2/Cfv:M2*/
		
		fInverse.Inverse(FiberStretch_v);
		double A1 = FiberStretch[0]*fInverse[0] + FiberStretch[4]*fInverse[4] + FiberStretch[5]*fInverse[5];
		double A2 = FiberStretch[0]*fInverse[4] + FiberStretch[4]*fInverse[2] + FiberStretch[5]*fInverse[3];
		double A3 = FiberStretch[0]*fInverse[5] + FiberStretch[4]*fInverse[3] + FiberStretch[5]*fInverse[1];
		double& I5 = fI[1] = (FiberStretch[0]*A1 + FiberStretch[4]*A2 + FiberStretch[5]*A3)/FiberStretch_v[0]; /* (C.Cv^-1.C):M1/(Cv:M1) */

		double B1 = FiberStretch[1]*fInverse[1] + FiberStretch[3]*fInverse[3] + FiberStretch[5]*fInverse[5];
		double B2 = FiberStretch[1]*fInverse[3] + FiberStretch[3]*fInverse[2] + FiberStretch[5]*fInverse[4];
		double B3 = FiberStretch[1]*fInverse[5] + FiberStretch[3]*fInverse[4] + FiberStretch[5]*fInverse[0];
		double& I7 = fI[3] = (FiberStretch[1]*B1+ FiberStretch[3]*B2 + FiberStretch[5]*B3)/FiberStretch_v[1]; /* (C.Cv^-1.C):M2/(Cv:M2)*/

/*		if (I4 < 1.0) I4 = 1.0;	
		if (I6 < 1.0) I6 = 1.0;
		if (I5 < 1.0) I5 = 1.0;
		if (I7 < 1.0) I7 = 1.0;
*/		
		/*calc 2dWeq/dI*/
		dWdI(fI, fdW, pindex);
		
		const double& s4 = fdW[0];
		const double& s5 = fdW[1];
		const double& s6 = fdW[2];
		const double& s7 = fdW[3];
		
		FiberStress[0] = s4/FiberStretch_v[0] + 2.0*s5*A1/FiberStretch_v[0]; 
		FiberStress[1] = s6/FiberStretch_v[1] + 2.0*s7*B1/FiberStretch_v[1];
		FiberStress[2] = 0.0;
		FiberStress[3] = s7*B2/FiberStretch_v[1];
		FiberStress[4] = s5*A2/FiberStretch_v[0];
		FiberStress[5] = s5*A3/FiberStretch_v[0] + s7*B3/FiberStretch_v[1];
	}
}
	
/*computes integrated fiber stress in local frame*/
void NLV_Ortho::ComputeFlowStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dSymMatrixT& FlowStress, const int pindex)
{

	/*calculate elastic structural invariants*/
	double& I4 = fI[0] = FiberStretch[0]/FiberStretch_v[0]; /*Cf:M1/Cfv:M1*/
	double& I6 = fI[2] = FiberStretch[1]/FiberStretch_v[1]; /*Cf:M2/Cfv:M2*/
		
	fInverse.Inverse(FiberStretch_v);
	double A1 = FiberStretch[0]*fInverse[0] + FiberStretch[4]*fInverse[4] + FiberStretch[5]*fInverse[5];
	double A2 = FiberStretch[0]*fInverse[4] + FiberStretch[4]*fInverse[2] + FiberStretch[5]*fInverse[3];
	double A3 = FiberStretch[0]*fInverse[5] + FiberStretch[4]*fInverse[3] + FiberStretch[5]*fInverse[1];
	double& I5 = fI[1] = (FiberStretch[0]*A1 + FiberStretch[4]*A2 + FiberStretch[5]*A3)/FiberStretch_v[0]; /* (C.Cv^-1.C):M1/(Cv:M1) */

	double B1 = FiberStretch[1]*fInverse[1] + FiberStretch[3]*fInverse[3] + FiberStretch[5]*fInverse[5];
	double B2 = FiberStretch[1]*fInverse[3] + FiberStretch[3]*fInverse[2] + FiberStretch[5]*fInverse[4];
	double B3 = FiberStretch[1]*fInverse[5] + FiberStretch[3]*fInverse[4] + FiberStretch[5]*fInverse[0];
	double& I7 = fI[3] = (FiberStretch[1]*B1+ FiberStretch[3]*B2 + FiberStretch[5]*B3)/FiberStretch_v[1]; /* (C.Cv^-1.C):M2/(Cv:M2)*/

/*	cout << "\nC: "<<FiberStretch;
	cout << "\nCv: "<<FiberStretch_v;
	cout << "\nI4: "<<I4;
	cout << "\nI6: "<<I6;
	cout << "\nI5: "<<I5;
	cout << "\nI7: "<<I7;
*/	
/*	if (I4 < 1.0) I4 = 1.0;	
	if (I6 < 1.0) I6 = 1.0;
	if (I5 < 1.0) I5 = 1.0;
	if (I7 < 1.0) I7 = 1.0;
*/		
	/*calc 2dW/dI*/
	dWdI(fI, fdW, pindex);
		
	const double& s4 = fdW[0];
	const double& s5 = fdW[1];
	const double& s6 = fdW[2];
	const double& s7 = fdW[3];

	FlowStress[0] = (s4*I4+s5*I5 + s5*A1*A1)/FiberStretch_v[0] + s7*B3*B3/FiberStretch_v[1]; 
	FlowStress[1] = (s6*I6+s7*I7 + s7*B1*B1)/FiberStretch_v[1] + s5*A3*A3/FiberStretch_v[0];
	FlowStress[2] = s5*A2*A2/FiberStretch_v[0] + s7*B2*B2/FiberStretch_v[1];
	
	FlowStress[3] = s5*A3*A2/FiberStretch_v[0] + s7*B1*B2/FiberStretch_v[1];
	FlowStress[4] = s5*A1*A2/FiberStretch_v[0] + s7*B3*B2/FiberStretch_v[1];
	FlowStress[5] = s5*A1*A3/FiberStretch_v[0] + s7*B3*B1/FiberStretch_v[1];
}

/*computes integrated moduli in local frame*/
void NLV_Ortho::ComputeFiberMod (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dSymMatrixT& FiberStress, dMatrixT& FiberMod,  const int pindex)
{
	FiberMod = 0.0;
	
	if (pindex < 0) /*eq*/
	{
		/*2del_SEQ/del_C + 2del_SNEQ/del_C*/
		/*calculate structural invariants*/
		double& I4 = fI[0] = FiberStretch[0]; /*Cf11*/
		double& I6 = fI[2] = FiberStretch[1]; /*Cf22*/
		double& I5 = fI[1] = FiberStretch[0]*FiberStretch[0] + FiberStretch[4]*FiberStretch[4] 
				+ FiberStretch[5]*FiberStretch[5]; /*I5=Cf11^2+Cf13^2+Cf12^2*/
		double& I7 = fI[3] = FiberStretch[1]*FiberStretch[1] + FiberStretch[3]*FiberStretch[3] 
				+ FiberStretch[5]*FiberStretch[5]; 	/*I7=Cf22^2+Cf23^2+Cf12^2*/
	
/*		if (I4 < 1.0) I4 = 1.0;	
		if (I6 < 1.0) I6 = 1.0;
		if (I5 < 1.0) I5 = 1.0;
		if (I7 < 1.0) I7 = 1.0;
*/	
		/*calc 2dW/dI and 2ddW/ddI*/
		ddWddI(fI, fdW, fddW, pindex);
		
		const double& s4 = fdW[0];
		const double& s5 = fdW[1];
		const double& s6 = fdW[2];
		const double& s7 = fdW[3];
	
		const double& d4 = fddW[0];
		const double& d5 = fddW[1];
		const double& d6 = fddW[2];
		const double& d7 = fddW[3];

		/* Mi x Mi*/
		FiberMod(0,0) += 2.0*d4;
		FiberMod(1,1) += 2.0*d6;
	
		/*(C.Mi + Mi.C) x (C.Mi + Mi.C)*/
		fMod1 = 0.0;
		fMod1[0] = 2.0*FiberStretch[0];
		fMod1[4] = FiberStretch[4];
		fMod1[5] = FiberStretch[5];
	
		fMod2 = 0.0;
		fMod2[1] = 2.0*FiberStretch[1];
		fMod2[3] = FiberStretch[3];
		fMod2[5] = FiberStretch[5];

		fMod3.DyadAB(fMod1,fMod1);
		FiberMod.AddScaled(2.0*d5, fMod3);
		fMod3.DyadAB(fMod2,fMod2);
		FiberMod.AddScaled(2.0*d7, fMod3);

		/* coef (1 o Mi + Mi o 1)*/
		/*4*dW/dIi*/
		double c5 = 2.0*s5;
		double c7 = 2.0*s7;

		/* fiber 1 */
		FiberMod(0,0) += 2.0*c5;
		FiberMod(4,4) += 0.5*c5;
		FiberMod(5,5) += 0.5*c5;
		/* fiber 2 */
		FiberMod(1,1) += 2.0*c7;
		FiberMod(3,3) += 0.5*c7;
		FiberMod(5,5) += 0.5*c7;
	}
	else /*neq*/
		ComputeCalg(FiberStretch, FiberStretch_v, FiberMod, pindex);
}

void NLV_Ortho::ComputeCalg(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,  dMatrixT& Calg, const int pindex)
{		
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
		
	/*Compute flow stress*/
	ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);

	/* FlowStress : (Cv.M1 + M1.Cv) */
	double sig1 = fFlowStress[0]*2.0*FiberStretch_v[0] + 2.0*fFlowStress[4]*FiberStretch_v[4] 
				+ 2.0*fFlowStress[5]*FiberStretch_v[5];
	double eta1 = fVisc1_f[pindex]->Function(sig1);
	double coeff1 = dt*sig1/(2.0*eta1);

	/* FlowStress : (Cv.M2 + M2.Cv) */
	double sig2 = fFlowStress[1]*2.0*FiberStretch_v[1] + 2.0*fFlowStress[3]*FiberStretch_v[3] 
				+ 2.0*fFlowStress[5]*FiberStretch_v[5];
	double eta2 = fVisc2_f[pindex]->Function(sig2);
	double coeff2 = dt*sig2/(2.0*eta2);

	/*deta_i/dsig_i*/
	double deta1 = fVisc1_f[pindex]->DFunction(sig1);
	double deta2 = fVisc2_f[pindex]->DFunction(sig2);
				
	double dcoeff1 = (1.0 - 1.0/eta1*deta1*sig1);
	double dcoeff2 = (1.0 - 1.0/eta2*deta2*sig2);

	/*compute -G= -dR/dCfv*/
	dFlowdC(FiberStretch, FiberStretch_v, fModMat,  pindex);	

	/*fiber 1*/
	/*(Cv.M1+M1.Cv)*/
	fMod1 = 0.0;
	fMod1[0] = 2.0*FiberStretch_v[0];
	fMod1[4] = FiberStretch_v[4];
	fMod1[5] = FiberStretch_v[5];
			
	/*(Cv.M1 + M1.Cv) : dFlowStress/dC*/
	fMod2[0] = 2.0*fModMat(0,0)*FiberStretch_v[0] + 2.0*fModMat(4,0)*FiberStretch_v[4] + 2.0*fModMat(5,0)*FiberStretch_v[5];
	fMod2[1] = 2.0*fModMat(0,1)*FiberStretch_v[0] + 2.0*fModMat(4,1)*FiberStretch_v[4] + 2.0*fModMat(5,1)*FiberStretch_v[5]; 
	fMod2[2] = 2.0*fModMat(0,2)*FiberStretch_v[0] + 2.0*fModMat(4,2)*FiberStretch_v[4] + 2.0*fModMat(5,2)*FiberStretch_v[5]; 
	fMod2[3] = 2.0*fModMat(0,3)*FiberStretch_v[0] + 2.0*fModMat(4,3)*FiberStretch_v[4] + 2.0*fModMat(5,3)*FiberStretch_v[5]; 
	fMod2[4] = 2.0*fModMat(0,4)*FiberStretch_v[0] + 2.0*fModMat(4,4)*FiberStretch_v[4] + 2.0*fModMat(5,4)*FiberStretch_v[5]; 
	fMod2[5] = 2.0*fModMat(0,5)*FiberStretch_v[0] + 2.0*fModMat(4,5)*FiberStretch_v[4] + 2.0*fModMat(5,5)*FiberStretch_v[5]; 

	/*(Cv.M1 + M1.Cv) x (Cv.Mi + Mi.Cv) : dFlowStress/dC*/
	fMod3.DyadAB(fMod1, fMod2);
	fG.SetToScaled(dt/(2.0*eta1)*dcoeff1, fMod3);	

	/*fiber 2*/
	/* 1/eta2 (Cv.M2+M2.Cv) */
	fMod1 = 0.0;
	fMod1[1] = 2.0*FiberStretch_v[1];
	fMod1[3] = FiberStretch_v[3];
	fMod1[5] = FiberStretch_v[5];

	fMod2[0] = 2.0*fModMat(1,0)*FiberStretch_v[1] + 2.0*fModMat(3,0)*FiberStretch_v[3] + 2.0*fModMat(5,0)*FiberStretch_v[5];
	fMod2[1] = 2.0*fModMat(1,1)*FiberStretch_v[1] + 2.0*fModMat(3,1)*FiberStretch_v[3] + 2.0*fModMat(5,1)*FiberStretch_v[5]; 
	fMod2[2] = 2.0*fModMat(1,2)*FiberStretch_v[1] + 2.0*fModMat(3,2)*FiberStretch_v[3] + 2.0*fModMat(5,2)*FiberStretch_v[5]; 
	fMod2[3] = 2.0*fModMat(1,3)*FiberStretch_v[1] + 2.0*fModMat(3,3)*FiberStretch_v[3] + 2.0*fModMat(5,3)*FiberStretch_v[5]; 
	fMod2[4] = 2.0*fModMat(1,4)*FiberStretch_v[1] + 2.0*fModMat(3,4)*FiberStretch_v[3] + 2.0*fModMat(5,4)*FiberStretch_v[5]; 
	fMod2[5] = 2.0*fModMat(1,5)*FiberStretch_v[1] + 2.0*fModMat(3,5)*FiberStretch_v[3] + 2.0*fModMat(5,5)*FiberStretch_v[5]; 

	/*(Cv.M2 + M2.Cv) x (Cv.M2 + M2.Cv) : dFlowStress/dC*/
	fMod3.DyadAB(fMod1, fMod2);
	fG.AddScaled(dt/(2.0*eta2)*dcoeff2, fMod3);
		
	/*convert to matrix system*/
	fG(0,3) *= 2.0;
	fG(0,4) *= 2.0;
	fG(0,5) *= 2.0;		
	fG(1,3) *= 2.0;
	fG(1,4) *= 2.0;
	fG(1,5) *= 2.0;
	fG(2,3) *= 2.0;
	fG(2,4) *= 2.0;
	fG(2,5) *= 2.0;
	fG(3,3) *= 2.0;
	fG(3,4) *= 2.0;
	fG(3,5) *= 2.0;
	fG(4,3) *= 2.0;
	fG(4,4) *= 2.0;
	fG(4,5) *= 2.0;
	fG(5,3) *= 2.0;
	fG(5,4) *= 2.0;
	fG(5,5) *= 2.0;	
	
//	cout << "\nFlowStress: "<<fFlowStress;
//	cout << "\nFiberStretch: "<<FiberStretch;
//	cout << "\nFiberStretch_v: "<<FiberStretch_v;	
//	cout << "\ndSig/dC: "<<fModMat;
//	cout << "\nfG: "<<fG;

	/*compute K = dR/dCv*/
	/*fiber family 1*/
	fiK_f.ReducedIndexI();
		

	/* coef (1 o Mi + Mi o 1)*/
	/* fiber 1 */
	fiK_f(0,0) -= 2.0*coeff1;
	fiK_f(4,4) -= 0.5*coeff1;
	fiK_f(5,5) -= 0.5*coeff1;
	/* fiber 2 */
	fiK_f(1,1) -= 2.0*coeff2;
	fiK_f(3,3) -= 0.5*coeff2;
	fiK_f(5,5) -= 0.5*coeff2;
	
	/*dFlowStress/dCv*/
	dFlowdCv(FiberStretch, FiberStretch_v, fModMat, pindex);
	
	/*fiber 1*/
	/* 1/(2*eta1) (Cv.M1+M1.Cv) x (S.M1+M1.S)*/
	fMod1 = 0.0;
	fMod1[0] = 2.0*FiberStretch_v[0];
	fMod1[4] = FiberStretch_v[4];
	fMod1[5] = FiberStretch_v[5];
		
	fMod2 = 0.0;
	fMod2[0] = 2.0*fFlowStress[0];
	fMod2[4] = fFlowStress[4];
	fMod2[5] = fFlowStress[5];
		
	fMod3.DyadAB(fMod1, fMod2);
	fiK_f.AddScaled(-dt/(2.0*eta1)*dcoeff1, fMod3);
		
	/*(Cv.M1 + M1.Cv) : dFlowStress/dCv*/
	fMod2[0] = 2.0*fModMat(0,0)*FiberStretch_v[0] + 2.0*fModMat(4,0)*FiberStretch_v[4] + 2.0*fModMat(5,0)*FiberStretch_v[5];
	fMod2[1] = 2.0*fModMat(0,1)*FiberStretch_v[0] + 2.0*fModMat(4,1)*FiberStretch_v[4] + 2.0*fModMat(5,1)*FiberStretch_v[5]; 
	fMod2[2] = 2.0*fModMat(0,2)*FiberStretch_v[0] + 2.0*fModMat(4,2)*FiberStretch_v[4] + 2.0*fModMat(5,2)*FiberStretch_v[5]; 
	fMod2[3] = 2.0*fModMat(0,3)*FiberStretch_v[0] + 2.0*fModMat(4,3)*FiberStretch_v[4] + 2.0*fModMat(5,3)*FiberStretch_v[5]; 
	fMod2[4] = 2.0*fModMat(0,4)*FiberStretch_v[0] + 2.0*fModMat(4,4)*FiberStretch_v[4] + 2.0*fModMat(5,4)*FiberStretch_v[5]; 
	fMod2[5] = 2.0*fModMat(0,5)*FiberStretch_v[0] + 2.0*fModMat(4,5)*FiberStretch_v[4] + 2.0*fModMat(5,5)*FiberStretch_v[5]; 
		
	/*(Cv.M1 + M1.Cv) x (Cv.Mi + Mi.Cv) : dFlowStress/dCv*/
	fMod3.DyadAB(fMod1, fMod2);
	fiK_f.AddScaled(-dt/(2.0*eta1)*dcoeff1, fMod3);		

	/*fiber 2*/
	/* 1/eta2 (Cv.M2+M2.Cv) x (S.M2+M2.S)*/
	fMod1 = 0.0;
	fMod1[1] = 2.0*FiberStretch_v[1];
	fMod1[3] = FiberStretch_v[3];
	fMod1[5] = FiberStretch_v[5];
		
	fMod2 = 0.0;
	fMod2[1] = 2.0*fFlowStress[1];
	fMod2[3] = fFlowStress[3];
	fMod2[5] = fFlowStress[5];
		
	fMod3.DyadAB(fMod1, fMod2);
	fiK_f.AddScaled(-dt/(2.0*eta2)*dcoeff2, fMod3);

	/*(Cv.M2 + M2.Cv) : dFlowStress/dCv*/
	fMod2[0] = 2.0*fModMat(1,0)*FiberStretch_v[1] + 2.0*fModMat(3,0)*FiberStretch_v[3] + 2.0*fModMat(5,0)*FiberStretch_v[5];
	fMod2[1] = 2.0*fModMat(1,1)*FiberStretch_v[1] + 2.0*fModMat(3,1)*FiberStretch_v[3] + 2.0*fModMat(5,1)*FiberStretch_v[5]; 
	fMod2[2] = 2.0*fModMat(1,2)*FiberStretch_v[1] + 2.0*fModMat(3,2)*FiberStretch_v[3] + 2.0*fModMat(5,2)*FiberStretch_v[5]; 
	fMod2[3] = 2.0*fModMat(1,3)*FiberStretch_v[1] + 2.0*fModMat(3,3)*FiberStretch_v[3] + 2.0*fModMat(5,3)*FiberStretch_v[5]; 
	fMod2[4] = 2.0*fModMat(1,4)*FiberStretch_v[1] + 2.0*fModMat(3,4)*FiberStretch_v[3] + 2.0*fModMat(5,4)*FiberStretch_v[5]; 
	fMod2[5] = 2.0*fModMat(1,5)*FiberStretch_v[1] + 2.0*fModMat(3,5)*FiberStretch_v[3] + 2.0*fModMat(5,5)*FiberStretch_v[5]; 
	/*(Cv.M2 + M2.Cv) x (Cv.M2 + M2.Cv) : dFlowStress/dCv*/
	fMod3.DyadAB(fMod1, fMod2);
	fiK_f.AddScaled(-dt/(2.0*eta2)*dcoeff2, fMod3);
	
	/*convert to matrix system*/
	fiK_f(0,3) *= 2.0;
	fiK_f(0,4) *= 2.0;
	fiK_f(0,5) *= 2.0;		
	fiK_f(1,3) *= 2.0;
	fiK_f(1,4) *= 2.0;
	fiK_f(1,5) *= 2.0;
	fiK_f(2,3) *= 2.0;
	fiK_f(2,4) *= 2.0;
	fiK_f(2,5) *= 2.0;
	fiK_f(3,3) *= 2.0;
	fiK_f(3,4) *= 2.0;
	fiK_f(3,5) *= 2.0;
	fiK_f(4,3) *= 2.0;
	fiK_f(4,4) *= 2.0;
	fiK_f(4,5) *= 2.0;
	fiK_f(5,3) *= 2.0;
	fiK_f(5,4) *= 2.0;
	fiK_f(5,5) *= 2.0;
	/*K^-1.(-R)*/
	fiK_f.Inverse();

//	cout << "\ndSig/dCv"<<fModMat;
//	cout << "\nfiK: "<<fiK_f;	
	
	/* 2dSNEQ/dCv*/
	fModMat = 0.0;
	/*calculate elastic structural invariants*/
	double& I4 = fI[0] = FiberStretch[0]/FiberStretch_v[0]; /*Cf:M1/Cfv:M1*/
	double& I6 = fI[2] = FiberStretch[1]/FiberStretch_v[1]; /*Cf:M2/Cfv:M2*/
		
	fInverse.Inverse(FiberStretch_v);
	double A1 = FiberStretch[0]*fInverse[0] + FiberStretch[4]*fInverse[4] + FiberStretch[5]*fInverse[5];
	double A2 = FiberStretch[0]*fInverse[4] + FiberStretch[4]*fInverse[2] + FiberStretch[5]*fInverse[3];
	double A3 = FiberStretch[0]*fInverse[5] + FiberStretch[4]*fInverse[3] + FiberStretch[5]*fInverse[1];
	double& I5 = fI[1] = (FiberStretch[0]*A1 + FiberStretch[4]*A2 + FiberStretch[5]*A3)/FiberStretch_v[0]; /* (C.Cv^-1.C):M1/(Cv:M1) */

	double B1 = FiberStretch[1]*fInverse[1] + FiberStretch[3]*fInverse[3] + FiberStretch[5]*fInverse[5];
	double B2 = FiberStretch[1]*fInverse[3] + FiberStretch[3]*fInverse[2] + FiberStretch[5]*fInverse[4];
	double B3 = FiberStretch[1]*fInverse[5] + FiberStretch[3]*fInverse[4] + FiberStretch[5]*fInverse[0];
	double& I7 = fI[3] = (FiberStretch[1]*B1+ FiberStretch[3]*B2 + FiberStretch[5]*B3)/FiberStretch_v[1]; /* (C.Cv^-1.C):M2/(Cv:M2)*/

/*	if (I4 < 1.0) I4 = 1.0;	
	if (I6 < 1.0) I6 = 1.0;
	if (I5 < 1.0) I5 = 1.0;
	if (I7 < 1.0) I7 = 1.0;
*/		
	/*calc 2dW/dI and 2ddW/ddI*/
	ddWddI(fI, fdW, fddW, pindex);
		
	const double& s4 = fdW[0];
	const double& s5 = fdW[1];
	const double& s6 = fdW[2];
	const double& s7 = fdW[3];
	
	const double& d4 = fddW[0];
	const double& d5 = fddW[1];
	const double& d6 = fddW[2];
	const double& d7 = fddW[3];
	
	double Iv4_2 = FiberStretch_v[0]*FiberStretch_v[0];
	double Iv6_2 = FiberStretch_v[1]*FiberStretch_v[1];
	
	/*Mi x Mi*/
	fModMat(0,0) -= 2.0*(s4 + d4*I4)/Iv4_2;
	fModMat(1,1) -= 2.0*(s6 + d6*I6)/Iv6_2;

	/*(Mi.C.Cv^-1 + Cv^-1.C.Mi) x Mi*/
	double c5 = 2.0*(s5 + d5*I5)/Iv4_2;
	double c7 = 2.0*(s7 + d7*I7)/Iv6_2;
	
	fModMat(0,0) -= c5 * 2.0*A1;
	fModMat(4,0) -= c5 * A2;
	fModMat(5,0) -= c5 * A3;

	fModMat(1,1) -= c7 * 2.0*B1;
	fModMat(3,1) -= c7 * B2;
	fModMat(5,1) -= c7 * B3;
		
	/* Mi.C.Cv^-1 o Cv^-1 + Cv^-1 o Mi.C.Cv^-1*/
	c5 = 2.0*s5/FiberStretch_v[0];
	c7 = 2.0*s7/FiberStretch_v[1];

	fModMat(0,0) -= 2.0*c5*A1*fInverse[0]; 
	fModMat(0,1) -= 2.0*c5*A3*fInverse[5] + 2.0*c7*B3*fInverse[5]; 
	fModMat(0,2) -= 2.0*c5*A2*fInverse[4];
	fModMat(0,3) -= c7*B3*fInverse[4] + c5*(A3*fInverse[4] + A2*fInverse[5]);
	fModMat(0,4) -= c5*(A2*fInverse[0] + A1*fInverse[4]);
	fModMat(0,5) -= c7*B3*fInverse[0] + c5*(A1*fInverse[5] + A3*fInverse[0]);
	
	fModMat(1,1) -= 2.0*c7*B1*fInverse[1]; 
	fModMat(1,3) -= c7*B1*fInverse[3];
	fModMat(1,5) -= c7*B1*fInverse[5];

	fModMat(2,1) -= 2.0*c7*B2*fInverse[3]; 
	fModMat(2,3) -= c7*B2*fInverse[2];
	fModMat(2,5) -= c7*B2*fInverse[4];
	
	fModMat(3,1) -= c7*(B2*fInverse[1] + B1*fInverse[3]); 
	fModMat(3,3) -= 0.5*c7*(B1*fInverse[2] + B2*fInverse[3]);
	fModMat(3,5) -= 0.5*c7*(B1*fInverse[4] + B2*fInverse[5]);

	fModMat(4,0) -= c5*A1*fInverse[4]; 
	fModMat(4,1) -= c5*A3*fInverse[3] + c7*(B3*fInverse[3] + B2*fInverse[5]); 
	fModMat(4,2) -= c5*A2*fInverse[2];
	fModMat(4,3) -= 0.5*c5*(A3*fInverse[2] + A2*fInverse[3]) + 0.5*c7*(B3*fInverse[2] + B2*fInverse[4]);
	fModMat(4,4) -= 0.5*c5*(A1*fInverse[2] + A2*fInverse[4]);
	fModMat(4,5) -= 0.5*c7*(B2*fInverse[0] + B3*fInverse[4]) + 0.5*c5*(A1*fInverse[3] + A3*fInverse[4]);

	fModMat(5,0) -= c5*A1*fInverse[5]; 
	fModMat(5,1) -= c7*(B3*fInverse[1] + B1*fInverse[5]) + c5*A3*fInverse[1];
	fModMat(5,2) -= c5*A2*fInverse[3]; 
	fModMat(5,3) -= 0.5*c5*(A2*fInverse[1] + A3*fInverse[3]) + 0.5*c7*(B3*fInverse[3] + B1*fInverse[4]);
	fModMat(5,4) -= 0.5*c5*(A1*fInverse[3] + A2*fInverse[5]);
	fModMat(5,5) -= 0.5*c5*(A1*fInverse[1] + A3*fInverse[5]) + 0.5*c7*(B1*fInverse[0] + B3*fInverse[5]);

	/* (Mi.C.Cv^-1 + Cv^-1.C.Mi) x Cv^-1.C.Mi.C.Cv^-1*/
	/*Cv^-1.C.M1.C.Cv^-1*/
	c5 = 2.0*d5/Iv4_2;
	c7 = 2.0*d7/Iv6_2;

	/*(M1.C.Cv^-1 + Cv^-1.C.M1)*/
	fMod1[0] = 2.0*A1; 
	fMod1[1] = 0.0;
	fMod1[2] = 0.0;
	fMod1[3] = 0.0;
	fMod1[4] = A2;
	fMod1[5] = A3;

	/*Cv^-1.C.M1.C.Cv^-1*/
 	fMod2[0] = A1*A1; 
	fMod2[1] = A3*A3;
	fMod2[2] = A2*A2;
	fMod2[3] = A3*A2;
	fMod2[4] = A1*A2;
	fMod2[5] = A1*A3;

	fMod3.DyadAB(fMod1,fMod2);
	fModMat.AddScaled(-c5, fMod3);

	/*(M2.C.Cv^-1 + Cv^-1.C.M2)*/
 	fMod1[0] = 0.0; 
	fMod1[1] = 2.0*B1;
	fMod1[2] = 0.0;
	fMod1[3] = B2;
	fMod1[4] = 0.0;
	fMod1[5] = B3;

	/*Cv^-1.C.M2.C.Cv^-1*/
 	fMod2[0] = B3*B3; 
	fMod2[1] = B1*B1;
	fMod2[2] = B2*B2;
	fMod2[3] = B1*B2;
	fMod2[4] = B3*B2;
	fMod2[5] = B3*B1;

	fMod3.DyadAB(fMod1,fMod2);
	fModMat.AddScaled(-c7, fMod3);

	fModMat(0,3) *= 2.0;
	fModMat(0,4) *= 2.0;
	fModMat(0,5) *= 2.0;		
	fModMat(1,3) *= 2.0;
	fModMat(1,4) *= 2.0;
	fModMat(1,5) *= 2.0;
	fModMat(2,3) *= 2.0;
	fModMat(2,4) *= 2.0;
	fModMat(2,5) *= 2.0;
	fModMat(3,3) *= 2.0;
	fModMat(3,4) *= 2.0;
	fModMat(3,5) *= 2.0;
	fModMat(4,3) *= 2.0;
	fModMat(4,4) *= 2.0;
	fModMat(4,5) *= 2.0;
	fModMat(5,3) *= 2.0;
	fModMat(5,4) *= 2.0;
	fModMat(5,5) *= 2.0;
	
	/*2.0 dSNEQ/dCv : Delta C/Delta Cv */
	/* Delta Cv = K^-1: (-G) Delta C */
	fMod3.MultAB(fiK_f, fG);
//	cout << "\nfiK_f: "<<fiK_f;
//	cout << "\nfG: "<<fG;
//	cout << "\niK.G: "<<fMod3;
//  cout << "\ndSNEQ/DCv: "<<fModMat;
	
	Calg.MultAB(fModMat, fMod3);

	/*convert back to tensor*/
	Calg(0,3) *= 0.5;
	Calg(0,4) *= 0.5;
	Calg(0,5) *= 0.5;		
	Calg(1,3) *= 0.5;
	Calg(1,4) *= 0.5;
	Calg(1,5) *= 0.5;
	Calg(2,3) *= 0.5;
	Calg(2,4) *= 0.5;
	Calg(2,5) *= 0.5;
	Calg(3,3) *= 0.5;
	Calg(3,4) *= 0.5;
	Calg(3,5) *= 0.5;
	Calg(4,3) *= 0.5;
	Calg(4,4) *= 0.5;
	Calg(4,5) *= 0.5;
	Calg(5,3) *= 0.5;
	Calg(5,4) *= 0.5;
	Calg(5,5) *= 0.5;
	
//	cout << "\nCalg: "<<Calg;
	
	/*+ 2.0 dSNEQ/dC*/		
	/*Mi x Mi*/
	Calg(0,0) += 2.0*d4/Iv4_2;
	Calg(1,1) += 2.0*d6/Iv6_2;
	
	/*(Cv^-1.C.Mi + Mi.C.Cv^-1) x (Cv^-1.C.Mi + Mi.C.Cv^-1)*/
	fMod1 = 0.0;
	fMod1[0] = 2.0*A1;
	fMod1[4] = A2;
	fMod1[5] = A3;
	
	fMod2 = 0.0;
	fMod2[1] = 2.0*B1;
	fMod2[3] = B2;
	fMod2[5] = B3;

	fMod3.DyadAB(fMod1,fMod1);
	Calg.AddScaled(2.0*d5/Iv4_2, fMod3);
	fMod3.DyadAB(fMod2,fMod2);
	Calg.AddScaled(2.0*d7/Iv6_2, fMod3);	
	
	/*Cv^-1 o Mi + Mi o Cv^-1*/
	c5 = 2.0*s5/FiberStretch_v[0];
	c7 = 2.0*s7/FiberStretch_v[1];
		
	Calg(0,0) += 2.0*c5*fInverse[0];
	Calg(0,4) += c5*fInverse[4];
	Calg(0,5) += c5*fInverse[5];
	
	Calg(1,1) += 2.0*c7*fInverse[1];
	Calg(1,3) += c7*fInverse[3];
	Calg(1,5) += c7*fInverse[5];
	
	Calg(3,1) += c7*fInverse[3];
	Calg(3,3) += 0.5*c7*fInverse[2];
	Calg(3,5) += 0.5*c7*fInverse[4];
	
	Calg(4,0) += c5*fInverse[4];
	Calg(4,4) += 0.5*c5*fInverse[2];
	Calg(4,5) += 0.5*c5*fInverse[3];
	
	Calg(5,0) += c5*fInverse[5];
	Calg(5,1) += c7*fInverse[5];
	Calg(5,3) += 0.5*c7*fInverse[4];
	Calg(5,4) += 0.5*c5*fInverse[3];
	Calg(5,5) += 0.5*(c7*fInverse[0]+c5*fInverse[1]);
}

/*local newton loop for viscous stretch tensor*/ 
void NLV_Ortho::Compute_Cv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_vn, 
	dSymMatrixT& FiberStretch_v, const int pindex)
{	
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
		
	/*Compute flow stress*/
	ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);

	/* FlowStress : (Cv.M1 + M1.Cv) */
	double sig1 = fFlowStress[0]*2.0*FiberStretch_v[0] + 2.0*fFlowStress[4]*FiberStretch_v[4] 
				+ 2.0*fFlowStress[5]*FiberStretch_v[5];
	double eta1 = fVisc1_f[pindex]->Function(sig1);
	double coeff1 = dt*sig1/(2.0*eta1);

	/* FlowStress : (Cv.M2 + M2.Cv) */
	double sig2 = fFlowStress[1]*2.0*FiberStretch_v[1] + 2.0*fFlowStress[3]*FiberStretch_v[3] 
				+ 2.0*fFlowStress[5]*FiberStretch_v[5];
	double eta2 = fVisc2_f[pindex]->Function(sig2);
	double coeff2 = dt*sig2/(2.0*eta2);

	fResidual[0] = FiberStretch_v[0] - FiberStretch_vn[0];
	fResidual[1] = FiberStretch_v[1] - FiberStretch_vn[1];
	fResidual[2] = FiberStretch_v[2] - FiberStretch_vn[2];
	fResidual[3] = FiberStretch_v[3] - FiberStretch_vn[3];
	fResidual[4] = FiberStretch_v[4] - FiberStretch_vn[4];
	fResidual[5] = FiberStretch_v[5] - FiberStretch_vn[5];
//	cout << "\nfResidual1: "<<fResidual;
		
	/*fiber family 1*/
	/* coeffi* (Cv.M1+M1.Cv)  */
	fResidual[0] -= coeff1 * 2.0*FiberStretch_v[0];
	fResidual[4] -= coeff1 * FiberStretch_v[4];
	fResidual[5] -= coeff1 * FiberStretch_v[5];
//	cout << "\nfResidual2: "<<fResidual;
		
	/*fiber family 2*/
	/* coeff* (Cv.M1+M1.Cv)  */
	fResidual[1] -= coeff2 * 2.0*FiberStretch_v[1];
	fResidual[3] -= coeff2 * FiberStretch_v[3];
	fResidual[5] -= coeff2 * FiberStretch_v[5];

//	cout << "\nfResidual3: "<<fResidual;
	double error = sqrt(fResidual[0]*fResidual[0] + fResidual[1]*fResidual[1] + fResidual[2]*fResidual[2] 
			+ 2.0*fResidual[3]*fResidual[3] + 2.0*fResidual[4]*fResidual[4] + 2.0*fResidual[5]*fResidual[5]);
	int iteration  = 0;	

	while (error > kSmall && iteration < 10)
	{		
//		cout << "\nfResidual: "<<fResidual;
		/*compute Kij = dR/dCv*/
		fiK_f.ReducedIndexI();
		
		/* coef (1 o Mi + Mi o 1)*/
		/* fiber 1 */
		fiK_f(0,0) -= 2.0*coeff1;
		fiK_f(4,4) -= 0.5*coeff1;
		fiK_f(5,5) -= 0.5*coeff1;
		/* fiber 2 */
		fiK_f(1,1) -= 2.0*coeff2;
		fiK_f(3,3) -= 0.5*coeff2;
		fiK_f(5,5) -= 0.5*coeff2;
	
		/*dFlowStress/dCv*/
		dFlowdCv(FiberStretch, FiberStretch_v, fModMat, pindex);

		/*deta_i/dsig_i*/
		double deta1 = fVisc1_f[pindex]->DFunction(sig1);
		double deta2 = fVisc2_f[pindex]->DFunction(sig2);
				
		double dcoeff1 = (1.0 - 1.0/eta1*deta1*sig1);
		double dcoeff2 = (1.0 - 1.0/eta2*deta2*sig2);

		/*fiber 1*/
		/* 1/(2*eta1) (Cv.M1+M1.Cv) x (S.M1+M1.S)*/
		fMod1 = 0.0;
		fMod1[0] = 2.0*FiberStretch_v[0];
		fMod1[4] = FiberStretch_v[4];
		fMod1[5] = FiberStretch_v[5];
		
		fMod2 = 0.0;
		fMod2[0] = 2.0*fFlowStress[0];
		fMod2[4] = fFlowStress[4];
		fMod2[5] = fFlowStress[5];
		
		fMod3.DyadAB(fMod1, fMod2);
		fiK_f.AddScaled(-dt/(2.0*eta1)*dcoeff1, fMod3);
		
		/*(Cv.M1 + M1.Cv) : dFlowStress/dCv*/
		fMod2[0] = 2.0*fModMat(0,0)*FiberStretch_v[0] + 2.0*fModMat(4,0)*FiberStretch_v[4] + 2.0*fModMat(5,0)*FiberStretch_v[5];
		fMod2[1] = 2.0*fModMat(0,1)*FiberStretch_v[0] + 2.0*fModMat(4,1)*FiberStretch_v[4] + 2.0*fModMat(5,1)*FiberStretch_v[5]; 
		fMod2[2] = 2.0*fModMat(0,2)*FiberStretch_v[0] + 2.0*fModMat(4,2)*FiberStretch_v[4] + 2.0*fModMat(5,2)*FiberStretch_v[5]; 
		fMod2[3] = 2.0*fModMat(0,3)*FiberStretch_v[0] + 2.0*fModMat(4,3)*FiberStretch_v[4] + 2.0*fModMat(5,3)*FiberStretch_v[5]; 
		fMod2[4] = 2.0*fModMat(0,4)*FiberStretch_v[0] + 2.0*fModMat(4,4)*FiberStretch_v[4] + 2.0*fModMat(5,4)*FiberStretch_v[5]; 
		fMod2[5] = 2.0*fModMat(0,5)*FiberStretch_v[0] + 2.0*fModMat(4,5)*FiberStretch_v[4] + 2.0*fModMat(5,5)*FiberStretch_v[5]; 
		
		/*(Cv.M1 + M1.Cv) x (Cv.Mi + Mi.Cv) : dFlowStress/dCv*/
		fMod3.DyadAB(fMod1, fMod2);
		fiK_f.AddScaled(-dt/(2.0*eta1)*dcoeff1, fMod3);		

		/*fiber 2*/
		/* 1/eta2 (Cv.M2+M2.Cv) x (S.M2+M2.S)*/
		fMod1 = 0.0;
		fMod1[1] = 2.0*FiberStretch_v[1];
		fMod1[3] = FiberStretch_v[3];
		fMod1[5] = FiberStretch_v[5];
		
		fMod2 = 0.0;
		fMod2[1] = 2.0*fFlowStress[1];
		fMod2[3] = fFlowStress[3];
		fMod2[5] = fFlowStress[5];
		
		fMod3.DyadAB(fMod1, fMod2);
		fiK_f.AddScaled(-dt/(2.0*eta2)*dcoeff2, fMod3);

		/*(Cv.M2 + M2.Cv) : dFlowStress/dCv*/
		fMod2[0] = 2.0*fModMat(1,0)*FiberStretch_v[1] + 2.0*fModMat(3,0)*FiberStretch_v[3] + 2.0*fModMat(5,0)*FiberStretch_v[5];
		fMod2[1] = 2.0*fModMat(1,1)*FiberStretch_v[1] + 2.0*fModMat(3,1)*FiberStretch_v[3] + 2.0*fModMat(5,1)*FiberStretch_v[5]; 
		fMod2[2] = 2.0*fModMat(1,2)*FiberStretch_v[1] + 2.0*fModMat(3,2)*FiberStretch_v[3] + 2.0*fModMat(5,2)*FiberStretch_v[5]; 
		fMod2[3] = 2.0*fModMat(1,3)*FiberStretch_v[1] + 2.0*fModMat(3,3)*FiberStretch_v[3] + 2.0*fModMat(5,3)*FiberStretch_v[5]; 
		fMod2[4] = 2.0*fModMat(1,4)*FiberStretch_v[1] + 2.0*fModMat(3,4)*FiberStretch_v[3] + 2.0*fModMat(5,4)*FiberStretch_v[5]; 
		fMod2[5] = 2.0*fModMat(1,5)*FiberStretch_v[1] + 2.0*fModMat(3,5)*FiberStretch_v[3] + 2.0*fModMat(5,5)*FiberStretch_v[5]; 
		/*(Cv.M2 + M2.Cv) x (Cv.M2 + M2.Cv) : dFlowStress/dCv*/
		fMod3.DyadAB(fMod1, fMod2);
		fiK_f.AddScaled(-dt/(2.0*eta2)*dcoeff2, fMod3);
	
		/*convert to matrix system*/
		fiK_f(0,3) *= 2.0;
		fiK_f(0,4) *= 2.0;
		fiK_f(0,5) *= 2.0;		
		fiK_f(1,3) *= 2.0;
		fiK_f(1,4) *= 2.0;
		fiK_f(1,5) *= 2.0;
		fiK_f(2,3) *= 2.0;
		fiK_f(2,4) *= 2.0;
		fiK_f(2,5) *= 2.0;
		fiK_f(3,3) *= 2.0;
		fiK_f(3,4) *= 2.0;
		fiK_f(3,5) *= 2.0;
		fiK_f(4,3) *= 2.0;
		fiK_f(4,4) *= 2.0;
		fiK_f(4,5) *= 2.0;
		fiK_f(5,3) *= 2.0;
		fiK_f(5,4) *= 2.0;
		fiK_f(5,5) *= 2.0;

//		cout << "\nK: "<<fiK_f;
		/*K^-1.(-R)*/
		fiK_f.Inverse();
		
		/*calculate update*/
		/*Cv_n+1 = Cvn + Delta Cvn*/
		FiberStretch_v[0] -= (fiK_f(0,0)*fResidual[0] + fiK_f(0,1)*fResidual[1] + fiK_f(0,2)*fResidual[2] 
			+ fiK_f(0,3)*fResidual[3] + fiK_f(0,4)*fResidual[4] + fiK_f(0,5)*fResidual[5]);
		FiberStretch_v[1] -= (fiK_f(1,0)*fResidual[0] + fiK_f(1,1)*fResidual[1] + fiK_f(1,2)*fResidual[2] 
			+ fiK_f(1,3)*fResidual[3] + fiK_f(1,4)*fResidual[4] + fiK_f(1,5)*fResidual[5]);
		FiberStretch_v[2] -= (fiK_f(2,0)*fResidual[0] + fiK_f(2,1)*fResidual[1] + fiK_f(2,2)*fResidual[2] 
			+ fiK_f(2,3)*fResidual[3] + fiK_f(2,4)*fResidual[4] + fiK_f(2,5)*fResidual[5]);
		FiberStretch_v[3] -= (fiK_f(3,0)*fResidual[0] + fiK_f(3,1)*fResidual[1] + fiK_f(3,2)*fResidual[2] 
			+ fiK_f(3,3)*fResidual[3] + fiK_f(3,4)*fResidual[4] + fiK_f(3,5)*fResidual[5]);
		FiberStretch_v[4] -= (fiK_f(4,0)*fResidual[0] + fiK_f(4,1)*fResidual[1] + fiK_f(4,2)*fResidual[2] 
			+ fiK_f(4,3)*fResidual[3] + fiK_f(4,4)*fResidual[4] + fiK_f(4,5)*fResidual[5]);
		FiberStretch_v[5] -= (fiK_f(5,0)*fResidual[0] + fiK_f(5,1)*fResidual[1] + fiK_f(5,2)*fResidual[2] 
			+ fiK_f(5,3)*fResidual[3] + fiK_f(5,4)*fResidual[4] + fiK_f(5,5)*fResidual[5]);
				
/*
		cout << "\niteration: "<<iteration;
		cout << "\nC: "<<FiberStretch;
		cout << "\nCv: "<<FiberStretch_v;
		cout << "\nSig: "<<fFlowStress;
		cout << "\nfResidual: "<<fResidual;
		cout << "\nDSig/DCv: "<<fModMat;
		cout << "\nfK: "<<fiK_f;
*/

		fResidual[0] = FiberStretch_v[0] - FiberStretch_vn[0];
		fResidual[1] = FiberStretch_v[1] - FiberStretch_vn[1];
		fResidual[2] = FiberStretch_v[2] - FiberStretch_vn[2];
		fResidual[3] = FiberStretch_v[3] - FiberStretch_vn[3];
		fResidual[4] = FiberStretch_v[4] - FiberStretch_vn[4];
		fResidual[5] = FiberStretch_v[5] - FiberStretch_vn[5];
		
		/*Compute flow stress*/		
		ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);
	
		/* FlowStress : (Cv.M1 + M1.Cv) */
		sig1 = fFlowStress[0]*2.0*FiberStretch_v[0] + 2.0*fFlowStress[4]*FiberStretch_v[4] 
			+ 2.0*fFlowStress[5]*FiberStretch_v[5];
		eta1 = fVisc1_f[pindex]->Function(sig1);
		coeff1 = dt*sig1/(2.0*eta1);

		/* FlowStress : (Cv.M2 + M2.Cv) */
		sig2 = fFlowStress[1]*2.0*FiberStretch_v[1] + 2.0*fFlowStress[3]*FiberStretch_v[3] 
			+ 2.0*fFlowStress[5]*FiberStretch_v[5];
		eta2 = fVisc2_f[pindex]->Function(sig2);
		coeff2 = dt*sig2/(2.0*eta2);

		/*fiber family 1*/
		/* coeff* (Cv.M1+M1.Cv)  */
		fResidual[0] -= coeff1 * 2.0*FiberStretch_v[0];
		fResidual[4] -= coeff1 * FiberStretch_v[4];
		fResidual[5] -= coeff1 * FiberStretch_v[5];
		
		/* coeff* (Cv.M1+M1.Cv)  */
		fResidual[1] -= coeff2 * 2.0*FiberStretch_v[1];
		fResidual[3] -= coeff2 * FiberStretch_v[3];
		fResidual[5] -= coeff2 * FiberStretch_v[5];
		
		error = sqrt(fResidual[0]*fResidual[0] + fResidual[1]*fResidual[1] + fResidual[2]*fResidual[2] 
			+ 2.0*fResidual[3]*fResidual[3] + 2.0*fResidual[4]*fResidual[4] + 2.0*fResidual[5]*fResidual[5]);
		iteration++;
		
/*		cout <<"\niteration: "<<iteration;
		cout<< "\nerror: "<<error;
*/
	}
	if (iteration >= 10) 
		ExceptionT::GeneralFail("NLV_Ortho::Compute_Cv", 
			"number of iteration exceeds maximum of 10");
//	cout << "\nC: "<<FiberStretch;
//	cout << "\nCv: "<<FiberStretch_v;
}


/*computes dFlowStress/dC.  Note for this model, dFlowStress/dC = - dSNEQ/dCv*/
void NLV_Ortho::dFlowdC (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dMatrixT& FiberMod,  const int pindex)
{	 
	/*calculate elastic structural invariants*/
	double& I4 = fI[0] = FiberStretch[0]/FiberStretch_v[0]; /*Cf:M1/Cfv:M1*/
	double& I6 = fI[2] = FiberStretch[1]/FiberStretch_v[1]; /*Cf:M2/Cfv:M2*/
		
	fInverse.Inverse(FiberStretch_v);
	double A1 = FiberStretch[0]*fInverse[0] + FiberStretch[4]*fInverse[4] + FiberStretch[5]*fInverse[5];
	double A2 = FiberStretch[0]*fInverse[4] + FiberStretch[4]*fInverse[2] + FiberStretch[5]*fInverse[3];
	double A3 = FiberStretch[0]*fInverse[5] + FiberStretch[4]*fInverse[3] + FiberStretch[5]*fInverse[1];
	double& I5 = fI[1] = (FiberStretch[0]*A1 + FiberStretch[4]*A2 + FiberStretch[5]*A3)/FiberStretch_v[0]; /* (C.Cv^-1.C):M1/(Cv:M1) */

	double B1 = FiberStretch[1]*fInverse[1] + FiberStretch[3]*fInverse[3] + FiberStretch[5]*fInverse[5];
	double B2 = FiberStretch[1]*fInverse[3] + FiberStretch[3]*fInverse[2] + FiberStretch[5]*fInverse[4];
	double B3 = FiberStretch[1]*fInverse[5] + FiberStretch[3]*fInverse[4] + FiberStretch[5]*fInverse[0];
	double& I7 = fI[3] = (FiberStretch[1]*B1+ FiberStretch[3]*B2 + FiberStretch[5]*B3)/FiberStretch_v[1]; /* (C.Cv^-1.C):M2/(Cv:M2)*/

/*	if (I4 < 1.0) I4 = 1.0;	
	if (I6 < 1.0) I6 = 1.0;
	if (I5 < 1.0) I5 = 1.0;
	if (I7 < 1.0) I7 = 1.0;
*/		
	/*calc 2dW/dI and 2ddW/ddI*/
	ddWddI(fI, fdW, fddW, pindex);
		
	const double& s4 = fdW[0];
	const double& s5 = fdW[1];
	const double& s6 = fdW[2];
	const double& s7 = fdW[3];
	
	const double& d4 = fddW[0];
	const double& d5 = fddW[1];
	const double& d6 = fddW[2];
	const double& d7 = fddW[3];
	
	double Iv4_2 = FiberStretch_v[0]*FiberStretch_v[0];
	double Iv6_2 = FiberStretch_v[1]*FiberStretch_v[1];
	
	FiberMod = 0.0;

	 /*Mi x Mi*/
	FiberMod(0,0) += (s4 + d4*I4)/Iv4_2;
	FiberMod(1,1) += (s6 + d6*I6)/Iv6_2;
	
	/*Mi x (Mi.C.Cv^-1 + Cv^-1.C.Mi)*/
	double c5 = (s5 + d5*I5)/Iv4_2;
	double c7 = (s7 + d7*I7)/Iv6_2;
	
	FiberMod(0,0) += c5 * 2.0*A1;
	FiberMod(0,4) += c5 * A2;
	FiberMod(0,5) += c5 * A3;

	FiberMod(1,1) += c7 * 2.0*B1;
	FiberMod(1,3) += c7 * B2;
	FiberMod(1,5) += c7 * B3;
	
	/* Cv^-1.C.Mi o Cv^-1 + Cv^-1 o Cv^-1.C.Mi*/
	c5 = s5/FiberStretch_v[0];
	c7 = s7/FiberStretch_v[1];

	FiberMod(0,0) += 2.0*c5*A1*fInverse[0]; 
	FiberMod(0,1) += 2.0*c7*B3*fInverse[5]; 
	FiberMod(0,3) += c7*B3*fInverse[4];
	FiberMod(0,4) += c5*A1*fInverse[4];
	FiberMod(0,5) += c7*B3*fInverse[0] + 2.0*c5*A1*fInverse[5];
	
	FiberMod(1,0) += 2.0*c5*A3*fInverse[5]; 
	FiberMod(1,1) += 2.0*c7*B1*fInverse[1]; 
	FiberMod(1,3) += c7*B1*fInverse[3];
	FiberMod(1,4) += c5*A3*fInverse[3];
	FiberMod(1,5) += c5*A3*fInverse[1] + 2.0*c7*B1*fInverse[5];

	FiberMod(2,0) += 2.0*c5*A2*fInverse[4]; 
	FiberMod(2,1) += 2.0*c7*B2*fInverse[3]; 
	FiberMod(2,3) += c7*B2*fInverse[2];
	FiberMod(2,4) += c5*A2*fInverse[2];
	FiberMod(2,5) += c5*A2*fInverse[3] + 2.0*c7*B2*fInverse[4];
	
	FiberMod(3,0) += c5*A3*fInverse[4] + 2.0*c5*A2*fInverse[5]; 
	FiberMod(3,1) += c7*B2*fInverse[1] + 2.0*c7*B1*fInverse[3]; 
	FiberMod(3,3) += 0.5*c7*(B1*fInverse[2] + B2*fInverse[3]);
	FiberMod(3,4) += 0.5*c5*(A3*fInverse[2] + A2*fInverse[3]);
	FiberMod(3,5) += 0.5*c5*(A2*fInverse[1] + A3*fInverse[3]) + 0.5*c7*(B1*fInverse[4] + B2*fInverse[5]);

	FiberMod(4,0) += c5*(A2*fInverse[0] + A1*fInverse[4]); 
	FiberMod(4,1) += c7*(B3*fInverse[3] + B2*fInverse[5]); 
	FiberMod(4,3) += 0.5*c7*(B3*fInverse[2] + B2*fInverse[4]);
	FiberMod(4,4) += 0.5*c5*(A1*fInverse[2] + A2*fInverse[4]);
	FiberMod(4,5) += 0.5*c7*(B2*fInverse[0] + B3*fInverse[4]) +c5*(A1*fInverse[3] + A2*fInverse[5]);

	FiberMod(5,0) += c5*(A3*fInverse[0] + A1*fInverse[5]); 
	FiberMod(5,1) += c7*(B3*fInverse[1] + B1*fInverse[5]); 
	FiberMod(5,3) += 0.5*c7*(B3*fInverse[3] + B1*fInverse[4]);
	FiberMod(5,4) += 0.5*c5*(A1*fInverse[3] + A3*fInverse[4]);
	FiberMod(5,5) += 0.5*c5*(A1*fInverse[1] + A3*fInverse[5]) + 0.5*c7*(B1*fInverse[0] + B3*fInverse[5]);

	/*Cv^-1.C.Mi.C.Cv^-1 x (Mi.C.Cv^-1 + Cv^-1.C.Mi) */
	c5 = d5/Iv4_2;
	c7 = d7/Iv6_2;

	/*(M1.C.Cv^-1 + Cv^-1.C.M1)*/
	fMod1[0] = 2.0*A1; 
	fMod1[1] = 0.0;
	fMod1[2] = 0.0;
	fMod1[3] = 0.0;
	fMod1[4] = A2;
	fMod1[5] = A3;

	/*Cv^-1.C.M1.C.Cv^-1*/
 	fMod2[0] = A1*A1; 
	fMod2[1] = A3*A3;
	fMod2[2] = A2*A2;
	fMod2[3] = A3*A2;
	fMod2[4] = A1*A2;
	fMod2[5] = A1*A3;

	fMod3.DyadAB(fMod2,fMod1);
	FiberMod.AddScaled(c5, fMod3);

	/*(M2.C.Cv^-1 + Cv^-1.C.M2)*/
 	fMod1[0] = 0.0; 
	fMod1[1] = 2.0*B1;
	fMod1[2] = 0.0;
	fMod1[3] = B2;
	fMod1[4] = 0.0;
	fMod1[5] = B3;

	/*Cv^-1.C.M2.C.Cv^-1*/
 	fMod2[0] = B3*B3; 
	fMod2[1] = B1*B1;
	fMod2[2] = B2*B2;
	fMod2[3] = B1*B2;
	fMod2[4] = B3*B2;
	fMod2[5] = B3*B1;

	fMod3.DyadAB(fMod2,fMod1);
	FiberMod.AddScaled(c7, fMod3);
}

/*computes dFlowStress/dCv*/
void NLV_Ortho::dFlowdCv (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
		dMatrixT& FiberMod,  const int pindex)
{
	/*calculate elastic structural invariants*/
	double& I4 = fI[0] = FiberStretch[0]/FiberStretch_v[0]; /*Cf:M1/Cfv:M1*/
	double& I6 = fI[2] = FiberStretch[1]/FiberStretch_v[1]; /*Cf:M2/Cfv:M2*/
		
	fInverse.Inverse(FiberStretch_v);
	double A1 = FiberStretch[0]*fInverse[0] + FiberStretch[4]*fInverse[4] + FiberStretch[5]*fInverse[5];
	double A2 = FiberStretch[0]*fInverse[4] + FiberStretch[4]*fInverse[2] + FiberStretch[5]*fInverse[3];
	double A3 = FiberStretch[0]*fInverse[5] + FiberStretch[4]*fInverse[3] + FiberStretch[5]*fInverse[1];
	double& I5 = fI[1] = (FiberStretch[0]*A1 + FiberStretch[4]*A2 + FiberStretch[5]*A3)/FiberStretch_v[0]; /* (C.Cv^-1.C):M1/(Cv:M1) */

	double B1 = FiberStretch[1]*fInverse[1] + FiberStretch[3]*fInverse[3] + FiberStretch[5]*fInverse[5];
	double B2 = FiberStretch[1]*fInverse[3] + FiberStretch[3]*fInverse[2] + FiberStretch[5]*fInverse[4];
	double B3 = FiberStretch[1]*fInverse[5] + FiberStretch[3]*fInverse[4] + FiberStretch[5]*fInverse[0];
	double& I7 = fI[3] = (FiberStretch[1]*B1+ FiberStretch[3]*B2 + FiberStretch[5]*B3)/FiberStretch_v[1]; /* (C.Cv^-1.C):M2/(Cv:M2)*/

/*	if (I4 < 1.0) I4 = 1.0;	
	if (I6 < 1.0) I6 = 1.0;
	if (I5 < 1.0) I5 = 1.0;
	if (I7 < 1.0) I7 = 1.0;
*/		
	/*calc 2dW/dI and 2ddW/ddI*/
	ddWddI(fI, fdW, fddW, pindex);
		
	const double& s4 = fdW[0];
	const double& s5 = fdW[1];
	const double& s6 = fdW[2];
	const double& s7 = fdW[3];
	
	const double& d4 = fddW[0];
	const double& d5 = fddW[1];
	const double& d6 = fddW[2];
	const double& d7 = fddW[3];
	
	double Iv4_2 = FiberStretch_v[0]*FiberStretch_v[0];
	double Iv6_2 = FiberStretch_v[1]*FiberStretch_v[1];
	
	/*Cv^-1.C.M1.C.Cv^-1*/
 	fMod1[0] = A1*A1; 
	fMod1[1] = A3*A3;
	fMod1[2] = A2*A2;
	fMod1[3] = A3*A2;
	fMod1[4] = A1*A2;
	fMod1[5] = A1*A3;

	/*Cv^-1.C.M2.C.Cv^-1*/
 	fMod2[0] = B3*B3; 
	fMod2[1] = B1*B1;
	fMod2[2] = B2*B2;
	fMod2[3] = B1*B2;
	fMod2[4] = B3*B2;
	fMod2[5] = B3*B1;

	FiberMod = 0.0;

	 /*Mi x Mi*/	
	FiberMod(0,0) -= (2.0*s4*I4 + d4*I4*I4 + 2.0*s5*I5 + d5*I5*I5)/Iv4_2;
	FiberMod(1,1) -= (2.0*s6*I6 + d6*I6*I6 + 2.0*s7*I7 + d7*I7*I7)/Iv6_2;

	/*Mi x Cv^-1.C.Mi.C.Cv^-1*/
	double c5 = (s5 + d5*I5)/Iv4_2;
	double c7 = (s7 + d7*I7)/Iv6_2;
	
	FiberMod(0,0) -= c5 * fMod1[0];
	FiberMod(0,1) -= c5 * fMod1[1];
	FiberMod(0,2) -= c5 * fMod1[2];
	FiberMod(0,3) -= c5 * fMod1[3];
	FiberMod(0,4) -= c5 * fMod1[4];
	FiberMod(0,5) -= c5 * fMod1[5];

	FiberMod(1,0) -= c7 * fMod2[0];
	FiberMod(1,1) -= c7 * fMod2[1];
	FiberMod(1,2) -= c7 * fMod2[2];
	FiberMod(1,3) -= c7 * fMod2[3];
	FiberMod(1,4) -= c7 * fMod2[4];
	FiberMod(1,5) -= c7 * fMod2[5];
	
	/*Cv^-1.C.Mi.C.Cv^-1 x Mi */
	FiberMod(0,0) -= c5 * fMod1[0];
	FiberMod(1,0) -= c5 * fMod1[1];
	FiberMod(2,0) -= c5 * fMod1[2];
	FiberMod(3,0) -= c5 * fMod1[3];
	FiberMod(4,0) -= c5 * fMod1[4];
	FiberMod(5,0) -= c5 * fMod1[5];

	FiberMod(0,1) -= c7 * fMod2[0];
	FiberMod(1,1) -= c7 * fMod2[1];
	FiberMod(2,1) -= c7 * fMod2[2];
	FiberMod(3,1) -= c7 * fMod2[3];
	FiberMod(4,1) -= c7 * fMod2[4];
	FiberMod(5,1) -= c7 * fMod2[5];
		
	/*Cv^-1.C.Mi.C.Cv^-1 x Cv^-1.C.Mi.C.Cv^-1 */
	c5 = d5/Iv4_2;
	c7 = d7/Iv6_2;

	fMod3.DyadAB(fMod1,fMod1);
	FiberMod.AddScaled(-c5, fMod3);
	fMod3.DyadAB(fMod2,fMod2);
	FiberMod.AddScaled(-c7, fMod3);
	
	/*Cv^-1.C.Mi.C.Cv^-1 o Cv^-1 + Cv^-1 o Cv^-1.C.Mi.C.Cv^-1*/
	c5 = s5/FiberStretch_v[0];
	c7 = s7/FiberStretch_v[1];

	fMod3.ReducedI_AB(fInverse,fMod1);
	FiberMod.AddScaled(-c5, fMod3);	

	fMod3.ReducedI_AB(fMod1,fInverse);
	FiberMod.AddScaled(-c5, fMod3);

	fMod3.ReducedI_AB(fInverse,fMod2);
	FiberMod.AddScaled(-c7, fMod3);
	fMod3.ReducedI_AB(fMod2,fInverse);
	FiberMod.AddScaled(-c7, fMod3);	
}

void NLV_Ortho::dWdI(const dArrayT& I, dArrayT& dW, const int pindex)
{
	double I4 = I[0];
	double I5 = I[1];
	double I6 = I[2];
	double I7 = I[3];


	dW[0] = fPot1_f4[pindex+1]->DFunction(I4);
	dW[2] = fPot2_f4[pindex+1]->DFunction(I6);

	dW[1] = fPot1_f5[pindex+1]->DFunction(I5);
	dW[3] = fPot2_f5[pindex+1]->DFunction(I7);
}

void NLV_Ortho::ddWddI(const dArrayT& I, dArrayT& dW, dArrayT& ddW, const int pindex)
{
	double I4 = I[0];
	double I5 = I[1];
	double I6 = I[2];
	double I7 = I[3];

	dW[0] = fPot1_f4[pindex+1]->DFunction(I4);
	dW[2] = fPot2_f4[pindex+1]->DFunction(I6);
	dW[1] = fPot1_f5[pindex+1]->DFunction(I5);
	dW[3] = fPot2_f5[pindex+1]->DFunction(I7);

	ddW[0] = fPot1_f4[pindex+1]->DDFunction(I4);
	ddW[2] = fPot2_f4[pindex+1]->DDFunction(I6);
	ddW[1] = fPot1_f5[pindex+1]->DDFunction(I5);
	ddW[3] = fPot2_f5[pindex+1]->DDFunction(I7);
}
