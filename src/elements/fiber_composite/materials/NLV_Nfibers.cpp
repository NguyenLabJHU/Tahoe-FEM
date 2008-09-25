/* $Id: NLV_Nfibers.cpp,v 1.1 2008-09-25 13:32:32 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "NLV_Nfibers.h"
#include <math.h>
#include "toolboxConstants.h"
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
NLV_Nfibers::NLV_Nfibers(void):
  FSFiberMatViscT(),
  ParameterInterfaceT("nonlinear_viscoelasticity_Nfibers")
{
	/*reset default*/
	fNumFibProcess = 0;
	fNumMatProcess = 0;
}

/* destructor */
NLV_Nfibers::~NLV_Nfibers(void) 
{ 
	/*allocated?*/
	double l1 = fPot_f.MinorDim();
	double l2 = fVisc_f.MinorDim();
	double nfibs = fPot_f.MajorDim();
	
	for (int j = 0; j<nfibs; j++)
	{
		for (int i = 0; i < l1; i++){
				delete fPot_f(j,i);
		}
		for (int i = 0; i < l2; i++){
				delete fVisc_f(j,i);
		}
	}
}


/* free energy density */
double NLV_Nfibers::StrainEnergyDensity(void)
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
		
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	int numfibers = Fibers.MajorDim();
	numfibers--;
		
	/*Get  orthogonal coordinate basis for fiber plane*/
	const dMatrixT& Q = GetRotation();
	const double* p1 = Q(0);
	const double* p2 = Q(1);		

	/*fiber contribution*/
	/*equilibrium contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	int fibnum = 0;
	for (int i = 0; i<numfibers; i++)
	{
		if(!fsame)
			fibnum = i;
			
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
	
		/*calculate I4*/
		double I4, Iv4, Ie4;
		I4 = fFiberStretch[0]*cost*cost + fFiberStretch[1]*sint*sint + 2.0*fFiberStretch[5]*sint*cost;	
		energy += fPot_f(fibnum,0)->Function(I4);
	}


	/*non-equilibrium*/
	int j = fNumMatProcess;
	for (int i = 0; i < fNumFibProcess; i++)
	{
		ComputeFiberStretch(fC_v[j], fFiberStretch_v);
		int fibnum = 0;
		for (int k = 0; k<numfibers; k++)
		{
			double cost = Fibers(k,0)*p1[0] + Fibers(k,1)*p1[1] + Fibers(k,2)*p1[2];
			double sint = Fibers(k,0)*p2[0] + Fibers(k,1)*p2[1] + Fibers(k,2)*p2[2];
		
			/*calculate I4*/
			double I4, Iv4, Ie4;
			I4 = fFiberStretch[0]*cost*cost + fFiberStretch[1]*sint*sint + 2.0*fFiberStretch[5]*sint*cost;	
			Iv4 = fFiberStretch_v[0]*cost*cost + fFiberStretch_v[1]*sint*sint + 2.0*fFiberStretch_v[5]*sint*cost;
			Ie4 = I4/Iv4;

			energy += fPot_f(fibnum,i+1)->Function(I4);
		}
		j++;
	}
	return energy;
}

int NLV_Nfibers::NumOutputVariables() const {
	return kNumOutputVar;
}

void NLV_Nfibers::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void NLV_Nfibers::ComputeOutput(dArrayT& output)
{
	/*calculates deformed fiber vectors*/
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();

	const double* p_nt = Fibers(0);
	const double* p_is = Fibers(1);
	
	const dMatrixT& F = F_mechanical();
	double* pb = output.Pointer();
	
	/*deformed P1 fiber orientation*/
	F.Multx(p_nt, pb);
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
		const dSymMatrixT& Cvf = fC_v[fNumMatProcess];
		*pb++ = Cvf[0]*p_nt[0]*p_nt[0] + Cvf[1]*p_nt[1]*p_nt[1]+Cvf[2]*p_nt[2]*p_nt[2]
				+2.0*Cvf[5]*p_nt[0]*p_nt[1] + 2.0*Cvf[4]*p_nt[0]*p_nt[2] + 2.0*Cvf[2]*p_nt[1]*p_nt[2];
		*pb++ = Cvf[0]*p_is[0]*p_is[0] + Cvf[1]*p_is[1]*p_is[1]+Cvf[2]*p_is[2]*p_is[2]
				+2.0*Cvf[5]*p_is[0]*p_is[1] + 2.0*Cvf[4]*p_is[0]*p_is[2] + 2.0*Cvf[2]*p_is[1]*p_is[2];
	}
	else {
		*pb++ = 1.0;
		*pb++ = 1.0;
	}
	output[9] = F.Det();
//	if (CurrElementNumber() == 0 && CurrIP() == 0)
//		cout <<F;
}


/* information about subordinate parameter lists */
void NLV_Nfibers::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSFiberMatViscT::DefineSubs(sub_list);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("fibers_params", ParameterListT::OnePlus);
}


/* return the description of the given inline subordinate parameter list */
void NLV_Nfibers::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if(name == "eq_fiber_pot_choice")
	{	
		order = ParameterListT::Choice;
		sub_lists.AddSub("eq_fiber_pot");	
	}
	else if(name == "neq_fiber_pot_choice")
	{		
		order = ParameterListT::Choice;
		sub_lists.AddSub("neq_fiber_pot");	
	}
	else if(name=="fiber_visc_potential_choice")
	{
		order = ParameterListT::Choice;
		sub_lists.AddSub("fiber_visc_potential");	
	}		
	else /* inherited */
		FSFiberMatViscT::DefineInlineSub(name, order, sub_lists);
}
/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* NLV_Nfibers::NewSub(const StringT& name) const
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

	/* inherited */
	if (name == "fibers_params")
	{
		ParameterContainerT* fiber = new ParameterContainerT(name);
//		fiber->SetListOrder(ParameterListT::Sequence);
		fiber->SetSubSource(this);
		fiber->AddSub("eq_fiber_pot_choice", ParameterListT::Once,true);
		fiber->AddSub("neq_fiber_pot_choice", ParameterListT::Any,true);
		fiber->AddSub("fiber_visc_potential_choice", ParameterListT::Any,true);

		return(fiber);		
	}
	else if(name == "eq_fiber_pot")
	{	
		ParameterContainerT* eq_choice = new ParameterContainerT(name);
		eq_choice->SetSubSource(this);
		eq_choice->SetListOrder(ParameterListT::Choice);
		eq_choice->SetSubSource(this);
		
		eq_choice->AddSub("parabola");
		eq_choice->AddSub("fung-type");
		eq_choice->AddSub("fung_type_square");
		
		return(eq_choice);
	}
	else if(name == "neq_fiber_pot")
	{		
		ParameterContainerT* neq_choice = new ParameterContainerT(name);
		neq_choice->SetSubSource(this);
		neq_choice->SetListOrder(ParameterListT::Choice);
		neq_choice->SetSubSource(this);
		
		neq_choice->AddSub("parabola");
		neq_choice->AddSub("fung-type");
		neq_choice->AddSub("fung_type_square");
		return(neq_choice);
	}
	else if(name=="fiber_visc_potential")
	{
		ParameterContainerT* visc_choice = new ParameterContainerT(name);
		visc_choice->SetSubSource(this);
		visc_choice->SetListOrder(ParameterListT::Choice);
		visc_choice->SetSubSource(this);
		
		visc_choice->AddSub("scaled-csch");
		visc_choice->AddSub("linear_exponential");
		return(visc_choice);
	}		
	else return(FSFiberMatViscT::NewSub(name));


}
void NLV_Nfibers::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSFiberMatViscT::DefineParameters(list);

	ParameterT same(fsame, "use_for_all_fibers");
	same.SetDefault(false);
	list.AddParameter(same);
}


/* accept parameter list */
void NLV_Nfibers::TakeParameterList(const ParameterListT& list)
{

	StringT caller = "NLV_Nfibers::TakeParameterList";
	
	
	/* inherited */
	FSFiberMatViscT::TakeParameterList(list);
	fNumFibStress = dSymMatrixT::NumValues(fNumSD);
	

	int num_fibers = list.NumLists("fibers_params");
	fsame = false;
	fsame = list.GetParameter("use_for_all_fibers");
	if (fsame & num_fibers > 1)
			ExceptionT::GeneralFail("NLV_Nfibers::TakeParameterList", 
				"more than one list of fiber parameters given with option use_for_all_fibers");
	
	for (int i = 0; i< num_fibers; i++)
	{
		const ParameterListT& fiberlist = list.GetList("fibers_params",i);
		int num_fib_neq =  fiberlist.NumLists("neq_fiber_pot");
		int num_fib_visc =  fiberlist.NumLists("fiber_visc_potential");
		
		if (num_fib_neq != num_fib_visc)
			ExceptionT::GeneralFail("NLV_Nfibers::TakeParameterList", 
				"number of fiber viscosity functions does not match number of fiber nonequilibrium potentials");
		fNumFibProcess = num_fib_neq;
		fVisc_f.Dimension(num_fibers, fNumFibProcess);
		fPot_f.Dimension(num_fibers, fNumFibProcess+1);	

		const ParameterListT& fiber_pot = fiberlist.GetListChoice(*this, "eq_fiber_pot");
		if (fiber_pot.Name() == "parabola")
			fPot_f(i,0) = new ParabolaT;
		else if (fiber_pot.Name() == "fung-type")
			fPot_f(i,0) = new FungType;
		else if (fiber_pot.Name() == "fung_type_square")
			fPot_f(i,0) = new FungType2;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot_f(i,0)) throw ExceptionT::kOutOfMemory;
		fPot_f(i,0)->TakeParameterList(fiber_pot);


		for (int j = 0; j < fNumFibProcess; j++)
		{
			const ParameterListT& fiber_neq = fiberlist.GetListChoice(*this, "neq_fiber_pot",j);
			if (fiber_neq.Name() == "parabola")
				fPot_f(i,j+1) = new ParabolaT;
			else if (fiber_neq.Name() == "fung-type")
				fPot_f(i,j+1) = new FungType;
			else if (fiber_neq.Name() == "fung_type_square")
				fPot_f(i,j+1) = new FungType2;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fPot_f(j,i+1)) throw ExceptionT::kOutOfMemory;
			fPot_f(i,j+1)->TakeParameterList(fiber_neq);

			const ParameterListT& fiber_visc = fiberlist.GetListChoice(*this, "fiber_visc_potential", j);
			if (fiber_visc.Name() == "linear_exponential")
				fVisc_f(j,i) = new LinearExponentialT;
			else if (fiber_visc.Name() == "scaled-csch")
				fVisc_f(j,i) = new ScaledCsch;
			else 
				ExceptionT::GeneralFail(caller, "no such potential");
			if (!fVisc_f(j,i)) throw ExceptionT::kOutOfMemory;
			fVisc_f(j,i)->TakeParameterList(fiber_visc);
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
	
}
	
/***********************************************************************
 * Protected
 ***********************************************************************/

/*computes integrated fiber stress in local frame*/
void NLV_Nfibers::ComputeFiberStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
			dSymMatrixT& FiberStress, const int pindex)
{
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	int numfibers = Fibers.MajorDim();
	numfibers--;
		
	FiberStress = 0.0;

	/*Get  orthogonal coordinate basis for fiber plane*/
	const dMatrixT& Q = GetRotation();
	const double* p1 = Q(0);
	const double* p2 = Q(1);

	int fibnum = 0;
	for (int i = 0; i<numfibers; i++)
	{
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
		
		if (!fsame)
			fibnum = i;
			
		if (pindex < 0)  /*eq*/
		{ 		
			/*calculate I4*/
			double I4, s4;
			I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		
			/*calc 2dW/dI*/
			s4 = fPot_f(fibnum,pindex+1)->DFunction(I4);
			FiberStress[0] += cost*cost*s4;
			FiberStress[1] += sint*sint*s4;
			FiberStress[5] += sint*cost*s4;				
		}
		else 
		{
			/*calculate elastic structural invariants*/
			double I4, Iv4, Ie4, s4;
			I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
			Iv4 = FiberStretch_v[0]*cost*cost + FiberStretch_v[1]*sint*sint + 2.0*FiberStretch_v[5]*sint*cost;
			Ie4 = I4/Iv4;
			/*calc 2dWeq/dI*/
			s4 = fPot_f(fibnum,pindex+1)->DFunction(Ie4);
		
			FiberStress[0] += cost*cost*s4/Iv4;
			FiberStress[1] += sint*sint*s4/Iv4;
			FiberStress[5] += sint*cost*s4/Iv4;				
		}
		
	}
}
	
/*computes integrated moduli in local frame*/
void NLV_Nfibers::ComputeFiberMod (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dSymMatrixT& FiberStress, dMatrixT& FiberMod,  const int pindex)
{

	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	int numfibers = Fibers.MajorDim();
	numfibers--;

	FiberStress = 0.0;
	FiberMod = 0.0;

	/*Get  orthogonal coordinate basis for fiber plane*/
	const dMatrixT& Q = GetRotation();
	const double* p1 = Q(0);
	const double* p2 = Q(1);

	int fibnum=0;
	for (int i = 0; i<numfibers; i++)
	{
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
			
		if (pindex < 0)  /*eq*/
		{ 		
			/*calculate I4*/
			double I4, s4, d4;
			I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		
			if (!fsame)
				fibnum = i;
				
			/*calc 2dW/dI*/
			s4 = fPot_f(fibnum,pindex+1)->DFunction(I4); 
			d4= fPot_f(fibnum,pindex+1)->DDFunction(I4);
		
			FiberStress[0] += cost*cost*s4;
			FiberStress[1] += sint*sint*s4;
			FiberStress[5] += sint*cost*s4;				
		
			/* Mi x Mi*/
			FiberMod(0,0) += cost*cost*cost*cost*2.0*d4;
			FiberMod(1,1) += sint*sint*sint*sint*2.0*d4;
			FiberMod(5,5) += sint*sint*cost*cost*2.0*d4;
		
			FiberMod(0,1) += sint*sint*cost*cost*2.0*d4;
			FiberMod(1,0) += sint*sint*cost*cost*2.0*d4;
		
			FiberMod(0,5) += sint*cost*cost*cost*2.0*d4;
			FiberMod(5,0) += sint*cost*cost*cost*2.0*d4;
		
			FiberMod(1,5) += sint*sint*sint*cost*2.0*d4;
			FiberMod(5,1) += sint*sint*sint*cost*2.0*d4;

		}
		else /*neq*/
		{
			ComputeFiberStress (FiberStretch, FiberStretch_v, FiberStress,  pindex);		
			ComputeCalg(FiberStretch, FiberStretch_v, FiberMod, pindex);
		}
/*		cout << "\nFiber: "<<i<<"\t"<<pindex;
		cout << "\nFiberStress: "<<FiberStress;
		cout << "\nFiberMod: "<<FiberMod;
*/
	}
}

void NLV_Nfibers::ComputeCalg(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,  dMatrixT& Calg, const int pindex)
{		
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
		
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	int numfibers = Fibers.MajorDim();
	numfibers--;
	
	/*Get  orthogonal coordinate basis for fiber plane*/
	const dMatrixT& Q = GetRotation();
	const double* p1 = Q(0);
	const double* p2 = Q(1);

	Calg = 0.0;
	int fibnum=0;
	for (int i = 0; i<numfibers; i++)
	{
		if (!fsame)
			fibnum = i;
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];

		/*calculate elastic structural invariants*/
		double I4, Iv4, Ie4;
		I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		Iv4 = FiberStretch_v[0]*cost*cost + FiberStretch_v[1]*sint*sint + 2.0*FiberStretch_v[5]*sint*cost;
		Ie4 = I4/Iv4;
		double lv1 = sqrt(Iv4);
		double s4, d4;
		s4 = fPot_f(fibnum,pindex+1)->DFunction(Ie4); d4= fPot_f(fibnum,pindex+1)->DDFunction(Ie4);		

		double sig1 = s4*Ie4;

		double eta1 = fVisc_f(fibnum,pindex)->Function(sig1);
	
		/*deta_i/dsig_i*/
		double deta1 = fVisc_f(fibnum,pindex)->DFunction(sig1);

		/*dsig_i/dlv_i*/
		double dsig1dlv1 = (-1.0/lv1)*(2.0*s4*Ie4 + 2.0*d4*Ie4*Ie4);

		double dsig1dl1 = (1.0/sqrt(I4))*(2.0*s4*Ie4 + 2.0*d4*Ie4*Ie4);

		/*stiffness*/
		double k1 = 1.0 - dt/eta1*((1.0 - (1.0/eta1)*deta1*sig1)*dsig1dlv1*lv1 + sig1);
	
		double g1 = - dt/eta1*((1.0 - (1.0/eta1)*deta1*sig1)*dsig1dl1*lv1);
		
		/*dSNEQ/dC + dSNEQ/dCv*/
		double c4 =  2.0*(s4*1.0/sqrt(Ie4)*g1/k1 + d4*(1.0+sqrt(Ie4)*g1/k1))/(Iv4*Iv4);  

		/* Mi x Mi*/
		Calg(0,0) += cost*cost*cost*cost*c4;
		Calg(1,1) += sint*sint*sint*sint*c4;
		Calg(5,5) += sint*sint*cost*cost*c4;
		
		Calg(0,1) += sint*sint*cost*cost*c4;
		Calg(1,0) += sint*sint*cost*cost*c4;
		
		Calg(0,5) += sint*cost*cost*cost*c4;
		Calg(5,0) += sint*cost*cost*cost*c4;
			
		Calg(1,5) += sint*sint*sint*cost*c4;
		Calg(5,1) += sint*sint*sint*cost*c4;
	}
}

/*local newton loop for viscous stretch tensor*/ 
void NLV_Nfibers::Compute_Cv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_vn, 
	dSymMatrixT& FiberStretch_v, const int pindex)
{	
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
				
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	int numfibers = Fibers.MajorDim();
	numfibers--;
	
	/*Get  orthogonal coordinate basis for fiber plane*/
	const dMatrixT& Q = GetRotation();
	const double* p1 = Q(0);
	const double* p2 = Q(1);

	double Cv11=0.0;
	double Cv22 =0.0;
	double Cv12 = 0.0;
	int fibnum = 0;
	for (int i = 0; i<numfibers; i++)
	{	
		if (!fsame)
			fibnum = i;
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];

		/*calculate elastic structural invariants*/
		double I4, Iv4, Iv4n, Ie4;
		I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		Iv4 = FiberStretch_v[0]*cost*cost + FiberStretch_v[1]*sint*sint + 2.0*FiberStretch_v[5]*sint*cost;
		double lv1 = sqrt(Iv4);
		Iv4n = FiberStretch_vn[0]*cost*cost + FiberStretch_vn[1]*sint*sint + 2.0*FiberStretch_vn[5]*sint*cost;
		double lv1n = sqrt(Iv4n);
		Ie4 = I4/Iv4;
		 
		double s4, d4;
		
		/*calc 2dW/dI and 2ddW/ddI*/
		s4 = fPot_f(fibnum,pindex+1)->DFunction(Ie4); d4= fPot_f(fibnum,pindex+1)->DDFunction(Ie4);		
	
		double sig1 = s4*Ie4;
		double eta1 = fVisc_f(fibnum,pindex)->Function(sig1);

		double r1 = lv1 - dt/eta1*sig1*lv1 - lv1n;

/*		cout <<"\ndt: "<<dt;
		cout << "\nI: "<<I4<<"\t"<<Iv4<<"\t"<<Ie4;
		cout << "\nsig1: "<<sig1;
		cout << "\nr1: "<<r1;
*/		
		double error = sqrt(r1*r1);
		int iteration = 0;
		
		while (error > kSmall && iteration < 10)
		{	
			/*deta_i/dsig_i*/
			double deta1 = fVisc_f(fibnum,pindex)->DFunction(sig1);

			/*dsig_i/dlv_i*/
			double dsig1 = (-1.0/lv1)*(2.0*s4*Ie4 + 2.0*d4*Ie4*Ie4);

			/*stiffness*/
			double k1 = 1.0 - dt/eta1*((1.0 - (1.0/eta1)*deta1*sig1)*dsig1*lv1 + sig1);

			/*calculate update to viscous fiber stretch*/
			lv1 -= r1/k1;
		
			Iv4 = lv1*lv1;
			Ie4 = I4/Iv4;
		
			/*calc 2dW/dI and 2ddW/ddI*/
			s4 = fPot_f(fibnum,pindex+1)->DFunction(Ie4); d4= fPot_f(fibnum,pindex+1)->DDFunction(Ie4);		
			sig1 = s4*Ie4;
								
			eta1 = fVisc_f(fibnum,pindex)->Function(sig1);
	
			r1 = lv1 - dt/eta1*sig1*lv1 - lv1n;
		
			error = sqrt(r1*r1);
			iteration++;
/*		7/2-7/12 <<"\niteration: "<<iteration;
		cout<< "\nerror: "<<error;
*/
		}
		if (iteration >= 10) 
			ExceptionT::GeneralFail("NLV_Nfibers::Compute_Cv", 
				"number of iteration exceeds maximum of 10");

		Cv11 += cost*cost*Iv4;
		Cv22 += sint*sint*Iv4;
		Cv12 += sint*cost*Iv4;
	}
	FiberStretch_v[0] = Cv11;
	FiberStretch_v[1] = Cv22;
	FiberStretch_v[5] = Cv12;				
}




