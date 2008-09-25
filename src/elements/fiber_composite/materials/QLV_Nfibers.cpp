/* $Id: QLV_Nfibers.cpp,v 1.1 2008-09-25 13:32:32 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "QLV_Nfibers.h"
#include <math.h>
#include "toolboxConstants.h"
#include "ParameterContainerT.h"

#include "ParabolaT.h"
#include "FungType.h"
#include "FungType2.h"

#include "MooneyRivlin.h"
#include "NeoHookean.h"
#include "VWPotentialT.h"

const double Pi = acos(-1.0);
static const double third = 1.0/3.0;
const int kNumOutputVar = 6;
static const char* Labels[kNumOutputVar] = {"p1_X", "p1_Y", "p1_Z","p2_X", "p2_Y", "p2_Z"};

using namespace Tahoe;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
QLV_Nfibers::QLV_Nfibers(void):
	fSpectralDecompSpat(3),
  ParameterInterfaceT("quasilinear_viscoelasticity_Nfibers")
{
	/*reset default*/
	fNumFibProcess = 0;
	fNumMatProcess = 0;
}

/* destructor */
QLV_Nfibers::~QLV_Nfibers(void) 
{ 
	/*allocated?*/
	double l1 = fPot_f.MajorDim();
	double k1 = fPot_f.MinorDim();
	for (int j = 0; j < l1; j++){
		for (int i = 0; i < k1; i++){
			delete fPot_f(j,i);
		}
	}
}

/* free energy density */
double QLV_Nfibers::StrainEnergyDensity(void)
{
	const dMatrixT& F = F_mechanical();
	fC.MultATA(F);
	fStress = S_IJ();
	double energy = 0.5*fStress.ScalarProduct(fC);
	return energy;
}

int QLV_Nfibers::NumOutputVariables() const {
	return kNumOutputVar;
}

void QLV_Nfibers::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void QLV_Nfibers::ComputeOutput(dArrayT& output)
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
}

/* modulus */
const dMatrixT& QLV_Nfibers::C_IJKL(void)
{
	int elem = CurrElementNumber();
	int ip = CurrIP();

	const dMatrixT& F = F_mechanical();  /*calculate inelastic deviatoric stress*/
	fb.MultAAT(F); 		/*b = FF^T*/
	fC.MultATA(F);		/*C = F^TF*/
	
	/*equilibrium contribution*/
	/*calculate eq. matrix contribution*/
	ComputeMatrixMod(fb, fStress, fModulus);

	/*	if (0)
	{
		cout << "\nelem: "<<elem<<"\tip: "<<ip;
		cout << "\nQ: "<<GetRotation();
		cout << "\nfC: "<<fC;
		cout<< "\nFiberStretch: "<<fFiberStretch;
		cout <<"\nMatrixMod: "<<fModulus;
	}
	*/

	/* eq. fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberMod(fFiberStretch, fFiberStress, fFiberMod);

	/* rotate and assemble eq. stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);
	
	/* rotate and assemble eq. modulus to lab coordinates */
	AssembleFiberModuli(fFiberMod, fModulus);
	
	/*
	if (elem == 0 || elem == 2209)
	{
		cout << "\nFiberMod: "<<fFiberMod;
		cout << "\nModulus: "<<fModulus;
	}
	*/

if (fNumFibProcess+fNumMatProcess > 0)
{
	/*calculate nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	for (int i = 0; i < fNumMatProcess; i++)
	{
		const dMatrixT& F_n = F_mechanical_last();  /*calculate inelastic deviatoric stress*/
		fb_n.MultAAT(F_n); 		/*b = FF^T*/

		/*calculate neq. matrix contribution*/
		ComputeMatrixAlgMod(fb, fb_n, fHm[i], fHm_n[i], fModulus, i, dSymMatrixT::kAccumulate);
		fStress.AddScaled(1.0, fHm[i]);
	}

	for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
	{
		const dMatrixT& F_n = F_mechanical_last();  /*calculate inelastic deviatoric stress*/
		fC_n.MultATA(F_n); 		/*b = FF^T*/
		/*fiber contribution*/
		ComputeFiberStretch(fC_n, fFiberStretch_n);

		/*calculate SNEQ and dSNEQ/dC*/
		ComputeFiberAlgMod(fFiberStretch, fFiberStretch_n, fFiberStress, fhf[i], fhf_n[i], fFiberMod, i);				

		/* rotate neq. stress to lab coordinates and assemble in fStress */
		AssembleFiberStress(fFiberStress, fStress);

		/* rotate and assemble neq. modulus to lab coordinates */
		AssembleFiberModuli(fFiberMod, fModulus);
	}
}
	return fModulus;
}
	
/* stress */
const dSymMatrixT& QLV_Nfibers::S_IJ(void)
{
		
	/* stretch */
	const dMatrixT& F = F_mechanical();  /*calculate inelastic deviatoric stress*/
	fb.MultAAT(F); 		/*b = FF^T*/
	fC.MultATA(F);		/*C = F^TF*/
	
	/*matrix contribution*/
	/*calculate matrix contribution*/
	ComputeMatrixStress(fb, fStress);

	/*fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberStress(fFiberStretch, fFiberStress);
	/* rotate and assemble stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);
	
	/*calculate nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
if (fNumMatProcess + fNumFibProcess > 0)
{

    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {	
		for (int i = 0; i < fNumMatProcess && fNumMatProcess > 0; i++)
		{
			const dMatrixT& F_n = F_mechanical_last();  /*calculate inelastic deviatoric stress*/
			fb_n.MultAAT(F_n); 		/*b = FF^T*/

			ComputeMatrixOverStress(fb, fb_n, fHm[i], fHm_n[i], i);
			fStress.AddScaled(1.0, fHm[i]);
		}
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{			
			const dMatrixT& F_n = F_mechanical_last();  /*calculate inelastic deviatoric stress*/
			fC_n.MultATA(F_n); 		/*b = FF^T*/
			/*fiber contribution*/
			ComputeFiberStretch(fC_n, fFiberStretch_n);

			/*compute fiber stress*/
			ComputeFiberOverStress(fFiberStretch, fFiberStretch_n, fFiberStress, fhf[i], fhf_n[i], i);

			/* rotate neq. stress to lab coordinates and assemble in fStress */
			AssembleFiberStress(fFiberStress, fStress);

		}
		Store(element, CurrIP());
	}
	else 
	{
		for (int i = 0; i < fNumMatProcess && fNumMatProcess > 0; i++)
		{
			const dMatrixT& F_n = F_mechanical_last();  /*calculate inelastic deviatoric stress*/
			fb_n.MultAAT(F_n); 		/*b = FF^T*/

			/*calculate neq. matrix contribution*/
			ComputeMatrixOverStress(fb, fb_n,  fHm[i], fHm_n[i],  i);
			fStress.AddScaled(1.0, fHm[i]);
		}
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{			
			/* neq. fiber contribution*/
			ComputeFiberOverStress(fFiberStretch, fFiberStretch_n, fFiberStress, fhf[i], fhf_n[i], i);
				
			/* rotate neq. stress to lab coordinates and assemble in fStress */
			AssembleFiberStress(fFiberStress, fStress);
		}
	}
}
	return(fStress);
}

/* material description */
const dMatrixT& QLV_Nfibers::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& QLV_Nfibers::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
	return fStress;
}


/*initializes history variable */
void  QLV_Nfibers::PointInitialize(void)
{
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0 && fNumMatProcess+fNumFibProcess > 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
		      
			  for (int i = 0; i < fNumMatProcess; i++)
			  {
				fHm_n[i] = 0.0;
				fHm[i] = 0.0;
			  }

			  for (int i = 0; i < fNumFibProcess; i++)
			  {
				fhf_n[i] = 0.0;
				fhf[i] = 0.0;
			  }
		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void QLV_Nfibers::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP() && fNumMatProcess+fNumFibProcess > 0; ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		for (int i = 0; i < fNumMatProcess; i++)
			fHm_n[i] = fHm[i];

		for (int i = 0; i < fNumFibProcess; i++)
			fhf_n[i] = fhf[i];

		/* write to storage */
		Store(element, ip);
	}
}

void QLV_Nfibers::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP() && fNumMatProcess+fNumFibProcess > 0; ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		/* assign "current" to "last" */	
		for (int i = 0; i < fNumMatProcess; i++)
			fHm[i] = fHm_n[i];
		
		for (int i = 0; i < fNumFibProcess; i++)
			fhf[i] = fhf_n[i];
		
		/* write to storage */
		Store(element, ip);
	}
}

void QLV_Nfibers::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}
void QLV_Nfibers::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;

}

/* information about subordinate parameter lists */
void QLV_Nfibers::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSFiberMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("eq_matrix_potential", ParameterListT::Once);
	sub_list.AddSub("neq_matrix_potential", ParameterListT::Any);
	sub_list.AddSub("matrix_relax_time", ParameterListT::Any);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("fibers_params", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
void QLV_Nfibers::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
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
	else if(name=="fiber_relaxtime")
	{
		sub_lists.AddSub("fiber_relax_time");	
	}		
	else /* inherited */
		FSFiberMatT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* QLV_Nfibers::NewSub(const StringT& name) const
{
	/* inherited */

	C1FunctionT* func = NULL;
	if (name == "fung-type")
		func = new FungType;
	else if (name == "parabola")
		func = new ParabolaT;
	else if (name == "fung_type_square")
		func = new FungType2;
	if (func)
		return func;

	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "mooney-rivlin")
		pot = new MooneyRivlin;
	else if (name == "veronda-westmann")
		pot = new VWPotentialT;
	if (pot)
		return pot;

	if (name == "eq_matrix_potential" || name == "neq_matrix_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
	
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		choice->AddSub("mooney-rivlin");
		choice->AddSub("veronda-westmann");
		return(choice);
	}
	else if (name == "matrix_relax_time")
	{
		ParameterContainerT* tau_m = new ParameterContainerT(name);
		ParameterT taum(ParameterT::Double, "tau_m");
		tau_m->AddParameter(taum);
		return(tau_m);
	}
	/* inherited */
	else if (name == "fibers_params")
	{
		ParameterContainerT* fiber = new ParameterContainerT(name);
//		fiber->SetListOrder(ParameterListT::Sequence);
		fiber->SetSubSource(this);
		fiber->AddSub("eq_fiber_pot_choice", ParameterListT::Once,true);
		fiber->AddSub("neq_fiber_pot_choice", ParameterListT::Any,true);
		fiber->AddSub("fiber_relaxtime", ParameterListT::Any,true);

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
	else if(name=="fiber_relax_time")
	{
		ParameterContainerT* tau_f = new ParameterContainerT(name);
		ParameterT tauf(ParameterT::Double, "tau_f");
		tau_f->AddParameter(tauf);
		return(tau_f);
	}		
	else return(FSFiberMatT::NewSub(name));
}

void QLV_Nfibers::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSFiberMatT::DefineParameters(list);

	ParameterT same(fsame, "use_for_all_fibers");
	same.SetDefault(false);
	list.AddParameter(same);
}

/* accept parameter list */
void QLV_Nfibers::TakeParameterList(const ParameterListT& list)
{

	StringT caller = "QLV_Nfibers::TakeParameterList";
		
	
	/* inherited */
	FSFiberMatT::TakeParameterList(list);
	fNumFibStress = dSymMatrixT::NumValues(fNumSD);
	

	int num_fibers = list.NumLists("fibers_params");
	fsame = false;
	fsame = list.GetParameter("use_for_all_fibers");
	if (fsame & num_fibers > 1)
			ExceptionT::GeneralFail("QLV_Nfibers::TakeParameterList", 
				"more than one list of fiber parameters given with option use_for_all_fibers");
	
	for (int i = 0; i< num_fibers; i++)
	{
		const ParameterListT& fiberlist = list.GetList("fibers_params",i);
		int num_fib_neq =  fiberlist.NumLists("neq_fiber_pot");
		int num_fib_visc =  fiberlist.NumLists("fiber_relax_time");
		
		if (num_fib_neq != num_fib_visc)
			ExceptionT::GeneralFail("QLV_Nfibers::TakeParameterList", 
				"number of fiber relaxation times does not match number of fiber nonequilibrium potentials");
		fNumFibProcess = num_fib_neq;
		ftime_f.Dimension(num_fibers, fNumFibProcess);
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

			if (!fPot_f(i,j+1)) throw ExceptionT::kOutOfMemory;
			fPot_f(i,j+1)->TakeParameterList(fiber_neq);

			const ParameterListT& fiber_time = fiberlist.GetList("fiber_relax_time",j);
			ftime_f(i,j) = fiber_time.GetParameter("tau_f");
		}
	}
	

	int num_mat_neq =  list.NumLists("neq_matrix_potential");
	int num_mat_tau =  list.NumLists("matrix_relax_time");

	if (num_mat_neq != num_mat_tau)
		ExceptionT::GeneralFail("QLV_Nfibers::TakeParameterList", 
			"number of matrix relaxation times does not match number of matrix nonequilibrium potentials");
	fNumMatProcess = num_mat_neq;
	fPot_m.Dimension(fNumMatProcess+1);	
	ftime_m.Dimension(fNumMatProcess);

	const ParameterListT& matrix_pot = list.GetListChoice(*this, "eq_matrix_potential");
	if(matrix_pot.Name() == "neo-hookean")
		fPot_m[0] = new NeoHookean;
	else if(matrix_pot.Name() == "mooney-rivlin")
		fPot_m[0] = new MooneyRivlin;
	else if(matrix_pot.Name() == "veronda-westmann")
		fPot_m[0] = new VWPotentialT;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot_m[0]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", matrix_pot.Name().Pointer());			
	fPot_m[0]->TakeParameterList(matrix_pot);
	
	for (int i = 0; i < fNumMatProcess; i++)
	{
		const ParameterListT& matrix_neq = list.GetListChoice(*this, "neq_matrix_potential",i);
		if(matrix_neq.Name() == "mooney-rivlin")
			fPot_m[i+1] = new MooneyRivlin;
		else if(matrix_neq.Name() == "neo-hookean")
			fPot_m[i+1] = new NeoHookean;
		else if(matrix_neq.Name() == "veronda-westmann")
			fPot_m[i+1] = new VWPotentialT;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fPot_m[i+1]) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", matrix_pot.Name().Pointer());			
		fPot_m[i+1]->TakeParameterList(matrix_neq);
		fPot_m[i+1]->SetKappa(0.0);

		const ParameterListT& matrix_time = list.GetList("matrix_relax_time",i);
		ftime_m[i] = matrix_time.GetParameter("tau_m");
	}

	/*dimension state variable storage arrays*/
	SetStateVariables(fNumMatProcess, fNumFibProcess);
	
	/*Dimension workspace*/
	fb.Dimension(fNumSD);
	fb_n.Dimension(fNumSD);
	fC_n.Dimension(fNumSD);
	fFiberStretch_n.Dimension(fNumSD);
	fEigs.Dimension(fNumSD);

	ftau.Dimension(fNumSD);
	ftau_n.Dimension(fNumSD);

	fdtau_dep.Dimension(fNumSD);
	fModMat.Dimension(fNumFibStress);
	fMod3.Dimension(fNumFibStress);



}
	
/* accept parameter list */
void QLV_Nfibers::SetStateVariables(const int num_mat_proc, const int num_fib_proc)
{
	/*dimension state variables*/
	fHm.Dimension(num_mat_proc);
	fHm_n.Dimension(num_mat_proc);

	fhf.Dimension(num_fib_proc);
	fhf_n.Dimension(num_fib_proc);
	
	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);
	int numfib = fPot_f.MajorDim();

	fnstatev = 0;
	fnstatev += numstress*num_mat_proc;   /*current H*/
	fnstatev += numstress*num_mat_proc;   /*last H_n*/
	
	fnstatev += numfib*num_fib_proc;	/*current h*/
	fnstatev += numfib*num_fib_proc;	/*current h_n*/
		
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	
	/* assign pointers to current and last blocks of state variable array */
	for (int i = 0; i < num_mat_proc; i++)
	{
		fHm[i].Set(ndof, pstatev);
		pstatev += numstress;
		fHm_n[i].Set(ndof, pstatev);
		pstatev += numstress;
	}
	for (int i = 0; i < num_fib_proc; i++)
	{
		fhf[i].Set(numfib, pstatev);
		pstatev += numfib;
		fhf_n[i].Set(numfib, pstatev);
		pstatev += numfib;
	}
}
	
/***********************************************************************
 * Protected
 ***********************************************************************/

void QLV_Nfibers::ComputeMatrixStress (const dSymMatrixT& Stretch,  dSymMatrixT& Stress)
{
	Stress = 0.0;

	/*compute eigenvalues of Cauchy-Green stretch tensor C (or Ce if neq)*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(Stretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	double J = sqrt(fEigs.Product());

	/*deviatoric part*/
	double J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[0]->DevStress(fEigs, ftau);
	/*mean stress*/
	ftau += fPot_m[0]->MeanStress(J);

	const dMatrixT& F_tot = F_total();
	Stress.AddScaled(1.0, PullBack(F_tot,fSpectralDecompSpat.EigsToRank2(ftau)));		
}

void QLV_Nfibers::ComputeMatrixOverStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_n, 
		dSymMatrixT& H, const dSymMatrixT& H_n, const int process_index, const int fillmode)
{
	/*compute invariant with flory decomposition*/
	/*Calculate deviatoric stress*/
	if (fillmode == dSymMatrixT::kOverwrite)
		H = 0.0;
		
	if (fillmode == dSymMatrixT::kOverwrite)
		H = 0.0;

	/*compute eigenvalues of Cauchy-Green stretch tensor C (or Ce if neq)*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(Stretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	double J = sqrt(fEigs.Product());

	/*deviatoric part*/
	double J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[process_index+1]->DevStress(fEigs, ftau);
	/*mean stress*/
	ftau += fPot_m[process_index+1]->MeanStress(J);

	/*stress here is the overstress in local fiber coords*/
	double dt = fFSMatSupport->TimeStep();
	double taudtS = dt/ftime_m[process_index];

	double alpha = exp(-0.5*taudtS);
	double beta = exp(-taudtS);
		
	/*compute eigenvalues of C_n */
	fSpectralDecompSpat.SpectralDecomp_Jacobi(Stretch_n, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	J = sqrt(fEigs.Product());

	/*deviatoric part*/
	J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[process_index+1]->DevStress(fEigs, ftau_n);
		
	/*calculate exp(-0.5*dt/tau_S)*(s_n+1 - s_n);*/
	ftau.AddScaled(-1.0, ftau_n);		
	const dMatrixT& F_tot = F_total();
	H.AddScaled(alpha, PullBack(F_tot,fSpectralDecompSpat.EigsToRank2(ftau)));

	/*Q_n+1 = exp(-dt/tau_S)*Q_n + exp(-0.5*dt/tau_S)*(s - s_n)*/
	H.AddScaled(beta, H_n);
	
}
void QLV_Nfibers::ComputeMatrixAlgMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_n,
		dSymMatrixT& H, const dSymMatrixT& H_n, dMatrixT& Mod, const int process_index, const int fillmode)
{
	if (fillmode == dSymMatrixT::kOverwrite)
	{
		H = 0.0;
		Mod = 0.0;
	}

	fSpectralDecompSpat.SpectralDecomp_Jacobi(Stretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	double l0 = fEigs[0];
	double l1 = fEigs[1];
	double l2 = fEigs[2];
	
	/*deviatoric part*/
	double ninth = third*third;
	double J = sqrt(fEigs.Product());
	double J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[process_index+1]->DevStress(fEigs, ftau);
	fPot_m[process_index+1]->DevMod(fEigs, fdtau_dep);

	/*stress here is the overstress in local fiber coords*/
	double dt = fFSMatSupport->TimeStep();
	double taudtS = dt/ftime_m[process_index];

	double alpha = exp(-0.5*taudtS);
	double beta = exp(-taudtS);
		
	/*compute eigenvalues of Cauchy-Green stretch tensor C (or Ce if neq)*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(Stretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	J = sqrt(fEigs.Product());

	/*deviatoric part*/
	J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[process_index+1]->DevStress(fEigs, ftau_n);
		
	/*calculate exp(-0.5*dt/tau_S)*(s_n+1 - s_n);*/
	ftau.AddScaled(-1.0, ftau_n);		
	const dMatrixT& F_tot = F_total();
	H.AddScaled(alpha, PullBack(F_tot,fSpectralDecompSpat.EigsToRank2(ftau)));
	/*Q_n+1 = exp(-dt/tau_S)*Q_n + exp(-0.5*dt/tau_S)*(s - s_n)*/
	H.AddScaled(beta, H_n);
		
		
	/*assemble eigenvalues to moduli tensor*/
	fdtau_dep[0] -= 2.0*ftau[0];
	fdtau_dep[1] -= 2.0*ftau[1];
	fdtau_dep[2] -= 2.0*ftau[2];
	fModMat = fSpectralDecompSpat.EigsToRank4(fdtau_dep);

	double dl, coeff;
	dl = l0 - l1;
	if (fabs(dl) > kSmall)
		coeff = (ftau[0]*l1 - ftau[1]*l0)/dl;
	else 
		coeff = 0.5*(fdtau_dep(0,0)-fdtau_dep(0,1))-ftau[0];
	MixedRank4_3D(eigenvectors[0], eigenvectors[1], fMod3);
	fModMat.AddScaled(2.0*coeff, fMod3);
    
	dl = l0 - l2;
	if (fabs(dl) > kSmall)
		coeff = (ftau[0]*l2 - ftau[2]*l0)/dl;
	else 
		coeff = 0.5*(fdtau_dep(0,0)-fdtau_dep(0,2))-ftau[2];	
	MixedRank4_3D(eigenvectors[0], eigenvectors[2], fMod3);
	fModMat.AddScaled(2.0*coeff, fMod3);

	dl = l1 - l2;
	if (fabs(dl) > kSmall)
		coeff  = (ftau[1]*l2 - ftau[2]*l1)/dl;
	else
		coeff = 0.5*(fdtau_dep(1,1)-fdtau_dep(1,2))-ftau[1];	
	MixedRank4_3D(eigenvectors[1], eigenvectors[2], fMod3);
	fModMat.AddScaled(2.0*coeff, fMod3);
	

	/*assemble eigenvalues to stress tensor*/
	/* transform spatial moduli and add to Mod*/
	Mod.AddScaled(alpha, PullBack(F_tot, fModMat));
}

void QLV_Nfibers::ComputeMatrixMod (const dSymMatrixT& Stretch,  dSymMatrixT& Stress, dMatrixT& Mod)
{
	Stress = 0.0;
	Mod = 0.0;

	fSpectralDecompSpat.SpectralDecomp_Jacobi(Stretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	double l0 = fEigs[0];
	double l1 = fEigs[1];
	double l2 = fEigs[2];
	
	/*deviatoric part*/
	double ninth = third*third;
	double J = sqrt(fEigs.Product());
	double J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[0]->DevStress(fEigs, ftau);
	fPot_m[0]->DevMod(fEigs, fdtau_dep);
	
	/*mean stress*/
	ftau += fPot_m[0]->MeanStress(J);
	/*add mean eigen mod*/
	fdtau_dep += fPot_m[0]->MeanMod(J);

	/*assemble eigenvalues to moduli tensor*/
	fdtau_dep[0] -= 2.0*ftau[0];
	fdtau_dep[1] -= 2.0*ftau[1];
	fdtau_dep[2] -= 2.0*ftau[2];
	fModMat = fSpectralDecompSpat.EigsToRank4(fdtau_dep);

	double dl, coeff;
	dl = l0 - l1;
	if (fabs(dl) > kSmall)
		coeff = (ftau[0]*l1 - ftau[1]*l0)/dl;
	else 
		coeff = 0.5*(fdtau_dep(0,0)-fdtau_dep(0,1))-ftau[0];
	MixedRank4_3D(eigenvectors[0], eigenvectors[1], fMod3);
	fModMat.AddScaled(2.0*coeff, fMod3);
    
	dl = l0 - l2;
	if (fabs(dl) > kSmall)
		coeff = (ftau[0]*l2 - ftau[2]*l0)/dl;
	else 
		coeff = 0.5*(fdtau_dep(0,0)-fdtau_dep(0,2))-ftau[2];	
	MixedRank4_3D(eigenvectors[0], eigenvectors[2], fMod3);
	fModMat.AddScaled(2.0*coeff, fMod3);
    
	dl = l1 - l2;
	if (fabs(dl) > kSmall)
		coeff  = (ftau[1]*l2 - ftau[2]*l1)/dl;
	else
		coeff = 0.5*(fdtau_dep(1,1)-fdtau_dep(1,2))-ftau[1];	
	MixedRank4_3D(eigenvectors[1], eigenvectors[2], fMod3);
	fModMat.AddScaled(2.0*coeff, fMod3);
	

	/*assemble eigenvalues to stress tensor*/
	const dMatrixT& F_tot = F_total();
	Stress.AddScaled(1.0, PullBack(F_tot,fSpectralDecompSpat.EigsToRank2(ftau)));	

	/* transform spatial moduli and add to Mod*/
	Mod.AddScaled(1.0, PullBack(F_tot, fModMat));
}

/*computes integrated fiber stress in local frame*/
void QLV_Nfibers::ComputeFiberOverStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_n,  dSymMatrixT& FiberStress, dArrayT& h, 
		const dArrayT& h_n,	const int pindex, const int fillmode)
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
		if (!fsame)
			fibnum = i;
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
			
		/*calculate I4*/
		double I4, s4;
		I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		/*calc 2dW/dI*/
		s4 = fPot_f(fibnum, pindex+1)->DFunction(I4); 

		/*stress here is the overstress in local fiber coords*/
		double dt = fFSMatSupport->TimeStep();
		
		double taudtS = dt/ftime_f(fibnum, pindex);

		double alpha = exp(-0.5*taudtS);
		double beta = exp(-taudtS);

		/*calculate elastic structural invariants*/
		double I4_n, s4_n;
		I4_n = FiberStretch_n[0]*cost*cost + FiberStretch_n[1]*sint*sint + 2.0*FiberStretch_n[5]*sint*cost;
		/*calc 2dWeq/dI*/
		s4_n = fPot_f(fibnum, pindex+1)->DFunction(I4_n);
		
		h[fibnum] = beta*h_n[fibnum] + alpha*(s4-s4_n);				
		FiberStress[0] += h[fibnum]*cost*cost;
		FiberStress[1] += h[fibnum]*sint*sint;
		FiberStress[5] += h[fibnum]*sint*cost;				
	}
}
	
/*computes integrated fiber stress in local frame*/
void QLV_Nfibers::ComputeFiberStress (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress)
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
		if (!fsame)
			fibnum = i;
			
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
			
		/*calculate I4*/
		double I4, s4;
		I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		/*calc 2dW/dI*/
		s4 = fPot_f(fibnum, 0)->DFunction(I4); 

		FiberStress[0] += cost*cost*s4;
		FiberStress[1] += sint*sint*s4;
		FiberStress[5] += sint*cost*s4;				
	}
}

/*computes integrated moduli in local frame*/
void QLV_Nfibers::ComputeFiberMod (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress, dMatrixT& FiberMod)
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
	
	int fibnum = 0;
	for (int i = 0; i<numfibers; i++)
	{
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
			
		if (!fsame)
			fibnum = i;

		/*calculate I4*/
		double I4, s4, d4;
		I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		/*calc 2dW/dI*/
		s4 = fPot_f(fibnum, 0)->DFunction(I4); 
		d4= fPot_f(fibnum, 0)->DDFunction(I4);

		/*calculate I4*/
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
}


void QLV_Nfibers::ComputeFiberAlgMod (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_n, dSymMatrixT& FiberStress, dArrayT& h, 
		const dArrayT& h_n, dMatrixT& FiberMod, const int pindex, const int fillmode)
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
	
	int fibnum = 0;
	for (int i = 0; i<numfibers; i++)
	{
		double cost = Fibers(i,0)*p1[0] + Fibers(i,1)*p1[1] + Fibers(i,2)*p1[2];
		double sint = Fibers(i,0)*p2[0] + Fibers(i,1)*p2[1] + Fibers(i,2)*p2[2];
			
		if (!fsame)
			fibnum = i;

		/*calculate I4*/
		double I4, s4, d4;
		I4 = FiberStretch[0]*cost*cost + FiberStretch[1]*sint*sint + 2.0*FiberStretch[5]*sint*cost;	
		/*calc 2dW/dI*/
		s4 = fPot_f(fibnum, pindex+1)->DFunction(I4); 
		d4= fPot_f(fibnum, pindex+1)->DDFunction(I4);

		/*stress here is the overstress in local fiber coords*/
		double dt = fFSMatSupport->TimeStep();
		double taudtS = dt/ftime_f(fibnum, pindex);

		double alpha = exp(-0.5*taudtS);
		double beta = exp(-taudtS);

		/*calculate elastic structural invariants*/
		double I4_n, s4_n;

		I4_n = FiberStretch_n[0]*cost*cost + FiberStretch_n[1]*sint*sint + 2.0*FiberStretch_n[5]*sint*cost;
		/*calc 2dWeq/dI*/
		s4_n = fPot_f(fibnum, pindex+1)->DFunction(I4_n);

		h[fibnum] = beta*h_n[fibnum] + alpha*(s4-s4_n);				
		FiberStress[0] += h[fibnum]*cost*cost;
		FiberStress[1] += h[fibnum]*sint*sint;
		FiberStress[5] += h[fibnum]*sint*cost;				

		FiberMod(0,0) += cost*cost*cost*cost*2.0*d4*alpha;
		FiberMod(1,1) += sint*sint*sint*sint*2.0*d4*alpha;
		FiberMod(5,5) += sint*sint*cost*cost*2.0*d4*alpha;
		
		FiberMod(0,1) += sint*sint*cost*cost*2.0*d4*alpha;
		FiberMod(1,0) += sint*sint*cost*cost*2.0*d4*alpha;
		
		FiberMod(0,5) += sint*cost*cost*cost*2.0*d4*alpha;
		FiberMod(5,0) += sint*cost*cost*cost*2.0*d4*alpha;
		
		FiberMod(1,5) += sint*sint*sint*cost*2.0*d4*alpha;
		FiberMod(5,1) += sint*sint*sint*cost*2.0*d4*alpha;
	}
}




