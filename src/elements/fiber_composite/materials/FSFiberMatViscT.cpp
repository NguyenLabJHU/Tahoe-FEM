/* $Id: FSFiberMatViscT.cpp,v 1.11 2008-10-25 20:44:53 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatViscT.h"
#include "FSFiberMatSupportT.h"
#include "ParameterContainerT.h"

#include "MooneyRivlin.h"
#include "NeoHookean.h"
#include "VWPotentialT.h"

#include "LinearExponentialT.h"
#include "ScaledCsch.h"

#include "iArray2DT.h"

static const double third = 1.0/3.0;

using namespace Tahoe;

/* constructor */
FSFiberMatViscT::FSFiberMatViscT(void):
	FSFiberMatT(),
	fSpectralDecompSpat(3),
	ParameterInterfaceT("vicoelastic_fiber_composite_material")
{
	/*Reset default*/
	fNumFibProcess = 0;
	fNumMatProcess = 0;
}


/* destructor */
FSFiberMatViscT::~FSFiberMatViscT(void) 
{ 
	/*allocated?*/
	double l1 = fPot_m.Length();
	double l2 = fVisc_m.Length();
	for (int i = 0; i < l1; i++){
			delete fPot_m[i];
	}
	for (int i = 0; i < l2; i++){
			delete fVisc_m[i];
	}
}

/* modulus */
const dMatrixT& FSFiberMatViscT::C_IJKL(void)
{
	int elem = CurrElementNumber();
	int ip = CurrIP();
	/* stretch */
	Compute_C(fC);
	
	/*equilibrium contribution*/
	/*calculate eq. matrix contribution*/
	ComputeMatrixMod(fC, fStress, fModulus);	

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
		/*calculate neq. matrix contribution*/
		ComputeMatrixMod(fC, fC_v[i], fStress, fModulus, i, dSymMatrixT::kAccumulate);
	}

	int j = fNumMatProcess;
	for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
	{
		/* neq. fiber contribution*/
		ComputeFiberStretch(fC_v[j], fFiberStretch_v);

		/*calculate SNEQ and dSNEQ/dC*/
		ComputeFiberMod(fFiberStretch, fFiberStretch_v, fFiberStress, fFiberMod, i);				

		/* rotate neq. stress to lab coordinates and assemble in fStress */
		AssembleFiberStress(fFiberStress, fStress);

		/* rotate and assemble neq. modulus to lab coordinates */
		AssembleFiberModuli(fFiberMod, fModulus);

		j++;
	}
}
	return fModulus;
}
	
/* stress */
const dSymMatrixT& FSFiberMatViscT::S_IJ(void)
{
		
	/* stretch */
	Compute_C(fC);

	/*matrix contribution*/
	/*calculate matrix contribution*/
	ComputeMatrixStress(fC, fStress);

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
			/*calculate neq. matrix contribution*/
			ComputeMatrixCv(fC, fC_vn[i], fC_v[i], i);
			ComputeMatrixStress(fC, fC_v[i], fStress, i, dSymMatrixT::kAccumulate);
		}
		int j = fNumMatProcess;
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{
			/* neq. fiber contribution*/
/*			if(Need_F_last())
			{
				const dMatrixT& F_last = F_mechanical_last();
				fC_n.MultATA(F_last);
			}
			Compute_Cv(fC_n, fC, fC_vn[j], fC_v[j], i);		
			*/
			/*rotate stretches in fiber plane*/
			ComputeFiberStretch(fC_vn[j], fFiberStretch_vn);
			ComputeFiberStretch(fC_v[j], fFiberStretch_v);
			
			Compute_Cv(fFiberStretch, fFiberStretch_vn, fFiberStretch_v, i);		
			/*Rotate fFiberStretch_v from local fiber coords to lab coordinates to get C_v*/		
			AssembleFiberStress(fFiberStretch_v, fC_v[j], dSymMatrixT::kOverwrite);

			/*compute fiber stress*/
			ComputeFiberStress(fFiberStretch, fFiberStretch_v, fFiberStress, i);

			/* rotate neq. stress to lab coordinates and assemble in fStress */
			AssembleFiberStress(fFiberStress, fStress);

			j++;
		}
		Store(element, CurrIP());
	}
	else 
	{
		for (int i = 0; i < fNumMatProcess && fNumMatProcess > 0; i++)
		{
			/*calculate neq. matrix contribution*/
			ComputeMatrixStress(fC, fC_v[i], fStress,  i, dSymMatrixT::kAccumulate);
		}
		int j = fNumMatProcess;
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{			
			/* neq. fiber contribution*/
			ComputeFiberStretch(fC_v[j], fFiberStretch_v);
			ComputeFiberStress(fFiberStretch, fFiberStretch_v, fFiberStress, i);
				
			/* rotate neq. stress to lab coordinates and assemble in fStress */
			AssembleFiberStress(fFiberStress, fStress);
			j++;
		}
	}
}
		if (0)
	{
		cout << "\nelem: "<<CurrElementNumber()<<"\tip: "<<CurrIP();
//		cout << "\nQ: "<<GetRotation();
		cout << "\nfC: "<<fC;
		cout << "\nS: "<<fStress;
//		cout<< "\nFiberStretch: "<<fFiberStretch;
//		cout <<"\nMatrixMod: "<<fModulus;
	}
	return(fStress);
}

/* material description */
const dMatrixT& FSFiberMatViscT::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& FSFiberMatViscT::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
	return fStress;
}

/*initializes history variable */
void  FSFiberMatViscT::PointInitialize(void)
{
	int numprocess = fNumFibProcess+fNumMatProcess;
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0 && numprocess > 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
		      
			  for (int i = 0; i < numprocess; i++)
			  {
				fC_vn[i].Identity();
				fC_v[i].Identity();
			  }

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void FSFiberMatViscT::UpdateHistory(void)
{
	int numprocess = fNumFibProcess+fNumMatProcess;
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP() && numprocess > 0; ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		for (int i = 0; i < numprocess; i++)
			fC_vn[i] = fC_v[i];

		/* write to storage */
		Store(element, ip);
	}
}

void FSFiberMatViscT::ResetHistory(void)
{
	int numprocess = fNumFibProcess+fNumMatProcess;
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP() && numprocess > 0; ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		/* assign "current" to "last" */	
		for (int i = 0; i < numprocess; i++)
			fC_v[i] = fC_vn[i];
		
		/* write to storage */
		Store(element, ip);
	}
}

void FSFiberMatViscT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}
void FSFiberMatViscT::Store(ElementCardT& element, int ip)
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
void FSFiberMatViscT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("eq_matrix_potential", ParameterListT::Once);
	sub_list.AddSub("neq_matrix_potential", ParameterListT::Any);

	/* choice of viscosity */
	sub_list.AddSub("matrix_visc_potential", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSFiberMatViscT::NewSub(const StringT& name) const
{
	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "mooney-rivlin")
		pot = new MooneyRivlin;
	else if (name == "veronda-westmann")
		pot = new VWPotentialT;
	if (pot)
		return pot;

	C1FunctionT* func = NULL;
	if (name == "scaled-csch")
		func = new ScaledCsch;
	else if (name == "linear_exponential")
		func = new LinearExponentialT;

	if (func)
		return func;

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
	else if (name == "matrix_visc_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		choice->AddSub("scaled-csch");
		choice->AddSub("linear_exponential");
		return(choice);
	}
	else return(FSFiberMatT::NewSub(name));
}

/* accept parameter list */
void FSFiberMatViscT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);
	StringT caller = "FSFiberMatViscT::TakeParameterList";
	int num_mat_neq =  list.NumLists("neq_matrix_potential");
	int num_mat_visc = list.NumLists("matrix_visc_potential");
	if (num_mat_neq != num_mat_visc)
		ExceptionT::GeneralFail("FSFiberMatViscT::TakeParameterList", 
			"number of matrix viscosity functions does not match number of matrix nonequilibrium potentials");
	fNumMatProcess = list.NumLists("matrix_visc_potential");


	fPot_m.Dimension(fNumMatProcess+1);
	fVisc_m.Dimension(fNumMatProcess);
		
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

		const ParameterListT& matrix_visc = list.GetListChoice(*this, "matrix_visc_potential", i);
		if (matrix_visc.Name() == "linear_exponential")
			fVisc_m[i] = new LinearExponentialT;
		else if (matrix_visc.Name() == "scaled-csch")
			fVisc_m[i] = new ScaledCsch;
		else 
			ExceptionT::GeneralFail(caller, "no such potential");
		if (!fVisc_m[i]) throw ExceptionT::kOutOfMemory;
		fVisc_m[i]->TakeParameterList(matrix_visc);
	}
				
	/*viscous fiber stretch at time step n, viscous stretch at time step n and vn*/
	fFiberStretch_v.Dimension(fNumSD);
	fFiberStretch_vn.Dimension(fNumSD);

	/*Dimension work spaces*/
	fCalg.Dimension(fNumFibStress);	
	/*Dimension work spaces for matrix calculation*/

	fInverse.Dimension(fNumSD);
	fModMat.Dimension(fNumFibStress);
	fMod3.Dimension(fNumFibStress);

	fb.Dimension(fNumSD);
	fEigs.Dimension(fNumSD);
	ftau.Dimension(fNumSD);
	fdtau_dep.Dimension(fNumSD);

  	fiK_m.Dimension(fNumSD);
	fCalg_m.Dimension(fNumSD);
}

/* accept parameter list */
void FSFiberMatViscT::SetStateVariables(const int numprocess)
{
	/*dimension state variables*/
	fC_v.Dimension(numprocess);
	fC_vn.Dimension(numprocess);

	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/

	fnstatev *= numprocess;
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	
	/* assign pointers to current and last blocks of state variable array */
	for (int i = 0; i < numprocess; i++)
	{
		fC_v[i].Set(ndof, pstatev);
		pstatev += numstress;
		fC_vn[i].Set(ndof, pstatev);
		pstatev += numstress;
	}
}

void FSFiberMatViscT::ComputeMatrixStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, 
			dSymMatrixT& Stress, const int process_index, const int fillmode)
{
	/*compute invariant with flory decomposition*/
	if (fillmode == dSymMatrixT::kOverwrite)
		Stress = 0.0;
		
	const dMatrixT& F = F_mechanical();
	if (process_index > -1)
	{
		/*be = F.Cv^-1.F^T*/
		fInverse.Inverse(Stretch_v);
		fb.MultQBQT(F, fInverse);
	}
	else
		fb.MultAAT(F); 		/*b = FF^T*/

	if (fillmode == dSymMatrixT::kOverwrite)
		Stress = 0.0;

	/*compute eigenvalues of Cauchy-Green stretch tensor C (or Ce if neq)*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	double J = sqrt(fEigs.Product());

	/*deviatoric part*/
	double J23 = pow(J,-2.0*third);
	fEigs *= J23;
		
	/*dev stress*/
	fPot_m[process_index+1]->DevStress(fEigs, ftau);
	/*mean stress*/
	ftau += fPot_m[process_index+1]->MeanStress(J);

	const dMatrixT& F_tot = F_total();
	Stress.AddScaled(1.0, PullBack(F_tot,fSpectralDecompSpat.EigsToRank2(ftau)));	
}

void FSFiberMatViscT::ComputeMatrixMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress,
			dMatrixT& Mod, const int process_index, const int fillmode)
{
	if (fillmode == dSymMatrixT::kOverwrite)
	{
		Stress = 0.0;
		Mod = 0.0;
	}

	/*compute invariant with flory decomposition*/
	const dMatrixT& F = F_mechanical();
	
	if (process_index < 0)  /*eq*/
	{
		/*b = FF^T*/
		fb.MultAAT(F);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
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
	}
	else  /*neq*/
	{
		/*b_tr = F.Cv^-1_n.F^T*/
		fInverse.Inverse(fC_vn[process_index]);
		fb.MultQBQT(F, fInverse);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
		const dArrayT& eigenvals_tr = fSpectralDecompSpat.Eigenvalues(); /*eigs_tr*/
		const double ltr0 = eigenvals_tr[0];
		const double ltr1 = eigenvals_tr[1];
		const double ltr2 = eigenvals_tr[2];

		/*b_e = F.Cv^-1_n.F^T*/
		fInverse.Inverse(Stretch_v);
		fb.MultQBQT(F, fInverse);

		fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
		fEigs = fSpectralDecompSpat.Eigenvalues();
		const ArrayT<dArrayT>& eigenvectors = fSpectralDecompSpat.Eigenvectors();
		double Je = sqrt(fEigs.Product());
		
		/*deviatoric part*/
		double Je23 = pow(Je,-2.0*third);
		fEigs *= Je23;
				
		/*dev stress*/
		fPot_m[process_index+1]->DevStress(fEigs, ftau);
		fPot_m[process_index+1]->DevMod(fEigs, fdtau_dep);
	
		/*viscosity*/
		double stress_mag = sqrt(ftau[0]*ftau[0]+ftau[1]*ftau[1]+ftau[2]*ftau[2]);
		double ietaS = 1.0/fVisc_m[process_index]->Function(stress_mag);

		/*calculates K_AB*/
		double dt = fFSMatSupport->TimeStep();
		fiK_m(0,0) = 1+0.5*ietaS*dt*fdtau_dep[0];
		fiK_m(1,1) = 1+0.5*ietaS*dt*fdtau_dep[1];
		fiK_m(2,2) = 1+0.5*ietaS*dt*fdtau_dep[2];

		fiK_m(1,2) = 0.5*ietaS*dt*fdtau_dep[3];
		fiK_m(0,2) = 0.5*ietaS*dt*fdtau_dep[4];
		fiK_m(0,1) = 0.5*ietaS*dt*fdtau_dep[5];
       
		fiK_m(2,1) = fiK_m(1,2);
		fiK_m(2,0) = fiK_m(0,2);
		fiK_m(1,0) = fiK_m(0,1);
		/*inverts KAB*/
		fiK_m.Inverse(fiK_m);

		/*calculate algorithmic moduli for matrix*/
		fCalg_m(0,0) = fdtau_dep(0,0)*fiK_m(0,0) + fdtau_dep(0,1)*fiK_m(1,0) + fdtau_dep(0,2)*fiK_m(2,0) 
			- 2.0*ftau[0];
		fCalg_m(1,0) = fdtau_dep(1,0)*fiK_m(0,0) + fdtau_dep(1,1)*fiK_m(1,0) + fdtau_dep(1,2)*fiK_m(2,0);
		fCalg_m(2,0) = fdtau_dep(2,0)*fiK_m(0,0) + fdtau_dep(2,1)*fiK_m(1,0) + fdtau_dep(2,2)*fiK_m(2,0);
		fCalg_m(0,1) = fdtau_dep(0,0)*fiK_m(0,1) + fdtau_dep(0,1)*fiK_m(1,1) + fdtau_dep(0,2)*fiK_m(2,1);
		fCalg_m(1,1) = fdtau_dep(1,0)*fiK_m(0,1) + fdtau_dep(1,1)*fiK_m(1,1) + fdtau_dep(1,2)*fiK_m(2,1) 
			- 2.0*ftau[1];
		fCalg_m(2,1) = fdtau_dep(2,0)*fiK_m(0,1) + fdtau_dep(2,1)*fiK_m(1,1) + fdtau_dep(2,2)*fiK_m(2,1);
		fCalg_m(0,2) = fdtau_dep(0,0)*fiK_m(0,2) + fdtau_dep(0,1)*fiK_m(1,2) + fdtau_dep(0,2)*fiK_m(2,2);
		fCalg_m(1,2) = fdtau_dep(1,0)*fiK_m(0,2) + fdtau_dep(1,1)*fiK_m(1,2) + fdtau_dep(1,2)*fiK_m(2,2);
		fCalg_m(2,2) = fdtau_dep(2,0)*fiK_m(0,2) + fdtau_dep(2,1)*fiK_m(1,2) + fdtau_dep(2,2)*fiK_m(2,2)
			- 2.0*ftau[2];
		
		fModMat = fSpectralDecompSpat.NonSymEigsToRank4(fCalg_m);

		double dl_tr, coeff;
	
		dl_tr = ltr0 - ltr1;
		if (fabs(dl_tr) > kSmall)
			coeff = (ftau[0]*ltr1 - ftau[1]*ltr0)/dl_tr;
		else 
			coeff = 0.5*(fCalg_m(0,0)-fCalg_m(0,1))-ftau[0];
		MixedRank4_3D(eigenvectors[0], eigenvectors[1], fMod3);
		fModMat.AddScaled(2.0*coeff, fMod3);
    
		dl_tr = ltr0 - ltr2;
		if (fabs(dl_tr) > kSmall)
			coeff =(ftau[0]*ltr2 - ftau[2]*ltr0)/dl_tr;
		else 
			coeff = 0.5*(fCalg_m(0,0)-fCalg_m(0,2))-ftau[2];	
		MixedRank4_3D(eigenvectors[0], eigenvectors[2], fMod3);
		fModMat.AddScaled(2.0*coeff, fMod3);
    
		dl_tr = ltr1 - ltr2;
		if (fabs(dl_tr) > kSmall)
			coeff  = (ftau[1]*ltr2 - ftau[2]*ltr1)/dl_tr;
		else
			coeff = 0.5*(fCalg_m(1,1)-fCalg_m(1,2))-ftau[1];	
		MixedRank4_3D(eigenvectors[1], eigenvectors[2], fMod3);
		fModMat.AddScaled(2.0*coeff, fMod3);
	}
	
	/*assemble eigenvalues to stress tensor*/
	const dMatrixT& F_tot = F_total();
	Stress.AddScaled(1.0, PullBack(F_tot,fSpectralDecompSpat.EigsToRank2(ftau)));	

    /* transform spatial moduli and add to Mod*/
    Mod.AddScaled(1.0, PullBack(F_tot, fModMat));
}

/*compute viscous stretch of matrix*/
void FSFiberMatViscT::ComputeMatrixCv(const dSymMatrixT& Stretch, const dSymMatrixT& Stretchv_last, 
	dSymMatrixT& Stretchv, const int process_index)
{
	/*calc trial state*/
	const dMatrixT& F = F_mechanical();
	
	fInverse.Inverse(Stretchv_last);
	
	fb.MultQBQT(F, fInverse);   /*btr = F.Cv_n^-1.F^T*/
	
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	const dArrayT& eigenvals_tr = fSpectralDecompSpat.Eigenvalues(); 
	
	/*trial eigen log strain*/
	const double ep_tr0 = 0.5*log(eigenvals_tr[0]);
	const double ep_tr1 = 0.5*log(eigenvals_tr[1]);
	const double ep_tr2 = 0.5*log(eigenvals_tr[2]);
	
	/*elastic eigen log strain*/
	dArrayT& le = fEigs;  /*using fEigs as work space*/
	le = eigenvals_tr;
	
	double ep_e0 = ep_tr0;
	double ep_e1 = ep_tr1;
	double ep_e2 = ep_tr2;

	/*kirchhoff eigen stress and dtau_dep*/
	double s0, s1, s2;
	double c00, c11, c22, c12, c02, c01;
	double ninth = third*third;
	
	double tol;
	
	/*initializes principle viscous stretch*/
	do 
	{
		/*compute eigen kirchhoff stress*/
		double J = sqrt(le.Product());

		/*deviatoric stress*/
		double J23 = pow(J,-2.0*third);

		/*temporarily stores lebar in le*/
		le *= J23;
		double l0_bar = le[0];
		double l1_bar = le[1];
		double l2_bar = le[2];
	
		fPot_m[process_index+1]->DevStress(le, ftau);
		fPot_m[process_index+1]->DevMod(le, fdtau_dep);
		le /= J23; /*reverting back to le*/
		
		s0 = ftau[0];
		s1 = ftau[1];
		s2 = ftau[2];
	    
		c00 = fdtau_dep(0,0);
		c11 = fdtau_dep(1,1);
		c22 = fdtau_dep(2,2);
		c12 = fdtau_dep(1,2);
		c02 = fdtau_dep(0,2);
		c01 = fdtau_dep(0,1);
	
	    /*calculate the residual*/
	    double dt = fFSMatSupport->TimeStep();
		const double stress_mag = sqrt(s0*s0+s1*s1+s2*s2);
		const double ietaS = 1.0/(fVisc_m[process_index]->Function(stress_mag));

	    double res0 = ep_e0 + dt*0.5*ietaS*s0  - ep_tr0;
	    double res1 = ep_e1 + dt*0.5*ietaS*s1  - ep_tr1;
	    double res2 = ep_e2 + dt*0.5*ietaS*s2  - ep_tr2;
		
		/*calc dres_A/dep_B*/
		fiK_m(0,0) = 1+0.5*ietaS*dt*c00;
		fiK_m(1,1) = 1+0.5*ietaS*dt*c11;
		fiK_m(2,2) = 1+0.5*ietaS*dt*c22;

		fiK_m(1,2) = 0.5*ietaS*dt*c12;
		fiK_m(0,2) = 0.5*ietaS*dt*c02;
		fiK_m(0,1) = 0.5*ietaS*dt*c01;
       
		fiK_m(2,1) = fiK_m(1,2);
		fiK_m(2,0) = fiK_m(0,2);
		fiK_m(1,0) = fiK_m(0,1);

		/*inverts KAB*/
		fiK_m.Inverse(fiK_m);

	    /*solve for the principal strain increments*/
	    double dep_e0 = -fiK_m(0,0)*res0 - fiK_m(0,1)*res1 - fiK_m(0,2)*res2;
	    double dep_e1 = -fiK_m(1,0)*res0 - fiK_m(1,1)*res1 - fiK_m(1,2)*res2;
	    double dep_e2 = -fiK_m(2,0)*res0 - fiK_m(2,1)*res1 - fiK_m(2,2)*res2;
	     
	    /*updates principal elastic stretches*/ 
	    ep_e0 += dep_e0;
	    ep_e1 += dep_e1;
	    ep_e2 += dep_e2;
	    
	    le[0] = exp(2.0*ep_e0);
	    le[1] = exp(2.0*ep_e1);
	    le[2] = exp(2.0*ep_e2);
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(res0*res0 + res1*res1+res2*res2);
		
	}while (tol > kSmall); 
	
	/*calculate be*/
	fb = fSpectralDecompSpat.EigsToRank2(le); /*be is colinear to b_tr*/
	/*Cv = F^T.be^-1.F*/
	Stretchv.MultQTBQ(F,fb.Inverse());
}

