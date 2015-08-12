/* $Id: AnisoScleraWaxs_expo.cpp,v 1.2 2015-08-12 14:55:29 tahoe.vickynguyen Exp $ */
/* created: TDN (01/22/2001) */

#include "AnisoScleraWaxs_expo.h"
#include <math.h>
#include "bessel.h"
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
//#include "FSFiberMatSupportT.h"


#ifdef VIB_MATERIAL

/* point generator */
#include "EvenSpacePtsT.h"

/*fiber potentials*/
#include "FungType.h"
#include "VWType.h"
#include "LanirFiber.h"
#include "LinearType.h"


/*viscosity functions*/
#include "ScaledCsch.h"

//#define MY_DEBUG

const double Pi = acos(-1.0);
const double third = 1.0/3.0;
const int kNumOutputVar = 19;
static const char* Labels[kNumOutputVar] = 
	{"Fiber_X", "Fiber_Y", "Fiber_Z","Orth_fiber_X", "Orth_fiber_Y", "Orth_fiber_Z","Normal_Sclera_X", "Normal_Sclera_Y", "Normal_Sclera_Z","J"};
static const int perm[3][3] = {0,1,2,1,2,0,2,0,1};
using namespace Tahoe;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
AnisoScleraWaxs_expo::AnisoScleraWaxs_expo(void):
  FSFiberMatViscT(),
	fSpectralDecompSpat(3),
	fCircle(NULL),
	fDType(kMeridional),
	ParameterInterfaceT("aniso_viscoelastic_sclera_waxs_expo")
{
	/*reset default*/
	fNumFibProcess = 0;
	fNumMatProcess = 0;

#ifndef VIB_MATERIAL
	ExceptionT::BadInputValue("AnisoScleraWaxs_expo::AnisoScleraWaxs_expo", 
		"VIB_MATERIAL must be enabled");
#endif
}

/* destructor */
AnisoScleraWaxs_expo::~AnisoScleraWaxs_expo(void) 
{ 
	delete fCircle; 
	/*allocated?*/
	if (fNumFibProcess > 0) 
	{	
		delete fPotential[0];	
		for (int i = 0; i < fNumFibProcess; i++)
		{
			delete fPotential[i+1];
			delete fViscosity[i];
		}
	}
}

/* free energy density */
double AnisoScleraWaxs_expo::StrainEnergyDensity(void)
{
	/*matrix contribution*/
	/* stretch */
	Compute_C(fC);

	double I3 = fC.Det();
	double I3rg = pow(I3, -fGamma);
	double I1 = fC[0]+fC[1]+fC[2];
	
	/*coupled compressible Neo-Hookean*/
	/* mu/2 (I1 -3) +mu/(2 gamma) (I3^-gamma -1)*/
	double energy = 0.5*fMu*(I1-3.0) + fMu/(2.0*fGamma)*(I3rg -1.0);
	/*fiber contribution*/
	/* stretched bonds */
	ComputeFiberStretch(fC, fFiberStretch);
	CompI4(fFiberStretch);

	/*equilibrium contribution*/
	/* update potential table */
	fPotential[0]->MapFunction(fI4,fU);

	/* sum contributions */
	double* pU = fU.Pointer();
	double* pj =  fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	for (int i = 0; i < fI4.Length(); i++) {
		energy += (*pU++)*(*pj++);
	}
	
		
	
	/*nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
  ElementCardT& element = CurrentElement();
  Load(element, CurrIP());

	for (int i = 0; i < fNumFibProcess; i++)
	{
		ComputeFiberStretch(fC_v[i], fFiberStretch_v);
		CompI4(fFiberStretch, fFiberStretch_v);
		
		dArrayT& I4e = fI4;
		I4e /= fI4v;
		
		fPotential[i+1]->MapFunction(I4e,fU);

		/* sum contributions */
		double* pU = fU.Pointer();
		double* pj =  fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

		for (int k = 0; k < I4e.Length(); k++)
			energy += (*pU++)*(*pj++);
	}
	return energy;
}

/* nonequilibrium free energy density */
double AnisoScleraWaxs_expo::NonequilibriumStrainEnergyDensity(void)
{
	/* no matrix contribution*/
	double energy = 0.0;
	/* stretch */
	Compute_C(fC);
	

	double I3 = fC.Det();
	double I3rg = pow(I3, -fGamma);
	double I1 = fC[0]+fC[1]+fC[2];
	
	/*fiber contribution*/
	/* stretched bonds */
	ComputeFiberStretch(fC, fFiberStretch);
	CompI4(fFiberStretch);

	/*nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	double* pj =  fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	
	for (int i = 0; i < fNumFibProcess; i++)
	{
		ComputeFiberStretch(fC_v[i], fFiberStretch_v);
		CompI4(fFiberStretch, fFiberStretch_v);
		
		dArrayT& I4e = fI4;
		I4e /= fI4v;
		
		fPotential[i+1]->MapFunction(I4e,fU);

		/* sum contributions */
		double* pU = fU.Pointer();
		for (int k = 0; k < I4e.Length(); k++)
			energy += (*pU++)*(*pj++);
	}
	return energy;
}

int AnisoScleraWaxs_expo::NumOutputVariables() const {
	return kNumOutputVar;
}

void AnisoScleraWaxs_expo::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void AnisoScleraWaxs_expo::ComputeOutput(dArrayT& output)
{
	/*calculates deformed fiber vectors*/
//	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	const dMatrixT& Q = GetRotation();
	const double* Fiber_dir = Q(0);
	const double* Fiber_orth = Q(1);
	const double* Fiber_norm = Q(2);
	
	const dMatrixT& F = F_mechanical();
	double* pb = output.Pointer();
	
	/*deformed NT fiber orientation*/
	F.Multx(Fiber_dir, pb);
	pb += NumSD();
	
	/*deformed IS fiber orientation*/
	F.Multx(Fiber_orth, pb);
	pb += NumSD();

	/*deformed OT fiber orientation*/
	F.Multx(Fiber_norm, pb);
	pb += NumSD();

	/*non-equilibrium strain energy density */
	*pb++ = F.Det();
	
	/*calculate the stress eigenvectors*/
	const dSymMatrixT& stress = s_ij();
	fSpectralDecompSpat.SpectralDecomp_Jacobi(stress, false);	
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();
	const double* v1 = eigenvectors[0].Pointer();
	*pb++ = *v1++;
	*pb++ = *v1++;
	*pb++ = *v1;
	const double* v2 = eigenvectors[1].Pointer();
	*pb++ = *v2++;
	*pb++ = *v2++;
	*pb++ = *v2;
	const double* v3 = eigenvectors[2].Pointer();
	*pb++ = *v3++;
	*pb++ = *v3++;
	*pb  =  *v3;
}

/* describe the parameters needed by the interface */
void AnisoScleraWaxs_expo::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);
	
	/* integration points */
	ParameterT points(ParameterT::Integer, "n_points");
	points.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(points);
	
	/*dispersion parameter*/
	ParameterT Dispersion(ParameterT::String,"Dispersion_file");
	Dispersion.SetDescription("Contains the dispersion parameter for each element");
	list.AddParameter(Dispersion);

	/* integration points */
	ParameterT vol_type(ParameterT::Enumeration, "pressure_model");
	vol_type.AddEnumeration("Ogden",kOgden);
	vol_type.AddEnumeration("Blatz",kBlatz);
	vol_type.SetDefault(kBlatz);
	list.AddParameter(vol_type);	
}

/* information about subordinate parameter lists */
void AnisoScleraWaxs_expo::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("matrix_material_params", ParameterListT::Once);
	
	/* choice of energy potential for fibrils */
	sub_list.AddSub("eq_fibril_potential_WAXS", ParameterListT::Once);
	sub_list.AddSub("neq_fibril_potential_WAXS", ParameterListT::Any);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("viscosity", ParameterListT::Any);

	}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AnisoScleraWaxs_expo::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = FSSolidMatT::NewSub(name);

	if (sub) 
	{
		return sub;
	}
	else if (name == "matrix_material_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* bound */
		LimitT lower(0.0, LimitT::Lower);
		
		/* exponential functions*/
		ParameterContainerT matrix("Neo-Hookean");
		ParameterT mu(ParameterT::Double, "shear_modulus");
		mu.AddLimit(lower);
		ParameterT gamma(ParameterT::Double, "bulk_modulus");
		gamma.AddLimit(lower);
		matrix.AddParameter(mu);
		matrix.AddParameter(gamma);
		choice->AddSub(matrix);
		return choice;
	}
	else if (name == "eq_fibril_potential_WAXS")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT fung("fung_type");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT alpha(ParameterT::Double, "alpha");
			ParameterT beta(ParameterT::Double, "beta");

			fung.AddParameter(alpha);
			fung.AddParameter(beta);
			alpha.AddLimit(lower);
			beta.AddLimit(lower);

			/* set the description */
			choice->AddSub(fung);
		}
		
		ParameterContainerT linearSE("linear_type_sclera");		
		{
			LimitT lower(0.0, LimitT::Lower);
			
			ParameterT theta(ParameterT::Double, "theta");
			
			linearSE.AddParameter(theta);
			theta.AddLimit(lower);
			
			/* set the description */
			choice->AddSub(linearSE);
		}

		
		ParameterContainerT fung0("vw_type");		
		{
			LimitT lower(0.0, LimitT::Lower);
			
			ParameterT alpha(ParameterT::Double, "alpha");
			ParameterT beta(ParameterT::Double, "beta");
			
			fung0.AddParameter(alpha);
			fung0.AddParameter(beta);
			alpha.AddLimit(lower);
			beta.AddLimit(lower);
			
			/* set the description */
			choice->AddSub(fung0);
		}
		
		ParameterContainerT lanir("lanir_fiber_recruit");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT K(ParameterT::Double, "fiber_stiffness_K");
			K.AddLimit(0,LimitT::Lower);
	
			ParameterT alpha(ParameterT::Double, "shape_param_alpha");
			alpha.AddLimit(1,LimitT::LowerInclusive);

			ParameterT beta(ParameterT::Double, "scale_param_beta");
			beta.AddLimit(0,LimitT::Lower);

			lanir.AddParameter(K);
			lanir.AddParameter(alpha);
			lanir.AddParameter(beta);

			/* set the description */
			choice->AddSub(lanir);
		}
		return(choice);
	}
	else if (name == "neq_fibril_potential_WAXS")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT fung("neq_fung_type");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT alpha(ParameterT::Double, "alpha");

			fung.AddParameter(alpha);
			alpha.AddLimit(lower);

			/* set the description */
			choice->AddSub(fung);
		}
		
		ParameterContainerT linearSE("neq_linear_type_sclera");		
		{
			LimitT lower(0.0, LimitT::Lower);
			
			ParameterT theta(ParameterT::Double, "theta");
			
			linearSE.AddParameter(theta);
			theta.AddLimit(lower);
			
			/* set the description */
			choice->AddSub(linearSE);
		}
		
		
		ParameterContainerT fung0("neq_vw_type");		
		{
			LimitT lower(0.0, LimitT::Lower);
			
			ParameterT alpha(ParameterT::Double, "alpha");
			
			fung0.AddParameter(alpha);
			alpha.AddLimit(lower);
			
			/* set the description */
			choice->AddSub(fung0);
		}
		
		ParameterContainerT lanir("neq_lanir_fiber_recruit");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT K(ParameterT::Double, "fiber_stiffness_K");
			K.AddLimit(0,LimitT::Lower);
	
			lanir.AddParameter(K);

			/* set the description */
			choice->AddSub(lanir);
		}
		return(choice);
	}
	else if (name == "viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT sinh("scaled_csch");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT eta0(ParameterT::Double, "eta0");
			ParameterT tau0(ParameterT::Double, "tau0");

			sinh.AddParameter(eta0);
			sinh.AddParameter(tau0);
			eta0.AddLimit(lower);
			tau0.AddLimit(0, LimitT::Lower);

			/* set the description */
			sinh.SetDescription("param oder, [eta0, tau0]");	
			choice->AddSub(sinh);
		}
		return(choice);
	}
}

/* accept parameter list */
void AnisoScleraWaxs_expo::TakeParameterList(const ParameterListT& list)
{
	char caller[] = "AnisoScleraWaxs_expo::TakeParameterList";
	/* inherited */
	FSFiberMatT::TakeParameterList(list); 

	/*initializing some parameters for cornea_mod fibril distribution*/
	int num_neq_pot = list.NumLists("neq_fibril_potential_WAXS");
	int num_visc = list.NumLists("viscosity");

	if (num_visc != num_neq_pot)
		ExceptionT::GeneralFail("AnisoScleraWaxs_expo::TakeParameterList", 
			"number of viscosity functions does not match number of nonequilibrium potentials");

	fNumFibProcess = num_neq_pot;
	fPotential.Dimension(fNumFibProcess+1);
	if(fNumFibProcess > 0)
		fViscosity.Dimension(fNumFibProcess);
		
	int b = list.GetParameter("pressure_model");
	fVolType = (b == kBlatz) ? kBlatz : kOgden;
	

	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_material_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fGamma = matrix.GetParameter("bulk_modulus");
	}

	const ParameterListT& potential = list.GetListChoice(*this, "eq_fibril_potential_WAXS");

	double temp1, temp2;
	if(potential.Name() == "fung_type")
	{
		double alpha_eq = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential[0] = new FungType(alpha_eq, beta);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
		
		temp1 = beta;
	}
	else if(potential.Name() == "vw_type")
	{
		double alpha_eq = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential[0] = new VWType(alpha_eq, beta);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
		
		temp1 = beta;
	}
	else if(potential.Name() == "linear_type_sclera")
	{
		double theta_eq = potential.GetParameter("theta");
		fPotential[0] = new LinearType(theta_eq);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
	}
	
	else if (potential.Name() == "lanir_fiber_recruit")
	{
		double K_eq= potential.GetParameter("fiber_stiffness_K");
		double alpha = potential.GetParameter("shape_param_alpha");
		double beta = potential.GetParameter("scale_param_beta");
		fPotential[0] = new LanirFiber(K_eq, alpha, beta);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
		
		temp1 = alpha;
		temp2 = beta;
	}
	else 
		ExceptionT::GeneralFail(caller, "no such potential");

	for (int i = 0; i < fNumFibProcess; i++)
	{
		const ParameterListT& neq_potential = list.GetListChoice(*this, "neq_fibril_potential_WAXS", i);
		if (neq_potential.Name() == "neq_fung_type")
		{
			double alpha_neq = neq_potential.GetParameter("alpha");

			fPotential[i+1] = new FungType(alpha_neq, temp1);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
			
		}
		if (neq_potential.Name() == "neq_linear_type_sclera")
		{
			double theta_neq = neq_potential.GetParameter("theta");
			
			fPotential[i+1] = new LinearType(theta_neq);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
			
		}
		else if (neq_potential.Name() == "neq_vw_type")
		{
			double alpha_neq = neq_potential.GetParameter("alpha");
			
			fPotential[i+1] = new VWType(alpha_neq, temp1);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
			
		}
		else if (neq_potential.Name() == "neq_lanir_fiber_recruit")
		{
			double K_neq = neq_potential.GetParameter("fiber_stiffness_K");

			fPotential[i+1] = new LanirFiber(K_neq, temp1, temp2);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
		}
		
		const ParameterListT& visc = list.GetListChoice(*this, "viscosity", i);
		if (visc.Name() == "scaled_csch")
		{
			double a = visc.GetParameter("eta0");
			double b = visc.GetParameter("tau0");
			fViscosity[i] = new ScaledCsch(a,b);
			if (!fViscosity[i]) throw ExceptionT::kOutOfMemory;		
		}		
	}
	
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fNumFibStress = dSymMatrixT::NumValues(fNumSD-1);
	fNumFibModuli = dSymMatrixT::NumValues(fNumFibStress);
	fFiberStretch.Dimension(fNumSD-1);
	fFiberStress.Dimension(fNumSD-1);
	fFiberMod.Dimension(fNumFibStress);

	/*viscous fiber stretch at time step n, viscous stretch at time step n and vn*/
	fFiberStretch_v.Dimension(fNumSD-1);
	fFiberStretch_vn.Dimension(fNumSD-1);	

	/*Dimension work spaces*/
	fInverse.Dimension(fNumSD);
	fCalg.Dimension(fNumFibStress);	
	/* allocate memory */
	/*dimension invserse viscosity matrix*/
	fiVisc.Dimension(fNumFibStress);

	/*Dimension work spaces*/
	fFlowStress.Dimension(fNumSD-1);
	fResidual.Dimension(fNumSD-1);
	fiK.Dimension(fNumFibStress);
	fG.Dimension(fNumFibStress);
	fMod1.Dimension(fNumFibStress);
	fMod2.Dimension(fNumFibStress);
	fMod3.Dimension(fNumFibStress);
	fVec.Dimension(fNumFibStress);

	/*dimension state variable storage arrays*/
	fNumMatProcess = 0;
	int numprocess = fNumFibProcess;
	SetStateVariables(numprocess);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	fUserFile = list.GetParameter("Dispersion_file");
	
	Construct();
}
	
/***********************************************************************
 * Protected
 ***********************************************************************/
void AnisoScleraWaxs_expo::ComputeMatrixStress(const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, 
				dSymMatrixT& Stress, const int process_index, const int fillmode)
{
	if (process_index > -1)
		ExceptionT::GeneralFail("AnisoScleraWaxs_expo::ComputeMatrixStress", 
			"Material model does not allow for neq matrix description");
	if (fillmode != dSymMatrixT::kOverwrite)
		ExceptionT::GeneralFail("AnisoScleraWaxs_expo::ComputeMatrixStress", 
			"Expects to overwrite Stress");
			

	if (fVolType == kOgden)
	{
		double I1 = Stretch[0]+Stretch[1]+Stretch[2];
		double I3 = Stretch.Det();

		double I3rthird = pow(I3, -third);

		fInverse.Inverse(Stretch);
	
		double coeff = -third*fMu*I3rthird*I1 +0.5*fGamma*(I3-1);
		Stress = fInverse;
		Stress *= coeff;
	
		Stress[0] += fMu*I3rthird;
		Stress[1] += fMu*I3rthird;
		Stress[2] += fMu*I3rthird;
	}
	else if (fVolType == kBlatz)
	{
		/*2pdf{W}{C_IJ} = mu ( del_IJ - I3^-gamma C^-1_IJ)*/
		double I3 = Stretch.Det();
		//cout << "\nstretch_ACV:"<< Stretch;
		//cout <<"\ndet:"<< I3;
		double I3rg = pow(I3, -fGamma);
		//cout<<"\ndetgamme:"<<I3rg;
		Stress.Inverse(Stretch);
		Stress *= -I3rg*fMu;
	
		Stress[0] += fMu;
		Stress[1] += fMu;
		Stress[2] += fMu;
		//cout <<"\nmatrix stress"<< Stress;
	}

}

void AnisoScleraWaxs_expo::ComputeMatrixMod(const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress,
				dMatrixT& Mod, const int process_index, const int fillmode)
{
	if (process_index > -1)
		ExceptionT::GeneralFail("AnisoScleraWaxs_expo::ComputeMatrixMod", 
			"Material model does not allow for neq matrix description");
	if (fillmode != dSymMatrixT::kOverwrite)
		ExceptionT::GeneralFail("AnisoScleraWaxs_expo::ComputeMatrixMod", 
			"Expects to overwrite Stress");

	if (fVolType == kOgden)
	{
		//	matrix stress
		double I1 = Stretch[0]+Stretch[1]+Stretch[2];
		double I3 = Stretch.Det();
		double I3rthird = pow(I3, -third);
		fInverse.Inverse(Stretch);
	
		double coeff = -third*fMu*I3rthird*I1 +0.5*fGamma*(I3-1);
		Stress = fInverse;
		Stress *= coeff;
	
		Stress[0] += fMu*I3rthird;
		Stress[1] += fMu*I3rthird;
		Stress[2] += fMu*I3rthird;
	
		//	modulus
		coeff = 2.0*third*fMu*I3rthird*I1 - fGamma*(I3-1.0);
	
		Mod.ReducedI_C(fInverse);
		Mod *= coeff;
	
		coeff = -2.0*third*fMu*I3rthird;
		Mod(0,0) += coeff*(fInverse[0] + fInverse[0]);
		Mod(0,1) += coeff*(fInverse[1] + fInverse[0]);
		Mod(0,2) += coeff*(fInverse[2] + fInverse[0]);
		Mod(0,3) += coeff*fInverse[3];
		Mod(0,4) += coeff*fInverse[4];
		Mod(0,5) += coeff*fInverse[5];
	
		Mod(1,0) += coeff*(fInverse[0] + fInverse[1]);
		Mod(1,1) += coeff*(fInverse[1] + fInverse[1]);
		Mod(1,2) += coeff*(fInverse[2] + fInverse[1]);
		Mod(1,3) += coeff*fInverse[3];
		Mod(1,4) += coeff*fInverse[4];
		Mod(1,5) += coeff*fInverse[5];
	
		Mod(2,0) += coeff*(fInverse[0] + fInverse[2]);
		Mod(2,1) += coeff*(fInverse[1] + fInverse[2]);
		Mod(2,2) += coeff*(fInverse[2] + fInverse[2]);
		Mod(2,3) += coeff*fInverse[3];
		Mod(2,4) += coeff*fInverse[4];
		Mod(2,5) += coeff*fInverse[5];

		Mod(3,0) += coeff*(fInverse[3]);
		Mod(3,1) += coeff*(fInverse[3]);
		Mod(3,2) += coeff*(fInverse[3]);
	
		Mod(4,0) += coeff*(fInverse[4]);
		Mod(4,1) += coeff*(fInverse[4]);
		Mod(4,2) += coeff*(fInverse[4]);
	
		Mod(5,0) += coeff*(fInverse[5]);
		Mod(5,1) += coeff*(fInverse[5]);
		Mod(5,2) += coeff*(fInverse[5]);
		
		coeff = 2.0*third*third*fMu*I3rthird*I1 + fGamma*I3;
		Mod(0,0) += coeff*(fInverse[0]*fInverse[0]);
		Mod(0,1) += coeff*(fInverse[0]*fInverse[1]);
		Mod(0,2) += coeff*(fInverse[0]*fInverse[2]);
		Mod(0,3) += coeff*(fInverse[0]*fInverse[3]);
		Mod(0,4) += coeff*(fInverse[0]*fInverse[4]);
		Mod(0,5) += coeff*(fInverse[0]*fInverse[5]);
		//
		Mod(1,0) += coeff*(fInverse[1]*fInverse[0]);
		Mod(1,1) += coeff*(fInverse[1]*fInverse[1]);
		Mod(1,2) += coeff*(fInverse[1]*fInverse[2]);
		Mod(1,3) += coeff*(fInverse[1]*fInverse[3]);
		Mod(1,4) += coeff*(fInverse[1]*fInverse[4]);
		Mod(1,5) += coeff*(fInverse[1]*fInverse[5]);
		//
		Mod(2,0) += coeff*(fInverse[2]*fInverse[0]);
		Mod(2,1) += coeff*(fInverse[2]*fInverse[1]);
		Mod(2,2) += coeff*(fInverse[2]*fInverse[2]);
		Mod(2,3) += coeff*(fInverse[2]*fInverse[3]);
		Mod(2,4) += coeff*(fInverse[2]*fInverse[4]);
		Mod(2,5) += coeff*(fInverse[2]*fInverse[5]);
		//
		Mod(3,0) += coeff*(fInverse[3]*fInverse[0]);
		Mod(3,1) += coeff*(fInverse[3]*fInverse[1]);
		Mod(3,2) += coeff*(fInverse[3]*fInverse[2]);
		Mod(3,3) += coeff*(fInverse[3]*fInverse[3]);
		Mod(3,4) += coeff*(fInverse[3]*fInverse[4]);
		Mod(3,5) += coeff*(fInverse[3]*fInverse[5]);
		//
		Mod(4,0) += coeff*(fInverse[4]*fInverse[0]);
		Mod(4,1) += coeff*(fInverse[4]*fInverse[1]);
		Mod(4,2) += coeff*(fInverse[4]*fInverse[2]);
		Mod(4,3) += coeff*(fInverse[4]*fInverse[3]);
		Mod(4,4) += coeff*(fInverse[4]*fInverse[4]);
		Mod(4,5) += coeff*(fInverse[4]*fInverse[5]);
		//
		Mod(5,0) += coeff*(fInverse[5]*fInverse[0]);
		Mod(5,1) += coeff*(fInverse[5]*fInverse[1]);
		Mod(5,2) += coeff*(fInverse[5]*fInverse[2]);
		Mod(5,3) += coeff*(fInverse[5]*fInverse[3]);
		Mod(5,4) += coeff*(fInverse[5]*fInverse[4]);
		Mod(5,5) += coeff*(fInverse[5]*fInverse[5]);
	}
	else if(fVolType==kBlatz)
	{
		//	matrix contribution
		double I3 = Stretch.Det();
		double I3rg = pow(I3, -fGamma);
		fInverse.Inverse(Stretch);
		fStress = fInverse;
		Stress *= -I3rg*fMu;
	
		Stress[0] += fMu;
		Stress[1] += fMu;
		Stress[2] += fMu;

	
		//	2pdf{S_IJ}{C_KL} = 2 mu I_3^-gamma (gamma C^-1_IJ C^-1_KL + 0.5(C^-1_IK C^-1_JL +C^-1_IL+C^-1_JK)
		//	modulus
		double coeff = 2.0*fMu*I3rg;
		Mod.ReducedI_C(fInverse);
		Mod *= coeff;
	
		coeff *= fGamma;
		Mod(0,0) += coeff*fInverse[0]*fInverse[0];
		Mod(0,1) += coeff*fInverse[0]*fInverse[1];
		Mod(0,2) += coeff*fInverse[0]*fInverse[2];
		Mod(0,3) += coeff*fInverse[0]*fInverse[3];
		Mod(0,4) += coeff*fInverse[0]*fInverse[4];
		Mod(0,5) += coeff*fInverse[0]*fInverse[5];
	
		Mod(1,0) += coeff*fInverse[1]*fInverse[0];
		Mod(1,1) += coeff*fInverse[1]*fInverse[1];
		Mod(1,2) += coeff*fInverse[1]*fInverse[2];
		Mod(1,3) += coeff*fInverse[1]*fInverse[3];
		Mod(1,4) += coeff*fInverse[1]*fInverse[4];
		Mod(1,5) += coeff*fInverse[1]*fInverse[5];

		Mod(2,0) += coeff*fInverse[2]*fInverse[0];
		Mod(2,1) += coeff*fInverse[2]*fInverse[1];
		Mod(2,2) += coeff*fInverse[2]*fInverse[2];
		Mod(2,3) += coeff*fInverse[2]*fInverse[3];
		Mod(2,4) += coeff*fInverse[2]*fInverse[4];
		Mod(2,5) += coeff*fInverse[2]*fInverse[5];

		Mod(3,0) += coeff*fInverse[3]*fInverse[0];
		Mod(3,1) += coeff*fInverse[3]*fInverse[1];
		Mod(3,2) += coeff*fInverse[3]*fInverse[2];
		Mod(3,3) += coeff*fInverse[3]*fInverse[3];
		Mod(3,4) += coeff*fInverse[3]*fInverse[4];
		Mod(3,5) += coeff*fInverse[3]*fInverse[5];

		Mod(4,0) += coeff*fInverse[4]*fInverse[0];
		Mod(4,1) += coeff*fInverse[4]*fInverse[1];
		Mod(4,2) += coeff*fInverse[4]*fInverse[2];
		Mod(4,3) += coeff*fInverse[4]*fInverse[3];
		Mod(4,4) += coeff*fInverse[4]*fInverse[4];
		Mod(4,5) += coeff*fInverse[4]*fInverse[5];

		Mod(5,0) += coeff*fInverse[5]*fInverse[0];
		Mod(5,1) += coeff*fInverse[5]*fInverse[1];
		Mod(5,2) += coeff*fInverse[5]*fInverse[2];
		Mod(5,3) += coeff*fInverse[5]*fInverse[3];
		Mod(5,4) += coeff*fInverse[5]*fInverse[4];
		Mod(5,5) += coeff*fInverse[5]*fInverse[5];
	}
 }

/*computes integrated fiber stress in local frame*/
void AnisoScleraWaxs_expo::ComputeFiberStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
			dSymMatrixT& FiberStress, const int pindex)
{

	/*initialize pointers*/
	/* PK2 values in local frame formed by NT and IS orientations*/	
	FiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_12*/

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);

	if (pindex == -1) { 
		CompI4(FiberStretch);

		/* derivatives of the potential */
//		cout << "\nip: "<<CurrIP();
//		cout << "\nfI4: "<<fI4;
		fPotential[pindex+1]->MapDFunction(fI4, fdU);
//		cout << "\nfdU: "<<fdU<<endl;
		
		/* initialize kernel pointers */
		double* pdU = fdU.Pointer();
		double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	
		/*integrate w.r.t in-plane orientation theta*/
		for (int i = 0; i < fI4.Length(); i++)
		{
			double factor = 2.0*(*pj++)*(*pdU++);
			s1 += factor*(*p1++);
			s2 += factor*(*p2++);
			s3 += factor*(*p3++);
		}
		
	}
	else {
		/*calculate I4e and I4v*/
		CompI4(FiberStretch, FiberStretch_v);
		dArrayT& I4e = fI4;
		I4e /= fI4v; /*fI4 = I4e*/

		/* derivatives of the potential */
		fPotential[pindex+1]->MapDFunction(I4e, fdU);

		/* initialize kernel pointers */
		double* pdU = fdU.Pointer();
		double* plv = fI4v.Pointer();
		double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	
		/*integrate w.r.t in-plane orientation theta*/
		for (int i = 0; i < I4e.Length(); i++)
		{
			double factor = 2.0*(*pj++)*(*pdU++)/(*plv++);
			s1 += factor*(*p1++);
			s2 += factor*(*p2++);
			s3 += factor*(*p3++);
		}
	}
}
	
/*computes integrated moduli in local frame*/
void AnisoScleraWaxs_expo::ComputeFiberMod (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dSymMatrixT& FiberStress, dMatrixT& FiberMod,  const int pindex)
{
	/*initialize pointers*/
	/* PK2 and Modulus values in local coordinate frame formed by a1 (fNT) and a2 (fIS) */	
	FiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_12*/
 
	FiberMod = 0.0;
	double& c11 = FiberMod[0]; /*cf_1111*/ 
	double& c22 = FiberMod[4]; /*cf_2222*/
	double& c33 = FiberMod[8]; /*cf_1212*/
	double& c23 = FiberMod[7]; /*cf_2212*/
	double& c13 = FiberMod[6]; /*cf_1112*/
	double& c12 = FiberMod[3]; /*cf_1122*/

	/*0  3  6
	  1  4  7
	  2  5  8*/

	/* stress */
	double* ps1  = fStressTable(0);
	double* ps2  = fStressTable(1);
	double* ps3  = fStressTable(2);

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);

	if (pindex == -1) 
	{ 
		CompI4(FiberStretch);

		/* derivatives of the potential */
		fPotential[pindex+1]->MapDFunction(fI4, fdU);
		fPotential[pindex+1]->MapDDFunction(fI4, fddU);	

		/* initialize kernel pointers */	
		double* pdU  = fdU.Pointer();
		double* pddU = fddU.Pointer();
		double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

		
		for (int i = 0; i < fI4.Length(); i++)
		{
			double sfactor =  2.0*(*pj)*(*pdU++);
			double cfactor = 4.0*(*pj++)*(*pddU++);

			s1 += sfactor*(*ps1++);
			s2 += sfactor*(*ps2++);
			s3 += sfactor*(*ps3++);

			c11 += cfactor*(*pc11++);
			c22 += cfactor*(*pc22++);
			c33 += cfactor*(*pc33++);

			c23 += cfactor*(*pc23++);
			c13 += cfactor*(*pc13++);
			c12 += cfactor*(*pc12++);
			
		}
		
		/*symmetric modulus*/
		FiberMod[1] = c12;
		FiberMod[2] = c13;
		FiberMod[5] = c23;
	}
	else 
	{
		/*calculate I4e and I4v*/
		CompI4(FiberStretch, FiberStretch_v);
		dArrayT& I4e = fI4;
		I4e /= fI4v;  /*fI4=I4e*/

		/* derivatives of the potential */
		fPotential[pindex+1]->MapDFunction(I4e, fdU);
		fPotential[pindex+1]->MapDDFunction(I4e, fddU);	

		/* initialize kernel pointers */	
		double* pdU  = fdU.Pointer();
		double* pddU = fddU.Pointer();
		double* plv = fI4v.Pointer();
		double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

		for (int i = 0; i < I4e.Length(); i++)
		{
			double sfactor =  2.0*(*pj)*(*pdU++)/(*plv);
			double cfactor = 4.0*(*pj++)*(*pddU++)/((*plv)*(*plv));
			plv ++;

			s1 += sfactor*(*ps1++);
			s2 += sfactor*(*ps2++);
			s3 += sfactor*(*ps3++);

			c11 += cfactor*(*pc11++);
			c22 += cfactor*(*pc22++);
			c33 += cfactor*(*pc33++);
			c23 += cfactor*(*pc23++);
			c13 += cfactor*(*pc13++);
			c12 += cfactor*(*pc12++);
		}
		/*symmetric modulus*/
		FiberMod[1] = c12;
		FiberMod[2] = c13;
		FiberMod[5] = c23;
		
		ComputeCalg(FiberStretch, FiberStretch_v, fCalg, pindex);
		FiberMod += fCalg;
	}
	
}

void AnisoScleraWaxs_expo::ComputeCalg(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,  dMatrixT& Calg, const int pindex)
{
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();

	/*Compute viscosity*/
	ComputeViscosity(FiberStretch, FiberStretch_v, fiVisc, pindex);
	fiVisc(0,2) *= 2.0;
	fiVisc(1,2) *= 2.0;
	fiVisc(2,2) *= 2.0;
	fiVisc.Inverse(); 
	/*compute Kij = del_ij - dt*V^-1_ik 2dS_k/dCfvt_j*/
	/*Cfvt = {Cfv11, Cfv22, 2 Cfv12}*/
	/*stress term*/
	dFlowdCv(FiberStretch, FiberStretch_v, fMod1, pindex);
	fMod1.ToMatrix(fMod2);
		
	fMod2(0,2) *= 2.0;
	fMod2(1,2) *= 2.0;
	fMod2(2,2) *= 2.0;

	fiK(0,0) = 1.0 - dt*(fiVisc(0,0)*fMod2(0,0)+fiVisc(0,1)*fMod2(1,0)+fiVisc(0,2)*fMod2(2,0));
	fiK(1,1) = 1.0 - dt*(fiVisc(1,0)*fMod2(0,1)+fiVisc(1,1)*fMod2(1,1)+fiVisc(1,2)*fMod2(2,1));
	fiK(2,2) = 1.0 - dt*(fiVisc(2,0)*fMod2(0,2)+fiVisc(2,1)*fMod2(1,2)+fiVisc(2,2)*fMod2(2,2));

	fiK(0,1) = -dt*(fiVisc(0,0)*fMod2(0,1)+fiVisc(0,1)*fMod2(1,1)+fiVisc(0,2)*fMod2(2,1));
	fiK(0,2) = -dt*(fiVisc(0,0)*fMod2(0,2)+fiVisc(0,1)*fMod2(1,2)+fiVisc(0,2)*fMod2(2,2));
	fiK(1,2) = -dt*(fiVisc(1,0)*fMod2(0,2)+fiVisc(1,1)*fMod2(1,2)+fiVisc(1,2)*fMod2(2,2));

	fiK(1,0) = -dt*(fiVisc(1,0)*fMod2(0,0)+fiVisc(1,1)*fMod2(1,0)+fiVisc(1,2)*fMod2(2,0));
	fiK(2,0) = -dt*(fiVisc(2,0)*fMod2(0,0)+fiVisc(2,1)*fMod2(1,0)+fiVisc(2,2)*fMod2(2,0));
	fiK(2,1) = -dt*(fiVisc(2,0)*fMod2(0,1)+fiVisc(2,1)*fMod2(1,1)+fiVisc(2,2)*fMod2(2,1));

	/*viscosity terms*/
	/*Compute flow stress*/
	ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);
	fVec[0] = (fiVisc(0,0)*fFlowStress[0] + fiVisc(0,1)*fFlowStress[1] + fiVisc(0,2)*fFlowStress[2]);
	fVec[1] = (fiVisc(1,0)*fFlowStress[0] + fiVisc(1,1)*fFlowStress[1] + fiVisc(1,2)*fFlowStress[2]);
	fVec[2] = (fiVisc(2,0)*fFlowStress[0] + fiVisc(2,1)*fFlowStress[1] + fiVisc(2,2)*fFlowStress[2]);

	ComputeDViscDCv(FiberStretch, FiberStretch_v, fVec, fMod2, pindex);
	fMod2(0,2) *= 2.0;
	fMod2(1,2) *= 2.0;
	fMod2(2,2) *= 2.0;

	fiK(0,0) += dt*(fiVisc(0,0)*fMod2(0,0)+fiVisc(0,1)*fMod2(1,0)+fiVisc(0,2)*fMod2(2,0));
	fiK(1,1) += dt*(fiVisc(1,0)*fMod2(0,1)+fiVisc(1,1)*fMod2(1,1)+fiVisc(1,2)*fMod2(2,1));
	fiK(2,2) += dt*(fiVisc(2,0)*fMod2(0,2)+fiVisc(2,1)*fMod2(1,2)+fiVisc(2,2)*fMod2(2,2));

	fiK(0,1) += dt*(fiVisc(0,0)*fMod2(0,1)+fiVisc(0,1)*fMod2(1,1)+fiVisc(0,2)*fMod2(2,1));
	fiK(0,2) += dt*(fiVisc(0,0)*fMod2(0,2)+fiVisc(0,1)*fMod2(1,2)+fiVisc(0,2)*fMod2(2,2));
	fiK(1,2) += dt*(fiVisc(1,0)*fMod2(0,2)+fiVisc(1,1)*fMod2(1,2)+fiVisc(1,2)*fMod2(2,2));

	fiK(1,0) += dt*(fiVisc(1,0)*fMod2(0,0)+fiVisc(1,1)*fMod2(1,0)+fiVisc(1,2)*fMod2(2,0));
	fiK(2,0) += dt*(fiVisc(2,0)*fMod2(0,0)+fiVisc(2,1)*fMod2(1,0)+fiVisc(2,2)*fMod2(2,0));
	fiK(2,1) += dt*(fiVisc(2,0)*fMod2(0,1)+fiVisc(2,1)*fMod2(1,1)+fiVisc(2,2)*fMod2(2,1));
	
	fiK.Inverse();

	/*compute G= dSNEQ/dCfv. Note for this model, dFlowStress/dC = - dSNEQ/dCfv*/
	/*stress term*/
	dFlowdC(FiberStretch, FiberStretch_v, fMod1,  pindex);
	fMod1.ToMatrix(fMod2);
	fMod2(0,2) *= 2.0;
	fMod2(1,2) *= 2.0;
	fMod2(2,2) *= 2.0;

	/*compute dt V^-1_KL 2dFlow_L/dCft_J*/
	fG(0,0) = dt*(fiVisc(0,0)*fMod2(0,0)+fiVisc(0,1)*fMod2(1,0)+fiVisc(0,2)*fMod2(2,0));
	fG(1,1) = dt*(fiVisc(1,0)*fMod2(0,1)+fiVisc(1,1)*fMod2(1,1)+fiVisc(1,2)*fMod2(2,1));
	fG(2,2) = dt*(fiVisc(2,0)*fMod2(0,2)+fiVisc(2,1)*fMod2(1,2)+fiVisc(2,2)*fMod2(2,2));

	fG(0,1) = dt*(fiVisc(0,0)*fMod2(0,1)+fiVisc(0,1)*fMod2(1,1)+fiVisc(0,2)*fMod2(2,1));
	fG(0,2) = dt*(fiVisc(0,0)*fMod2(0,2)+fiVisc(0,1)*fMod2(1,2)+fiVisc(0,2)*fMod2(2,2));
	fG(1,2) = dt*(fiVisc(1,0)*fMod2(0,2)+fiVisc(1,1)*fMod2(1,2)+fiVisc(1,2)*fMod2(2,2));

	fG(1,0) = dt*(fiVisc(1,0)*fMod2(0,0)+fiVisc(1,1)*fMod2(1,0)+fiVisc(1,2)*fMod2(2,0));
	fG(2,0) = dt*(fiVisc(2,0)*fMod2(0,0)+fiVisc(2,1)*fMod2(1,0)+fiVisc(2,2)*fMod2(2,0));
	fG(2,1) = dt*(fiVisc(2,0)*fMod2(0,1)+fiVisc(2,1)*fMod2(1,1)+fiVisc(2,2)*fMod2(2,1));

	/*viscosity terms*/
	ComputeDViscDC(FiberStretch, FiberStretch_v, fVec, fMod2, pindex);
	fMod2(0,2) *= 2.0;
	fMod2(1,2) *= 2.0;
	fMod2(2,2) *= 2.0;
	
	
	fG(0,0) += -dt*(fiVisc(0,0)*fMod2(0,0)+fiVisc(0,1)*fMod2(1,0)+fiVisc(0,2)*fMod2(2,0));
	fG(1,1) += -dt*(fiVisc(1,0)*fMod2(0,1)+fiVisc(1,1)*fMod2(1,1)+fiVisc(1,2)*fMod2(2,1));
	fG(2,2) += -dt*(fiVisc(2,0)*fMod2(0,2)+fiVisc(2,1)*fMod2(1,2)+fiVisc(2,2)*fMod2(2,2));

	fG(0,1) += -dt*(fiVisc(0,0)*fMod2(0,1)+fiVisc(0,1)*fMod2(1,1)+fiVisc(0,2)*fMod2(2,1));
	fG(0,2) += -dt*(fiVisc(0,0)*fMod2(0,2)+fiVisc(0,1)*fMod2(1,2)+fiVisc(0,2)*fMod2(2,2));
	fG(1,2) += -dt*(fiVisc(1,0)*fMod2(0,2)+fiVisc(1,1)*fMod2(1,2)+fiVisc(1,2)*fMod2(2,2));

	fG(1,0) += -dt*(fiVisc(1,0)*fMod2(0,0)+fiVisc(1,1)*fMod2(1,0)+fiVisc(1,2)*fMod2(2,0));
	fG(2,0) += -dt*(fiVisc(2,0)*fMod2(0,0)+fiVisc(2,1)*fMod2(1,0)+fiVisc(2,2)*fMod2(2,0));
	fG(2,1) += -dt*(fiVisc(2,0)*fMod2(0,1)+fiVisc(2,1)*fMod2(1,1)+fiVisc(2,2)*fMod2(2,1));
		
	/*compute Calg = dSNEQ_I/dCfv.(K^-1.G) */
	dFlowdC(FiberStretch, FiberStretch_v, fMod1,  pindex);

	fMod1.ToMatrix(fMod3);
	fMod3(0,2) *= 2.0;
	fMod3(1,2) *= 2.0;
	fMod3(2,2) *= 2.0;

	fMod2(0,0) = -(fMod3(0,0)*fiK(0,0)+fMod3(0,1)*fiK(1,0)+fMod3(0,2)*fiK(2,0));
	fMod2(1,1) = -(fMod3(1,0)*fiK(0,1)+fMod3(1,1)*fiK(1,1)+fMod3(1,2)*fiK(2,1));
	fMod2(2,2) = -(fMod3(2,0)*fiK(0,2)+fMod3(2,1)*fiK(1,2)+fMod3(2,2)*fiK(2,2));

	fMod2(0,1) = -(fMod3(0,0)*fiK(0,1)+fMod3(0,1)*fiK(1,1)+fMod3(0,2)*fiK(2,1));
	fMod2(0,2) = -(fMod3(0,0)*fiK(0,2)+fMod3(0,1)*fiK(1,2)+fMod3(0,2)*fiK(2,2));
	fMod2(1,2) = -(fMod3(1,0)*fiK(0,2)+fMod3(1,1)*fiK(1,2)+fMod3(1,2)*fiK(2,2));

	fMod2(1,0) = -(fMod3(1,0)*fiK(0,0)+fMod3(1,1)*fiK(1,0)+fMod3(1,2)*fiK(2,0));
	fMod2(2,0) = -(fMod3(2,0)*fiK(0,0)+fMod3(2,1)*fiK(1,0)+fMod3(2,2)*fiK(2,0));
	fMod2(2,1) = -(fMod3(2,0)*fiK(0,1)+fMod3(2,1)*fiK(1,1)+fMod3(2,2)*fiK(2,1));
		
	Calg.MultAB(fMod2, fG);
	/*convert Calg back to matrix system*/
	Calg(0,2) *= 0.5;
	Calg(1,2) *= 0.5;
	Calg(2,2) *= 0.5; 

}

/*local newton loop for viscous stretch tensor*/ 
void AnisoScleraWaxs_expo::Compute_Cv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_vn, 
	dSymMatrixT& FiberStretch_v, const int pindex)
{
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
			
	double error;
	int iteration  = 0;

	/*compute  viscosity  based on currn values of stretch matrices, C and Cv*/
	ComputeViscosity(FiberStretch, FiberStretch_v, fiVisc, pindex);
	fiVisc(0,2) *= 2.0;
	fiVisc(1,2) *= 2.0;
	fiVisc(2,2) *= 2.0;	
	fiVisc.Inverse(); 

	/*Compute flow stress*/
	ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);

	fResidual[0] = FiberStretch_v[0] - 2.0*dt*(fiVisc(0,0)*fFlowStress[0] + fiVisc(0,1)*fFlowStress[1] 
			+ fiVisc(0,2)*fFlowStress[2]) - FiberStretch_vn[0];
	fResidual[1] = FiberStretch_v[1] - 2.0*dt*(fiVisc(1,0)*fFlowStress[0] + fiVisc(1,1)*fFlowStress[1] 
			+ fiVisc(1,2)*fFlowStress[2]) - FiberStretch_vn[1];
	fResidual[2] = FiberStretch_v[2] - 2.0*dt*(fiVisc(2,0)*fFlowStress[0] + fiVisc(2,1)*fFlowStress[1] 
			+ fiVisc(2,2)*fFlowStress[2]) - FiberStretch_vn[2];
	
	do 
	{
		/*compute Kij = del_ij - dt*V^-1_ik 2dS_k/dCfvt_j*/
		/*Cfvt = {Cfv11, Cfv22, 2 Cfv12}*/
		/*stress term*/
		dFlowdCv(FiberStretch, FiberStretch_v, fMod1, pindex);
		fMod1.ToMatrix(fMod2);
		
		fMod2(0,2) *= 2.0;
		fMod2(1,2) *= 2.0;
		fMod2(2,2) *= 2.0;
		
		fiK(0,0) = 1.0 - dt*(fiVisc(0,0)*fMod2(0,0)+fiVisc(0,1)*fMod2(1,0)+fiVisc(0,2)*fMod2(2,0));
		fiK(1,1) = 1.0 - dt*(fiVisc(1,0)*fMod2(0,1)+fiVisc(1,1)*fMod2(1,1)+fiVisc(1,2)*fMod2(2,1));
		fiK(2,2) = 1.0 - dt*(fiVisc(2,0)*fMod2(0,2)+fiVisc(2,1)*fMod2(1,2)+fiVisc(2,2)*fMod2(2,2));

		fiK(0,1) = -dt*(fiVisc(0,0)*fMod2(0,1)+fiVisc(0,1)*fMod2(1,1)+fiVisc(0,2)*fMod2(2,1));
		fiK(0,2) = -dt*(fiVisc(0,0)*fMod2(0,2)+fiVisc(0,1)*fMod2(1,2)+fiVisc(0,2)*fMod2(2,2));
		fiK(1,2) = -dt*(fiVisc(1,0)*fMod2(0,2)+fiVisc(1,1)*fMod2(1,2)+fiVisc(1,2)*fMod2(2,2));

		fiK(1,0) = -dt*(fiVisc(1,0)*fMod2(0,0)+fiVisc(1,1)*fMod2(1,0)+fiVisc(1,2)*fMod2(2,0));
		fiK(2,0) = -dt*(fiVisc(2,0)*fMod2(0,0)+fiVisc(2,1)*fMod2(1,0)+fiVisc(2,2)*fMod2(2,0));
		fiK(2,1) = -dt*(fiVisc(2,0)*fMod2(0,1)+fiVisc(2,1)*fMod2(1,1)+fiVisc(2,2)*fMod2(2,1));


		/*viscosity terms*/
		/*V^-1.Sig*/
		fVec[0] = (fiVisc(0,0)*fFlowStress[0] + fiVisc(0,1)*fFlowStress[1] + fiVisc(0,2)*fFlowStress[2]);
		fVec[1] = (fiVisc(1,0)*fFlowStress[0] + fiVisc(1,1)*fFlowStress[1] + fiVisc(1,2)*fFlowStress[2]);
		fVec[2] = (fiVisc(2,0)*fFlowStress[0] + fiVisc(2,1)*fFlowStress[1] + fiVisc(2,2)*fFlowStress[2]);

		/*2 d(V)/dC.V^-1.Sig*/
		ComputeDViscDCv(FiberStretch, FiberStretch_v, fVec, fMod2, pindex);
		fMod2(0,2) *= 2.0;
		fMod2(1,2) *= 2.0;
		fMod2(2,2) *= 2.0;

		fiK(0,0) += dt*(fiVisc(0,0)*fMod2(0,0)+fiVisc(0,1)*fMod2(1,0)+fiVisc(0,2)*fMod2(2,0));
		fiK(1,1) += dt*(fiVisc(1,0)*fMod2(0,1)+fiVisc(1,1)*fMod2(1,1)+fiVisc(1,2)*fMod2(2,1));
		fiK(2,2) += dt*(fiVisc(2,0)*fMod2(0,2)+fiVisc(2,1)*fMod2(1,2)+fiVisc(2,2)*fMod2(2,2));

		fiK(0,1) += dt*(fiVisc(0,0)*fMod2(0,1)+fiVisc(0,1)*fMod2(1,1)+fiVisc(0,2)*fMod2(2,1));
		fiK(0,2) += dt*(fiVisc(0,0)*fMod2(0,2)+fiVisc(0,1)*fMod2(1,2)+fiVisc(0,2)*fMod2(2,2));
		fiK(1,2) += dt*(fiVisc(1,0)*fMod2(0,2)+fiVisc(1,1)*fMod2(1,2)+fiVisc(1,2)*fMod2(2,2));

		fiK(1,0) += dt*(fiVisc(1,0)*fMod2(0,0)+fiVisc(1,1)*fMod2(1,0)+fiVisc(1,2)*fMod2(2,0));
		fiK(2,0) += dt*(fiVisc(2,0)*fMod2(0,0)+fiVisc(2,1)*fMod2(1,0)+fiVisc(2,2)*fMod2(2,0));
		fiK(2,1) += dt*(fiVisc(2,0)*fMod2(0,1)+fiVisc(2,1)*fMod2(1,1)+fiVisc(2,2)*fMod2(2,1));
	
		fiK.Inverse();
				
		/*calculate increment of Cfvt*/
		double dCfvt0 = -(fiK(0,0)*fResidual[0] + fiK(0,1)*fResidual[1] + fiK(0,2)*fResidual[2]);
		double dCfvt1 = -(fiK(1,0)*fResidual[0] + fiK(1,1)*fResidual[1] + fiK(1,2)*fResidual[2]);
		double dCfvt2 = -(fiK(2,0)*fResidual[0] + fiK(2,1)*fResidual[1] + fiK(2,2)*fResidual[2]);
		
		/*calculate update*/
		FiberStretch_v[0] += dCfvt0;
		FiberStretch_v[1] += dCfvt1;
		FiberStretch_v[2] += dCfvt2;
		
//		cout << "\niteration: "<<iteration;
//		cout <<setprecision(12)<< "\nCf: "<<Cf;			
//		cout << "\nCfv: "<<Cfv;
//		cout << "\nflowstress: " << S;
//		cout << "\nresidual: " << R;
//		cout << "\n2dFdCfvt: "<<fMod1;
//		cout << "\nfiK: "<<fiK;
//		cout << "\ndCfv: "<<dCfv0<<"\t"<<dCfv1<<"\t"<<dCfv2;
//		cout << "\nCfv: "<<Cfv;

		/*Compute viscosity*/
		ComputeViscosity(FiberStretch, FiberStretch_v, fiVisc, pindex);
		fiVisc(0,2) *= 2.0;
		fiVisc(1,2) *= 2.0;
		fiVisc(2,2) *= 2.0;
		fiVisc.Inverse(); 

		/*compute flow stress based on current values of stretch matrices, C and Cv*/
		ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);
		
		/*residual at the zeroth iteration*/
		fResidual[0] = FiberStretch_v[0] - 2.0*dt*(fiVisc(0,0)*fFlowStress[0] + fiVisc(0,1)*fFlowStress[1] 
			+ fiVisc(0,2)*fFlowStress[2]) - FiberStretch_vn[0];
		fResidual[1] = FiberStretch_v[1] - 2.0*dt*(fiVisc(1,0)*fFlowStress[0] + fiVisc(1,1)*fFlowStress[1] 
			+ fiVisc(1,2)*fFlowStress[2]) - FiberStretch_vn[1];
		fResidual[2] = FiberStretch_v[2] - 2.0*dt*(fiVisc(2,0)*fFlowStress[0] + fiVisc(2,1)*fFlowStress[1] 
			+ fiVisc(2,2)*fFlowStress[2]) - FiberStretch_vn[2];
		
		error = sqrt(fResidual[0]*fResidual[0] + fResidual[1]*fResidual[1] + fResidual[2]*fResidual[2]);

//		cout <<"\niteration: "<<iteration;
//		cout<< "\nerror: "<<error;
	
		iteration++;
	}while (error > kSmall && iteration < 6);
	if (iteration >= 10) 
		ExceptionT::GeneralFail("AnisoScleraWaxs_expo::Compute_Cv", 
			"number of iteration exceeds maximum of 10");
}


/*computes integrated fiber stress in local frame*/
void AnisoScleraWaxs_expo::ComputeFlowStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dSymMatrixT& FlowStress, const int pindex)
{
	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT& I4e = fI4;
	I4e /= fI4v;
	/* derivatives of the potential */
	fPotential[pindex+1]->MapDFunction(I4e, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* ple  = I4e.Pointer();
	double* plv = fI4v.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);

	/* PK2 values in local frame formed by NT and IS orientations*/	
	FlowStress = 0.0;
	double& s1 = FlowStress[0]; /*sf_11*/
	double& s2 = FlowStress[1]; /*sf_22*/
	double& s3 = FlowStress[2]; /*sf_12*/
	
	/*integrate w.r.t in-plane orientation theta*/
	for (int i = 0; i < I4e.Length(); i++)
	{
		double factor = 2.0*(*pj++)*(*pdU++)*(*ple++)/(*plv++);
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}
}

/*computes dFlowStress/dC.  Note for this model, dFlowStress/dC = - dSNEQ/dCv*/
void AnisoScleraWaxs_expo::dFlowdC (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
		dSymMatrixT& FiberMod,  const int pindex)
{
	/*initialize pointers*/
	/* PK2 and Modulus values in local coordinate frame formed by a1 (fNT) and a2 (fIS) */	
	 
	FiberMod = 0.0;
	double& c11 = FiberMod[0]; /*cf_1111*/
	double& c22 = FiberMod[1]; /*cf_2222*/
	double& c33 = FiberMod[2]; /*cf_1212*/
	double& c23 = FiberMod[3]; /*cf_2212*/
	double& c13 = FiberMod[4]; /*cf_1112*/
	double& c12 = FiberMod[5]; /*cf_1122*/

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);

	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT I4e = fI4;
	I4e /= fI4v;  /*fI4=I4e*/

	/* derivatives of the potential */
	fPotential[pindex+1]->MapDFunction(I4e, fdU);
	fPotential[pindex+1]->MapDDFunction(I4e, fddU);	

	/* initialize kernel pointers */	
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* ple = I4e.Pointer();
	double* plv = fI4v.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	for (int i = 0; i < I4e.Length(); i++)
	{
		double cfactor = 4.0*(*pj++)*( (*pddU++)*(*ple)/((*plv)*(*plv)) + (*pdU++)/((*plv)*(*plv)));
		plv ++;
		ple ++;

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c33 += cfactor*(*pc33++);
		c23 += cfactor*(*pc23++);
		c13 += cfactor*(*pc13++);
		c12 += cfactor*(*pc12++);
	}
}

/*computes dFlowStress/dCv*/
void AnisoScleraWaxs_expo::dFlowdCv (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
		dSymMatrixT& FiberMod,  const int pindex)
{
	/*initialize pointers*/
	/* PK2 and Modulus values in local coordinate frame formed by a1 (fNT) and a2 (fIS) */	
	 
	FiberMod = 0.0;
	double& c11 = FiberMod[0]; /*cf_1111*/
	double& c22 = FiberMod[1]; /*cf_2222*/
	double& c33 = FiberMod[2]; /*cf_1212*/
	double& c23 = FiberMod[3]; /*cf_2212*/
	double& c13 = FiberMod[4]; /*cf_1112*/
	double& c12 = FiberMod[5]; /*cf_1122*/

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);

	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT& I4e = fI4;
	I4e /= fI4v;  /*fI4=I4e*/

	/* derivatives of the potential */
	fPotential[pindex+1]->MapDFunction(I4e, fdU);
	fPotential[pindex+1]->MapDDFunction(I4e, fddU);	

	/* initialize kernel pointers */	
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* ple = I4e.Pointer();
	double* plv = fI4v.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();


	for (int i = 0; i < I4e.Length(); i++)
	{
		double cfactor = -4.0*(*pj++)*( (*pddU++)*(*ple)*(*ple)/((*plv)*(*plv)) 
											+ 2.0*(*pdU++)*(*ple)/((*plv)*(*plv)));
		plv ++;
		ple ++;
	
		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c33 += cfactor*(*pc33++);
		c23 += cfactor*(*pc23++);
		c13 += cfactor*(*pc13++);
		c12 += cfactor*(*pc12++);
	}
}

void AnisoScleraWaxs_expo::ComputeViscosity(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,  
			dMatrixT& Visc,  const int pindex)
{
	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT& I4e = fI4;
	I4e /= fI4v;
	
	double* ple   = I4e.Pointer();
	double* plv  = fI4v.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);
	
	Visc = 0.0;
	double& v11 = Visc(0,0); /*cf_1111*/
	double& v22 = Visc(1,1); /*cf_2222*/
	double& v33 = Visc(2,2); /*cf_1212*/
	double& v23 = Visc(1,2); /*cf_2212*/
	double& v13 = Visc(0,2); /*cf_1112*/
	double& v12 = Visc(0,1); /*cf_1122*/

	for (int i = 0; i < I4e.Length(); i++)
	{
//		cout << "\ncalculating viscosity: ";
//		cout << "\nI4e: "<<*pl;
		double tau = fPotential[pindex+1]->DFunction(*ple);
		tau *= 2.0*(*ple++);
		double visc = fViscosity[pindex]->Function(tau);
//		cout << "\ntau: "<<tau;
//		cout << "\nvisc: "<<visc;
		double vfactor =  (*pj++)*visc/((*plv)*(*plv++));
	
		v11 += vfactor*(*pc11++);
		v22 += vfactor*(*pc22++);
		v33 += vfactor*(*pc33++);
		v23 += vfactor*(*pc23++);
		v13 += vfactor*(*pc13++);
		v12 += vfactor*(*pc12++);
	}
	Visc(2,1) = Visc(1,2);
	Visc(2,0) = Visc(0,2);
	Visc(1,0) = Visc(0,1);
}
				

/*returns viscosity tensor for given C and Cv in local frame, dV^-1_IK/dCv_J Sig_K*/
void  AnisoScleraWaxs_expo::ComputeDViscDCv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
		const dArrayT& Vec, dMatrixT& DVisc, const int pindex)
{
	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT& I4e = fI4;
	I4e /= fI4v;

/*	cout << "\nDCv";
	cout << "\nFiberStretch: "<< FiberStretch;
	cout << "\nFiberStretch_v: "<< FiberStretch_v;
	cout << "\nI4e: "<< I4e;
*/
	/* derivatives of the potential */
	fPotential[pindex+1]->MapDFunction(I4e, fdU);
	fPotential[pindex+1]->MapDDFunction(I4e, fddU);	

	/* initialize kernel pointers */	
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* ple   = I4e.Pointer();
	double* plv  = fI4v.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	/*stress*/
	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);
	
	DVisc = 0.0;
	double& dv11 = DVisc(0,0); /*cf_1111*/
	double& dv22 = DVisc(1,1); /*cf_2222*/
	double& dv33 = DVisc(2,2); /*cf_1212*/
	double& dv23 = DVisc(1,2); /*cf_2212*/
	double& dv13 = DVisc(0,2); /*cf_1112*/
	double& dv12 = DVisc(0,1); /*cf_1122*/

	for (int i = 0; i < I4e.Length(); i++)
	{
		double tau = 2.0*(*pdU)*(*ple);
		double dtau = 4.0*(*pddU++)*(*ple)+4.0*(*pdU++);
		double visc = fViscosity[pindex]->Function(tau);
		double dvisc = fViscosity[pindex]->DFunction(tau);

		double coeff1 = - (dvisc*dtau*(*ple) + 4.0*visc)/((*plv)*(*plv)*(*plv));
		double coeff2 =  coeff1*((*p1++)*Vec[0]+(*p2++)*Vec[1]+2.0*(*p3++)*Vec[2]);
		double coeff3 =  coeff2*(*pj++);
		
		ple++;
		plv++;

/*		cout << "\ntau: "<<tau;
		cout << "\ndtau: "<<dtau;
		cout << "\nvisc: "<<visc;
		cout << "\ndvisc: "<<dvisc;
		cout << "\ncoeff1: "<<coeff1;
		cout << "\ncoeff2: "<<coeff2;
		cout << "\ncoeff3: "<<coeff3;
*/
	
		dv11 += coeff3*(*pc11++);
		dv22 += coeff3*(*pc22++);
		dv33 += coeff3*(*pc33++);
		dv23 += coeff3*(*pc23++);
		dv13 += coeff3*(*pc13++);
		dv12 += coeff3*(*pc12++);
	}
	DVisc(2,1) = DVisc(1,2);
	DVisc(2,0) = DVisc(0,2);
	DVisc(1,0) = DVisc(0,1);
}
				
/*returns viscosity tensor for given C and Cv in local frame, dV^-1_IK/dC_J Sig_K/*/
void  AnisoScleraWaxs_expo::ComputeDViscDC(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
		const dArrayT& Vec, dMatrixT& DVisc, const int pindex)
{
	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT& I4e = fI4;
	I4e /= fI4v;
	
/*	cout << "\nDC";
	cout << "\nFiberStretch: "<< FiberStretch;
	cout << "\nFiberStretch_v: "<< FiberStretch_v;
	cout << "\nI4e: "<< I4e;
*/	
	/* derivatives of the potential */
	fPotential[pindex+1]->MapDFunction(I4e, fdU);
	fPotential[pindex+1]->MapDDFunction(I4e, fddU);	

	/* initialize kernel pointers */	
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* ple   = I4e.Pointer();
	double* plv  = fI4v.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	/*stress*/
	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);
	
	DVisc = 0.0;
	double& dv11 = DVisc(0,0); /*cf_1111*/
	double& dv22 = DVisc(1,1); /*cf_2222*/
	double& dv33 = DVisc(2,2); /*cf_1212*/
	double& dv23 = DVisc(1,2); /*cf_2212*/
	double& dv13 = DVisc(0,2); /*cf_1112*/
	double& dv12 = DVisc(0,1); /*cf_1122*/

	for (int i = 0; i < I4e.Length(); i++)
	{
		double tau = 2.0*(*pdU)*(*ple);
		double dtau = 4.0*(*pddU++)*(*ple)+4.0*(*pdU++);
		double visc = fViscosity[pindex]->Function(tau);
		double dvisc = fViscosity[pindex]->DFunction(tau);

		double coeff1 = (dvisc*dtau)/((*plv)*(*plv)*(*plv));
		double coeff2 =  coeff1*((*p1++)*Vec[0]+(*p2++)*Vec[1]+2.0*(*p3++)*Vec[2]);
		double coeff3 =  coeff2*(*pj++);
		
		ple++;
		plv++;

/*		cout << "\ntau: "<<tau;
		cout << "\ndtau: "<<dtau;
		cout << "\nvisc: "<<visc;
		cout << "\ndvisc: "<<dvisc;
		cout << "\ncoeff1: "<<coeff1;
		cout << "\ncoeff2: "<<coeff2;
		cout << "\ncoeff3: "<<coeff3;
*/		
	
		dv11 += coeff3*(*pc11++);
		dv22 += coeff3*(*pc22++);
		dv33 += coeff3*(*pc33++);
		dv23 += coeff3*(*pc23++);
		dv13 += coeff3*(*pc13++);
		dv12 += coeff3*(*pc12++);
	}
	DVisc(2,1) = DVisc(1,2);
	DVisc(2,0) = DVisc(0,2);
	DVisc(1,0) = DVisc(0,1);
}

/*calculates  I4=C:M*/
void AnisoScleraWaxs_expo::CompI4(const dSymMatrixT& FiberStretch)
{	
	/*calculate fibril lengths                                                        *
	 *I4 = C*:M where M = cos^2 a1 x a1 + sin^2 a2 x a2 + sin cos (a1 x a2 + a2 x a1) */

	const double& C11 = FiberStretch[0];
	const double& C22 = FiberStretch[1];
	const double& C12 = FiberStretch[2];

	/* initialize kernel pointers */
	double* pl = fI4.Pointer();
	double* s1 = fStressTable(0);
	double* s2 = fStressTable(1);
	double* s3 = fStressTable(2);
		
	for (int i = 0; i < fI4.Length(); i++)
		*pl++ = C11*(*s1++) + C22*(*s2++) + 2.0*C12*(*s3++);		
}

/*calculates both I4=C:M and I4v=Cv:M*/
void AnisoScleraWaxs_expo::CompI4(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v)
{	
	/*calculate fibril lengths                                                        *
	 *I4 = Ce*:\widtilde{M} = C:M/(Cv:M) where M = cos^2 a1 x a1 + sin^2 a2 x a2 + sin cos (a1 x a2 + a2 x a1) */
	/*rotate viscous stretch onto plane of fibers*/
	CompI4(FiberStretch_v);
	fI4v = fI4;

	/*calculate total I4*/
	CompI4(FiberStretch);
}

/* Initialize angle tables */
void AnisoScleraWaxs_expo::Construct(void)
{
	
	// open the dispersion file
        ifstreamT fDataInput(fUserFile);
  //	fDataInput.open(fUserFile);
	if (fDataInput.is_open()) 
	{
		cout << "\n WAXS file successfully open \n";
		int num_elem=NumElements();
		fFiberDispersion_list.Dimension(num_elem);
	
		//cout << "num_elem"<<num_elem<< endl;
		// number of division of 2pi
		int deg;
		fDataInput >> deg;
		//cout << "angle"<<deg<<endl;
		for (int k=0; k<num_elem;  k++)
		{
			dArrayT& Disp_scl = fFiberDispersion_list[k];
			Disp_scl.Dimension(deg);
			for (int j=0; j<deg; j++)
			{
				double disp;
				fDataInput >> disp;
				Disp_scl[j]=disp;
				//cout << "j=="<<j<< " and waxs==%"<<Disp_scl[j]<<"\n"<<endl;
			}	
		}
	}
	 else
	{
		cout << "Error opening file"<<endl;
	}
	fDataInput.close();

	
	//dArrayT& Disp_scl_temp = fFiberDispersion_list[9];
	
	
	//for (int j=0; j<64; j++) {
//		cout << Disp_scl_temp[j]<< "  ";
//	}
//	cout<< "------------"<<endl;
	/*
	for (int k=0; k<NumElements();  k++) //change to while file is still good;
		{
		cout<< "\n element == "<< k <<"  dipersion =="<<fFiberDispersion_list[k]<<endl;
		}
	
	*/
	cout<< "dispersion file close"<<endl;
	const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/

	/* fetch points */
	int numbonds = points.MajorDim();
		
	int nel = NumElements();     
	fjacobians.Dimension(nel);

	/* length table */
	fI4.Dimension(numbonds);
	fI4v.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumFibStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumFibModuli, numbonds);	

	/* set pointers */
	double *s1 = fStressTable(0);
	double *s2 = fStressTable(1);
	double *s3 = fStressTable(2);

	double *c11 = fModuliTable(0);
	double *c22 = fModuliTable(1);
	double *c33 = fModuliTable(2);
	double *c23 = fModuliTable(3);
	double *c13 = fModuliTable(4);
	double *c12 = fModuliTable(5);

//	StringT name("distributions.dat");
//	ofstreamT dist_out(name);
			 		
	fFSFiberMatSupport->TopElement();
	/*	
		for (int i=0; i<20; i++) {
		cout<<"dispersion\n"<<FiberMatSupportT().Fiber_Disp(i);
		}
		*/
	
	/*integrate dtemp - min to normalize distribution function*/
	const dArrayT& angles = fCircle->CircleAngles(0.0);
	const dArrayT& weights = fCircle->Jacobians();
		while (fFSFiberMatSupport->NextElement()) 
		{
			int ielm = fFSFiberMatSupport->CurrElementNumber();
			dArrayT& disp_elem=fFiberDispersion_list[ielm];
			
				//cout<<disp_elem<<endl;
//				
//			cout<<"-----------------------------"<<endl;
//			
			dArrayT& jac = fjacobians[ielm];
			jac.Dimension(numbonds);
			//cout<<numbonds<<endl;
			const dArrayT& angles = fCircle->CircleAngles(0.0);
			//cout<<"angle\n:"<<angles<<endl;
			
			const dArrayT& weights = fCircle->Jacobians();
//			cout<<"weight\n:"<<weights<<endl;
			int sum=0;
			for (int i = 0; i < numbonds; i++)
			{
				double f = disp_elem[i];
				jac[i] = f*weights[i]; 
			}
			
			double energy_anis=0;
			double* pl =  jac.Pointer();
			for (int i = 0; i < numbonds; i++) {
				energy_anis += (*pl++);
			}
			//make sure the xray scatter is normalized
			for (int i = 0; i < numbonds; i++)
			{
				jac[i] = jac[i]/energy_anis; 
			}
			
			//cout<< jac <<endl;
			//cout << angles <<endl;
//			cout<<"-------------"<<endl;
//			cout <<disp_elem<<endl;
			/*if (ielm ==10) {
				cout<< "Distribution"<<endl;
			cout<< jac <<endl;
				cout << "Angles" << endl;
				cout << angles << endl;
			}
			 */
			
			//if (ielm==2) {
			//			cout << jac <<endl;
//				cout << "---------------------"<<endl;
//				cout << angles <<endl;
//				cout << "---------------------"<<endl;
//				cout << numbonds <<endl;
//				cout << "---------------------"<<endl;
//				cout << weights <<endl;
//
//			}
			
		
		}
		
		
		for (int i = 0; i < numbonds; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s1[i] = cosi*cosi;      
		s2[i] = sini*sini;
		s3[i] = sini*cosi;
	
		/* moduli angle tables */
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];
		c33[i] = s3[i]*s3[i];
		c23[i] = s2[i]*s3[i];
		c13[i] = s1[i]*s3[i];
		c12[i] = s2[i]*s1[i];
	}
	
	
	//check that the distribution function integrates to one
	/*
	fFSFiberMatSupport->TopElement();
	cout<< "ENERGY "<<endl;
	while (fFSFiberMatSupport->NextElement()) 
	{
		double energy_anis_1=0;
		
		double* pl =  fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
		for (int i = 0; i < fI4.Length(); i++) {
			energy_anis_1 += (*pl++);
		}
		cout << energy_anis_1 <<endl;
	}
	cout << "----------------------" << endl;
	*/	
	 
}


#endif




