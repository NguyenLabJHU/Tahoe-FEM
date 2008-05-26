/* $Id: AnisoCorneaVisco.cpp,v 1.11 2008-05-26 21:30:12 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "AnisoCorneaVisco.h"
#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "ofstreamT.h"

#ifdef VIB_MATERIAL

/* point generator */
#include "EvenSpacePtsT.h"

#include "FungType.h"
#include "PowerTrig.h"
#include "ScaledCsch.h"
#include "LinearExponentialT.h"

#define MY_DEBUG

const double Pi = acos(-1.0);
const int kNumOutputVar = 7;
static const char* Labels[kNumOutputVar] = 
	{"NT_X", "NT_Y", "NT_Z","IS_X", "IS_Y", "IS_Z","phiNEQ"};

static const int perm[3][3] = {0,1,2,1,2,0,2,0,1};

using namespace Tahoe;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
AnisoCorneaVisco::AnisoCorneaVisco(void):
  FSFiberMatViscT(),
  ParameterInterfaceT("aniso_viscoelastic_cornea"),
	fCircle(NULL),
	fDistribution(NULL),
	finhomogeneous(kHomogeneous),
	r1(0.0),r2(0.0),r3(0.0),r4(0.0)
{
	/*reset default*/
	fNumFibProcess = 0;
	fNumMatProcess = 0;

#ifndef VIB_MATERIAL
	ExceptionT::BadInputValue("AnisoCorneaVisco::AnisoCorneaVisco", 
		"VIB_MATERIAL must be enabled");
#endif
}

/* destructor */
AnisoCorneaVisco::~AnisoCorneaVisco(void) 
{ 
	delete fCircle; 
	delete fDistribution;
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
double AnisoCorneaVisco::StrainEnergyDensity(void)
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
	double* pj = fjacobian.Pointer();	
	if (finhomogeneous) {
 		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}
	for (int i = 0; i < fI4.Length(); i++)
		energy += (*pU++)*(*pj++);
		
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
		double* pj = fjacobian.Pointer();	
		if (finhomogeneous) {
			pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
		}
		for (int k = 0; k < I4e.Length(); k++)
			energy += (*pU++)*(*pj++);
	}
	return energy;
}

/* nonequilibrium free energy density */
double AnisoCorneaVisco::NonequilibriumStrainEnergyDensity(void)
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

  double* pj = fjacobian.Pointer();	
  if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
  }
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

int AnisoCorneaVisco::NumOutputVariables() const {
	return kNumOutputVar;
}

void AnisoCorneaVisco::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void AnisoCorneaVisco::ComputeOutput(dArrayT& output)
{
	/*calculates deformed fiber vectors*/
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	
	const double* p_nt = Fibers(0);
	const double* p_is = Fibers(1);
	
	const dMatrixT& F = F_mechanical();
	double* pb = output.Pointer();
	
	/*deformed NT fiber orientation*/
	F.Multx(p_nt, pb);
	pb += NumSD();
	
	/*deformed IS fiber orientation*/
	F.Multx(p_is, pb);
	pb += NumSD();

	/*non-equilibrium strain energy density */
	*pb = NonequilibriumStrainEnergyDensity();
}

/* describe the parameters needed by the interface */
void AnisoCorneaVisco::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);
	
	/* 2D option must be plain stress !!!!NOTE!!!! */ 
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);

	/* integration points */
	ParameterT points(ParameterT::Integer, "n_points");
	points.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(points);
}

/* information about subordinate parameter lists */
void AnisoCorneaVisco::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("matrix_material_params", ParameterListT::Once);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("eq_fibril_potential", ParameterListT::Once);
	sub_list.AddSub("neq_fibril_potential", ParameterListT::Any);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("viscosity", ParameterListT::Any);

	/* choice of fibril distribution funcion */
	sub_list.AddSub("fibril_distribution", ParameterListT::Once);
	
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AnisoCorneaVisco::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = FSSolidMatT::NewSub(name);
//	ParameterInterfaceT* sub = FSFiberMatViscT::NewSub(name);
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
	else if (name == "eq_fibril_potential" || name == "neq_fibril_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT fung("fung_type");		
		LimitT lower(0.0, LimitT::Lower);

		ParameterT alpha(ParameterT::Double, "alpha");
		ParameterT beta(ParameterT::Double, "beta");

		fung.AddParameter(alpha);
		fung.AddParameter(beta);
		alpha.AddLimit(lower);
		beta.AddLimit(lower);

		/* set the description */
		fung.SetDescription("f(I) = alpha*(exp(beta*(I - 1.0)) + beta/I)");	
		choice->AddSub(fung);
		return(choice);
	}
	else if (name == "viscosity")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT sinh("scaled_csch");		
		LimitT lower(0.0, LimitT::Lower);

		ParameterT eta0(ParameterT::Double, "eta0");
		ParameterT tau0(ParameterT::Double, "tau0");

		sinh.AddParameter(eta0);
		sinh.AddParameter(tau0);
		eta0.AddLimit(lower);
		tau0.AddLimit(lower);

		/* set the description */
		sinh.SetDescription("f(tau) = eta0* tau/tau0 / Sinh(tau/tau0)");	
		choice->AddSub(sinh);
		
		ParameterContainerT exp("linear-exponential");
 
		ParameterT a(ParameterT::Double, "a");
		ParameterT b(ParameterT::Double, "b");
		ParameterT c(ParameterT::Double, "c");
		ParameterT d(ParameterT::Double, "d");

		exp.AddParameter(a);
		exp.AddParameter(b);
		exp.AddParameter(c);
		exp.AddParameter(d);
		/* set the description */
		exp.SetDescription("f(x) = a + b x + c (1 - exp[-x/d])");	
		choice->AddSub(exp);

		return(choice);
	}
	else if (name == "fibril_distribution")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		LimitT lower(0.0, LimitT::Lower);
		LimitT upper(1.0, LimitT::Upper);


		/* sine & cosine raised to a power */
		ParameterContainerT powertrig("power_trig");		

		ParameterT a(ParameterT::Double, "a");
		ParameterT b(ParameterT::Double, "b");
		ParameterT c(ParameterT::Double, "c");
		ParameterT n(ParameterT::Double, "n");
		ParameterT phi(ParameterT::Double, "phi");

		powertrig.AddParameter(a);
		powertrig.AddParameter(b);
		powertrig.AddParameter(c);
		powertrig.AddParameter(n);
		powertrig.AddParameter(phi);
		c.AddLimit(lower);
		c.AddLimit(upper);
		n.AddLimit(lower);

		/* set the description */
		powertrig.SetDescription("f(theta) = a*cos(theta+phi)^n + b*sin(theta+phi)^n + c");	
		choice->AddSub(powertrig);

		/* spatially interpolated blending of two distributions */
		ParameterContainerT blend_powertrig("blend_power_trig");		
		{
		ParameterT a1(ParameterT::Double, "a1");
		ParameterT b1(ParameterT::Double, "b1");
		ParameterT c1(ParameterT::Double, "c1");
		ParameterT n1(ParameterT::Double, "n1");

		blend_powertrig.AddParameter(a1);
		blend_powertrig.AddParameter(b1);
		blend_powertrig.AddParameter(c1);
		blend_powertrig.AddParameter(n1);
		c1.AddLimit(lower);
		c1.AddLimit(upper);
		n1.AddLimit(lower);

		ParameterT a2(ParameterT::Double, "a2");
		ParameterT b2(ParameterT::Double, "b2");
		ParameterT c2(ParameterT::Double, "c2");
		ParameterT n2(ParameterT::Double, "n2");

		blend_powertrig.AddParameter(a2);
		blend_powertrig.AddParameter(b2);
		blend_powertrig.AddParameter(c2);
		blend_powertrig.AddParameter(n2);
		c2.AddLimit(lower);
		c2.AddLimit(upper);
		n2.AddLimit(lower);

		ParameterT r1(ParameterT::Double, "r1");
		ParameterT r2(ParameterT::Double, "r2");
		blend_powertrig.AddParameter(r1);
		blend_powertrig.AddParameter(r2);
		r1.AddLimit(lower);
		r2.AddLimit(lower);
		}
		/* set the description */
		blend_powertrig.SetDescription("f(theta) = dist1(theta,phi=0)*(r-r2)/(r1-r2) + dist1(theta,phi(x))*(r-r1)/(r2-r1) ");	
		choice->AddSub(blend_powertrig);

		/* a full map of a cornea's distributions */
		ParameterContainerT cornea("cornea");		
		{
		ParameterT a1(ParameterT::Double, "a1");
		ParameterT b1(ParameterT::Double, "b1");
		ParameterT c1(ParameterT::Double, "c1");
		ParameterT n1(ParameterT::Double, "n1");

		cornea.AddParameter(a1);
		cornea.AddParameter(b1);
		cornea.AddParameter(c1);
		cornea.AddParameter(n1);
		c1.AddLimit(lower);
		c1.AddLimit(upper);
		n1.AddLimit(lower);

		ParameterT a2(ParameterT::Double, "a2");
		ParameterT b2(ParameterT::Double, "b2");
		ParameterT c2(ParameterT::Double, "c2");
		ParameterT n2(ParameterT::Double, "n2");

		cornea.AddParameter(a2);
		cornea.AddParameter(b2);
		cornea.AddParameter(c2);
		cornea.AddParameter(n2);
		c2.AddLimit(lower);
		c2.AddLimit(upper);
		n2.AddLimit(lower);

		ParameterT c3(ParameterT::Double, "c3");
		cornea.AddParameter(c3);
		c3.AddLimit(lower);

		ParameterT r1(ParameterT::Double, "r1"); // NT-IS end
		cornea.AddParameter(r1);
		r1.AddLimit(lower);
		ParameterT r2(ParameterT::Double, "r2"); // circumferential begin
		cornea.AddParameter(r2);
		r2.AddLimit(lower);
		ParameterT r3(ParameterT::Double, "r3"); // circumferential end
		cornea.AddParameter(r3);
		r3.AddLimit(lower);
		ParameterT r4(ParameterT::Double, "r4"); // uniform begin
		cornea.AddParameter(r4);
		r4.AddLimit(lower);
		}
		/* set the description */
		cornea.SetDescription("f(theta) = distNT-IS(0:r1) + distCIRC(r2:r3) + distISO(r4:)  ");	
		choice->AddSub(cornea);

		/* element distributions read from a file */
		ParameterContainerT inhomo_dist("inhomogeneous_distribution");		
		ParameterT dist_file(ParameterT::String, "element_distributions_file");
		inhomo_dist.AddParameter(dist_file);
		inhomo_dist.SetDescription("read element distributions from file");	
		choice->AddSub(inhomo_dist);

	  /* a full map of a cornea's distributions */
		ParameterContainerT cornea_mod("cornea_mod");
		{
			ParameterT a1(ParameterT::Double, "a1");
			ParameterT b1(ParameterT::Double, "b1");
			ParameterT c1(ParameterT::Double, "c1");
			ParameterT n1(ParameterT::Double, "n1");

			cornea_mod.AddParameter(a1);
			cornea_mod.AddParameter(b1);
			cornea_mod.AddParameter(c1);
			cornea_mod.AddParameter(n1);
			c1.AddLimit(lower);
			c1.AddLimit(upper);
			n1.AddLimit(lower);

			ParameterT a2(ParameterT::Double, "a2");
			ParameterT b2(ParameterT::Double, "b2");
			ParameterT c2(ParameterT::Double, "c2");
			ParameterT n2(ParameterT::Double, "n2");

			cornea_mod.AddParameter(a2);
			cornea_mod.AddParameter(b2);
			cornea_mod.AddParameter(c2);
			cornea_mod.AddParameter(n2);
			c2.AddLimit(lower);
			c2.AddLimit(upper);
			n2.AddLimit(lower);

			ParameterT r1(ParameterT::Double, "r1"); // NT-IS end
			cornea_mod.AddParameter(r1);
			r1.AddLimit(lower);
			ParameterT r2(ParameterT::Double, "r2"); // circumferential begin
			cornea_mod.AddParameter(r2);
			r2.AddLimit(lower);
			ParameterT r3(ParameterT::Double, "R_IS"); // circumferential end
			cornea_mod.AddParameter(r3);
			r3.AddLimit(lower);
			ParameterT r4(ParameterT::Double, "R_NT"); // circumferential end
			cornea_mod.AddParameter(r4);
			r4.AddLimit(lower);
		}
	
    /* set the description */
    cornea_mod.SetDescription("f(theta) = distNT-IS(0:r1) + distCIRC(r2:r3) ");
    choice->AddSub(cornea_mod);

	return(choice);
	}
}

/* accept parameter list */
void AnisoCorneaVisco::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);

	/*initializing some parameters for cornea_mod fibril distribution*/
	int num_neq_pot = list.NumLists("neq_fibril_potential");
	int num_visc = list.NumLists("viscosity");

	if (num_visc != num_neq_pot)
		ExceptionT::GeneralFail("AnisoCorneaVisco::TakeParameterList", 
			"number of viscosity functions does not match number of nonequilibrium potentials");

	fNumFibProcess = num_neq_pot;
	fPotential.Dimension(fNumFibProcess+1);
	if(fNumFibProcess > 0)
		fViscosity.Dimension(fNumFibProcess);
		
	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_material_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fGamma = matrix.GetParameter("bulk_modulus");
	}

	const ParameterListT& potential = list.GetListChoice(*this, "eq_fibril_potential");
	if (potential.Name() == "fung_type")
	{
		double alpha = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential[0] = new FungType(alpha, beta);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
	}

	for (int i = 0; i < fNumFibProcess; i++)
	{
		const ParameterListT& neq_potential = list.GetListChoice(*this, "neq_fibril_potential", i);
		if (neq_potential.Name() == "fung_type")
		{
			double alpha = neq_potential.GetParameter("alpha");
			double beta = neq_potential.GetParameter("beta");
			fPotential[i+1] = new FungType(alpha,beta);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
		}
		
		const ParameterListT& visc = list.GetListChoice(*this, "viscosity", i);
		if (visc.Name() == "linear-exponential")
		{
			double a = visc.GetParameter("a");
			double b = visc.GetParameter("b");
			double c = visc.GetParameter("c");
			double d = visc.GetParameter("d");			
			fViscosity[i] = new LinearExponentialT(a,b,c,d);
			if (!fViscosity[i]) throw ExceptionT::kOutOfMemory;
		}
		else if (visc.Name() == "scaled_csch")
		{
			double a = visc.GetParameter("eta0");
			double b = visc.GetParameter("tau0");
			fViscosity[i] = new ScaledCsch(a,b);
			if (!fViscosity[i]) throw ExceptionT::kOutOfMemory;
		}
	}
	
	const ParameterListT& distr = list.GetListChoice(*this, "fibril_distribution");
	if (distr.Name() == "power_trig")
	{
		double a = distr.GetParameter("a");
		double b = distr.GetParameter("b");
		double c = distr.GetParameter("c");
		double n = distr.GetParameter("n");
		double phi = distr.GetParameter("phi");
		fDistribution = new PowerTrig(a,b,c,n,phi); 
		if (!fDistribution) throw ExceptionT::kOutOfMemory;
	}
	else if (distr.Name() == "blend_power_trig")
	{
		double a1 = distr.GetParameter("a1");
		double b1 = distr.GetParameter("b1");
		double c1 = distr.GetParameter("c1");
		double n1 = distr.GetParameter("n1");
		fDistribution = new PowerTrig(a1,b1,c1,n1,0.0); 
		if (!fDistribution) throw ExceptionT::kOutOfMemory;

		a2 = distr.GetParameter("a2");
		b2 = distr.GetParameter("b2");
		c2 = distr.GetParameter("c2");
		n2 = distr.GetParameter("n2");

		r1 = distr.GetParameter("r1");
		r2 = distr.GetParameter("r2");
		finhomogeneous = kBlend;
	}
	else if (distr.Name() == "cornea")
	{
		double a1 = distr.GetParameter("a1");
		double b1 = distr.GetParameter("b1");
		double c1 = distr.GetParameter("c1");
		double n1 = distr.GetParameter("n1");
		fDistribution = new PowerTrig(a1,b1,c1,n1,0.0); 
		if (!fDistribution) throw ExceptionT::kOutOfMemory;

		a2 = distr.GetParameter("a2");
		b2 = distr.GetParameter("b2");
		c2 = distr.GetParameter("c2");
		n2 = distr.GetParameter("n2");

		c3 = distr.GetParameter("c3");

		r1 = distr.GetParameter("r1");
		r2 = distr.GetParameter("r2");
		r3 = distr.GetParameter("r3");
		r4 = distr.GetParameter("r4");
		finhomogeneous = kBlend;
	}
  else if (distr.Name() == "cornea_mod")
  {
    double a1 = distr.GetParameter("a1");
    double b1 = distr.GetParameter("b1");
    double c1 = distr.GetParameter("c1");
    double n1 = distr.GetParameter("n1");
    fDistribution = new PowerTrig(a1,b1,c1,n1,0.0);
    if (!fDistribution) throw ExceptionT::kOutOfMemory;

    a2 = distr.GetParameter("a2");
    b2 = distr.GetParameter("b2");
    c2 = distr.GetParameter("c2");
    n2 = distr.GetParameter("n2");

    r1 = distr.GetParameter("r1");
    r2 = distr.GetParameter("r2");
    r3 = distr.GetParameter("R_IS");
    r4 = distr.GetParameter("R_NT");
	
	if (r2 < r1) 
			ExceptionT::GeneralFail("AnisoCorneaVisco::TakeParameterList", 
			"radius of limbal annulus has to be greater than radius of central region");
	if (r3 < r2) 
			ExceptionT::GeneralFail("AnisoCorneaVisco::TakeParameterList", 
			"IS meridan has to be greater than radius of limbal annulus");
	if (r4 < r3) 
			ExceptionT::GeneralFail("AnisoCorneaVisco::TakeParameterList", 
			"NT meridan has to be greater than IS meridian");

    finhomogeneous = kBlend_Mod;
  }
	else if (distr.Name() == "inhomogeneous_distribution") {
    ParameterContainerT inhomo_dist("inhomogeneous_distribution");
    StringT dist_file =  distr.GetParameter("element_distributions_file");
		finhomogeneous = kCornea;
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
	Construct();
}
	
/***********************************************************************
 * Protected
 ***********************************************************************/
void AnisoCorneaVisco::ComputeMatrixStress(const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, 
				dSymMatrixT& Stress, const int process_index, const int fillmode)
{
	if (process_index > -1)
		ExceptionT::GeneralFail("AnisoCorneaVisco::ComputeMatrixStress", 
			"Material model does not allow for neq matrix description");
	if (fillmode != dSymMatrixT::kOverwrite)
		ExceptionT::GeneralFail("AnisoCorneaVisco::ComputeMatrixStress", 
			"Expects to overwrite Stress");
			
	
	/*2pdf{W}{C_IJ} = mu ( del_IJ - I3^-gamma C^-1_IJ)*/
	double I3 = Stretch.Det();
	double I3rg = pow(I3, -fGamma);
	Stress.Inverse(Stretch);
	Stress *= -I3rg*fMu;
	
	Stress[0] += fMu;
	Stress[1] += fMu;
	Stress[2] += fMu;
}

void AnisoCorneaVisco::ComputeMatrixMod(const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress,
				dMatrixT& Mod, const int process_index, const int fillmode)
{
	if (process_index > -1)
		ExceptionT::GeneralFail("AnisoCorneaVisco::ComputeMatrixMod", 
			"Material model does not allow for neq matrix description");
	if (fillmode != dSymMatrixT::kOverwrite)
		ExceptionT::GeneralFail("AnisoCorneaVisco::ComputeMatrixMod", 
			"Expects to overwrite Stress");

	/*matrix contribution*/
	double I3 = fC.Det();
	double I3rg = pow(I3,-fGamma);
	Stress.Inverse(fC);
	
	/*2pdf{S_IJ}{C_KL} = 2 mu I_3^-gamma (gamma C^-1_IJ C^-1_KL + 0.5(C^-1_IK C^-1_JL +C^-1_IL+C^-1_JK)*/
	/*modulus*/
	double coeff = 2.0*fMu*I3rg;
	Mod.ReducedI_C(Stress);
	Mod *= coeff;
	
	coeff *= fGamma;
	Mod(0,0) += coeff*Stress[0]*Stress[0];
	Mod(0,1) += coeff*Stress[0]*Stress[1];
	Mod(0,2) += coeff*Stress[0]*Stress[2];
	Mod(0,3) += coeff*Stress[0]*Stress[3];
	Mod(0,4) += coeff*Stress[0]*Stress[4];
	Mod(0,5) += coeff*Stress[0]*Stress[5];
	
	Mod(1,0) += coeff*Stress[1]*Stress[0];
	Mod(1,1) += coeff*Stress[1]*Stress[1];
	Mod(1,2) += coeff*Stress[1]*Stress[2];
	Mod(1,3) += coeff*Stress[1]*Stress[3];
	Mod(1,4) += coeff*Stress[1]*Stress[4];
	Mod(1,5) += coeff*Stress[1]*Stress[5];

	Mod(2,0) += coeff*Stress[2]*Stress[0];
	Mod(2,1) += coeff*Stress[2]*Stress[1];
	Mod(2,2) += coeff*Stress[2]*Stress[2];
	Mod(2,3) += coeff*Stress[2]*Stress[3];
	Mod(2,4) += coeff*Stress[2]*Stress[4];
	Mod(2,5) += coeff*Stress[2]*Stress[5];

	Mod(3,0) += coeff*Stress[3]*Stress[0];
	Mod(3,1) += coeff*Stress[3]*Stress[1];
	Mod(3,2) += coeff*Stress[3]*Stress[2];
	Mod(3,3) += coeff*Stress[3]*Stress[3];
	Mod(3,4) += coeff*Stress[3]*Stress[4];
	Mod(3,5) += coeff*Stress[3]*Stress[5];

	Mod(4,0) += coeff*Stress[4]*Stress[0];
	Mod(4,1) += coeff*Stress[4]*Stress[1];
	Mod(4,2) += coeff*Stress[4]*Stress[2];
	Mod(4,3) += coeff*Stress[4]*Stress[3];
	Mod(4,4) += coeff*Stress[4]*Stress[4];
	Mod(4,5) += coeff*Stress[4]*Stress[5];

	Mod(5,0) += coeff*Stress[5]*Stress[0];
	Mod(5,1) += coeff*Stress[5]*Stress[1];
	Mod(5,2) += coeff*Stress[5]*Stress[2];
	Mod(5,3) += coeff*Stress[5]*Stress[3];
	Mod(5,4) += coeff*Stress[5]*Stress[4];
	Mod(5,5) += coeff*Stress[5]*Stress[5];
}

/*computes integrated fiber stress in local frame*/
void AnisoCorneaVisco::ComputeFiberStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
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
		fPotential[pindex+1]->MapDFunction(fI4, fdU);

		/* initialize kernel pointers */
		double* pdU = fdU.Pointer();
		double* pj  = fjacobian.Pointer();
		if (finhomogeneous) {
			pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
		}
	
		/*integrate w.r.t in-plane orientation theta*/
		for (int i = 0; i < fI4.Length(); i++)
		{
			double factor = (*pj++)*(*pdU++)/Pi;
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
		double* pj  = fjacobian.Pointer();
		if (finhomogeneous) {
			pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
		}
	
		/*integrate w.r.t in-plane orientation theta*/
		for (int i = 0; i < I4e.Length(); i++)
		{
			double factor = (*pj++)*(*pdU++)/(*plv++)/Pi;
			s1 += factor*(*p1++);
			s2 += factor*(*p2++);
			s3 += factor*(*p3++);
		}
	}
}
	
/*computes integrated moduli in local frame*/
void AnisoCorneaVisco::ComputeFiberMod (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
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
		double* pj   = fjacobian.Pointer();
		if (finhomogeneous) {
			pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
		}

		for (int i = 0; i < fI4.Length(); i++)
		{
			double sfactor =  (*pj)*(*pdU++)/Pi;
			double cfactor = 2.0/Pi*(*pj++)*(*pddU++);

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
		double* pj   = fjacobian.Pointer();
		if (finhomogeneous) {
			pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
		}

		for (int i = 0; i < I4e.Length(); i++)
		{
			double sfactor =  (*pj)*(*pdU++)/(*plv)/Pi;
			double cfactor = 2.0/Pi*(*pj++)*(*pddU++)/((*plv)*(*plv));
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

void AnisoCorneaVisco::ComputeCalg(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,  dMatrixT& Calg, const int pindex)
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
	
//	cout << "\ncomputeCalg: ";
//	cout << "\nfMod2: "<<fMod2;
	
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
void AnisoCorneaVisco::Compute_Cv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_vn, 
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
		ExceptionT::GeneralFail("AnisoCorneaVisco::Compute_Cv", 
			"number of iteration exceeds maximum of 10");
}


/*computes integrated fiber stress in local frame*/
void AnisoCorneaVisco::ComputeFlowStress (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
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
	double* pj  = fjacobian.Pointer();
	if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}

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
		double factor = (*pj++)*(*pdU++)*(*ple++)/(*plv++)/Pi;
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}
}

/*computes dFlowStress/dC.  Note for this model, dFlowStress/dC = - dSNEQ/dCv*/
void AnisoCorneaVisco::dFlowdC (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v, 
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
	double* pj   = fjacobian.Pointer();
	if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}

	for (int i = 0; i < I4e.Length(); i++)
	{
		double cfactor = 2.0/Pi*(*pj++)*( (*pddU++)*(*ple)/((*plv)*(*plv)) + (*pdU++)/((*plv)*(*plv)));
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
void AnisoCorneaVisco::dFlowdCv (const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
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
	double* pj   = fjacobian.Pointer();
	if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}

	for (int i = 0; i < I4e.Length(); i++)
	{
		double cfactor = -2.0/Pi*(*pj++)*( (*pddU++)*(*ple)*(*ple)/((*plv)*(*plv)) 
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

void AnisoCorneaVisco::ComputeViscosity(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,  
			dMatrixT& Visc,  const int pindex)
{
	/*calculate I4e and I4v*/
	CompI4(FiberStretch, FiberStretch_v);
	dArrayT& I4e = fI4;
	I4e /= fI4v;
	
	double* ple   = I4e.Pointer();
	double* plv  = fI4v.Pointer();
	double* pj   = fjacobian.Pointer();
	if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}

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
		double vfactor =  (*pj++)*visc/((*plv)*(*plv++)) * 0.5/Pi;
	
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
void  AnisoCorneaVisco::ComputeDViscDCv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
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
	double* pj   = fjacobian.Pointer();
	if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}

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
		double coeff3 =  coeff2*(*pj++)*0.5/Pi;
		
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
void  AnisoCorneaVisco::ComputeDViscDC(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v,
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
	double* pj   = fjacobian.Pointer();
	if (finhomogeneous) {
		pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	}

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
		double coeff3 =  coeff2*(*pj++)*0.5/Pi;
		
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
void AnisoCorneaVisco::CompI4(const dSymMatrixT& FiberStretch)
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
void AnisoCorneaVisco::CompI4(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_v)
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
void AnisoCorneaVisco::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numbonds = points.MajorDim();
		
	/* length table */
	fI4.Dimension(numbonds);
	fI4v.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* jacobian table */
	fjacobian.Dimension(numbonds);

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

	/* fetch jacobians with fibril distribution function */
	fjacobian = fCircle->Jacobians(0.0, fDistribution);
#ifdef MY_DEBUG
  /* output stream */
  ofstreamT& out   = fFSFiberMatSupport->Output();
  bool print_input = fFSFiberMatSupport->PrintInput();
{
	out << "#DISTRIBUTION: " ;
  for (int i = 0; i < numbonds; i++) {
		out << fjacobian[i] << " ";
	}
	out << "\n";
}
#endif

	if (finhomogeneous == kBlend || finhomogeneous == kBlend_Mod) 
	{
		const dArray2DT& coordinates = fFSFiberMatSupport->InitialCoordinates();
		int nelm = fFSFiberMatSupport->NumElements(); 
		fjacobians.Dimension(nelm);
		dArrayT jac(numbonds);
		StringT name("distributions.dat");
		ofstreamT dist_out(name);
		dArrayT xc(3);
		fFSFiberMatSupport->TopElement();

		double inv_d_wg = numbonds/(2.0*Pi);
		while (fFSFiberMatSupport->NextElement()) 
		{
			int ielm = fFSFiberMatSupport->CurrElementNumber();
			fjacobians[ielm].Dimension(numbonds);
			fjacobians[ielm] = 0.0;
			
			/* calculate element centroid */
			iArrayT nodes = fFSFiberMatSupport->CurrentElement()->NodesX();
			int nen = NumElementNodes();       
			xc = 0.0;
			for (int i = 0; i < nen; i++) 
			{
				for (int j = 0; j < coordinates.MinorDim(); j++)
				xc[j] += coordinates(nodes[i],j);
			}
			xc /= nen;

			int inormal = 2; // HACK
			double x1 = xc[perm[inormal][1]];
			double x2 = xc[perm[inormal][2]];

			// projected polar coordinates
			double r = sqrt(x1*x1 + x2*x2);
			double phi = atan2(x2,x1);
			double wg = 1.0;

			if (finhomogeneous == kBlend)
			{
				C1FunctionT* distribution, * iso_distribution;
				iso_distribution = new PowerTrig(0.0,0.0,c3,0,0.0);

				// blending function : assuming linear radial weighting function
				if (r3 > 0.0 && r > r3) {  /*the point is in the sclera*/
					// distribution three : uniform
					wg = (r-r3)/(r4-r3); 
					wg = min(max(wg,0.0),1.0);
					jac = fCircle->Jacobians(0.0, iso_distribution);
				} 
				else {
					// distribution one : aligned w/ "NT-IS" system
					wg = (r-r2)/(r1-r2); // HACK : use ellipsoidal coord system
					wg = min(max(wg,0.0),1.0);
					jac = fCircle->Jacobians(0.0, fDistribution);  /*Cos[theta]^8 + Sin[theta]^8+c1*/
				}
				fjacobians[ielm].AddScaled(wg,jac);

				// distribution two : rotated wrt/ "NT-IS" system
				wg = 1.0 - wg;
				distribution = new PowerTrig(a2,b2,c2,n2,-phi); // assuming power trig
				jac = fCircle->Jacobians(0.0, distribution);
				delete distribution;
				fjacobians[ielm].AddScaled(wg,jac);
			}
			if (finhomogeneous == kBlend_Mod)
			{
				// projected polar coordinates
				double ecc = sqrt(1.0 - r3*r3/(r4*r4)); //eccentricity: sqrt( 1 - RIS^2/RNT^2 )
				r = sqrt( (1.0 - ecc*ecc) * x1*x1 + x2*x2); /*map r to ellipse*/
				phi = atan2(x2,x1);

				C1FunctionT* Dlimbus, *Dcentral;
				Dlimbus = new PowerTrig(a2,b2,c2,n2,-phi);
				Dcentral = fDistribution;
				// blending function : assuming linear radial weighting function
				if (r < r1){  /*central region*/
					fjacobians[ielm] = fCircle->Jacobians(0.0, Dcentral);  
				}
				else if (r > r2 && (r2-r1) > kSmall) { /*limbus*/
					// distribution one : aligned w/ "NT-IS" system
					fjacobians[ielm] = fCircle->Jacobians(0.0, Dlimbus); 
				}				
				else { /*peripheral region: linear blending*/
					wg = (r2-r)/(r2-r1); 
					jac = fCircle->Jacobians(0.0, Dcentral);  
					fjacobians[ielm].AddScaled(wg,jac);

					wg = (r-r1)/(r2-r1); 
					jac = fCircle->Jacobians(0.0, Dlimbus);
					fjacobians[ielm].AddScaled(wg,jac);
				}
			}
#ifdef MY_DEBUG
{
			out << "element: " << ielm << " x: " << x1 << " " << x2 << " r: " << r << " phi: " << phi << " wg: " << wg << "\n";
			out << "#distribution: " << ielm << " ";
      for (int i = 0; i < numbonds; i++) {
				out << fjacobians[ielm][i] << " ";
			}
			out << "\n";
}
#endif
			dist_out << "element: " << ielm+1 << " x: " << x1 << " " << x2 
           << " r: " << r << " phi: " << phi << " wg: " << wg << "\n";
			dist_out << "#distribution: " << ielm << " ";
      for (int i = 0; i < numbonds; i++) {
				dist_out << inv_d_wg*fjacobians[ielm][i] << " ";
			}
			dist_out << inv_d_wg*fjacobians[ielm][0] << " ";
			dist_out << "\n";
		}
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
}
#endif /*VIB_MATERIAL*/
