/* $Id: AnisoCornea.cpp,v 1.11 2008-06-01 01:05:34 thao Exp $ */
/* created: paklein (11/08/1997) */

#include "AnisoCornea.h"
#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"
#include "ofstreamT.h"

#if defined(VIB_MATERIAL)

/* point generator */
#include "EvenSpacePtsT.h"

#include "FungType.h"
#include "PowerTrig.h"

const double Pi = acos(-1.0);

const int kNumOutputVar = 6;
static const int perm[3][3] = {0,1,2,1,2,0,2,0,1};
static const char* Labels[kNumOutputVar] = {"NT_X", "NT_Y", "NT_Z","IS_X", "IS_Y", "IS_Z"};

using namespace Tahoe;

/* constructors */
AnisoCornea::AnisoCornea(void):
    FSFiberMatT(),
	ParameterInterfaceT("aniso_cornea"),
	fCircle(NULL),
	fPotential(NULL),
	fDistribution(NULL)
{
#ifndef VIB_MATERIAL
	ExceptionT::BadInputValue("AnisoCornea::AnisoCornea", 
		"VIB_MATERIAL must be enabled");
#endif

}

/* destructor */
AnisoCornea::~AnisoCornea(void) 
{ 
	delete fCircle; 
	delete fPotential;
	delete fDistribution;
}

/* strain energy density */
double AnisoCornea::StrainEnergyDensity(void)
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
	ComputeLengths(fFiberStretch);

	/* update potential table */
	fPotential->MapFunction(fI4,fU);

	/* sum contributions */
	double* pU = fU.Pointer();
	double* pj = fjacobian.Pointer();	
	for (int i = 0; i < fI4.Length(); i++)
		energy += (*pU++)*(*pj++);
	
	return energy;
}

int AnisoCornea::NumOutputVariables() const {
	return kNumOutputVar;
}

void AnisoCornea::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void AnisoCornea::ComputeOutput(dArrayT& output)
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
}

/* describe the parameters needed by the interface */
void AnisoCornea::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSFiberMatT::DefineParameters(list);
	
	/* integration points */
	ParameterT points(ParameterT::Integer, "n_points");
	points.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(points);
}

/* information about subordinate parameter lists */
void AnisoCornea::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("matrix_material_params", ParameterListT::Once);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("fibril_potential", ParameterListT::Once);

	/* choice of fibril distribution funcion */
	sub_list.AddSub("fibril_distribution", ParameterListT::Once);
}

/* return the description of the given inline subordinate parameter list */
/*void AnisoCornea::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "fibril_potential")
	{
		order = ParameterListT::Choice;
		sub_lists.AddSub("fung_type", ParameterListT::Once);
	}
	else if (name == "fibril_distribution")
	{
		order = ParameterListT::Choice;
		sub_lists.AddSub("power_trig", ParameterListT::Once);
	}
	else 
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);

}
*/
/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AnisoCornea::NewSub(const StringT& name) const
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
	else if (name == "fibril_potential")
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
		ParameterT phi1(ParameterT::Double, "phi1");

		blend_powertrig.AddParameter(a1);
		blend_powertrig.AddParameter(b1);
		blend_powertrig.AddParameter(c1);
		blend_powertrig.AddParameter(n1);
		blend_powertrig.AddParameter(phi1);
		
		c1.AddLimit(lower);
		c1.AddLimit(upper);
		n1.AddLimit(lower);

		ParameterT a2(ParameterT::Double, "a2");
		ParameterT b2(ParameterT::Double, "b2");
		ParameterT c2(ParameterT::Double, "c2");
		ParameterT n2(ParameterT::Double, "n2");
		ParameterT phi2(ParameterT::Double, "phi2");

		blend_powertrig.AddParameter(a2);
		blend_powertrig.AddParameter(b2);
		blend_powertrig.AddParameter(c2);
		blend_powertrig.AddParameter(n2);
		blend_powertrig.AddParameter(phi2);

		c2.AddLimit(lower);
		c2.AddLimit(upper);
		n2.AddLimit(lower);

		ParameterT r1(ParameterT::Double, "r1");
		ParameterT r2(ParameterT::Double, "r2");
		blend_powertrig.AddParameter(r1);
		blend_powertrig.AddParameter(r2);
		r1.AddLimit(lower);
		r2.AddLimit(lower);

		ParameterT r3(ParameterT::Double, "R_IS"); // circumferential end
		blend_powertrig.AddParameter(r3);
		r3.AddLimit(lower);
		ParameterT r4(ParameterT::Double, "R_NT"); // circumferential end
		blend_powertrig.AddParameter(r4);
		r4.AddLimit(lower);
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
void AnisoCornea::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);
		
	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_material_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fGamma = matrix.GetParameter("bulk_modulus");
	}

	const ParameterListT& potential = list.GetListChoice(*this, "fibril_potential");
	if (potential.Name() == "fung_type")
	{
		double alpha = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential = new FungType(alpha, beta);
		if (!fPotential) throw ExceptionT::kOutOfMemory;
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
		double phi1 = distr.GetParameter("phi1");
		fDistribution = new PowerTrig(a1,b1,c1,n1,phi1); 
		if (!fDistribution) throw ExceptionT::kOutOfMemory;

		a2 = distr.GetParameter("a2");
		b2 = distr.GetParameter("b2");
		c2 = distr.GetParameter("c2");
		n2 = distr.GetParameter("n2");
		xi = distr.GetParameter("phi2");

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
		finhomogeneous = kCornea;
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

    finhomogeneous = kCornea_Mod;
  }
	else if (distr.Name() == "inhomogeneous_distribution") {
    ParameterContainerT inhomo_dist("inhomogeneous_distribution");
    StringT dist_file =  distr.GetParameter("element_distributions_file");
		finhomogeneous = kFile;
	}
		
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fNumFibStress = dSymMatrixT::NumValues(fNumSD-1);
	fNumFibModuli = dSymMatrixT::NumValues(fNumFibStress);
	fFiberStretch.Dimension(fNumSD-1);
	fFiberStress.Dimension(fNumSD-1);
	fFiberMod.Dimension(fNumFibStress);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	Construct();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void AnisoCornea::ComputeMatrixStress(const dSymMatrixT& C, dSymMatrixT& Stress)
{
	/*2pdf{W}{C_IJ} = mu ( del_IJ - I3^-gamma C^-1_IJ)*/
	double I3 = C.Det();
	double I3rg = pow(I3, -fGamma);
	Stress.Inverse(C);
	Stress *= -I3rg*fMu;
	
	Stress[0] += fMu;
	Stress[1] += fMu;
	Stress[2] += fMu;
}

void AnisoCornea::ComputeMatrixMod(const dSymMatrixT& C, dSymMatrixT& Stress, dMatrixT& Mod)
{
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
void AnisoCornea::ComputeFiberStress (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress)
{
	/* stretched bonds */
	ComputeLengths(FiberStretch);

	/* derivatives of the potential */
	fPotential->MapDFunction(fI4, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fI4.Pointer();
	double* pj  = fjacobian.Pointer();

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);

	/* PK2 values in local frame formed by NT and IS orientations*/	
	fFiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_12*/
	
	/*integrate w.r.t in-plane orientation theta*/
	for (int i = 0; i < fI4.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/Pi;
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}
}
	
/*computes integrated moduli in local frame*/
void AnisoCornea::ComputeFiberMod (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress, dMatrixT& FiberMod)
{
	/* stretched bonds */
	ComputeLengths(FiberStretch);
	
	/* derivatives of the potential */
	fPotential->MapDFunction(fI4, fdU);
	fPotential->MapDDFunction(fI4, fddU);	

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pl   = fI4.Pointer();
	double* pj   = fjacobian.Pointer();

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
	
	/* PK2 and Modulus values in local coordinate frame formed by a1 (fNT) and a2 (fIS) */	
	fFiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_12*/
 
	fFiberMod = 0.0;
	double& c11 = FiberMod[0]; /*cf_1111*/ 
	double& c22 = FiberMod[4]; /*cf_2222*/
	double& c33 = FiberMod[8]; /*cf_1212*/
	double& c23 = FiberMod[7]; /*cf_2212*/
	double& c13 = FiberMod[6]; /*cf_1112*/
	double& c12 = FiberMod[3]; /*cf_1122*/

	/*0  3  6
	  1  4  7
	  2  5  8*/
	  
	for (int i = 0; i < fI4.Length(); i++)
	{
		double sfactor =  (*pj)*(*pdU++)/Pi;
		double cfactor = 2.0/Pi*(*pj++)*(*pddU++);
		pl++;

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

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void AnisoCornea::ComputeLengths(const dSymMatrixT& FiberStretch)
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

/***********************************************************************
* Private
***********************************************************************/
/* Initialize angle tables */
void AnisoCornea::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numbonds = points.MajorDim();
		
	/* length table */
	fI4.Dimension(numbonds);
	
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

	if (finhomogeneous == kBlend || finhomogeneous == kCornea || finhomogeneous == kCornea_Mod) 
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

			if (finhomogeneous == kCornea)
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
			if (finhomogeneous == kCornea_Mod)
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
			if (finhomogeneous == kBlend)
			{
				// projected polar coordinates
				double ecc = sqrt(1.0 - r3*r3/(r4*r4)); //eccentricity: sqrt( 1 - RIS^2/RNT^2 )
				r = sqrt( (1.0 - ecc*ecc) * x1*x1 + x2*x2); /*map r to ellipse*/

				C1FunctionT* Dlimbus, *Dcentral;
				Dlimbus = new PowerTrig(a2,b2,c2,n2,xi);
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
			dist_out << "element: " << ielm+1 << " x: " << x1 << " " << x2<<"\n"; 
//           << " r: " << r << " phi: " << phi << " wg: " << wg << "\n";
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

