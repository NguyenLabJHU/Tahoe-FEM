/* $Id: AnisoFiber3D.cpp,v 1.1 2007-12-19 23:35:54 thao Exp $ */
/* created: paklein (11/08/1997) */

#include "AnisoFiber3D.h"
#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"

#include "FungType.h"

#ifdef VIB_MATERIAL
/* point generators */
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"
#include "FCCPtsT.h"

#include "PowerTrig.h"

const double Pi = acos(-1.0);

const int kNumOutputVar = 6;
static const char* Labels[kNumOutputVar] = {"NT_X", "NT_Y", "NT_Z","IS_X", "IS_Y", "IS_Z"};

using namespace Tahoe;

/* constructors */
AnisoFiber3D::AnisoFiber3D(void):
    FSFiberMatT(),
	ParameterInterfaceT("aniso_fiber_3D"),
	fSphere(NULL),
//	fDistribution(NULL),
	fPotential(NULL)
{
#ifndef VIB_MATERIAL
	ExceptionT::BadInputValue("AnisoFiber3D::AnisoFiber3D", 
		"VIB_MATERIAL must be enabled");
#endif

}

/* destructor */
AnisoFiber3D::~AnisoFiber3D(void) 
{ 
	delete fSphere; 
	delete fPotential; 
//	delete fDistribution;
}

/* strain energy density */
double AnisoFiber3D::StrainEnergyDensity(void)
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

int AnisoFiber3D::NumOutputVariables() const {
	return kNumOutputVar;
}

void AnisoFiber3D::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void AnisoFiber3D::ComputeOutput(dArrayT& output)
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
void AnisoFiber3D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSFiberMatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void AnisoFiber3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("matrix_material_params", ParameterListT::Once);

	sub_list.AddSub("sphere_integration_choice", ParameterListT::Once, true);

	sub_list.AddSub("3D_fiber_potential", ParameterListT::Once, true);

	/* choice of fibril distribution funcion */
//	sub_list.AddSub("3D_fiber_distribution", ParameterListT::Once);
}

/* return the description of the given inline subordinate parameter list */
/*void AnisoFiber3D::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
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
ParameterInterfaceT* AnisoFiber3D::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = FSFiberMatT::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "sphere_integration_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* bound */
		LimitT lower(1, LimitT::LowerInclusive);
  
		/* like the grid on a sphere */
		ParameterContainerT lat_long("latitude_longitude");
		ParameterT n_phi(ParameterT::Integer, "n_phi");
		n_phi.AddLimit(lower);
		ParameterT n_theta(ParameterT::Integer, "n_theta");
		n_theta.AddLimit(lower);
		lat_long.AddParameter(n_phi);
		lat_long.AddParameter(n_theta);
		choice->AddSub(lat_long);
    
		/* icosahedral points */
		ParameterContainerT ico("icosahedral");
		ParameterT ico_points(ParameterT::Integer, "points");
		ico_points.AddLimit(6, LimitT::Only);
		ico_points.AddLimit(10, LimitT::Only);
		ico_points.AddLimit(40, LimitT::Only);
		ico_points.AddLimit(160, LimitT::Only);
		ico.AddParameter(ico_points);
		choice->AddSub(ico);
    
		/* FCC point arrangement */
		ParameterContainerT fcc("fcc_points");
		ParameterT n_shells(ParameterT::Integer, "shells");
		n_shells.AddLimit(lower);
		fcc.AddParameter(n_shells);
		fcc.AddParameter(ParameterT::Double, "nearest_neighbor_distance");
		choice->AddSub(fcc);
		return choice;
	}
	else if (name == "3D_fiber_potential")
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
  // NOTE : this definition _masks_ the definition in AnisoViscoCornea
/*	else if (name == "3D_fiber_distribution")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT powertrig("power_trig");		
		LimitT lower(0.0, LimitT::Lower);
		LimitT upper(1.0, LimitT::Upper);

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

		powertrig.SetDescription("f(theta) = a*cos(theta+phi)^n + b*sin(theta+phi)^n + c");	
		choice->AddSub(powertrig);

		return(choice);
	}
*/
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
	else /* inherited */
	{
		/*Need to create new sub*/
		return ParameterInterfaceT::NewSub(name);
	}
}

/* accept parameter list */
void AnisoFiber3D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);

	fNumFibStress = dSymMatrixT::NumValues(fNumSD);
	fNumFibModuli = fNumFibStress*fNumFibStress;
	
	/* allocate memory */
	fFiberStretch.Dimension(fNumSD);
	fFiberStress.Dimension(fNumSD);
	fFiberMod.Dimension(fNumFibStress);
	
	/*initializing some parameters for cornea_mod fibril distribution*/
	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_material_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fGamma = matrix.GetParameter("bulk_modulus");
	}

	const ParameterListT& potential = list.GetListChoice(*this, "3D_fiber_potential");
	if (potential.Name() == "fung_type")
	{
		double alpha = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential = new FungType(alpha,beta);
		if (!fPotential) throw ExceptionT::kOutOfMemory;
	}
	
/*	const ParameterListT& distr = list.GetListChoice(*this, "3D_fiber_distribution");
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
  else
    ExceptionT::GeneralFail("AnisoFiber3D::TakeParameterList", "unrecognized potential \"%s\"", potential.Name().Pointer());
*/
	const ParameterListT& points = list.GetListChoice(*this, "sphere_integration_choice");
	if (points.Name() == "latitude_longitude")
	{
		int n_phi = points.GetParameter("n_phi");
		int n_theta = points.GetParameter("n_theta");
		fSphere = new LatLongPtsT(n_phi, n_theta);
    }
	else if (points.Name() == "icosahedral")
	{
		int np = points.GetParameter("points");
		fSphere = new IcosahedralPtsT(np);
    }
	else if (points.Name() == "fcc_points")
	{
		int num_shells = points.GetParameter("shells");
		double bond_length = points.GetParameter("nearest_neighbor_distance");
		fSphere = new FCCPtsT(num_shells, bond_length);
    }
	else
		ExceptionT::GeneralFail("AnisoFiber3D::TakeParameterList", "unrecognized point scheme \"%s\"", points.Name().Pointer());

  /* set tables */
  Construct();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void AnisoFiber3D::ComputeMatrixStress(const dSymMatrixT& C, dSymMatrixT& Stress)
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

void AnisoFiber3D::ComputeMatrixMod(const dSymMatrixT& C, dSymMatrixT& Stress, dMatrixT& Mod)
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
void AnisoFiber3D::ComputeFiberStress (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress)
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
	double *p4  = fStressTable(3);
	double *p5  = fStressTable(4);
	double *p6  = fStressTable(5);

	/* PK2 values in local frame formed by NT and IS orientations*/	
	fFiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_33*/
	double& s4 = FiberStress[3]; /*sf_33*/
	double& s5 = FiberStress[4]; /*sf_33*/
	double& s6 = FiberStress[5]; /*sf_33*/
	
	/*integrate w.r.t in-plane orientation theta*/
	for (int i = 0; i < fI4.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/Pi;
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
		s4 += factor*(*p4++);
		s5 += factor*(*p5++);
		s6 += factor*(*p6++);
	}
}
	
/*computes integrated moduli in local frame*/
void AnisoFiber3D::ComputeFiberMod (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress, dMatrixT& FiberMod)
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

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);
	double *p4  = fStressTable(3);
	double *p5  = fStressTable(4);
	double *p6  = fStressTable(5);

	/* PK2 values in local frame formed by NT and IS orientations*/	
	fFiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_33*/
	double& s4 = FiberStress[3]; /*sf_23*/
	double& s5 = FiberStress[4]; /*sf_13*/
	double& s6 = FiberStress[5]; /*sf_12*/

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc44  = fModuliTable(3);
	double* pc55  = fModuliTable(4);
	double* pc66  = fModuliTable(5);

	double* pc12 = fModuliTable(6);
	double* pc13 = fModuliTable(7);
	double* pc14 = fModuliTable(8);
	double* pc15 = fModuliTable(9);
	double* pc16 = fModuliTable(10);

	double* pc23 = fModuliTable(11);
	double* pc24 = fModuliTable(12);
	double* pc25 = fModuliTable(13);
	double* pc26 = fModuliTable(14);

	double* pc34 = fModuliTable(15);
	double* pc35 = fModuliTable(16);
	double* pc36 = fModuliTable(17);

	double* pc45 = fModuliTable(18);
	double* pc46 = fModuliTable(19);


	double* pc56 = fModuliTable(20);
	
	/* PK2 and Modulus values in local coordinate frame formed by a1 (fNT) and a2 (fIS) */	 
	fFiberMod = 0.0;
	double& c11 = FiberMod(0,0); /*cf_1111*/ 
	double& c22 = FiberMod(1,1); /*cf_2222*/
	double& c33 = FiberMod(2,2); /*cf_3333*/
	double& c44 = FiberMod(3,3); /*cf_2323*/
	double& c55 = FiberMod(4,4); /*cf_1313*/
	double& c66 = FiberMod(5,5); /*cf_1212*/
	
	double& c12 = FiberMod(0,1); /*cf_1122*/
	double& c13 = FiberMod(0,2); /*cf_1133*/
	double& c14 = FiberMod(0,3); /*cf_1123*/
	double& c15 = FiberMod(0,4); /*cf_1113*/
	double& c16 = FiberMod(0,5); /*cf_1112*/

	double& c23 = FiberMod(1,2); /*cf_233*/
	double& c24 = FiberMod(1,3); /*cf_2223*/
	double& c25 = FiberMod(1,4); /*cf_2213*/
	double& c26 = FiberMod(1,5); /*cf_2212*/
	  
	double& c34 = FiberMod(2,3); /*cf_3323*/
	double& c35 = FiberMod(2,4); /*cf_3313*/
	double& c36 = FiberMod(2,5); /*cf_3312*/

	double& c45 = FiberMod(3,4); /*cf_2313*/
	double& c46 = FiberMod(3,5); /*cf_2312*/

	double& c56 = FiberMod(4,5); /*cf_1312*/

	for (int i = 0; i < fI4.Length(); i++)
	{
		double sfactor =  (*pj)*(*pdU++)/Pi;
		double cfactor = 2.0/Pi*(*pj++)*(*pddU++);
		pl++;

		s1 += sfactor*(*p1++);
		s2 += sfactor*(*p2++);
		s3 += sfactor*(*p3++);
		s4 += sfactor*(*p4++);
		s5 += sfactor*(*p5++);
		s6 += sfactor*(*p6++);

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c33 += cfactor*(*pc33++);
		c44 += cfactor*(*pc44++);
		c55 += cfactor*(*pc55++);
		c66 += cfactor*(*pc66++);

		c12 += cfactor*(*pc12++);
		c13 += cfactor*(*pc13++);
		c14 += cfactor*(*pc14++);
		c15 += cfactor*(*pc15++);
		c16 += cfactor*(*pc16++);

		c23 += cfactor*(*pc23++);
		c24 += cfactor*(*pc24++);
		c25 += cfactor*(*pc25++);
		c26 += cfactor*(*pc26++);

		c34 += cfactor*(*pc34++);
		c35 += cfactor*(*pc35++);
		c36 += cfactor*(*pc36++);

		c45 += cfactor*(*pc45++);
		c46 += cfactor*(*pc46++);

		c56 += cfactor*(*pc56++);

	}
	/*symmetric modulus*/
	FiberMod.CopySymmetric();
}

void AnisoFiber3D::ComputeLengths(const dSymMatrixT& FiberStretch)
{	
	/*calculate fibril lengths                                                        *
	 *I4 = C*:M where M = cos^2 a1 x a1 + sin^2 a2 x a2 + sin cos (a1 x a2 + a2 x a1) */

	const double& C11 = FiberStretch[0];
	const double& C22 = FiberStretch[1];
	const double& C33 = FiberStretch[2];
	const double& C23 = FiberStretch[3];
	const double& C13 = FiberStretch[4];
	const double& C12 = FiberStretch[5];

	/* initialize kernel pointers */
	double* pl = fI4.Pointer();

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);
	double *p4  = fStressTable(3);
	double *p5  = fStressTable(4);
	double *p6  = fStressTable(5);
		
	for (int i = 0; i < fI4.Length(); i++)
		*pl++ = C11*(*p1++) + C22*(*p2++) + C33*(*p3++) + 2.0*C23*(*p4++) + 2.0*C13*(*p5++) + 2.0*C12*(*p6++);
		
}

/***********************************************************************
* Private
***********************************************************************/
/* Initialize angle tables */
void AnisoFiber3D::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fSphere->SpherePoints(0.0, 0.0);
	int numbonds = points.MajorDim();
	
	/* allocate memory */
	/*3D fiber stress and modulus*/
	fFiberStress.Dimension(fNumSD);
	fFiberMod.Dimension(fNumFibStress);
	
	/* length table */
	fI4.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* jacobian table */
	fjacobian.Dimension(numbonds);
	/* angles table */
//	fangles.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumFibStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(21, numbonds);	

	/* set pointers */
	double *s1 = fStressTable(0);
	double *s2 = fStressTable(1);
	double *s3 = fStressTable(2);
	double *s4 = fStressTable(3);
	double *s5 = fStressTable(4);
	double *s6 = fStressTable(5);

	double *c11 = fModuliTable(0);
	double *c22 = fModuliTable(1);
	double *c33 = fModuliTable(2);
	double *c44 = fModuliTable(3);
	double *c55 = fModuliTable(4);
	double *c66 = fModuliTable(5);

	double *c12 = fModuliTable(6);
	double *c13 = fModuliTable(7);
	double *c14 = fModuliTable(8);
	double *c15 = fModuliTable(9);
	double *c16 = fModuliTable(10);

	double *c23 = fModuliTable(11);
	double *c24 = fModuliTable(12);
	double *c25 = fModuliTable(13);
	double *c26 = fModuliTable(14);

	double *c34 = fModuliTable(15);
	double *c35 = fModuliTable(16);
	double *c36 = fModuliTable(17);

	double *c45 = fModuliTable(18);
	double *c46 = fModuliTable(19);

	double *c56 = fModuliTable(20);

//	fjacobian = fSphere->Jacobians(0.0, fDistribution);
	fjacobian = fSphere->Jacobians();

	for (int i = 0; i < numbonds; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);
		double xsi1 = xsi[0];
		double xsi2 = xsi[1];
		double xsi3 = xsi[2];
		
		/* stress angle tables */
		s1[i] = xsi1*xsi1;      
		s2[i] = xsi2*xsi2;
		s3[i] = xsi3*xsi3;
		s4[i] = xsi2*xsi3;
		s5[i] = xsi1*xsi3;
		s6[i] = xsi1*xsi2;
	
		/* moduli angle tables */
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];
		c33[i] = s3[i]*s3[i];
		c44[i] = s4[i]*s4[i];
		c55[i] = s5[i]*s5[i];
		c66[i] = s6[i]*s6[i];
		
		c12[i] = s1[i]*s2[i];
		c13[i] = s1[i]*s3[i];
		c14[i] = s1[i]*s4[i];
		c15[i] = s1[i]*s5[i];
		c16[i] = s1[i]*s6[i];

		c23[i] = s2[i]*s3[i];
		c24[i] = s2[i]*s4[i];
		c25[i] = s2[i]*s5[i];
		c26[i] = s2[i]*s6[i];

		c34[i] = s3[i]*s4[i];
		c35[i] = s3[i]*s5[i];
		c36[i] = s3[i]*s6[i];

		c45[i] = s4[i]*s5[i];
		c46[i] = s4[i]*s6[i];

		c56[i] = s5[i]*s6[i];

	}
}

#endif /*VIB_MATERIAL*/

