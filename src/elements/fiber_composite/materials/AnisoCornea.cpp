/* $Id: AnisoCornea.cpp,v 1.1 2006-08-03 01:10:41 thao Exp $ */
/* created: paklein (11/08/1997) */
#include "AnisoCornea.h"

#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"

/* point generator */
#include "EvenSpacePtsT.h"

#include "FungType.h"
#include "PowerTrig.h"

const double Pi = acos(-1.0);
const int kNumOutputVar = 6;
static const char* Labels[kNumOutputVar] = {"NT_X", "NT_Y", "NT_Z","IS_X", "IS_Y", "IS_Z"};

using namespace Tahoe;

/* constructors */
AnisoCornea::AnisoCornea(void):
	ParameterInterfaceT("aniso_cornea"),
	fCircle(NULL),
	fPotential(NULL),
	fDistribution(NULL)
{

}

/* destructor */
AnisoCornea::~AnisoCornea(void) 
{ 
	delete fCircle; 
	delete fPotential;
	delete fDistribution;
}

/* modulus */
const dMatrixT& AnisoCornea::C_IJKL(void)
{
	/* stretch */
	Compute_C(fC);
	
	
	/*matrix contribution*/
	double I3 = fC.Det();
	double I3rg = pow(I3,-fGamma);
	fMatStress.Inverse(fC);
	
	/*2pdf{S_IJ}{C_KL} = 2 mu I_3^-gamma (gamma C^-1_IJ C^-1_KL + 0.5(C^-1_IK C^-1_JL +C^-1_IL+C^-1_JK)*/
	/*modulus*/
	double coeff = 2.0*fMu*I3rg;
	fMatMod.ReducedI_C(fMatStress);
	fMatMod *= coeff;
	
	coeff *= fGamma;
	fMatMod(0,0) += coeff*fMatStress[0]*fMatStress[0];
	fMatMod(0,1) += coeff*fMatStress[0]*fMatStress[1];
	fMatMod(0,2) += coeff*fMatStress[0]*fMatStress[2];
	fMatMod(0,3) += coeff*fMatStress[0]*fMatStress[3];
	fMatMod(0,4) += coeff*fMatStress[0]*fMatStress[4];
	fMatMod(0,5) += coeff*fMatStress[0]*fMatStress[5];
	
	fMatMod(1,0) += coeff*fMatStress[1]*fMatStress[0];
	fMatMod(1,1) += coeff*fMatStress[1]*fMatStress[1];
	fMatMod(1,2) += coeff*fMatStress[1]*fMatStress[2];
	fMatMod(1,3) += coeff*fMatStress[1]*fMatStress[3];
	fMatMod(1,4) += coeff*fMatStress[1]*fMatStress[4];
	fMatMod(1,5) += coeff*fMatStress[1]*fMatStress[5];

	fMatMod(2,0) += coeff*fMatStress[2]*fMatStress[0];
	fMatMod(2,1) += coeff*fMatStress[2]*fMatStress[1];
	fMatMod(2,2) += coeff*fMatStress[2]*fMatStress[2];
	fMatMod(2,3) += coeff*fMatStress[2]*fMatStress[3];
	fMatMod(2,4) += coeff*fMatStress[2]*fMatStress[4];
	fMatMod(2,5) += coeff*fMatStress[2]*fMatStress[5];

	fMatMod(3,0) += coeff*fMatStress[3]*fMatStress[0];
	fMatMod(3,1) += coeff*fMatStress[3]*fMatStress[1];
	fMatMod(3,2) += coeff*fMatStress[3]*fMatStress[2];
	fMatMod(3,3) += coeff*fMatStress[3]*fMatStress[3];
	fMatMod(3,4) += coeff*fMatStress[3]*fMatStress[4];
	fMatMod(3,5) += coeff*fMatStress[3]*fMatStress[5];

	fMatMod(4,0) += coeff*fMatStress[4]*fMatStress[0];
	fMatMod(4,1) += coeff*fMatStress[4]*fMatStress[1];
	fMatMod(4,2) += coeff*fMatStress[4]*fMatStress[2];
	fMatMod(4,3) += coeff*fMatStress[4]*fMatStress[3];
	fMatMod(4,4) += coeff*fMatStress[4]*fMatStress[4];
	fMatMod(4,5) += coeff*fMatStress[4]*fMatStress[5];

	fMatMod(5,0) += coeff*fMatStress[5]*fMatStress[0];
	fMatMod(5,1) += coeff*fMatStress[5]*fMatStress[1];
	fMatMod(5,2) += coeff*fMatStress[5]*fMatStress[2];
	fMatMod(5,3) += coeff*fMatStress[5]*fMatStress[3];
	fMatMod(5,4) += coeff*fMatStress[5]*fMatStress[4];
	fMatMod(5,5) += coeff*fMatStress[5]*fMatStress[5];

	/*stress*/
	/*2pdf{W}{C_IJ} = mu ( del_IJ - I3^-gamma C^-1_IJ)*/
	fMatStress *= -I3rg*fMu;
	
	fMatStress[0] += fMu;
	fMatStress[1] += fMu;
	fMatStress[2] += fMu;

//	cout <<setprecision(8)<< "\nStretch: "<<fC;
//	cout <<setprecision(8)<< "\nMatrix Stress: "<<fMatStress;
//	cout <<setprecision(8)<< "\nMatrix Mod: "<<fMatMod;

	/*fiber contribution*/
	/* stretched bonds */
	ComputeLengths(fC);

//	cout << "\nfLengths: "<<fLengths;
	
	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
	fPotential->MapDDFunction(fLengths, fddU);	

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pl   = fLengths.Pointer();
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
	double& s1 = fFiberStress[0]; /*sf_11*/
	double& s2 = fFiberStress[1]; /*sf_22*/
	double& s3 = fFiberStress[2]; /*sf_12*/
 
	fFiberMod = 0.0;
	double& c11 = fFiberMod[0]; /*cf_1111*/
	double& c22 = fFiberMod[1]; /*cf_2222*/
	double& c33 = fFiberMod[2]; /*cf_1212*/
	double& c23 = fFiberMod[3]; /*cf_2212*/
	double& c13 = fFiberMod[4]; /*cf_1112*/
	double& c12 = fFiberMod[5]; /*cf_1122*/

	for (int i = 0; i < fLengths.Length(); i++)
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
		
	/* rotate stress to lab coordinates */
	/* sig = s11 a1 x a1 + s22 a2 x a2 + s12 (a1 x a2 + a2 x a1) */
	AssembleFiberStress(fFiberStress, fStress);

	//cout  <<setprecision(16)<< "\nunass stress: "<<fFiberStress;
	//cout  <<setprecision(16)<< "\nassemble stress: "<<fStress;

	/* rotate modulus to lab coordinates */
	/* c = c1111 a1 x a1 x a1 x a1 + c1122 a1 x a1 x a2 x a2 + c1112 a1 x a1 x a1 x a2
		 + c2211 a2 x a2 x a1 x a1 + c2222 a2 x a2 x a2 x a2 + c2212 a2 x a2 x a1 x a2
		 + c1211 a1 x a2 x a1 x a1 + c1222 a1 x a2 x a2 x a2 + c1212 a1 x a2 x a1 x a2 */
	AssembleFiberModuli(fFiberMod, fModulus);

//	cout <<setprecision(16)<< "\nunass mod: "<<fFiberMod;
//	cout <<setprecision(16)<< "\nassemble mod: "<<fModulus;

//	cout <<setprecision(8)<< "\nFiber Stress: "<<fStress;
//	cout <<setprecision(8)<< "\nMat Stress: "<<fMatrix;
//	cout <<setprecision(8)<< "\nFiber Mod: "<<fModulus;

	/*add matrix and fiber components*/
	fStress += fMatStress;
	fModulus += fMatMod;

/*	cout << "\nStretch: "<<fC;
	cout << "\nfModulus: "<<fModulus;
	cout << "\nStress: "<<fStress; 
*/	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& AnisoCornea::S_IJ(void)
{
	
	/* stretch */
	Compute_C(fC);
	
	/*matrix contribution*/
	/*2pdf{W}{C_IJ} = mu ( del_IJ - I3^-gamma C^-1_IJ)*/
	/*Note: pdf{I_3}{C} = I_3 C^-1*/
	double I3 = fC.Det();
	double I3rg = pow(I3, -fGamma);
	fMatStress.Inverse(fC);
	fMatStress *= -I3rg*fMu;
	
	fMatStress[0] += fMu;
	fMatStress[1] += fMu;
	fMatStress[2] += fMu;

//	cout <<setprecision(12)<< "\nStretch: "<<fC;
//	cout <<setprecision(12)<< "\nMatrix Stress: "<<fMatStress;
		
	/*fiber contribution*/
	/* stretched bonds */
	ComputeLengths(fC);

//	cout << "\nfLength: "<<fLengths;
	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);
//	cout << "\ndU: "<<fdU;
	/* PK2 values in local frame formed by NT and IS orientations*/	
	fFiberStress = 0.0;
	double& s1 = fFiberStress[0]; /*sf_11*/
	double& s2 = fFiberStress[1]; /*sf_22*/
	double& s3 = fFiberStress[2]; /*sf_12*/
	
	/*integrate w.r.t in-plane orientation theta*/
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/Pi;
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}

	/* rotate stress to lab coordinates */
	/* sig = s1 a1 x a1 + s2 a2 x a2 + s3 (a1 x a2 + a2 x a1) */
	AssembleFiberStress(fFiberStress, fStress);

//	cout << "\n stretch: "<<fC;
//	cout << "\n unassemble Stress: "<<fFiberStress;
//	cout << "\nMat Stress: "<<fMatStress;
//	cout << "\nFiber Stress: "<<fStress;
	
	fStress += fMatStress;
	return(fStress);
}

/* material description */
const dMatrixT& AnisoCornea::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& AnisoCornea::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
	return fStress;
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
	ComputeLengths(fC);

	/* update potential table */
	fPotential->MapFunction(fLengths,fU);

	/* sum contributions */
	double* pU = fU.Pointer();
	double* pj = fjacobian.Pointer();	
	for (int i = 0; i < fLengths.Length(); i++)
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
	FSSolidMatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);

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
	sub_list.AddSub("fibril_potential_choice", ParameterListT::Once);

	/* choice of fibril distribution funcion */
	sub_list.AddSub("fibril_distribution", ParameterListT::Once);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AnisoCornea::NewSub(const StringT& name) const
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
	else if (name == "fibril_potential_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* bound */
		LimitT lower(0.0, LimitT::Lower);
		
		/* exponential functions*/
		ParameterContainerT fungtype("exponential");
		ParameterT alpha(ParameterT::Double, "multiplier_alpha");
		alpha.AddLimit(lower);
		ParameterT beta(ParameterT::Double, "exponent_beta");
		beta.AddLimit(lower);
		fungtype.AddParameter(alpha);
		fungtype.AddParameter(beta);
		choice->AddSub(fungtype);
		return choice;
	}
	else if (name == "fibril_distribution")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* bound */
		LimitT lower(0.0, LimitT::Lower);
		
		/* exponential functions*/
		ParameterContainerT distribution("power_trig");
		ParameterT A(ParameterT::Double, "A");
		ParameterT B(ParameterT::Double, "B");
		ParameterT C(ParameterT::Double, "C");
		ParameterT n(ParameterT::Double, "n");
		ParameterT phi(ParameterT::Double, "p");
		distribution.AddParameter(A);
		distribution.AddParameter(B);
		distribution.AddParameter(C);
		distribution.AddParameter(n);
		distribution.AddParameter(phi);
		choice->AddSub(distribution);
		return choice;
	}
}

/* accept parameter list */
void AnisoCornea::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);

	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_material_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fGamma = matrix.GetParameter("bulk_modulus");
	}

	const ParameterListT& potential = list.GetListChoice(*this, "fibril_potential_choice");
	if (potential.Name() == "exponential")
	{
		double alpha = potential.GetParameter("multiplier_alpha");
		double beta = potential.GetParameter("exponent_beta");
		fPotential = new FungType(alpha,beta);
		if (!fPotential) throw ExceptionT::kOutOfMemory;
	}

	const ParameterListT& distr = list.GetListChoice(*this, "fibril_distribution");
	if (distr.Name() == "power_trig")
	{
		double A = distr.GetParameter("A");
		double B = distr.GetParameter("B");
		double C = distr.GetParameter("C");
		double n = distr.GetParameter("n");
		double phi = distr.GetParameter("p");
		fDistribution = new PowerTrig(A,B,C,n,phi); 
		if (!fDistribution) throw ExceptionT::kOutOfMemory;
	}

	/* dimension work space */
	fNumSD = 3;
	fC.Dimension(fNumSD);
	fModulus.Dimension(dSymMatrixT::NumValues(fNumSD));
	fStress.Dimension(fNumSD);

	fMatStress.Dimension(fNumSD);
	fMatMod.Dimension(dSymMatrixT::NumValues(fNumSD));
	
/*	fIS.Dimension(fNumSD);
	fNT.Dimension(fNumSD);
	fOP.Dimension(fNumSD);
	fQ.Dimension(fNumSD);
*/	
	
	/*assuming in-plane fibers*/
	fNumFibStress = dSymMatrixT::NumValues(2);
	fNumFibModuli = dSymMatrixT::NumValues(fNumFibStress);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	Construct();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void AnisoCornea::ComputeLengths(const dSymMatrixT& C)
{	
	/*Rotate stretch to plane of fibril families
	  C* = QT C Q 
	  Qij = ei aj is rotation matrix, a1 = fNT, a2 = fIS, a3 = fOP
	  fQ(0,0) = fNT[0];
	  fQ(0,1) = fIS[0];
	  fQ(0,2) = fOP[0];

	 fQ(1,0) = fNT[1];
	 fQ(1,1) = fIS[1];
	 fQ(1,2) = fOP[1];

	 fQ(2,0) = fNT[2];
	 fQ(2,1) = fIS[2];
	 fQ(2,2) = fOP[2];
	*/
	
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
//	cout << "\nAnisoCornea::ComputeLengths Fibers: "<<Fibers;

	const double x1 = Fibers(0,0);
	const double x2 = Fibers(0,1);
	const double x3 = Fibers(0,2);
	
	const double y1 = Fibers(1,0);
	const double y2 = Fibers(1,1);
	const double y3 = Fibers(1,2);

	double C11, C22, C12;

	C11 = C[0]*x1*x1 + C[1]*x2*x2 + C[2]*x3*x3 + 2.0*(C[3]*x2*x3 + C[4]*x1*x3 + C[5]*x1*x2);
	C22 = C[0]*y1*y1 + C[1]*y2*y2 + C[2]*y3*y3 + 2.0*(C[3]*y2*y3 + C[4]*y1*y3 + C[5]*y1*y2);
	C12 = C[0]*x1*y1 + C[1]*x2*y2 + C[2]*x3*y3
		+ C[3]*(x2*y3 + y2*x3) + C[4]*(x1*y3 + y1*x3)   + C[5]*(x1*y2 + y1*x2);

	/*calculate fibril lengths                                                        *
	 *I4 = C*:M where M = cos^2 a1 x a1 + sin^2 a2 x a2 + sin cos (a1 x a2 + a2 x a1) */
	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s1 = fStressTable(0);
	double* s2 = fStressTable(1);
	double* s3 = fStressTable(2);
		
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = C11*(*s1++) + C22*(*s2++) + 2.0*C12*(*s3++);
}

/***********************************************************************
* Private
***********************************************************************/
void AnisoCornea::AssembleFiberStress(const dArrayT& sigf, dSymMatrixT& sig)
{
	const double& s1 = sigf[0]; /*sf11*/
	const double& s2 = sigf[1]; /*sf22*/
	const double& s3 = sigf[2]; /*sf12*/
	
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
//	cout << "\nFibers: "<<Fibers;

	const double x1 = Fibers(0,0);
	const double x2 = Fibers(0,1);
	const double x3 = Fibers(0,2);
	
	const double y1 = Fibers(1,0);
	const double y2 = Fibers(1,1);
	const double y3 = Fibers(1,2);
//	const double y3 = 0.0;
	
	/* sig = s11 a1 x a1 + s22 a2 x a2 + s12 (a1 x a2 + a2 x a1) 
	   sig = Q.sf.QT                                            */
	sig[0] = s1*x1*x1 + s2*y1*y1 + 2.0*s3*x1*y1;
	sig[1] = s1*x2*x2 + s2*y2*y2 + 2.0*s3*x2*y2;
	sig[2] = s1*x3*x3 + s2*y3*y3 + 2.0*s3*x3*y3;
	sig[3] = s1*x2*x3 + s2*y2*y3 + s3*(x2*y3 + x3*y2);
	sig[4] = s1*x1*x3 + s2*y1*y3 + s3*(x1*y3 + x3*y1);
	sig[5] = s1*x1*x2 + s2*y1*y2 + s3*(x1*y2 + x2*y1);
}

void AnisoCornea::AssembleFiberModuli(const dSymMatrixT& cf, dMatrixT& mod)
{
	const double& c11 = cf[0];
	const double& c22 = cf[1];
	const double& c33 = cf[2];
	const double& c23 = cf[3];
	const double& c13 = cf[4];
	const double& c12 = cf[5];
	
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
//	cout << "\nFibers: "<<Fibers;
	const double x1 = Fibers(0,0);
	const double x2 = Fibers(0,1);
	const double x3 = Fibers(0,2);
//	const double x3 = 0.0;
	
	const double y1 = Fibers(1,0);
	const double y2 = Fibers(1,1);
	const double y3 = Fibers(1,2);
//	const double y3 = 0.0;
	
		/* c = c1111 a1 x a1 x a1 x a1 + c1122 a1 x a1 x a2 x a2 + c1112 a1 x a1 x a1 x a2
		 + c2211 a2 x a2 x a1 x a1 + c2222 a2 x a2 x a2 x a2 + c2212 a2 x a2 x a1 x a2
		 + c1211 a1 x a2 x a1 x a1 + c1222 a1 x a2 x a2 x a2 + c1212 a1 x a2 x a1 x a2 */

	/*diagonal terms*/
	/*1111*/
	mod(0,0) = c11*x1*x1*x1*x1 + 4.0*c13*x1*x1*x1*y1 + 2.0*c12*x1*x1*y1*y1 
			+ 4.0*c33*x1*x1*y1*y1 + 4.0*c23*x1*y1*y1*y1 + c22*y1*y1*y1*y1;
	/*2222*/
	mod(1,1) = c11*x2*x2*x2*x2 + 4.0*c13*x2*x2*x2*y2 + 2.0*c12*x2*x2*y2*y2 
			+ 4.0*c33*x2*x2*y2*y2 + 4.0*c23*x2*y2*y2*y2 + c22*y2*y2*y2*y2;
	/*3333*/
	mod(2,2) = c11*x3*x3*x3*x3 + 4.0*c13*x3*x3*x3*y3 + 2.0*c12*x3*x3*y3*y3 
			+ 4.0*c33*x3*x3*y3*y3 + 4.0*c23*x3*y3*y3*y3 + c22*y3*y3*y3*y3;
	/*2323*/
	mod(3,3) = c11*x2*x2*x3*x3 + 2.0*c13*x2*x3*x3*y2 + c33*x3*x3*y2*y2 
			+ 2.0*c13*x2*x2*x3*y3 + 2.0*c12*x2*x3*y2*y3 
			+ 2.0*c33*x2*x3*y2*y3 + 2.0*c23*x3*y2*y2*y3 + c33*x2*x2*y3*y3 
			+ 2.0*c23*x2*y2*y3*y3 + c22*y2*y2*y3*y3;
	/*1313*/
	mod(4,4) = c11*x1*x1*x3*x3 + 2.0*c13*x1*x3*x3*y1 + c33*x3*x3*y1*y1 
			+ 2.0*c13*x1*x1*x3*y3 + 2.0*c12*x1*x3*y1*y3  
			+ 2.0*c33*x1*x3*y1*y3 + 2.0*c23*x3*y1*y1*y3 + c33*x1*x1*y3*y3 
			+ 2.0*c23*x1*y1*y3*y3 + c22*y1*y1*y3*y3;
	/*1212*/
	mod(5,5) = c11*x1*x1*x2*x2 + 2.0*c13*x1*x2*x2*y1 + c33*x2*x2*y1*y1 
			+ 2.0*c13*x1*x1*x2*y2 + 2.0*c12*x1*x2*y1*y2 
			+ 2.0*c33*x1*x2*y1*y2 + 2.0*c23*x2*y1*y1*y2 + c33*x1*x1*y2*y2 
			+ 2.0*c23*x1*y1*y2*y2 + c22*y1*y1*y2*y2;
	
	/*upper diagonal terms*/
	/*1122*/
	mod(0,1) = c11*x1*x1*x2*x2 + 2.0*c13*x1*x2*x2*y1 + c12*x2*x2*y1*y1 
			+ 2.0*c13*x1*x1*x2*y2 + 4.0*c33*x1*x2*y1*y2 + 2.0*c23*x2*y1*y1*y2 
			+ c12*x1*x1*y2*y2 + 2.0*c23*x1*y1*y2*y2 + c22*y1*y1*y2*y2;
	/*1133*/
	mod(0,2) = c11*x1*x1*x3*x3 + 2.0*c13*x1*x3*x3*y1 + c12*x3*x3*y1*y1 
			+ 2.0*c13*x1*x1*x3*y3 + 4.0*c33*x1*x3*y1*y3 
			+ 2.0*c23*x3*y1*y1*y3 + c12*x1*x1*y3*y3 + 2.0*c23*x1*y1*y3*y3 + c22*y1*y1*y3*y3;
	/*1123*/	
	mod(0,3) = c11*x1*x1*x2*x3 + 2.0*c13*x1*x2*x3*y1 + c12*x2*x3*y1*y1 
			+ c13*x1*x1*x3*y2 + 2.0*c33*x1*x3*y1*y2 + c23*x3*y1*y1*y2 + c13*x1*x1*x2*y3 
			+ 2.0*c33*x1*x2*y1*y3 + c23*x2*y1*y1*y3 + c12*x1*x1*y2*y3 
			+ 2.0*c23*x1*y1*y2*y3 + c22*y1*y1*y2*y3;
	/*1113*/
	mod(0,4) = c11*x1*x1*x1*x3 + 3.0*c13*x1*x1*x3*y1 + c12*x1*x3*y1*y1 + 2.0*c33*x1*x3*y1*y1 
			+ c23*x3*y1*y1*y1 + c13*x1*x1*x1*y3 + c12*x1*x1*y1*y3 
			+ 2.0*c33*x1*x1*y1*y3 + 3.0*c23*x1*y1*y1*y3 + c22*y1*y1*y1*y3;
	/*1112*/
	mod(0,5) = c11*x1*x1*x1*x2 + 3.0*c13*x1*x1*x2*y1 + c12*x1*x2*y1*y1 + 2.0*c33*x1*x2*y1*y1 
			+ c23*x2*y1*y1*y1 + c13*x1*x1*x1*y2 + c12*x1*x1*y1*y2 + 2.0*c33*x1*x1*y1*y2 
			+ 3.0*c23*x1*y1*y1*y2 + c22*y1*y1*y1*y2;
   
	/*2233*/
	mod(1,2) = c11*x2*x2*x3*x3 + 2*c13*x2*x3*x3*y2 + c12*x3*x3*y2*y2 + 2.0*c13*x2*x2*x3*y3 
			+ 4.0*c33*x2*x3*y2*y3 + 2.0*c23*x3*y2*y2*y3 + c12*x2*x2*y3*y3 
			+ 2.0*c23*x2*y2*y3*y3 + c22*y2*y2*y3*y3;
	/*2223*/
	mod(1,3) = c11*x2*x2*x2*x3 + 3.0*c13*x2*x2*x3*y2 + c12*x2*x3*y2*y2 + 2.0*c33*x2*x3*y2*y2 
			+ c23*x3*y2*y2*y2 + c13*x2*x2*x2*y3 + c12*x2*x2*y2*y3 + 2.0*c33*x2*x2*y2*y3 
			+ 3.0*c23*x2*y2*y2*y3 + c22*y2*y2*y2*y3;
	/*2213*/
	mod(1,4) = c11*x1*x2*x2*x3 + c13*x2*x2*x3*y1 + 2.0*c13*x1*x2*x3*y2 + 2.0*c33*x2*x3*y1*y2 + c12*x1*x3*y2*y2 + c23*x3*y1*y2*y2 + 
   c13*x1*x2*x2*y3 + c12*x2*x2*y1*y3 + 2*c33*x1*x2*y2*y3 + 2*c23*x2*y1*y2*y3 + c23*x1*y2*y2*y3 + c22*y1*y2*y2*y3;
	/*2212*/
	mod(1,5) = c11*x1*x2*x2*x2 + c13*x2*x2*x2*y1 + 3.0*c13*x1*x2*x2*y2 + c12*x2*x2*y1*y2 
			+ 2.0*c33*x2*x2*y1*y2 + c12*x1*x2*y2*y2 + 2.0*c33*x1*x2*y2*y2 + 3.0*c23*x2*y1*y2*y2 + c23*x1*y2*y2*y2 + c22*y1*y2*y2*y2;
   
	/*3323*/
	mod(2,3) = c11*x2*x3*x3*x3 + c13*x3*x3*x3*y2 + 3.0*c13*x2*x3*x3*y3 + c12*x3*x3*y2*y3 
			+ 2.0*c33*x3*x3*y2*y3 + c12*x2*x3*y3*y3 + 2.0*c33*x2*x3*y3*y3 
			+ 3.0*c23*x3*y2*y3*y3 + c23*x2*y3*y3*y3 + c22*y2*y3*y3*y3;
	/*3313*/
	mod(2,4) = c11*x1*x3*x3*x3 + c13*x3*x3*x3*y1 + 3.0*c13*x1*x3*x3*y3 + c12*x3*x3*y1*y3 
		+ 2.0*c33*x3*x3*y1*y3 + c12*x1*x3*y3*y3 + 2.0*c33*x1*x3*y3*y3 + 3.0*c23*x3*y1*y3*y3 
		+ c23*x1*y3*y3*y3 + c22*y1*y3*y3*y3;
	/*3312*/
	mod(2,5) = c11*x1*x2*x3*x3 + c13*x2*x3*x3*y1 + c13*x1*x3*x3*y2 + c12*x3*x3*y1*y2 
			+ 2.0*c13*x1*x2*x3*y3 + 2.0*c33*x2*x3*y1*y3 + 2.0*c33*x1*x3*y2*y3 
			+ 2.0*c23*x3*y1*y2*y3 + c12*x1*x2*y3*y3 + c23*x2*y1*y3*y3 + c23*x1*y2*y3*y3 + c22*y1*y2*y3*y3;

	/*2313*/
	mod(3,4) = c11*x1*x2*x3*x3 + c13*x2*x3*x3*y1 + c13*x1*x3*x3*y2 + c33*x3*x3*y1*y2 
			+ 2.0*c13*x1*x2*x3*y3 + c12*x2*x3*y1*y3 + c33*x2*x3*y1*y3 + c12*x1*x3*y2*y3 
			+ c33*x1*x3*y2*y3 + 2.0*c23*x3*y1*y2*y3 + c33*x1*x2*y3*y3 + c23*x2*y1*y3*y3 + c23*x1*y2*y3*y3 + 
   c22*y1*y2*y3*y3; 
	/*2312*/
	mod(3,5) = c11*x1*x2*x2*x3 + c13*x2*x2*x3*y1 + 2.0*c13*x1*x2*x3*y2 + c12*x2*x3*y1*y2 
			+ c33*x2*x3*y1*y2 + c33*x1*x3*y2*y2 + c23*x3*y1*y2*y2 + c13*x1*x2*x2*y3 
			+ c33*x2*x2*y1*y3 + c12*x1*x2*y2*y3 + c33*x1*x2*y2*y3 + 2.0*c23*x2*y1*y2*y3 
			+  c23*x1*y2*y2*y3 + c22*y1*y2*y2*y3; 
   
	/*1312*/
	mod(4,5) = c11*x1*x1*x2*x3 + 2.0*c13*x1*x2*x3*y1 + c33*x2*x3*y1*y1 + c13*x1*x1*x3*y2 
			+ c12*x1*x3*y1*y2 + c33*x1*x3*y1*y2 + c23*x3*y1*y1*y2 + c13*x1*x1*x2*y3 + c12*x1*x2*y1*y3 
			+ c33*x1*x2*y1*y3 + c23*x2*y1*y1*y3 + c33*x1*x1*y2*y3 + 2.0*c23*x1*y1*y2*y3 + c22*y1*y1*y2*y3; 
	
	/*lower diagonal terms*/
	mod(1,0) = mod(0,1);
	mod(2,0) = mod(0,2);
	mod(2,1) = mod(1,2);
	mod(3,0) = mod(0,3);
	mod(3,1) = mod(1,3);
	mod(3,2) = mod(2,3);
	mod(4,0) = mod(0,4);
	mod(4,1) = mod(1,4);
	mod(4,2) = mod(2,4);
	mod(4,3) = mod(3,4);
	mod(5,0) = mod(0,5);
	mod(5,1) = mod(1,5);
	mod(5,2) = mod(2,5);
	mod(5,3) = mod(3,5);
	mod(5,4) = mod(4,5);
}

/* Initialize angle tables */
void AnisoCornea::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numbonds = points.MajorDim();
	
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fFiberStress.Dimension(fNumFibStress);
	fFiberMod.Dimension(fNumFibStress);
	
	/* length table */
	fLengths.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* jacobian table */
	fjacobian.Dimension(numbonds);
	/* angles table */
	fangles.Dimension(numbonds);

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

	/* fetch jacobians */
	fjacobian = fCircle->Jacobians();
	/*multiply by distribution function*/
	fangles = fCircle->CircleAngles(0.0);
//	cout << "\nfangles: "<<fangles;
	fDistribution->MapFunction(fangles,fU);
	const double* pd = fU.Pointer();
	double* pj = fjacobian.Pointer();	
	for (int i = 0; i < fangles.Length(); i++)
		pj[i] = fU[i]*pj[i];

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
