/* $Id: FSFiberMatT.cpp,v 1.2 2006-08-10 01:46:53 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatT.h"
#include "FSFiberMatSupportT.h"
#include "iArray2DT.h"

const int kNumOutputVar = 6;
static const char* Labels[kNumOutputVar] = {"NT_X", "NT_Y", "NT_Z","IS_X", "IS_Y", "IS_Z"};

using namespace Tahoe;

/* constructor */
FSFiberMatT::FSFiberMatT(void):
	ParameterInterfaceT("fiber_composite_material"),
	fFSFiberMatSupport(NULL)
{

}

/* set the material support or pass NULL to clear */
void FSFiberMatT::SetFSFiberMatSupport(const FSFiberMatSupportT* support)
{
	/* set inherited material support */
	FSSolidMatT::SetFSMatSupport(support);

	fFSFiberMatSupport = support;

}

/* modulus */
const dMatrixT& FSFiberMatT::C_IJKL(void)
{
	/* stretch */
	Compute_C(fC);
	
	/*calculate matrix contribution*/
	ComputeMatrixMod(fC, fStress, fModulus);
	
	/*fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberMod(fFiberStretch, fFiberStress, fFiberMod);
				
	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);

	/* rotate modulus to lab coordinates */
	AssembleFiberModuli(fFiberMod, fModulus);

	return fModulus;
}
	
/* stress */
const dSymMatrixT& FSFiberMatT::S_IJ(void)
{
	
	/* stretch */
	Compute_C(fC);
	
	/*matrix contribution*/
	/*calculate matrix contribution*/
	ComputeMatrixStress(fC, fStress);

	/*fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberStress(fFiberStretch, fFiberStress);
	
	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);

	return(fStress);
}

/* material description */
const dMatrixT& FSFiberMatT::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& FSFiberMatT::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
	return fStress;
}

int FSFiberMatT::NumOutputVariables() const {
	return kNumOutputVar;
}

void FSFiberMatT::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}
	
void FSFiberMatT::ComputeOutput(dArrayT& output)
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
void FSFiberMatT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void FSFiberMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);

	/* dimension work space */
	fNumSD = 3;
	fC.Dimension(fNumSD);
	fModulus.Dimension(dSymMatrixT::NumValues(fNumSD));
	fStress.Dimension(fNumSD);

	/*assuming in-plane fibers*/
	fNumFibStress = dSymMatrixT::NumValues(fNumSD-1);
	fNumFibModuli = dSymMatrixT::NumValues(fNumFibStress);
	
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fFiberStretch.Dimension(fNumSD-1);
	fFiberStress.Dimension(fNumSD-1);
	fFiberMod.Dimension(fNumFibStress);
}


/*********************************************************************************************
 *protected                                                                                  *
 *********************************************************************************************/
void FSFiberMatT::ComputeFiberStretch(const dSymMatrixT& C, dSymMatrixT& Cf)
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
//	cout << "\nFSFiberMatT::AssembleFiberStress: "<<Fibers;
	const double x1 = Fibers(0,0);
	const double x2 = Fibers(0,1);
	const double x3 = Fibers(0,2);
	
	const double y1 = Fibers(1,0);
	const double y2 = Fibers(1,1);
	const double y3 = Fibers(1,2);


	Cf[0] = C[0]*x1*x1 + C[1]*x2*x2 + C[2]*x3*x3 + 2.0*(C[3]*x2*x3 + C[4]*x1*x3 + C[5]*x1*x2);
	Cf[1] = C[0]*y1*y1 + C[1]*y2*y2 + C[2]*y3*y3 + 2.0*(C[3]*y2*y3 + C[4]*y1*y3 + C[5]*y1*y2);
	Cf[2] = C[0]*x1*y1 + C[1]*x2*y2 + C[2]*x3*y3
		+ C[3]*(x2*y3 + y2*x3) + C[4]*(x1*y3 + y1*x3)   + C[5]*(x1*y2 + y1*x2);
}

void FSFiberMatT::AssembleFiberStress(const dSymMatrixT& sigf, dSymMatrixT& sig, const int fillmode)
{
	if (fillmode == dSymMatrixT::kOverwrite)
		sig = 0.0;
	
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
	sig[0] += s1*x1*x1 + s2*y1*y1 + 2.0*s3*x1*y1;
	sig[1] += s1*x2*x2 + s2*y2*y2 + 2.0*s3*x2*y2;
	sig[2] += s1*x3*x3 + s2*y3*y3 + 2.0*s3*x3*y3;
	sig[3] += s1*x2*x3 + s2*y2*y3 + s3*(x2*y3 + x3*y2);
	sig[4] += s1*x1*x3 + s2*y1*y3 + s3*(x1*y3 + x3*y1);
	sig[5] += s1*x1*x2 + s2*y1*y2 + s3*(x1*y2 + x2*y1);
}

void FSFiberMatT::AssembleFiberModuli(const dSymMatrixT& cf, dMatrixT& mod, const int fillmode)
{
	const double& c11 = cf[0];
	const double& c22 = cf[1];
	const double& c33 = cf[2];
	const double& c23 = cf[3];
	const double& c13 = cf[4];
	const double& c12 = cf[5];
	
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();

	const double x1 = Fibers(0,0);
	const double x2 = Fibers(0,1);
	const double x3 = Fibers(0,2);
	
	const double y1 = Fibers(1,0);
	const double y2 = Fibers(1,1);
	const double y3 = Fibers(1,2);
	
	/*Rotate moduli from local frame (defined by fibrils) to global cartesian frame*/
	/*C_IJKL = QIa QJb QKc QLd Cf_abcd*/

	if (fillmode == dSymMatrixT::kOverwrite)
		mod = 0.0;

	/*diagonal terms*/
	/*1111*/
	mod(0,0) += c11*x1*x1*x1*x1 + 4.0*c13*x1*x1*x1*y1 + 2.0*c12*x1*x1*y1*y1 
			+ 4.0*c33*x1*x1*y1*y1 + 4.0*c23*x1*y1*y1*y1 + c22*y1*y1*y1*y1;
	/*2222*/
	mod(1,1) += c11*x2*x2*x2*x2 + 4.0*c13*x2*x2*x2*y2 + 2.0*c12*x2*x2*y2*y2 
			+ 4.0*c33*x2*x2*y2*y2 + 4.0*c23*x2*y2*y2*y2 + c22*y2*y2*y2*y2;
	/*3333*/
	mod(2,2) += c11*x3*x3*x3*x3 + 4.0*c13*x3*x3*x3*y3 + 2.0*c12*x3*x3*y3*y3 
			+ 4.0*c33*x3*x3*y3*y3 + 4.0*c23*x3*y3*y3*y3 + c22*y3*y3*y3*y3;
	/*2323*/
	mod(3,3) += c11*x2*x2*x3*x3 + 2.0*c13*x2*x3*x3*y2 + c33*x3*x3*y2*y2 
			+ 2.0*c13*x2*x2*x3*y3 + 2.0*c12*x2*x3*y2*y3 
			+ 2.0*c33*x2*x3*y2*y3 + 2.0*c23*x3*y2*y2*y3 + c33*x2*x2*y3*y3 
			+ 2.0*c23*x2*y2*y3*y3 + c22*y2*y2*y3*y3;
	/*1313*/
	mod(4,4) += c11*x1*x1*x3*x3 + 2.0*c13*x1*x3*x3*y1 + c33*x3*x3*y1*y1 
			+ 2.0*c13*x1*x1*x3*y3 + 2.0*c12*x1*x3*y1*y3  
			+ 2.0*c33*x1*x3*y1*y3 + 2.0*c23*x3*y1*y1*y3 + c33*x1*x1*y3*y3 
			+ 2.0*c23*x1*y1*y3*y3 + c22*y1*y1*y3*y3;
	/*1212*/
	mod(5,5) += c11*x1*x1*x2*x2 + 2.0*c13*x1*x2*x2*y1 + c33*x2*x2*y1*y1 
			+ 2.0*c13*x1*x1*x2*y2 + 2.0*c12*x1*x2*y1*y2 
			+ 2.0*c33*x1*x2*y1*y2 + 2.0*c23*x2*y1*y1*y2 + c33*x1*x1*y2*y2 
			+ 2.0*c23*x1*y1*y2*y2 + c22*y1*y1*y2*y2;
	
	/*off diagonal terms*/
	/*1122*/
	double C01 = c11*x1*x1*x2*x2 + 2.0*c13*x1*x2*x2*y1 + c12*x2*x2*y1*y1 
			+ 2.0*c13*x1*x1*x2*y2 + 4.0*c33*x1*x2*y1*y2 + 2.0*c23*x2*y1*y1*y2 
			+ c12*x1*x1*y2*y2 + 2.0*c23*x1*y1*y2*y2 + c22*y1*y1*y2*y2;
	mod(1,0) += C01;
	mod(0,1) += C01;

	/*1133*/
	double C02 = c11*x1*x1*x3*x3 + 2.0*c13*x1*x3*x3*y1 + c12*x3*x3*y1*y1 
			+ 2.0*c13*x1*x1*x3*y3 + 4.0*c33*x1*x3*y1*y3 
			+ 2.0*c23*x3*y1*y1*y3 + c12*x1*x1*y3*y3 + 2.0*c23*x1*y1*y3*y3 + c22*y1*y1*y3*y3;
	mod(0,2) += C02;
	mod(2,0) += C02;
	
	/*1123*/	
	double C03 = c11*x1*x1*x2*x3 + 2.0*c13*x1*x2*x3*y1 + c12*x2*x3*y1*y1 
			+ c13*x1*x1*x3*y2 + 2.0*c33*x1*x3*y1*y2 + c23*x3*y1*y1*y2 + c13*x1*x1*x2*y3 
			+ 2.0*c33*x1*x2*y1*y3 + c23*x2*y1*y1*y3 + c12*x1*x1*y2*y3 
			+ 2.0*c23*x1*y1*y2*y3 + c22*y1*y1*y2*y3;
	mod(0,3) += C03;
	mod(3,0) += C03;
	
	/*1113*/
	double C04 = c11*x1*x1*x1*x3 + 3.0*c13*x1*x1*x3*y1 + c12*x1*x3*y1*y1 + 2.0*c33*x1*x3*y1*y1 
			+ c23*x3*y1*y1*y1 + c13*x1*x1*x1*y3 + c12*x1*x1*y1*y3 
			+ 2.0*c33*x1*x1*y1*y3 + 3.0*c23*x1*y1*y1*y3 + c22*y1*y1*y1*y3;
	mod(0,4) += C04;
	mod(4,0) += C04;
	
	/*1112*/
	double C05 = c11*x1*x1*x1*x2 + 3.0*c13*x1*x1*x2*y1 + c12*x1*x2*y1*y1 + 2.0*c33*x1*x2*y1*y1 
			+ c23*x2*y1*y1*y1 + c13*x1*x1*x1*y2 + c12*x1*x1*y1*y2 + 2.0*c33*x1*x1*y1*y2 
			+ 3.0*c23*x1*y1*y1*y2 + c22*y1*y1*y1*y2;
	mod(0,5) += C05;
	mod(5,0) += C05;
	
	/*2233*/
	double C12 = c11*x2*x2*x3*x3 + 2*c13*x2*x3*x3*y2 + c12*x3*x3*y2*y2 + 2.0*c13*x2*x2*x3*y3 
			+ 4.0*c33*x2*x3*y2*y3 + 2.0*c23*x3*y2*y2*y3 + c12*x2*x2*y3*y3 
			+ 2.0*c23*x2*y2*y3*y3 + c22*y2*y2*y3*y3;
	mod(1,2) += C12;
	mod(2,1) += C12;
	
	/*2223*/
	double C13 = c11*x2*x2*x2*x3 + 3.0*c13*x2*x2*x3*y2 + c12*x2*x3*y2*y2 + 2.0*c33*x2*x3*y2*y2 
			+ c23*x3*y2*y2*y2 + c13*x2*x2*x2*y3 + c12*x2*x2*y2*y3 + 2.0*c33*x2*x2*y2*y3 
			+ 3.0*c23*x2*y2*y2*y3 + c22*y2*y2*y2*y3;
	mod(1,3) += C13;
	mod(3,1) += C13;
	
	/*2213*/
	double C14 = c11*x1*x2*x2*x3 + c13*x2*x2*x3*y1 + 2.0*c13*x1*x2*x3*y2 + 2.0*c33*x2*x3*y1*y2 
			+ c12*x1*x3*y2*y2 + c23*x3*y1*y2*y2 + c13*x1*x2*x2*y3 + c12*x2*x2*y1*y3 
			+ 2.0*c33*x1*x2*y2*y3 + 2.0*c23*x2*y1*y2*y3 + c23*x1*y2*y2*y3 + c22*y1*y2*y2*y3;
	mod(1,4) += C14;
	mod(4,1) += C14;
	
	/*2212*/
	double C15 = c11*x1*x2*x2*x2 + c13*x2*x2*x2*y1 + 3.0*c13*x1*x2*x2*y2 + c12*x2*x2*y1*y2 
			+ 2.0*c33*x2*x2*y1*y2 + c12*x1*x2*y2*y2 + 2.0*c33*x1*x2*y2*y2 + 3.0*c23*x2*y1*y2*y2 
			+ c23*x1*y2*y2*y2 + c22*y1*y2*y2*y2;
    mod(1,5) += C15;
	mod(5,1) += C15;
	
	/*3323*/
	double C23 = c11*x2*x3*x3*x3 + c13*x3*x3*x3*y2 + 3.0*c13*x2*x3*x3*y3 + c12*x3*x3*y2*y3 
			+ 2.0*c33*x3*x3*y2*y3 + c12*x2*x3*y3*y3 + 2.0*c33*x2*x3*y3*y3 
			+ 3.0*c23*x3*y2*y3*y3 + c23*x2*y3*y3*y3 + c22*y2*y3*y3*y3;
	mod(2,3) += C23;
	mod(3,2) += C23;
	
	/*3313*/
	double C24 = c11*x1*x3*x3*x3 + c13*x3*x3*x3*y1 + 3.0*c13*x1*x3*x3*y3 + c12*x3*x3*y1*y3 
		+ 2.0*c33*x3*x3*y1*y3 + c12*x1*x3*y3*y3 + 2.0*c33*x1*x3*y3*y3 + 3.0*c23*x3*y1*y3*y3 
		+ c23*x1*y3*y3*y3 + c22*y1*y3*y3*y3;
	mod(2,4) += C24;
	mod(4,2) += C24;
	
	/*3312*/
	double C25 = c11*x1*x2*x3*x3 + c13*x2*x3*x3*y1 + c13*x1*x3*x3*y2 + c12*x3*x3*y1*y2 
			+ 2.0*c13*x1*x2*x3*y3 + 2.0*c33*x2*x3*y1*y3 + 2.0*c33*x1*x3*y2*y3 
			+ 2.0*c23*x3*y1*y2*y3 + c12*x1*x2*y3*y3 + c23*x2*y1*y3*y3 + c23*x1*y2*y3*y3 + c22*y1*y2*y3*y3;
	mod(2,5) += C25;
	mod(5,2) += C25;
	
	/*2313*/
	double C34 = c11*x1*x2*x3*x3 + c13*x2*x3*x3*y1 + c13*x1*x3*x3*y2 + c33*x3*x3*y1*y2 
			+ 2.0*c13*x1*x2*x3*y3 + c12*x2*x3*y1*y3 + c33*x2*x3*y1*y3 + c12*x1*x3*y2*y3 
			+ c33*x1*x3*y2*y3 + 2.0*c23*x3*y1*y2*y3 + c33*x1*x2*y3*y3 + c23*x2*y1*y3*y3 
			+ c23*x1*y2*y3*y3 + c22*y1*y2*y3*y3; 
	mod(3,4) += C34;
	mod(4,3) += C34;
	
	/*2312*/
	double C35 = c11*x1*x2*x2*x3 + c13*x2*x2*x3*y1 + 2.0*c13*x1*x2*x3*y2 + c12*x2*x3*y1*y2 
			+ c33*x2*x3*y1*y2 + c33*x1*x3*y2*y2 + c23*x3*y1*y2*y2 + c13*x1*x2*x2*y3 
			+ c33*x2*x2*y1*y3 + c12*x1*x2*y2*y3 + c33*x1*x2*y2*y3 + 2.0*c23*x2*y1*y2*y3 
			+  c23*x1*y2*y2*y3 + c22*y1*y2*y2*y3; 
    mod(3,5) += C35;
	mod(5,3) += C35;
	
	/*1312*/
	double C45 = c11*x1*x1*x2*x3 + 2.0*c13*x1*x2*x3*y1 + c33*x2*x3*y1*y1 + c13*x1*x1*x3*y2 
			+ c12*x1*x3*y1*y2 + c33*x1*x3*y1*y2 + c23*x3*y1*y1*y2 + c13*x1*x1*x2*y3 + c12*x1*x2*y1*y3 
			+ c33*x1*x2*y1*y3 + c23*x2*y1*y1*y3 + c33*x1*x1*y2*y3 + 2.0*c23*x1*y1*y2*y3 + c22*y1*y1*y2*y3; 
	mod(4,5) += C45;
	mod(5,4) += C45;
	
}
