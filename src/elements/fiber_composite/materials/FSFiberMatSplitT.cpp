/* $Id: FSFiberMatSplitT.cpp,v 1.1 2010-06-24 14:35:49 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatSplitT.h"
#include "FSFiberMatSupportT.h"
#include "iArray2DT.h"

const double third = 1.0/3.0;
using namespace Tahoe;

/* constructor */
FSFiberMatSplitT::FSFiberMatSplitT(void):
	FSFiberMatT(),
	ParameterInterfaceT("fiber_composite_material_split")
{
}

/* modulus */
const dMatrixT& FSFiberMatSplitT::C_IJKL(void)
{
	/* stretch */
	Compute_C(fC);
	fI3 = fC.Det();
	double I3rthird = pow(fI3,-third);

	/*deviatoric part*/
	fCbar = fC;
	fCbar *= I3rthird;
	
	fInverse.Inverse(fC);

	/*calculate matrix contribution*/
	fSbar = 0.0;
	fModbar = 0.0;  /*initialize*/
	ComputeDevMatrixMod(fCbar,  fSbar, fModbar);
	
	/*fiber contribution*/
	ComputeFiberStretch(fCbar, fFiberStretch);
	ComputeFiberMod(fFiberStretch, fFiberStress, fFiberMod);

	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fSbar);
	/* rotate modulus to lab coordinates */
	AssembleFiberModuli(fFiberMod, fModbar);
	
	/*Deviatoric stress*/
	/*Sdev = I3^(-1/3)*(Sbar - 1/3 (Sbar:C) C^-1)*/
	fStress = fSbar;
	
	double coeff = fSbar.ScalarProduct(fC);
	coeff *= third;
	fStress[0] -= coeff*fInverse[0];
	fStress[1] -= coeff*fInverse[1];
	fStress[2] -= coeff*fInverse[2];
	fStress[3] -= coeff*fInverse[3];
	fStress[4] -= coeff*fInverse[4];
	fStress[5] -= coeff*fInverse[5];

	fStress *= I3rthird;

	double p =  ComputeVolMatrixStress(fI3);
		
	/*modulus*/
	
	/*I3^-2/3 Cbar_dev*/
	fModulus = fModbar; 
	fModulus *= I3rthird*I3rthird;

	/*2/3 I3^-1/3(Sbar:C) (-pdf(C^-1)(C))*/
	coeff = fSbar.ScalarProduct(fC);
	coeff *= 2.0*third*I3rthird;

	fMat.ReducedI_C(fInverse);
	coeff -= 2.0*p;
	
	
	fMat *= coeff;
	fModulus(0,0) += fMat(0,0);
	fModulus(0,1) += fMat(0,1);
	fModulus(0,2) += fMat(0,2);
	fModulus(0,3) += fMat(0,3);
	fModulus(0,4) += fMat(0,4);
	fModulus(0,5) += fMat(0,5);
	//
	fModulus(1,0) += fMat(1,0);
	fModulus(1,1) += fMat(1,1);
	fModulus(1,2) += fMat(1,2);
	fModulus(1,3) += fMat(1,3);
	fModulus(1,4) += fMat(1,4);
	fModulus(1,5) += fMat(1,5);
//
	fModulus(2,0) += fMat(2,0);
	fModulus(2,1) += fMat(2,1);
	fModulus(2,2) += fMat(2,2);
	fModulus(2,3) += fMat(2,3);
	fModulus(2,4) += fMat(2,4);
	fModulus(2,5) += fMat(2,5);
//
	fModulus(3,0) += fMat(3,0);
	fModulus(3,1) += fMat(3,1);
	fModulus(3,2) += fMat(3,2);
	fModulus(3,3) += fMat(3,3);
	fModulus(3,4) += fMat(3,4);
	fModulus(3,5) += fMat(3,5);
//
	fModulus(4,0) += fMat(4,0);
	fModulus(4,1) += fMat(4,1);
	fModulus(4,2) += fMat(4,2);
	fModulus(4,3) += fMat(4,3);
	fModulus(4,4) += fMat(4,4);
	fModulus(4,5) += fMat(4,5);
//
	fModulus(5,0) += fMat(5,0);
	fModulus(5,1) += fMat(5,1);
	fModulus(5,2) += fMat(5,2);
	fModulus(5,3) += fMat(5,3);
	fModulus(5,4) += fMat(5,4);
	fModulus(5,5) += fMat(5,5);
	
	/*Add volumetric stress*/
	fStress[0] += p*fInverse[0];
	fStress[1] += p*fInverse[1];
	fStress[2] += p*fInverse[2];
	fStress[3] += p*fInverse[3];
	fStress[4] += p*fInverse[4];
	fStress[5] += p*fInverse[5];

	
	/*-1/3*I3^-2/3 ( (Cbar:C) otimes C^-1 + C^-1 otimes (C:Cbar))*/
	fSymMat1.A_ijkl_B_kl(fModbar,fC);
	fSymMat2.A_ijkl_B_ij(fModbar,fC);
	coeff = -third*I3rthird*I3rthird;
	
	fModulus(0,0) += coeff*(fInverse[0]*fSymMat2[0] + fInverse[0]*fSymMat1[0]);
	fModulus(0,1) += coeff*(fInverse[0]*fSymMat2[1] + fInverse[1]*fSymMat1[0]);
	fModulus(0,2) += coeff*(fInverse[0]*fSymMat2[2] + fInverse[2]*fSymMat1[0]);
	fModulus(0,3) += coeff*(fInverse[0]*fSymMat2[3] + fInverse[3]*fSymMat1[0]);
	fModulus(0,4) += coeff*(fInverse[0]*fSymMat2[4] + fInverse[4]*fSymMat1[0]);
	fModulus(0,5) += coeff*(fInverse[0]*fSymMat2[5] + fInverse[5]*fSymMat1[0]);
	//
	fModulus(1,0) += coeff*(fInverse[1]*fSymMat2[0] + fInverse[0]*fSymMat1[1]);
	fModulus(1,1) += coeff*(fInverse[1]*fSymMat2[1] + fInverse[1]*fSymMat1[1]);
	fModulus(1,2) += coeff*(fInverse[1]*fSymMat2[2] + fInverse[2]*fSymMat1[1]);
	fModulus(1,3) += coeff*(fInverse[1]*fSymMat2[3] + fInverse[3]*fSymMat1[1]);
	fModulus(1,4) += coeff*(fInverse[1]*fSymMat2[4] + fInverse[4]*fSymMat1[1]);
	fModulus(1,5) += coeff*(fInverse[1]*fSymMat2[5] + fInverse[5]*fSymMat1[1]);
	//
	fModulus(2,0) += coeff*(fInverse[2]*fSymMat2[0] + fInverse[0]*fSymMat1[2]);
	fModulus(2,1) += coeff*(fInverse[2]*fSymMat2[1] + fInverse[1]*fSymMat1[2]);
	fModulus(2,2) += coeff*(fInverse[2]*fSymMat2[2] + fInverse[2]*fSymMat1[2]);
	fModulus(2,3) += coeff*(fInverse[2]*fSymMat2[3] + fInverse[3]*fSymMat1[2]);
	fModulus(2,4) += coeff*(fInverse[2]*fSymMat2[4] + fInverse[4]*fSymMat1[2]);
	fModulus(2,5) += coeff*(fInverse[2]*fSymMat2[5] + fInverse[5]*fSymMat1[2]);
	//
	fModulus(3,0) += coeff*(fInverse[3]*fSymMat2[0] + fInverse[0]*fSymMat1[3]);
	fModulus(3,1) += coeff*(fInverse[3]*fSymMat2[1] + fInverse[1]*fSymMat1[3]);
	fModulus(3,2) += coeff*(fInverse[3]*fSymMat2[2] + fInverse[2]*fSymMat1[3]);
	fModulus(3,3) += coeff*(fInverse[3]*fSymMat2[3] + fInverse[3]*fSymMat1[3]);
	fModulus(3,4) += coeff*(fInverse[3]*fSymMat2[4] + fInverse[4]*fSymMat1[3]);
	fModulus(3,5) += coeff*(fInverse[3]*fSymMat2[5] + fInverse[5]*fSymMat1[3]);
	//
	fModulus(4,0) += coeff*(fInverse[4]*fSymMat2[0] + fInverse[0]*fSymMat1[4]);
	fModulus(4,1) += coeff*(fInverse[4]*fSymMat2[1] + fInverse[1]*fSymMat1[4]);
	fModulus(4,2) += coeff*(fInverse[4]*fSymMat2[2] + fInverse[2]*fSymMat1[4]);
	fModulus(4,3) += coeff*(fInverse[4]*fSymMat2[3] + fInverse[3]*fSymMat1[4]);
	fModulus(4,4) += coeff*(fInverse[4]*fSymMat2[4] + fInverse[4]*fSymMat1[4]);
	fModulus(4,5) += coeff*(fInverse[4]*fSymMat2[5] + fInverse[5]*fSymMat1[4]);
	//
	fModulus(5,0) += coeff*(fInverse[5]*fSymMat2[0] + fInverse[0]*fSymMat1[5]);
	fModulus(5,1) += coeff*(fInverse[5]*fSymMat2[1] + fInverse[1]*fSymMat1[5]);
	fModulus(5,2) += coeff*(fInverse[5]*fSymMat2[2] + fInverse[2]*fSymMat1[5]);
	fModulus(5,3) += coeff*(fInverse[5]*fSymMat2[3] + fInverse[3]*fSymMat1[5]);
	fModulus(5,4) += coeff*(fInverse[5]*fSymMat2[4] + fInverse[4]*fSymMat1[5]);
	fModulus(5,5) += coeff*(fInverse[5]*fSymMat2[5] + fInverse[5]*fSymMat1[5]);
	
	/*-2/3 I3^-1/3(Sbar otimes C^-1 + C^-1 otimes Sbar)*/
	coeff = -2.0*third*I3rthird;
	fModulus(0,0) += 2.0*coeff*(fSbar[0]*fInverse[0]);
	fModulus(0,1) += coeff*(fSbar[0]*fInverse[1] + fSbar[1]*fInverse[0]);
	fModulus(0,2) += coeff*(fSbar[0]*fInverse[2] + fSbar[2]*fInverse[0]);
	fModulus(0,3) += coeff*(fSbar[0]*fInverse[3] + fSbar[3]*fInverse[0]);
	fModulus(0,4) += coeff*(fSbar[0]*fInverse[4] + fSbar[4]*fInverse[0]);
	fModulus(0,5) += coeff*(fSbar[0]*fInverse[5] + fSbar[5]*fInverse[0]);
	//
	fModulus(1,0) += coeff*(fSbar[1]*fInverse[0] + fSbar[0]*fInverse[1]);
	fModulus(1,1) += 2.0*coeff*(fSbar[1]*fInverse[1]);
	fModulus(1,2) += coeff*(fSbar[1]*fInverse[2] + fSbar[2]*fInverse[1]);
	fModulus(1,3) += coeff*(fSbar[1]*fInverse[3] + fSbar[3]*fInverse[1]);
	fModulus(1,4) += coeff*(fSbar[1]*fInverse[4] + fSbar[4]*fInverse[1]);
	fModulus(1,5) += coeff*(fSbar[1]*fInverse[5] + fSbar[5]*fInverse[1]);
	//
	fModulus(2,0) += coeff*(fSbar[2]*fInverse[0] + fSbar[0]*fInverse[2]);
	fModulus(2,1) += coeff*(fSbar[2]*fInverse[1] + fSbar[1]*fInverse[2]);
	fModulus(2,2) += 2.0*coeff*(fSbar[2]*fInverse[2]);
	fModulus(2,3) += coeff*(fSbar[2]*fInverse[3] + fSbar[3]*fInverse[2]);
	fModulus(2,4) += coeff*(fSbar[2]*fInverse[4] + fSbar[4]*fInverse[2]);
	fModulus(2,5) += coeff*(fSbar[2]*fInverse[5] + fSbar[5]*fInverse[2]);
	//
	fModulus(3,0) += coeff*(fSbar[3]*fInverse[0] + fSbar[0]*fInverse[3]);
	fModulus(3,1) += coeff*(fSbar[3]*fInverse[1] + fSbar[1]*fInverse[3]);
	fModulus(3,2) += coeff*(fSbar[3]*fInverse[2] + fSbar[2]*fInverse[3]);
	fModulus(3,3) += 2.0*coeff*(fSbar[3]*fInverse[3]);
	fModulus(3,4) += coeff*(fSbar[3]*fInverse[4] + fSbar[4]*fInverse[3]);
	fModulus(3,5) += coeff*(fSbar[3]*fInverse[5] + fSbar[5]*fInverse[3]);
	//
	fModulus(4,0) += coeff*(fSbar[4]*fInverse[0] + fSbar[0]*fInverse[4]);
	fModulus(4,1) += coeff*(fSbar[4]*fInverse[1] + fSbar[1]*fInverse[4]);
	fModulus(4,2) += coeff*(fSbar[4]*fInverse[2] + fSbar[2]*fInverse[4]);
	fModulus(4,3) += coeff*(fSbar[4]*fInverse[3] + fSbar[3]*fInverse[4]);
	fModulus(4,4) += 2.0*coeff*(fSbar[4]*fInverse[4]);
	fModulus(4,5) += coeff*(fSbar[4]*fInverse[5] + fSbar[5]*fInverse[4]);
	//
	fModulus(5,0) += coeff*(fSbar[5]*fInverse[0] + fSbar[0]*fInverse[5]);
	fModulus(5,1) += coeff*(fSbar[5]*fInverse[1] + fSbar[1]*fInverse[5]);
	fModulus(5,2) += coeff*(fSbar[5]*fInverse[2] + fSbar[2]*fInverse[5]);
	fModulus(5,3) += coeff*(fSbar[5]*fInverse[3] + fSbar[3]*fInverse[5]);
	fModulus(5,4) += coeff*(fSbar[5]*fInverse[4] + fSbar[4]*fInverse[5]);
	fModulus(5,5) += 2.0*coeff*(fSbar[5]*fInverse[5]);

	/* ( 1/9 I3^-2/3 C:Modbar:C - 2/9 I3^-1/3 Sbar:C) C^-1 otimes C^-1;*/
	double ninth = third*third;
	coeff = ninth*I3rthird*I3rthird*(fC.B_ij_A_ijkl_B_kl(fModbar));
	coeff += 2.0*ninth*I3rthird*(fSbar.ScalarProduct(fC));
	
	/*add volumetric*/
	coeff += ComputeVolMatrixMod(fI3);
	
	fModulus(0,0) += coeff*(fInverse[0]*fInverse[0]);
	fModulus(0,1) += coeff*(fInverse[0]*fInverse[1]);
	fModulus(0,2) += coeff*(fInverse[0]*fInverse[2]);
	fModulus(0,3) += coeff*(fInverse[0]*fInverse[3]);
	fModulus(0,4) += coeff*(fInverse[0]*fInverse[4]);
	fModulus(0,5) += coeff*(fInverse[0]*fInverse[5]);
	//
	fModulus(1,0) += coeff*(fInverse[1]*fInverse[0]);
	fModulus(1,1) += coeff*(fInverse[1]*fInverse[1]);
	fModulus(1,2) += coeff*(fInverse[1]*fInverse[2]);
	fModulus(1,3) += coeff*(fInverse[1]*fInverse[3]);
	fModulus(1,4) += coeff*(fInverse[1]*fInverse[4]);
	fModulus(1,5) += coeff*(fInverse[1]*fInverse[5]);
	//
	fModulus(2,0) += coeff*(fInverse[2]*fInverse[0]);
	fModulus(2,1) += coeff*(fInverse[2]*fInverse[1]);
	fModulus(2,2) += coeff*(fInverse[2]*fInverse[2]);
	fModulus(2,3) += coeff*(fInverse[2]*fInverse[3]);
	fModulus(2,4) += coeff*(fInverse[2]*fInverse[4]);
	fModulus(2,5) += coeff*(fInverse[2]*fInverse[5]);
	//
	fModulus(3,0) += coeff*(fInverse[3]*fInverse[0]);
	fModulus(3,1) += coeff*(fInverse[3]*fInverse[1]);
	fModulus(3,2) += coeff*(fInverse[3]*fInverse[2]);
	fModulus(3,3) += coeff*(fInverse[3]*fInverse[3]);
	fModulus(3,4) += coeff*(fInverse[3]*fInverse[4]);
	fModulus(3,5) += coeff*(fInverse[3]*fInverse[5]);
	//
	fModulus(4,0) += coeff*(fInverse[4]*fInverse[0]);
	fModulus(4,1) += coeff*(fInverse[4]*fInverse[1]);
	fModulus(4,2) += coeff*(fInverse[4]*fInverse[2]);
	fModulus(4,3) += coeff*(fInverse[4]*fInverse[3]);
	fModulus(4,4) += coeff*(fInverse[4]*fInverse[4]);
	fModulus(4,5) += coeff*(fInverse[4]*fInverse[5]);
	//
	fModulus(5,0) += coeff*(fInverse[5]*fInverse[0]);
	fModulus(5,1) += coeff*(fInverse[5]*fInverse[1]);
	fModulus(5,2) += coeff*(fInverse[5]*fInverse[2]);
	fModulus(5,3) += coeff*(fInverse[5]*fInverse[3]);
	fModulus(5,4) += coeff*(fInverse[5]*fInverse[4]);
	fModulus(5,5) += coeff*(fInverse[5]*fInverse[5]);
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& FSFiberMatSplitT::S_IJ(void)
{
	
	/* stretch */
	Compute_C(fC);
	fI3 = fC.Det();
	double I3rthird = pow(fI3,-third);

	/*deviatoric part*/
	fCbar = fC;
	fCbar *= I3rthird;

	fInverse = fC;
	fInverse.Inverse();
 
	/*calculate deviatoric matrix contribution*/	
	/*calculate matrix contribution*/
	fSbar = 0.0; /*initialize*/
	ComputeDevMatrixStress(fCbar, fSbar);

	/*fiber contribution*/
	ComputeFiberStretch(fCbar, fFiberStretch);
	ComputeFiberStress(fFiberStretch, fFiberStress);
	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fSbar);
	
	/*Dot with projection tensor*/
	/*Sdev = I3^(-1/3)*(Sbar - 1/3 (Sbar:C) C^-1)*/
	fStress = fSbar;
	
	double coeff = fSbar.ScalarProduct(fC);
	coeff *= -third;
	fStress[0] += coeff*fInverse[0];
	fStress[1] += coeff*fInverse[1];
	fStress[2] += coeff*fInverse[2];
	fStress[3] += coeff*fInverse[3];
	fStress[4] += coeff*fInverse[4];
	fStress[5] += coeff*fInverse[5];
	
	fStress *= I3rthird;
	
	/*Add volumetric*/
	double p =  ComputeVolMatrixStress(fI3);
	fStress[0] += p*fInverse[0];
	fStress[1] += p*fInverse[1];
	fStress[2] += p*fInverse[2];
	fStress[3] += p*fInverse[3];
	fStress[4] += p*fInverse[4];
	fStress[5] += p*fInverse[5];
	
	return(fStress);
}

/* accept parameter list */
void FSFiberMatSplitT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);

	/* dimension work space */
	fCbar.Dimension(fNumSD);
	fInverse.Dimension(fNumSD);
	
	fModbar.Dimension(dSymMatrixT::NumValues(fNumSD));
	
	fSbar.Dimension(fNumSD);
	
	fSymMat1.Dimension(fNumSD);
	fSymMat2.Dimension(fNumSD);
	fMat.Dimension(dSymMatrixT::NumValues(fNumSD));

}


/* stress */
void FSFiberMatSplitT::ComputeMatrixStress(const dSymMatrixT& Stretch, dSymMatrixT& Stress)
{
	
	/* stretch */
	double I3 = Stretch.Det();
	double I3rthird = pow(I3,-third);

	/*deviatoric part*/
	fCbar = Stretch;
	fCbar *= I3rthird;

	fInverse = Stretch;
	fInverse.Inverse();
 
	/*calculate deviatoric matrix contribution*/	
	/*calculate matrix contribution*/
	fSbar = 0.0; /*initialize*/
	ComputeDevMatrixStress(fCbar, fSbar);
		
	/*Dot with projection tensor*/
	/*Sdev = I3^(-1/3)*(Sbar - 1/3 (Sbar:C) C^-1)*/
	Stress = fSbar;
	
	double coeff = fSbar.ScalarProduct(Stretch);
	coeff *= -third;
	Stress[0] += coeff*fInverse[0];
	Stress[1] += coeff*fInverse[1];
	Stress[2] += coeff*fInverse[2];
	Stress[3] += coeff*fInverse[3];
	Stress[4] += coeff*fInverse[4];
	Stress[5] += coeff*fInverse[5];
	
	Stress *= I3rthird;
	
	/*Add volumetric*/
	double p =  ComputeVolMatrixStress(I3);
	Stress[0] += p*fInverse[0];
	Stress[1] += p*fInverse[1];
	Stress[2] += p*fInverse[2];
	Stress[3] += p*fInverse[3];
	Stress[4] += p*fInverse[4];
	Stress[5] += p*fInverse[5];
	
}


/* modulus */
void FSFiberMatSplitT::ComputeMatrixMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod)
{
	/* stretch */
	double I3;
	I3 = Stretch.Det();
	double I3rthird = pow(I3,-third);

	/*deviatoric part*/
	fCbar = Stretch;
	fCbar *= I3rthird;
	
	fInverse = Stretch;
	fInverse.Inverse();

	/*calculate matrix contribution*/
	fSbar = 0.0;
	fModbar = 0.0;  /*initialize*/
	ComputeDevMatrixMod(fCbar,  fSbar, fModbar);
	
	/*fiber contribution*/
	ComputeFiberStretch(fCbar, fFiberStretch);
	ComputeFiberMod(fFiberStretch, fFiberStress, fFiberMod);

	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fSbar);
	/* rotate modulus to lab coordinates */
	AssembleFiberModuli(fFiberMod, fModbar);
	
	/*Deviatoric stress*/
	/*Sdev = I3^(-1/3)*(Sbar - 1/3 (Sbar:C) C^-1)*/
	Stress = fSbar;
	
	double coeff = fSbar.ScalarProduct(Stretch);
	coeff *= third;
	Stress[0] -= coeff*fInverse[0];
	Stress[1] -= coeff*fInverse[1];
	Stress[2] -= coeff*fInverse[2];
	Stress[3] -= coeff*fInverse[3];
	Stress[4] -= coeff*fInverse[4];
	Stress[5] -= coeff*fInverse[5];

	fStress *= I3rthird;
		
	/*modulus*/
	
	/*I3^-2/3 Cbar_dev*/
	Mod = fModbar; 
	Mod *= I3rthird*I3rthird;
	
	/*-1/3*I3^-2/3 ( (Cbar:C) otimes C^-1 + C^-1 otimes (C:Cbar))*/
	fSymMat1.A_ijkl_B_kl(fModbar,Stretch);
	fSymMat2.A_ijkl_B_ij(fModbar,Stretch);
	coeff = -third*I3rthird*I3rthird;
	
	Mod(0,0) += coeff*(fInverse[0]*fSymMat2[0] + fInverse[0]*fSymMat1[0]);
	Mod(0,1) += coeff*(fInverse[0]*fSymMat2[1] + fInverse[1]*fSymMat1[0]);
	Mod(0,2) += coeff*(fInverse[0]*fSymMat2[2] + fInverse[2]*fSymMat1[0]);
	Mod(0,3) += coeff*(fInverse[0]*fSymMat2[3] + fInverse[3]*fSymMat1[0]);
	Mod(0,4) += coeff*(fInverse[0]*fSymMat2[4] + fInverse[4]*fSymMat1[0]);
	Mod(0,5) += coeff*(fInverse[0]*fSymMat2[5] + fInverse[5]*fSymMat1[0]);
	//
	Mod(1,0) += coeff*(fInverse[1]*fSymMat2[0] + fInverse[0]*fSymMat1[1]);
	Mod(1,1) += coeff*(fInverse[1]*fSymMat2[1] + fInverse[1]*fSymMat1[1]);
	Mod(1,2) += coeff*(fInverse[1]*fSymMat2[2] + fInverse[2]*fSymMat1[1]);
	Mod(1,3) += coeff*(fInverse[1]*fSymMat2[3] + fInverse[3]*fSymMat1[1]);
	Mod(1,4) += coeff*(fInverse[1]*fSymMat2[4] + fInverse[4]*fSymMat1[1]);
	Mod(1,5) += coeff*(fInverse[1]*fSymMat2[5] + fInverse[5]*fSymMat1[1]);
	//
	Mod(2,0) += coeff*(fInverse[2]*fSymMat2[0] + fInverse[0]*fSymMat1[2]);
	Mod(2,1) += coeff*(fInverse[2]*fSymMat2[1] + fInverse[1]*fSymMat1[2]);
	Mod(2,2) += coeff*(fInverse[2]*fSymMat2[2] + fInverse[2]*fSymMat1[2]);
	Mod(2,3) += coeff*(fInverse[2]*fSymMat2[3] + fInverse[3]*fSymMat1[2]);
	Mod(2,4) += coeff*(fInverse[2]*fSymMat2[4] + fInverse[4]*fSymMat1[2]);
	Mod(2,5) += coeff*(fInverse[2]*fSymMat2[5] + fInverse[5]*fSymMat1[2]);
	//
	Mod(3,0) += coeff*(fInverse[3]*fSymMat2[0] + fInverse[0]*fSymMat1[3]);
	Mod(3,1) += coeff*(fInverse[3]*fSymMat2[1] + fInverse[1]*fSymMat1[3]);
	Mod(3,2) += coeff*(fInverse[3]*fSymMat2[2] + fInverse[2]*fSymMat1[3]);
	Mod(3,3) += coeff*(fInverse[3]*fSymMat2[3] + fInverse[3]*fSymMat1[3]);
	Mod(3,4) += coeff*(fInverse[3]*fSymMat2[4] + fInverse[4]*fSymMat1[3]);
	Mod(3,5) += coeff*(fInverse[3]*fSymMat2[5] + fInverse[5]*fSymMat1[3]);
	//
	Mod(4,0) += coeff*(fInverse[4]*fSymMat2[0] + fInverse[0]*fSymMat1[4]);
	Mod(4,1) += coeff*(fInverse[4]*fSymMat2[1] + fInverse[1]*fSymMat1[4]);
	Mod(4,2) += coeff*(fInverse[4]*fSymMat2[2] + fInverse[2]*fSymMat1[4]);
	Mod(4,3) += coeff*(fInverse[4]*fSymMat2[3] + fInverse[3]*fSymMat1[4]);
	Mod(4,4) += coeff*(fSymMat1[4]*fSymMat2[4] + fInverse[4]*fSymMat1[4]);
	Mod(4,5) += coeff*(fInverse[4]*fSymMat2[5] + fInverse[5]*fSymMat1[4]);
	//
	Mod(5,0) += coeff*(fInverse[5]*fSymMat2[0] + fInverse[0]*fSymMat1[5]);
	Mod(5,1) += coeff*(fInverse[5]*fSymMat2[1] + fInverse[1]*fSymMat1[5]);
	Mod(5,2) += coeff*(fInverse[5]*fSymMat2[2] + fInverse[2]*fSymMat1[5]);
	Mod(5,3) += coeff*(fInverse[5]*fSymMat2[3] + fInverse[3]*fSymMat1[5]);
	Mod(5,4) += coeff*(fInverse[5]*fSymMat2[4] + fInverse[4]*fSymMat1[5]);
	Mod(5,5) += coeff*(fInverse[5]*fSymMat2[5] + fInverse[5]*fSymMat1[5]);

	/*2/3 (Sdev otimes C^-1 + C^-1 otimes Sdev)*/
	coeff = -2.0*third;
	Mod(0,0) += 2.0*coeff*(Stress[0]*fInverse[0]);
	Mod(0,1) += coeff*(Stress[0]*fInverse[1] + Stress[1]*fInverse[0]);
	Mod(0,2) += coeff*(Stress[0]*fInverse[2] + Stress[2]*fInverse[0]);
	Mod(0,3) += coeff*(Stress[0]*fInverse[3] + Stress[3]*fInverse[0]);
	Mod(0,4) += coeff*(Stress[0]*fInverse[4] + Stress[4]*fInverse[0]);
	Mod(0,5) += coeff*(Stress[0]*fInverse[5] + Stress[5]*fInverse[0]);
	//
	Mod(1,0) += coeff*(Stress[1]*fInverse[0] + Stress[0]*fInverse[1]);
	Mod(1,1) += 2.0*coeff*(Stress[1]*fInverse[1]);
	Mod(1,2) += coeff*(Stress[1]*fInverse[2] + Stress[2]*fInverse[1]);
	Mod(1,3) += coeff*(Stress[1]*fInverse[3] + Stress[3]*fInverse[1]);
	Mod(1,4) += coeff*(Stress[1]*fInverse[4] + Stress[4]*fInverse[1]);
	Mod(1,5) += coeff*(Stress[1]*fInverse[5] + Stress[5]*fInverse[1]);
	//
	Mod(2,0) += coeff*(Stress[2]*fInverse[0] + Stress[0]*fInverse[2]);
	Mod(2,1) += coeff*(Stress[2]*fInverse[1] + Stress[1]*fInverse[2]);
	Mod(2,2) += 2.0*coeff*(Stress[2]*fInverse[2]);
	Mod(2,3) += coeff*(Stress[2]*fInverse[3] + Stress[3]*fInverse[2]);
	Mod(2,4) += coeff*(Stress[2]*fInverse[4] + Stress[4]*fInverse[2]);
	Mod(2,5) += coeff*(Stress[2]*fInverse[5] + Stress[5]*fInverse[2]);
	//
	Mod(3,0) += coeff*(Stress[3]*fInverse[0] + Stress[0]*fInverse[3]);
	Mod(3,1) += coeff*(Stress[3]*fInverse[1] + Stress[1]*fInverse[3]);
	Mod(3,2) += coeff*(Stress[3]*fInverse[2] + Stress[2]*fInverse[3]);
	Mod(3,3) += 2.0*coeff*(Stress[3]*fInverse[3]);
	Mod(3,4) += coeff*(Stress[3]*fInverse[4] + Stress[4]*fInverse[3]);
	Mod(3,5) += coeff*(Stress[3]*fInverse[5] + Stress[5]*fInverse[3]);
	//
	Mod(4,0) += coeff*(Stress[4]*fInverse[0] + Stress[0]*fInverse[4]);
	Mod(4,1) += coeff*(Stress[4]*fInverse[1] + Stress[1]*fInverse[4]);
	Mod(4,2) += coeff*(Stress[4]*fInverse[2] + Stress[2]*fInverse[4]);
	Mod(4,3) += coeff*(Stress[4]*fInverse[3] + Stress[3]*fInverse[4]);
	Mod(4,4) += 2.0*coeff*(Stress[4]*fInverse[4]);
	Mod(4,5) += coeff*(Stress[4]*fInverse[5] + Stress[5]*fInverse[4]);
	//
	Mod(5,0) += coeff*(Stress[5]*fInverse[0] + Stress[0]*fInverse[5]);
	Mod(5,1) += coeff*(Stress[5]*fInverse[1] + Stress[1]*fInverse[5]);
	Mod(5,2) += coeff*(Stress[5]*fInverse[2] + Stress[2]*fInverse[5]);
	Mod(5,3) += coeff*(Stress[5]*fInverse[3] + Stress[3]*fInverse[5]);
	Mod(5,4) += coeff*(Stress[5]*fInverse[4] + Stress[4]*fInverse[5]);
	Mod(5,5) += 2.0*coeff*(Stress[5]*fInverse[5]);

	/* ( 1/9 I3^-2/3 C:Modbar:C - 2/9 I3^-1/3 Sbar:C) C^-1 otimes C^-1;*/
	double ninth = third*third;
	coeff = ninth*I3rthird*I3rthird*(Stretch.B_ij_A_ijkl_B_kl(fModbar));
	coeff -= 2.0*ninth*I3rthird*(fSbar.ScalarProduct(Stretch));
	
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

	/*2/3 I3^-1/3(Sbar:C) (-pdf(C^-1)(C))*/
	coeff = fSbar.ScalarProduct(Stretch);
	coeff *= 2.0*third*I3rthird;

	fMat.ReducedI_C(fInverse);
	fMat *= coeff;
	Mod(0,0) += fMat(0,0);
	Mod(0,1) += fMat(0,1);
	Mod(0,2) += fMat(0,2);
	Mod(0,3) += fMat(0,3);
	Mod(0,4) += fMat(0,4);
	Mod(0,5) += fMat(0,5);
	//
	Mod(1,0) += fMat(1,0);
	Mod(1,1) += fMat(1,1);
	Mod(1,2) += fMat(1,2);
	Mod(1,3) += fMat(1,3);
	Mod(1,4) += fMat(1,4);
	Mod(1,5) += fMat(1,5);
//
	Mod(2,0) += fMat(2,0);
	Mod(2,1) += fMat(2,1);
	Mod(2,2) += fMat(2,2);
	Mod(2,3) += fMat(2,3);
	Mod(2,4) += fMat(2,4);
	Mod(2,5) += fMat(2,5);
//
	Mod(3,0) += fMat(3,0);
	Mod(3,1) += fMat(3,1);
	Mod(3,2) += fMat(3,2);
	Mod(3,3) += fMat(3,3);
	Mod(3,4) += fMat(3,4);
	Mod(3,5) += fMat(3,5);
//
	Mod(4,0) += fMat(4,0);
	Mod(4,1) += fMat(4,1);
	Mod(4,2) += fMat(4,2);
	Mod(4,3) += fMat(4,3);
	Mod(4,4) += fMat(4,4);
	Mod(4,5) += fMat(4,5);
//
	Mod(5,0) += fMat(5,0);
	Mod(5,1) += fMat(5,1);
	Mod(5,2) += fMat(5,2);
	Mod(5,3) += fMat(5,3);
	Mod(5,4) += fMat(5,4);
	Mod(5,5) += fMat(5,5);
}
