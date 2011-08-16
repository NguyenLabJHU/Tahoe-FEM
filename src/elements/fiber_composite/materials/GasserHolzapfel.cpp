/* $Id: GasserHolzapfel.cpp,v 1.2 2011-08-16 13:00:32 thao Exp $ */

#include "GasserHolzapfel.h"
#include <math.h>
#include "ParameterContainerT.h"

const double third = 1.0/3.0;
const int kNumOutputVar = 4;
static const char* Labels[kNumOutputVar] = 
	{"Fib_X", "Fib_Y", "Fib_Z","J"};
static const int perm[3][3] = {0,1,2,1,2,0,2,0,1};

using namespace Tahoe;

/* constructors */
GasserHolzapfel::GasserHolzapfel(void):
    FSFiberMatSplitT(),
	ParameterInterfaceT("gasser_holzapfel")			
{
	
}													

/* destructor */
GasserHolzapfel::~GasserHolzapfel(void) 
{ 
}

/* strain energy density */
double GasserHolzapfel::StrainEnergyDensity(void)
{

	/*matrix contribution*/
	/* stretch */
	Compute_C(fC);
	double I3 = fC.Det();
	
	double I3rthird = pow(I3,-third);
	dSymMatrixT fCbar = fC;								
	fCbar *= I3rthird;
	ComputeFiberStretch(fCbar, fFiberStretch);
		
	double I1bar = (fCbar[0]+fCbar[1]+fCbar[2]);
	double I4bar = fFiberStretch[0];
		
	/*coupled compressible Neo-Hookean*/
	/* mu/2 (I1 -3) + 0.25*fKappa*(I3-1-log(I3)))*/
	double energyMatrix = 0.5*fMu*(I1bar-3.0) + 0.25*fBulkMod*((I3-1.0) + (exp(fk2*(I3-1.0))-1.0) -log(I3) );
	
	/*fiber contribution*/
	//TO DO::  Add fiber contribution
	
	double beta = fKappa*I1bar+(1-3*fKappa)*I4bar-1;
	
	double energyFiber = 0.5* fk1/fk2*exp(fk2*pow(beta,2)-1);
	
	double energy = energyMatrix + energyFiber;
	
	return energy;
}

int GasserHolzapfel::NumOutputVariables() const {
	return kNumOutputVar;
}

void GasserHolzapfel::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void GasserHolzapfel::ComputeOutput(dArrayT& output)
{
	/*calculates deformed fiber vectors*/
//	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	const dMatrixT& Q = GetRotation();
	const double* p_fib = Q(0);
	
	const dMatrixT& F = F_mechanical();
	double* pb = output.Pointer();
	
	/*deformed NT fiber orientation*/
	F.Multx(p_fib, pb);
	pb += NumSD();
	
	/*non-equilibrium strain energy density */
//	*pb = NonequilibriumStrainEnergyDensity();

	*pb = F.Det();
}

/* describe the parameters needed by the interface */
void GasserHolzapfel::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSFiberMatSplitT::DefineParameters(list);
	
	/* matrix properties */
	ParameterT mu(ParameterT::Double, "matrix_shear_mod");
	mu.AddLimit(0, LimitT::Lower);
	list.AddParameter(mu);

	ParameterT bulk(ParameterT::Double, "matrix_bulk_mod");
	bulk.AddLimit(0, LimitT::Lower);
	list.AddParameter(bulk);

	/* fiber properties */
	ParameterT kappa(ParameterT::Double, "fiber_dispersion");
	kappa.AddLimit(0, LimitT::LowerInclusive);
	kappa.AddLimit(0.33, LimitT::UpperInclusive);
	list.AddParameter(kappa);
	
	ParameterT k1(ParameterT::Double, "fiber_modulus_k1");
	k1.AddLimit(0, LimitT::Lower);
	list.AddParameter(k1);

	ParameterT k2(ParameterT::Double, "fiber_stiffening_k2");
	k2.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(k2);
	
}

/* accept parameter list */
void GasserHolzapfel::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatSplitT::TakeParameterList(list);

	/*matrix params*/
	fMu = list.GetParameter("matrix_shear_mod");
	fBulkMod = list.GetParameter("matrix_bulk_mod");
	fKappa = list.GetParameter("fiber_dispersion");
	fk1 = list.GetParameter("fiber_modulus_k1");
	fk2 = list.GetParameter("fiber_stiffening_k2");
	fb2 = fk2;
	
	/* allocate memory */
	/*3D fiber stress and modulus*/
	fNumFibStress = dSymMatrixT::NumValues(fNumSD);
	fNumFibModuli = fNumFibStress*fNumFibStress;
	fFiberStretch.Dimension(fNumSD);
	fFiberStress.Dimension(fNumSD);
	fFiberMod.Dimension(fNumFibStress);
	
	/*allocate workspace*/
	fMat.Dimension(fNumSD);

}

/***********************************************************************
 * Protected
 ***********************************************************************/

void GasserHolzapfel::ComputeDevMatrixStress(const dSymMatrixT& Cbar, dSymMatrixT& Sbar)
{
	/*TO DO:  Verify that this is right*/
	/*Sbar = 2pdf{W_matrix}{Cbar_IJ} defined in eq. 6.91 Holzapfel textbook*/
	/*     = mu  for NeoHookean */  
	
	Sbar[0] += fMu;
	Sbar[1] += fMu;
	Sbar[2] += fMu;
	
	/*Sbar[0] is equivalent to Sbar(0,0)*/
	/*Sbar[3] is equivalent to Sbar(2,3) equivalent to Sbar(3,2)*/

}

void GasserHolzapfel::ComputeDevMatrixMod(const dSymMatrixT& Cbar, dSymMatrixT& Stress, dMatrixT& Mod)
{
	/*TO DO:  Verify that this is right*/
	/*2pdf{Sbar_matrix_IJ}{Cbar_KL} = 0 like in 6.169 but without the J^-4/3 term*/
	
	Stress[0] += fMu;
	Stress[1] += fMu;
	Stress[2] += fMu;

	Mod += 0.0;

}

double GasserHolzapfel::ComputeVolMatrixStress(const double I3)
{
	/*p = 2pdf(Wvol)(I3)I3 = kappa/2(I3-1) C^-1*/
	double x = I3-1.0;
	return( 0.5*fBulkMod*x*(1.0 + 2.0* fb2* exp(fb2*x*x) *I3) )  ;
}

double GasserHolzapfel::ComputeVolMatrixMod(const double I3)
{
	/*2pdf(p)(I3)*I3*/
	double x = I3-1.0;
	double beta = exp(fb2*x*x);
	return( 2.0*fBulkMod* I3* (0.5 + 2.0*fb2*fb2*beta*x*x*I3 + fb2*beta* (-1.0 + 2.0*I3)) );
}
	
/*computes integrated fiber stress in local frame*/
void GasserHolzapfel::ComputeFiberStress (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress)  // FiberStretch = deviatoric component of stretch tensor in fiber coordinants
{
	
	//Define Invariants
	
	Compute_C(fC);
	double I3 = fC.Det();
	double I3rthird = pow(I3,-third);
	fMat = fC;
	fMat *= I3rthird;
	double I1bar = (fMat[0]+fMat[1]+fMat[2]);
		
	double I4bar=FiberStretch[0];
	
	/*initialize pointers*/
	/*Fiberstretch oriented in local fiber coordinates*/
	/* Sbar_fiber_IJ = 2pdf(W_fiber)(Cbar_IJ) values in local frame formed by fiber (1) and perpendicular (2) orientations*/	
	FiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_33*/

	double beta = fKappa*I1bar+(1.0-3.0*fKappa)*I4bar-1.0;
	s1 = 2.0*fk1*beta*(fKappa+(1.0-3.0*fKappa))*exp(fk2*pow(beta,2.0));       /* Fiber direction is in 11  */
	s2 = 2.0*fk1*beta*fKappa*exp(fk2*pow(beta,2.0));
	s3 = 2.0*fk1*beta*fKappa*exp(fk2*pow(beta,2.0));
/* all shear components, sf_23, sf_13, sf_12, are zero*/
	
}
	
/*computes integrated moduli in local frame*/
void GasserHolzapfel::ComputeFiberMod (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress, dMatrixT& FiberMod)
{
	
	Compute_C(fC);
	double I3 = fC.Det();
	double I3rthird = pow(I3,-third);
	fMat = fC;									//fMat is a temp value for fCbar
	fMat *= I3rthird;
	double I1bar = (fMat[0]+fMat[1]+fMat[2]);
	
	double I4bar=FiberStretch[0];
	
	/*initialize pointers*/
	/*Fiberstretch oriented in local fiber coordinates*/
	/* Cdev_bar_fiber_IJKL = 2pdf(Sbar_fiber_IJ)(Cbar_KL) Modulus values in local coordinate frame formed by a1 (fiber) and a2 (perp) */	
	FiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_33*/
 
	double beta = fKappa*I1bar+(1.0-3.0*fKappa)*I4bar-1.0;
	s1 = 2.0*fk1*beta*(fKappa+(1.0-3.0*fKappa))*exp(fk2*pow(beta,2.0));       /* Fiber direction is in 11  */
	s2 = 2.0*fk1*beta*fKappa*exp(fk2*pow(beta,2.0));
	s3 = 2.0*fk1*beta*fKappa*exp(fk2*pow(beta,2.0));
	/* all shear components, sf_23, sf_13, sf_12, are zero*/
	
	FiberMod = 0.0;
	// re-define invarients inside this function?
	// Define FiberMod (automatically updated as c11 etc?)
	FiberMod(0,0) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*(fKappa+(1.0-3.0*fKappa))*(fKappa+(1.0-3.0*fKappa));
	FiberMod(0,1) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*(fKappa+(1.0-3.0*fKappa))*fKappa;
	FiberMod(0,2) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*(fKappa+(1.0-3.0*fKappa))*fKappa;
		
	
	FiberMod(1,0) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*(fKappa+(1.0-3.0*fKappa))*fKappa;
	FiberMod(1,1) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*fKappa*fKappa;
	FiberMod(1,2) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*fKappa*fKappa;
	
	FiberMod(2,0) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*fKappa*(fKappa+(1-3*fKappa));
	FiberMod(2,1) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*fKappa*fKappa;
	FiberMod(2,2) = 4.0*fk1*(1.0+2.0*fk2*pow(beta,2.0))*exp(fk2*pow(beta,2.0))*fKappa*fKappa;
	
}


