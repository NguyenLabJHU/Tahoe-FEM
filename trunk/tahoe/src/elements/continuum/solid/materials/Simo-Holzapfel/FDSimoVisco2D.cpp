/* $Id: FDSimoVisco2D.cpp,v 1.7 2003-01-29 07:34:49 paklein Exp $ */
/* created:   TDN (5/31/2001) */
#include "FDSimoVisco2D.h"

#include <math.h>
#include <iostream.h>
#include "fstreamT.h"
#include "ExceptionT.h"
#include "FSMatSupportT.h"

using namespace Tahoe;

FDSimoVisco2D::FDSimoVisco2D(ifstreamT& in, const FSMatSupportT& support):
        Material2DT(in),
	FDSimoViscoBaseT(in, support),
	fStress(2),
	fModulus(3),
	fFbar_E(2),
	fFbar_I(2),
	fF(2),
	fthird(1.0/3.0)
{
	in >> ftauS;
	in >> ftauB;
}	

void FDSimoVisco2D::Print(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::Print(out);
	out<<"Relaxation Times \n";
	out<<"     Shear relaxation time: "<<ftauS<<'\n';
	out<<"     Bulk relaxation time: "<<ftauB<<'\n';
}

void FDSimoVisco2D::PrintName(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::PrintName(out);
}

const dMatrixT& FDSimoVisco2D::c_ijkl(void)
{
        /*equilibrium component*/

        /*deviatoric part*/
        Devmod(fFbar_E, fJ_E,fModulus, kEquilibrium);

        /*volumetric part*/
	dMatrixT bulkmod(3);
	dMatrixT IxI(3);

	IxI.ReducedIndexII();
	bulkmod.ReducedIndexI();
	bulkmod *= 2.0;
	bulkmod -= IxI;
	bulkmod *= -dUdJ(fJ_E, kEquilibrium);

	IxI *= fJ_E*ddUddJ(fJ_E, kEquilibrium);
	bulkmod += IxI;

	fModulus += bulkmod;
	/*non-equilibrium component*/

        /*deviatoric part*/
	dMatrixT devmod(3);
        Devmod(fFbar_I,fJ_I,devmod, kNonEquilibrium);
	devmod *= falphaS;

        /*volumetric part*/

	IxI.ReducedIndexII();
	bulkmod.ReducedIndexI();
	bulkmod *= 2.0;
	bulkmod -= IxI;
	bulkmod *= -dUdJ(fJ_I,kNonEquilibrium)*falphaB;

	IxI *= fJ_I*ddUddJ(fJ_I,kNonEquilibrium)*falphaB;
	bulkmod += IxI;

	fModulus += devmod;
	fModulus += bulkmod;

	return(fModulus);
}
       
const dSymMatrixT& FDSimoVisco2D::s_ij(void)
{
	double dt = fFSMatSupport.TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);
	fbetaS = exp(-taudtS);
	fbetaB = exp(-taudtB);

	fF = F_mechanical();
	fJ = fF.Det();
	double Js = pow(fJ, -fthird);
	fFbar_E = fF;
	fFbar_E *= Js;
	fJ_E = fJ;
	OutOfPlaneStretch(fFbar_E, fJ_E, kEquilibrium);

	/*equilibrium components*/

	/*deviatoric part*/
	Sig_dev(fFbar_E, fJ_E, fStress, kEquilibrium);
	
	/*volumetric part*/
	fStress.PlusIdentity(dUdJ(fJ_E, kEquilibrium));

	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	/*overstress*/
	if(fFSMatSupport.RunState() == GlobalT::kFormRHS)
	{
	        fJ_I = fJ;
		fFbar_I = fF;
		fFbar_I *= Js;
		OutOfPlaneStretch(fFbar_I, fJ_I, kNonEquilibrium);

	        /*deviatoric part*/
		Sig_dev(fFbar_I, fJ_I, fDevInStress, kNonEquilibrium);
	        fDevOverStress.SetToCombination(fbetaS, fDevOverStress_n, 
						falphaS, fDevInStress,
						-falphaS, fDevInStress_n);
	
		/*volumetric part*/
		fVolInStress[0] = dUdJ(fJ_I, kNonEquilibrium);
		fVolOverStress[0] = fbetaB*fVolOverStress_n[0] + falphaB* 
		                   (fVolInStress[0]-fVolInStress_n[0]);
		Store(element,CurrIP());
	}
	fStress += fDevOverStress;
	fStress.PlusIdentity(fVolOverStress[0]);

	return(fStress);
}

const dMatrixT& FDSimoVisco2D::C_IJKL(void)
{
        /* deformation gradient */
        const dMatrixT& Fmat = F();
        
        /* transform */
        fModulus.SetToScaled(fJ, PullBack(Fmat, c_ijkl()));
        return fModulus;        
}

const dSymMatrixT& FDSimoVisco2D::S_IJ(void)
{
        /* deformation gradient */
        const dMatrixT& Fmat = F();
        
        /* transform */
        fStress.SetToScaled(fJ, PullBack(Fmat, s_ij()));
        return fStress;
}
