/* $Id: FDSimoVisco3D.cpp,v 1.6 2002-11-14 17:06:12 paklein Exp $ */
/* created:   TDN (5/31/2001) */
#include "FDSimoVisco3D.h"
#include "FDMatSupportT.h"

#include <math.h>
#include <iostream.h>
#include "fstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

FDSimoVisco3D::FDSimoVisco3D(ifstreamT& in, const FDMatSupportT& support):
	FDSimoViscoBaseT(in, support),
	fStress(3),
	fModulus(6),
	fFbar_E(3),
	fFbar_I(3),
	fF(3),
	fthird(1.0/3.0)
{
	in >> ftauS;
	in >> ftauB;
}	

void FDSimoVisco3D::Print(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::Print(out);
	out<<"Relaxation Times \n";
	out<<"     Shear relaxation time: "<<ftauS<<'\n';
	out<<"     Bulk relaxation time: "<<ftauB<<'\n';
}

void FDSimoVisco3D::PrintName(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::PrintName(out);
}

const dMatrixT& FDSimoVisco3D::c_ijkl(void)
{
        /*equilibrium component*/

        /*deviatoric part*/
        Devmod(fFbar_E, fJ_E,fModulus, kEquilibrium);

        /*volumetric part*/
	dMatrixT bulkmod(6);
	dMatrixT IxI(6);

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
	dMatrixT devmod(6);
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
       
const dSymMatrixT& FDSimoVisco3D::s_ij(void)
{
	double dt = fFDMatSupport.TimeStep();
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

	/*equilibrium components*/

	/*deviatoric part*/
	Sig_dev(fFbar_E, fJ_E, fStress, kEquilibrium);
	
	/*volumetric part*/
	fStress.PlusIdentity(dUdJ(fJ_E, kEquilibrium));

	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	/*overstress*/
	if(fFDMatSupport.RunState() == GlobalT::kFormRHS)
	{
	        fJ_I = fJ;
		fFbar_I = fF;
		fFbar_I *= Js;
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

const dMatrixT& FDSimoVisco3D::C_IJKL(void)
{
        /* deformation gradient */
        const dMatrixT& Fmat = F();
        
        /* transform */
        fModulus.SetToScaled(fJ, PullBack(Fmat, c_ijkl()));
        return fModulus;        
}

const dSymMatrixT& FDSimoVisco3D::S_IJ(void)
{
        /* deformation gradient */
        const dMatrixT& Fmat = F();
        
        /* transform */
        fStress.SetToScaled(fJ, PullBack(Fmat, s_ij()));
        return fStress;
}
