/* $Id: FDSV_KStV3D.cpp,v 1.7 2003-01-29 07:34:50 paklein Exp $ */
/* created:   TDN (5/31/2001) */

#include "FDSV_KStV3D.h"

#include <math.h>
#include <iostream.h>
#include "fstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

FDSV_KStV3D::FDSV_KStV3D(ifstreamT& in, const FSMatSupportT& support):
	FDSimoViscoBaseT(in, support),
	fStress(3),
	fModulus(6),
	fE(3),
	fMu(2),
	fKappa(2),
	fthird(1.0/3.0)
{
	in >> ftauS;
	in >> ftauB;

        double& mu_EQ = fMu[kEquilibrium];
	double& mu_NEQ = fMu[kNonEquilibrium]; 
	double& kappa_EQ = fKappa[kEquilibrium]; 
	double& kappa_NEQ = fKappa[kNonEquilibrium];

	in >> mu_EQ;
	in >> kappa_EQ;

	in >> mu_NEQ;
	in >> kappa_NEQ;
}	

void FDSV_KStV3D::Print(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::Print(out);
	out << "Equilibrium Potential";
	out << "     Shear Modulus: "<<fMu[0]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[0]<<'\n';
	out << "Non-Equilibrium Potential";
	out << "     Shear Modulus: "<<fMu[1]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[1]<<'\n';
	out<<"Relaxation times\n";
	out<<"     Shear relaxation time: "<<ftauS<<'\n';
	out<<"     Bulk relaxation time: "<<ftauB<<'\n';
}

void FDSV_KStV3D::PrintName(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::PrintName(out);
	out << "Equilibrium/Non-Equilibrium Potential:\n";
	out << "Kirchoff St. Venant\n";
	out << "Kirchoff St. Venant\n";
}

const dMatrixT& FDSV_KStV3D::c_ijkl(void) 
{ 
        const dMatrixT& F_mech = F_mechanical(); 
        fModulus = PushForward(F_mech,C_IJKL());
        fModulus /= F_mech.Det();
        return fModulus;
}

const dSymMatrixT& FDSV_KStV3D::s_ij(void)
{
        const dMatrixT& F_mech = F_mechanical(); 
        fStress = PushForward(F_mech,S_IJ());
        fStress /= F_mech.Det();
        return fStress;
}
       
const dMatrixT& FDSV_KStV3D::C_IJKL(void)
{        
        /*equilibrium component*/
	double& mu = fMu[kEquilibrium];
	double& kappa = fKappa[kEquilibrium];

        /*deviatoric part*/
	dMatrixT IxI(6);
	IxI.ReducedIndexII();
	IxI *= fthird;
	fModulus.ReducedIndexI();
	fModulus -= IxI;
	fModulus *= 2.0*mu;

        /*volumetric part*/
	IxI *= 3.0*kappa;
	fModulus += IxI;

	/*non-equilibrium component*/
	mu = fMu[kNonEquilibrium];
	kappa = fKappa[kNonEquilibrium];

        /*deviatoric part*/
	dMatrixT mod(6);
	IxI.ReducedIndexII();
	IxI *= fthird;
	mod.ReducedIndexI();
	mod -= IxI;
	mod *= 2.0*mu*falphaS;
	fModulus += mod;
	
        /*volumetric part*/
	mod = IxI;
	mod *= 3.0*kappa*falphaB;
	fModulus += mod;

	return(fModulus);
}

const dSymMatrixT& FDSV_KStV3D::S_IJ(void)
{
	double dt = fFSMatSupport.TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);
	fbetaS = exp(-taudtS);
	fbetaB = exp(-taudtB);

        Compute_C(fE);
	
	fE.PlusIdentity(-1.0);
	fE *= 0.5;

	/*equilibrium components*/
	double& mu = fMu[kEquilibrium];
	double& kappa = fKappa[kEquilibrium];

	double I1 = fE[0]+fE[1]+fE[2]; 
	dSymMatrixT devE = fE;
	devE[0] -= fthird*I1;
	devE[1] -= fthird*I1;
	devE[2] -= fthird*I1;

	/*deviatoric part*/
	fStress = devE;
	fStress *= 2.0*mu;

	/*volumetric part*/
	fStress.PlusIdentity(kappa*I1);

	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	if(fFSMatSupport.RunState() == GlobalT::kFormRHS)
	{
	        double& mu = fMu[kNonEquilibrium];
		double& kappa = fKappa[kEquilibrium];

		/*deviatoric part*/       
		fDevInStress = devE;
		fDevInStress *= 2.0*mu;
	        fDevOverStress.SetToCombination(fbetaS, fDevOverStress_n, 
						falphaS, fDevInStress,
						-falphaS, fDevInStress_n);

		/*volumetric part*/
		fVolInStress[0] = kappa*I1;
		fVolOverStress[0] = fbetaB*fVolOverStress_n[0] + falphaB* 
		                   (fVolInStress[0]-fVolInStress_n[0]);

		Store(element,CurrIP());
	}
	fStress += fDevOverStress;
	fStress.PlusIdentity(fVolOverStress[0]);

	return(fStress);
}
