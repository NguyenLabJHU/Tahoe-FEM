/* $Id: FDSV_KStV2D.cpp,v 1.1 2003-03-19 19:03:19 thao Exp $ */
/* created:   TDN (5/31/2001) */

#include "FDSV_KStV2D.h"

#include <math.h>
#include <iostream.h>
#include "fstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

FDSV_KStV2D::FDSV_KStV2D(ifstreamT& in, const FSMatSupportT& support):
        Material2DT(in),
	FDSimoViscoBaseT(in, support),
	fStress(2),
	fModulus(3),
	fE(2),
	fMu(2),
	fKappa(2),
	fthird(1.0/3.0)
{
	if (fConstraintOption == Material2DT::kPlaneStress)
	{
	        cout << "Plane Stress formulation is not implemented\n";
		throw ExceptionT::kBadInputValue;
	}

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

void FDSV_KStV2D::Print(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::Print(out);
	out << "Equilibrium Potential";
	out << "     Shear Modulus: "<<fMu[0]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[0]<<'\n';
	out << "Non-Equilibrium Potential";
	out << "     Shear Modulus: "<<fMu[1]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[1]<<'\n';
	out<<"Constant relaxation time\n";
	out<<"     Shear relaxation time: "<<ftauS<<'\n';
	out<<"     Bulk relaxation time: "<<ftauB<<'\n';
}

void FDSV_KStV2D::PrintName(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::PrintName(out);
	out << "Equilibrium/Non-Equilibrium Potential:\n";
	out << "Kirchoff St. Venant\n";
	out << "Kirchoff St. Venant\n";
}

const dMatrixT& FDSV_KStV2D::c_ijkl(void) 
{ 
        const dMatrixT& F_mech = F_mechanical(); 
	/*put in plane stress correction*/
        fModulus = PushForward(F_mech,C_IJKL());
        fModulus /= F_mech.Det();
        return fModulus;
}

const dSymMatrixT& FDSV_KStV2D::s_ij(void)
{
        const dMatrixT& F_mech = F_mechanical(); 
	/*put in plane stress correction*/
        fStress = PushForward(F_mech,S_IJ());
        fStress /= F_mech.Det();
        return fStress;
}

const dMatrixT& FDSV_KStV2D::C_IJKL(void)
{        
        /*equilibrium component*/
	double& mu = fMu[kEquilibrium];
	double& kappa = fKappa[kEquilibrium];

        /*deviatoric part*/
	dMatrixT IxI(3);
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
	dMatrixT mod(3);
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

const dSymMatrixT& FDSV_KStV2D::S_IJ(void)
{
	double dt = fFSMatSupport.TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);
	fbetaS = exp(-taudtS);
	fbetaB = exp(-taudtB);

        Compute_C(fE);
	
	fE[0] -= 1.0;
	fE[1] -= 1.0;
	fE *= 0.5;

	/*equilibrium components*/
	double& mu = fMu[kEquilibrium];
	double& kappa = fKappa[kEquilibrium];

	double I1 = fE[0]+fE[1]; 
	dSymMatrixT devE = fE;

	devE[0] -= fthird*I1;
	devE[1] -= fthird*I1;
	
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
	        mu = fMu[kNonEquilibrium];
		kappa = fKappa[kEquilibrium];

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
