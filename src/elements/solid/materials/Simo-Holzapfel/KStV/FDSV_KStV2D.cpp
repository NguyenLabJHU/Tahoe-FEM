/* $Id: FDSV_KStV2D.cpp,v 1.4 2004-07-15 08:28:49 paklein Exp $ */
/* created:   TDN (5/31/2001) */

#include "FDSV_KStV2D.h"

#include <math.h>
#include <iostream.h>
#include "ifstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

FDSV_KStV2D::FDSV_KStV2D(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("FDSV_KStV2D"),
//        Material2DT(in),
	FDSimoViscoBaseT(in, support),
	fStress(2),
	fModulus(3),
	fModMat(3),
	fE(2),
	fMu(2),
	fKappa(2),
	fthird(1.0/3.0)
{
ExceptionT::GeneralFail("", "out of date");
#if 0
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
#endif
}	

#if 0
void FDSV_KStV2D::Print(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::Print(out);
	out << "Equilibrium Potential\n";
	out << "     Shear Modulus: "<<fMu[0]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[0]<<'\n';
	out << "Non-Equilibrium Potential\n";
	out << "     Shear Modulus: "<<fMu[1]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[1]<<'\n';
	out << "Constant relaxation time\n";
	out << "     Shear relaxation time: "<<ftauS<<'\n';
	out << "     Bulk relaxation time: "<<ftauB<<'\n';
}

void FDSV_KStV2D::PrintName(ostream& out) const
{
	/* inherited */
	FDSimoViscoBaseT::PrintName(out);
	out<<"2D plane strain formulation \n";
	out << "Equilibrium/Non-Equilibrium Potential:\n";
	out << "Kirchoff St. Venant\n";
	out << "Kirchoff St. Venant\n";
}
#endif

const dMatrixT& FDSV_KStV2D::c_ijkl(void) 
{ 
        const dMatrixT& F_mech = F_mechanical(); 
        fModulus = PushForward(F_mech,C_IJKL());
        fModulus /= F_mech.Det();
        return fModulus;
}

const dSymMatrixT& FDSV_KStV2D::s_ij(void)
{
        const dMatrixT& F_mech = F_mechanical(); 
        fStress = PushForward(F_mech,S_IJ());
        fStress /= F_mech.Det();
        return fStress;
}

const dMatrixT& FDSV_KStV2D::C_IJKL(void)
{        
        /*equilibrium component*/
 	double dt = fFSMatSupport->TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);

    /*equilibrium component*/
	double& mu = fMu[kEquilibrium];
	double& kappa = fKappa[kEquilibrium];

    /*deviatoric part*/
	fModulus = 0.0;
	fModulus(0,0) = fModulus(1,1) = 2.0*mu*(1.0 - fthird);
	fModulus(2,2) = mu;
	fModulus(0,1) =	fModulus(1,0) = -2.0*mu*fthird;

    /*volumetric part*/
	fModulus(0,0) += kappa; fModulus(1,1) += kappa; 
	fModulus(0,1) += kappa; fModulus(1,0) += kappa; 
	
	/*non-equilibrium component*/
	mu = fMu[kNonEquilibrium];
	kappa = fKappa[kNonEquilibrium];

	/*deviatoric part*/
	fModMat = 0.0;
	fModMat(0,0) = fModMat(1,1) = 2.0*mu*falphaS*(1.0 - fthird);
	fModMat(2,2) = mu*falphaS;
	fModMat(0,1) = fModMat(1,0) = -2.0*mu*falphaS*fthird;
	
    /*volumetric part*/
	fModMat(0,0) += kappa*falphaB; fModMat(1,1) += kappa*falphaB; 
	fModMat(0,1) += kappa*falphaB; fModMat(1,0) += kappa*falphaB; 
	
    fModulus += fModMat;
 
 	return(fModulus);
}

const dSymMatrixT& FDSV_KStV2D::S_IJ(void)
{
	double dt = fFSMatSupport->TimeStep();
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

	fE[0] -= fthird*I1;
	fE[1] -= fthird*I1;
	
	/*deviatoric part*/
	fStress = fE;
	fStress *= 2.0*mu;

	/*volumetric part*/
	fStress[0] += kappa*I1;
	fStress[1] += kappa*I1;

	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	if(fFSMatSupport->RunState() == GlobalT::kFormRHS)
	{
		mu = fMu[kNonEquilibrium];
		kappa = fKappa[kEquilibrium];

		/*deviatoric part*/       
		fdevSin = fE;
		fdevSin *= 2.0*mu;
		
		fdevQ[0] = fbetaS*fdevQ_n[0] + falphaS*(fdevSin[0]-fdevSin_n[0]);
		fdevQ[1] = fbetaS*fdevQ_n[1] + falphaS*(fdevSin[1]-fdevSin_n[1]);
		fdevQ[2] = fbetaS*fdevQ_n[2] + falphaS*(fdevSin[2]-fdevSin_n[2]);
		
		/*volumetric part*/
		fmeanSin[0] = kappa*I1;
		fmeanQ[0] = fbetaB*fmeanQ_n[0] + falphaB * (fmeanSin[0]-fmeanSin_n[0]);

		Store(element,CurrIP());
	}
	fStress += fdevQ;

	fStress[0] += fmeanQ[0];
	fStress[1] += fmeanQ[0];

	return(fStress);
}
