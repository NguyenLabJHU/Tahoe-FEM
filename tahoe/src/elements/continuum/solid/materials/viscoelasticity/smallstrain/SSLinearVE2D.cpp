/* $Id: SSLinearVE2D.cpp,v 1.4.4.1 2005-02-24 01:14:20 thao Exp $ */
/* created: TDN (5/31/2001) */
#include "SSLinearVE2D.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>
#include "ifstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

const int kNumOutputVar = 1;
static const char* Labels[kNumOutputVar] = {"Dvisc"};

SSLinearVE2D::SSLinearVE2D(ifstreamT& in, const SSMatSupportT& support):
	SSViscoelasticityT(in, support),
	fStress(2),
	fModulus(3),
	fModMat(3),
	fStrain3D(3),
 	fStress3D(3),
	fMu(2),
	fKappa(2),
	fthird(1.0/3.0)
{
        double& mu_EQ = fMu[kEquilibrium];
	double& kappa_EQ = fKappa[kEquilibrium]; 

	double& mu_NEQ = fMu[kNonEquilibrium]; 
	double& kappa_NEQ = fKappa[kNonEquilibrium];

	in >> ftauS;
	in >> ftauB;

	in >> mu_EQ;
	in >> kappa_EQ;

	in >> mu_NEQ;
	in >> kappa_NEQ;
}	

void SSLinearVE2D::Print(ostream& out) const
{
	/* inherited */
	SSViscoelasticityT::Print(out);
	out << "Equilibrium:\n";
	out << "     Shear Modulus: "<<fMu[0]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[0]<<'\n';
	out << "Non-Equilibrium:\n";
	out << "     Shear Modulus: "<<fMu[1]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[1]<<'\n';
	out << "Relaxation time: \n";
	out << "     Shear relaxation time: "<<ftauS<<'\n';
	out << "     Bulk relaxation time: "<<ftauB<<'\n';
}

void SSLinearVE2D::PrintName(ostream& out) const
{
	/* inherited */
	SSViscoelasticityT::PrintName(out);
	out << "2D linear, plane strain formulation\n";
}

double SSLinearVE2D::StrainEnergyDensity(void)
{
        const dSymMatrixT& strain = e();
	fStrain3D = 0;
	fStrain3D[0] = strain[0];
	fStrain3D[1] = strain[1];
	fStrain3D[5] = strain[2];
	
	/*equilibrium components*/
	double mu = fMu[kEquilibrium];
	double kappa = fKappa[kEquilibrium];

	double I1 = strain[0]+strain[1]; 

	fStrain3D[0] -= fthird*I1;
	fStrain3D[1] -= fthird*I1;
	fStrain3D[2] -= fthird*I1;
	
	/*deviatoric part*/
	fStress3D = fStrain3D;
	fStress3D *= 2.0*mu;

	/*volumetric part*/
	fStress3D[0] += kappa*I1;
	fStress3D[1] += kappa*I1;
	fStress3D[2] += kappa*I1;
	
	/*reduce to 2D*/
	fStress[0] = fStress3D[0];
	fStress[1] = fStress3D[1];
	fStress[2] = fStress3D[5];

	double energy = 0.5*fStress.ScalarProduct(e());
	
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	mu = fMu[kNonEquilibrium];
	kappa = fKappa[kNonEquilibrium];
	
	fStrain3D = fdevQ;
	fStrain3D /= 2.0*mu;
	fStrain3D[0] += fmeanQ[0]/kappa*fthird;
	fStrain3D[1] += fmeanQ[0]/kappa*fthird;
	fStrain3D[2] += fmeanQ[0]/kappa*fthird;

	fStress3D = fdevQ;
	fStress3D[0] += fmeanQ[0];
	fStress3D[1] += fmeanQ[0];
	fStress3D[2] += fmeanQ[0];

	energy += 0.5*fStress3D.ScalarProduct(fStrain3D);
	return(energy);
}

const dMatrixT& SSLinearVE2D::C_IJKL(void) 
{ 
	return(c_ijkl());
}

const dSymMatrixT& SSLinearVE2D::S_IJ(void)
{
	return(s_ij());
}

const dMatrixT& SSLinearVE2D::c_ijkl(void)
{        
 	double dt = fSSMatSupport->TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);

	/*equilibrium component*/
	double mu = fMu[kEquilibrium];
	double kappa = fKappa[kEquilibrium];
	
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

const dSymMatrixT& SSLinearVE2D::s_ij(void)
{
	double dt = fSSMatSupport->TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);
	fbetaS = exp(-taudtS);
	fbetaB = exp(-taudtB);

	const dSymMatrixT& strain = e();
	fStrain3D = 0;
	fStrain3D[0] = strain[0];
	fStrain3D[1] = strain[1];
	fStrain3D[5] = strain[2];
	
	/*equilibrium components*/
	double mu = fMu[kEquilibrium];
	double kappa = fKappa[kEquilibrium];

	double I1 = strain[0]+strain[1]; 

	fStrain3D[0] -= fthird*I1;
	fStrain3D[1] -= fthird*I1;
	fStrain3D[2] -= fthird*I1;
	
	/*deviatoric part*/
	fStress3D = fStrain3D;
	fStress3D *= 2.0*mu;

	/*volumetric part*/
	fStress3D[0] += kappa*I1;
	fStress3D[1] += kappa*I1;
	fStress3D[2] += kappa*I1;
	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	if(fSSMatSupport->RunState() == GlobalT::kFormRHS)
	{
		mu = fMu[kNonEquilibrium];
		kappa = fKappa[kNonEquilibrium];

		/*deviatoric part*/       
		fdevSin = fStrain3D;
		fdevSin *= 2.0*mu;
		
		fdevQ[0] = fbetaS*fdevQ_n[0] + falphaS*(fdevSin[0]-fdevSin_n[0]);
		fdevQ[1] = fbetaS*fdevQ_n[1] + falphaS*(fdevSin[1]-fdevSin_n[1]);
		fdevQ[2] = fbetaS*fdevQ_n[2] + falphaS*(fdevSin[2]-fdevSin_n[2]);
		fdevQ[3] = fbetaS*fdevQ_n[3] + falphaS*(fdevSin[3]-fdevSin_n[3]);
		fdevQ[4] = fbetaS*fdevQ_n[4] + falphaS*(fdevSin[4]-fdevSin_n[4]);
		fdevQ[5] = fbetaS*fdevQ_n[5] + falphaS*(fdevSin[5]-fdevSin_n[5]);
		
		/*volumetric part*/
		fmeanSin[0] = kappa*I1;
		fmeanQ[0] = fbetaB*fmeanQ_n[0] + falphaB * (fmeanSin[0]-fmeanSin_n[0]);
        
		Store(element,CurrIP());
	}
	fStress3D += fdevQ;

	fStress3D[0] += fmeanQ[0];
	fStress3D[1] += fmeanQ[0];
	fStress3D[2] += fmeanQ[0];

	fStress[0] = fStress3D[0];
	fStress[1] = fStress3D[1];
	fStress[2] = fStress3D[5];
	return(fStress);
}
int SSLinearVE2D::NumOutputVariables() const {return kNumOutputVar;}

void SSLinearVE2D::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}
	
void SSLinearVE2D::ComputeOutput(dArrayT& output)
{
	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	double etaS = fMu[kNonEquilibrium]*ftauS;
	double etaB = fKappa[kNonEquilibrium]*ftauB;
	
	output[0] = 0.5*(0.5/etaS*fdevQ.ScalarProduct() + 1.0/etaB*fmeanQ[0]*fmeanQ[0]); 
}	
