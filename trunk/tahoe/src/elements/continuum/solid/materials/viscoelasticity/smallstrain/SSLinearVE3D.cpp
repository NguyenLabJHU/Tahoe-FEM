/* $Id: SSLinearVE3D.cpp,v 1.4 2004-06-17 07:41:22 paklein Exp $ */
/* created: TDN (5/31/2001) */
#include "SSLinearVE3D.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>
#include "ifstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

const int kNumOutputVar = 1;
static const char* Labels[kNumOutputVar] = {"Dvisc"};

SSLinearVE3D::SSLinearVE3D(ifstreamT& in, const SSMatSupportT& support):
	SSViscoelasticityT(in, support),
	fe(3),
	fStress(3),
	fModulus(6),
	fModMat(6),
	fMu(2),
	fKappa(2),
	fthird(1.0/3.0)
{
	in >> ftauS;
	in >> ftauB;

        double& mu_EQ = fMu[kEquilibrium];
	double& kappa_EQ = fKappa[kEquilibrium]; 

	double& mu_NEQ = fMu[kNonEquilibrium]; 
	double& kappa_NEQ = fKappa[kNonEquilibrium];

	in >> mu_EQ;
	in >> kappa_EQ;

	in >> mu_NEQ;
	in >> kappa_NEQ;
}	

void SSLinearVE3D::Print(ostream& out) const
{
	/* inherited */
	SSViscoelasticityT::Print(out);
	out << "Equilibrium Potential";
	out << "     Shear Modulus: "<<fMu[0]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[0]<<'\n';
	out << "Non-Equilibrium Potential";
	out << "     Shear Modulus: "<<fMu[1]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[1]<<'\n';
	out<<"Relaxation times \n";
	out<<"     Shear relaxation time: "<<ftauS<<'\n';
	out<<"     Bulk relaxation time: "<<ftauB<<'\n';
}

void SSLinearVE3D::PrintName(ostream& out) const
{
	/* inherited */
	SSViscoelasticityT::PrintName(out);
	out << "Equilibrium/Non-Equilibrium Potential:\n";
	out << "Kirchoff St. Venant\n";
	out << "Kirchoff St. Venant\n";
}

double SSLinearVE3D::StrainEnergyDensity(void)
{
        /*get strains*/
        fe = e();
	
	/*equilibrium component of energy*/
	double mu = fMu[kEquilibrium];
	double kappa = fKappa[kNonEquilibrium];

	double I1 = fe[0]+fe[1]+fe[2]; 

	fe[0] -= fthird*I1;
	fe[1] -= fthird*I1;
	fe[2] -= fthird*I1;
	
	/*deviatoric part*/
	fStress = fe;
	fStress *= 2.0*mu;

	/*volumetric part*/
	fStress[0] += kappa*I1;
	fStress[1] += kappa*I1;
	fStress[2] += kappa*I1;

        double energy = 0.5*fStress.ScalarProduct(e());

	/*non-equilibrium component of energy*/
	/*equilibrium component of energy*/
	mu = fMu[kNonEquilibrium];
        kappa = fKappa[kNonEquilibrium];

	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

        fe = fdevQ;
        fe /= 2.0*mu;
        fe[0] += fmeanQ[0]/kappa*fthird;
        fe[1] += fmeanQ[0]/kappa*fthird;
        fe[2] += fmeanQ[0]/kappa*fthird;
    
	fStress = fdevQ;
	fStress[0] += fmeanQ[0];
	fStress[1] += fmeanQ[0];
	fStress[2] += fmeanQ[0];
    
        energy += 0.5*fStress.ScalarProduct(fe);
    
        return(energy);
}

const dMatrixT& SSLinearVE3D::C_IJKL(void) 
{ 
	return(c_ijkl());
}

const dSymMatrixT& SSLinearVE3D::S_IJ(void)
{
	return(s_ij());
}

const dMatrixT& SSLinearVE3D::c_ijkl(void)
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
	fModulus(0,0) = fModulus(1,1) =  fModulus(2,2) = 2.0*mu*(1.0 - fthird);
	fModulus(3,3) = fModulus(4,4) =  fModulus(5,5) = mu;
	fModulus(0,1) = fModulus(0,2) =  fModulus(1,2) = -2.0*mu*fthird;
	fModulus(1,0) = fModulus(2,0) =  fModulus(2,1) = -2.0*mu*fthird;

        /*volumetric part*/
	fModulus(0,0) += kappa; fModulus(1,1) += kappa; fModulus(2,2) += kappa;
	fModulus(0,1) += kappa; fModulus(0,2) += kappa; fModulus(1,2) += kappa;
	fModulus(1,0) += kappa; fModulus(2,0) += kappa; fModulus(2,1) += kappa;

	/*non-equilibrium component*/
	mu = fMu[kNonEquilibrium];
	kappa = fKappa[kNonEquilibrium];

	/*deviatoric part*/
	fModMat = 0.0;
	fModMat(0,0) = fModMat(1,1) =  fModMat(2,2) = 2.0*mu*falphaS*(1.0 - fthird);
	fModMat(3,3) = fModMat(4,4) =  fModMat(5,5) = mu*falphaS;
	fModMat(0,1) = fModMat(0,2) =  fModMat(1,2) = -2.0*mu*falphaS*fthird;
	fModMat(1,0) = fModMat(2,0) =  fModMat(2,1) = -2.0*mu*falphaS*fthird;
	
	/*volumetric part*/
	fModMat(0,0) += kappa*falphaB; fModMat(1,1) += kappa*falphaB; fModMat(2,2) += kappa*falphaB;
	fModMat(0,1) += kappa*falphaB; fModMat(0,2) += kappa*falphaB; fModMat(1,2) += kappa*falphaB;
	fModMat(1,0) += kappa*falphaB; fModMat(2,0) += kappa*falphaB; fModMat(2,1) += kappa*falphaB;

	fModulus += fModMat;
    
	return(fModulus);
}

const dSymMatrixT& SSLinearVE3D::s_ij(void)
{
	double dt = fSSMatSupport->TimeStep();
	double taudtS = dt/ftauS;
	double taudtB = dt/ftauB;

	falphaS = exp(-0.5*taudtS);
	falphaB = exp(-0.5*taudtB);
	fbetaS = exp(-taudtS);
	fbetaB = exp(-taudtB);

	fe = e();
	
	/*equilibrium components*/
	double mu = fMu[kEquilibrium];
	double kappa = fKappa[kEquilibrium];

	double I1 = fe[0]+fe[1]+fe[2]; 

	fe[0] -= fthird*I1;
	fe[1] -= fthird*I1;
	fe[2] -= fthird*I1;
	
	/*deviatoric part*/
	fStress = fe;
	fStress *= 2.0*mu;

	/*volumetric part*/
	fStress[0] += kappa*I1;
	fStress[1] += kappa*I1;
	fStress[2] += kappa*I1;

	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	if(fSSMatSupport->RunState() == GlobalT::kFormRHS)
	{
		mu = fMu[kNonEquilibrium];
		kappa = fKappa[kNonEquilibrium];

		/*deviatoric part*/       
		fdevSin = fe;
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
	fStress += fdevQ;

	fStress[0] += fmeanQ[0];
	fStress[1] += fmeanQ[0];
	fStress[2] += fmeanQ[0];

	return(fStress);
}

int SSLinearVE3D::NumOutputVariables() const {return kNumOutputVar;}

void SSLinearVE3D::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}
	
void SSLinearVE3D::ComputeOutput(dArrayT& output)
{
	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	double etaS = fMu[kNonEquilibrium]*ftauS;
	double etaB = fKappa[kNonEquilibrium]*ftauB;
	
	output[0] = 0.5*(0.5/etaS*fdevQ.ScalarProduct() + 1.0/etaB*fmeanQ[0]*fmeanQ[0]); 
}	
