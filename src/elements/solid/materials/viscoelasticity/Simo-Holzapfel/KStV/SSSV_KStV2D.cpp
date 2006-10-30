/* $Id: SSSV_KStV2D.cpp,v 1.1 2006-10-30 23:32:06 thao Exp $ */
/* created: TDN (5/31/2001) */
#include "SSSV_KStV2D.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>
#include "ifstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

const int kNumOutputVar = 4;
static const char* Labels[kNumOutputVar] = {"Dvisc","Er", "Iep_v", "IIe_v"};

SSSV_KStV2D::SSSV_KStV2D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("SSSV_KStV2D"),
	SSSimoViscoT(in, support),
        fStress(2),
	fModulus(3),
	fModMat(3),
	fStrain3D(3),
 	fStress3D(3),
	fMu(2),
	fKappa(2),
	fthird(1.0/3.0)
{
ExceptionT::GeneralFail("SSSV_KStV2D::SSSV_KStV2D", "out of date");
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

        IsotropicT::Set_mu_kappa(mu_EQ, kappa_EQ);
}	

#if 0
void SSSV_KStV2D::Print(ostream& out) const
{
	/* inherited */
	SSSimoViscoT::Print(out);
	out << "Equilibrium Potential:\n";
	out << "     Shear Modulus: "<<fMu[0]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[0]<<'\n';
	out << "Non-Equilibrium Potential:\n";
	out << "     Shear Modulus: "<<fMu[1]<<'\n';
	out << "     Bulk Modulus: "<<fKappa[1]<<'\n';
	out << "Relaxation time: \n";
	out << "     Shear relaxation time: "<<ftauS<<'\n';
	out << "     Bulk relaxation time: "<<ftauB<<'\n';
}

void SSSV_KStV2D::PrintName(ostream& out) const
{
	/* inherited */
	SSSimoViscoT::PrintName(out);
	out << "2D plane strain formulation\n";
	out << "Equilibrium/Non-Equilibrium Potential:\n";
	out << "Kirchoff St. Venant\n";
	out << "Kirchoff St. Venant\n";
}
#endif

double SSSV_KStV2D::StrainEnergyDensity(void)
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

const dMatrixT& SSSV_KStV2D::C_IJKL(void) 
{ 
	return(c_ijkl());
}

const dSymMatrixT& SSSV_KStV2D::S_IJ(void)
{
	return(s_ij());
}

const dMatrixT& SSSV_KStV2D::c_ijkl(void)
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

const dSymMatrixT& SSSV_KStV2D::s_ij(void)
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

	if (fSSMatSupport->RunState() == GlobalT::kFormRHS)
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

		/*evaluate viscous strains*/
		fViscStrain = 0.0;
		fViscStrain -= fdevQ;
		fViscStrain /= 2.0*mu;
		fViscStrain[0] -= fmeanQ[0]*fthird/kappa;
		fViscStrain[1] -= fmeanQ[0]*fthird/kappa;
		fViscStrain[2] -= fmeanQ[0]*fthird/kappa;

		fViscStrain[0] += strain[0];
		fViscStrain[1] += strain[1];
		fViscStrain[5] += strain[2];

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

/*Note to be called only during post processing*/
const dArrayT& SSSV_KStV2D::InternalStrainVars(void)
{
        double mu = fMu[kNonEquilibrium];
	double kappa = fKappa[kNonEquilibrium];
	const dSymMatrixT& strain = e();

	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
        /*evaluate viscous strains*/
        fViscStrain = 0.0;
	fViscStrain -= fdevQ;
	fViscStrain /= 2.0*mu;
        fViscStrain[0] -= fmeanQ[0]*fthird/kappa;
        fViscStrain[1] -= fmeanQ[0]*fthird/kappa;
        fViscStrain[2] -= fmeanQ[0]*fthird/kappa;

        fViscStrain[0] += strain[0];
	fViscStrain[1] += strain[1];
	fViscStrain[5] += strain[2];

	return(fViscStrain);
}

const dArrayT& SSSV_KStV2D::InternalStressVars(void)
{
	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
        /*evaluate viscous stresses*/
	fViscStress = fdevQ;
	fViscStress[0] += fmeanQ[0];
	fViscStress[1] += fmeanQ[0];
	fViscStress[2] += fmeanQ[0];
	
	return(fViscStress);
}

int SSSV_KStV2D::NumOutputVariables() const {return kNumOutputVar;}

void SSSV_KStV2D::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}
	
void SSSV_KStV2D::ComputeOutput(dArrayT& output)
{
        double mu = fMu[kNonEquilibrium];
	double kappa = fKappa[kNonEquilibrium];

	/*non-equilibrium components*/
        ElementCardT& element = CurrentElement(); 
	Load(element, CurrIP());

	const dArrayT& viscstrain = InternalStrainVars();
        
	double etaS = fMu[kNonEquilibrium]*ftauS;
	double etaB = fKappa[kNonEquilibrium]*ftauB;

	output[0] = 0.5*(0.5/etaS*fdevQ.ScalarProduct() + 1.0/etaB*fmeanQ[0]*fmeanQ[0]);
	double I1 = output[1] = fViscStrain[0]+fViscStrain[1]+fViscStrain[2];
	I1 *= fthird;
	output[2] = sqrt(2.0*fthird*((fViscStrain[0]-I1)*(fViscStrain[0]-I1)
                  +(fViscStrain[1]-I1)*(fViscStrain[1]-I1) 
                  + (fViscStrain[2]-I1)*(fViscStrain[2]-I1)
                  +2.0*fViscStrain[3]*fViscStrain[3]
                  +2.0*fViscStrain[4]*fViscStrain[4]
                  +2.0*fViscStrain[5]*fViscStrain[5]));
}	

/* information about subordinate parameter lists */
void SSSV_KStV2D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	IsotropicT::DefineSubs(sub_list);
	SSSimoViscoT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SSSV_KStV2D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = IsotropicT::NewSub(name);
	if (sub)
		return sub;
	else
		return SSSimoViscoT::NewSub(name);
}

/* accept parameter list */
void SSSV_KStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list);
	SSSimoViscoT::TakeParameterList(list);
}
