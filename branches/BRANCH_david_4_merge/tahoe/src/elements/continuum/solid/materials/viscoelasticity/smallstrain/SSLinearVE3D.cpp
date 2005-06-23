/* $Id: SSLinearVE3D.cpp,v 1.5 2004-07-15 08:29:34 paklein Exp $ */
/* created: TDN (5/31/2001) */
#include "SSLinearVE3D.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>

#include "ExceptionT.h"

using namespace Tahoe;

const double third = 1.0/3.0;
const int kNumOutputVar = 1;
static const char* Labels[kNumOutputVar] = {"Dvisc"};

SSLinearVE3D::SSLinearVE3D(void):
	ParameterInterfaceT("linear_viscoelastic")
{

}	

double SSLinearVE3D::StrainEnergyDensity(void)
{
        /*get strains*/
        fe = e();
	
	/*equilibrium component of energy*/
	double mu = fMu[kEquilibrium];
	double kappa = fKappa[kNonEquilibrium];

	double I1 = fe[0]+fe[1]+fe[2]; 

	fe[0] -= third*I1;
	fe[1] -= third*I1;
	fe[2] -= third*I1;
	
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
        fe[0] += fmeanQ[0]/kappa*third;
        fe[1] += fmeanQ[0]/kappa*third;
        fe[2] += fmeanQ[0]/kappa*third;
    
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
	fModulus(0,0) = fModulus(1,1) =  fModulus(2,2) = 2.0*mu*(1.0 - third);
	fModulus(3,3) = fModulus(4,4) =  fModulus(5,5) = mu;
	fModulus(0,1) = fModulus(0,2) =  fModulus(1,2) = -2.0*mu*third;
	fModulus(1,0) = fModulus(2,0) =  fModulus(2,1) = -2.0*mu*third;

        /*volumetric part*/
	fModulus(0,0) += kappa; fModulus(1,1) += kappa; fModulus(2,2) += kappa;
	fModulus(0,1) += kappa; fModulus(0,2) += kappa; fModulus(1,2) += kappa;
	fModulus(1,0) += kappa; fModulus(2,0) += kappa; fModulus(2,1) += kappa;

	/*non-equilibrium component*/
	mu = fMu[kNonEquilibrium];
	kappa = fKappa[kNonEquilibrium];

	/*deviatoric part*/
	fModMat = 0.0;
	fModMat(0,0) = fModMat(1,1) =  fModMat(2,2) = 2.0*mu*falphaS*(1.0 - third);
	fModMat(3,3) = fModMat(4,4) =  fModMat(5,5) = mu*falphaS;
	fModMat(0,1) = fModMat(0,2) =  fModMat(1,2) = -2.0*mu*falphaS*third;
	fModMat(1,0) = fModMat(2,0) =  fModMat(2,1) = -2.0*mu*falphaS*third;
	
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

	fe[0] -= third*I1;
	fe[1] -= third*I1;
	fe[2] -= third*I1;
	
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

/* describe the parameters needed by the interface */
void SSLinearVE3D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSViscoelasticityT::DefineParameters(list);

	/* common limit */
	LimitT positive(0.0, LimitT::Lower);

	/* relaxation times */
	ParameterT tau_shear(ParameterT::Double, "tau_shear");
	ParameterT tau_bulk(ParameterT::Double, "tau_bulk");
	tau_shear.AddLimit(positive);
	tau_bulk.AddLimit(positive);
	list.AddParameter(tau_shear);
	list.AddParameter(tau_bulk);

	/* elastic properties */
	ParameterT mu_EQ(ParameterT::Double, "mu_EQ");
	ParameterT kappa_EQ(ParameterT::Double, "kappa_EQ");
	ParameterT mu_NEQ(ParameterT::Double, "mu_NEQ");
	ParameterT kappa_NEQ(ParameterT::Double, "kappa_NEQ");
	mu_EQ.AddLimit(positive);
	kappa_EQ.AddLimit(positive);
	mu_NEQ.AddLimit(positive);
	kappa_NEQ.AddLimit(positive);
	list.AddParameter(mu_EQ);
	list.AddParameter(kappa_EQ);
	list.AddParameter(mu_NEQ);
	list.AddParameter(kappa_NEQ);
}

/* accept parameter list */
void SSLinearVE3D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSViscoelasticityT::TakeParameterList(list);

	/* dimension work space */
	fe.Dimension(3);
	fStress.Dimension(3);
	fModulus.Dimension(6);
	fModMat.Dimension(6);
	fMu.Dimension(2);
	fKappa.Dimension(2);

	/* relaxation times */
	ftauS = list.GetParameter("tau_shear");
	ftauB = list.GetParameter("tau_bulk");

	/* elastic properties */
	fMu[kEquilibrium] = list.GetParameter("mu_EQ");
	fKappa[kEquilibrium] = list.GetParameter("kappa_EQ");
	fMu[kNonEquilibrium] = list.GetParameter("mu_NEQ");
	fKappa[kNonEquilibrium] = list.GetParameter("kappa_NEQ");
}
