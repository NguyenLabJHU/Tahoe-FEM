 /* $Id: SSViscoelasticityT.cpp,v 1.3 2008-08-09 15:00:37 tdnguye Exp $ */
#include "SSViscoelasticityT.h"
#include "ExceptionT.h"
#include "ParameterContainerT.h"
#include "SSMatSupportT.h"

#include <math.h>
#include <iostream.h>

#include "ExceptionT.h"

using namespace Tahoe;

const double third = 1.0/3.0;
const int kNumOutputVar = 1;
static const char* Labels[kNumOutputVar] = {"Dvisc"};

using namespace Tahoe;

SSViscoelasticityT::SSViscoelasticityT(void):
	ParameterInterfaceT("linear_prony_series_viscoelasticity")
{
	/*set default*/
	fnumprocess = 1;
	if (NumSD() == 2)
		fConstraint = kPlaneStrain;
}	

void SSViscoelasticityT::PointInitialize(void)
{
	int ip = CurrIP();
	/* allocate element storage */
	ElementCardT& element = CurrentElement();

	if (element.Flag() == ElementCardT::kON){
		if (ip == 0) {
			element.Dimension(0, fnstatev*NumIP());
			element.DoubleData() = 0.0;
		}
	if (ip == NumIP() - 1)
		UpdateHistory();   /*set current to last*/
	}

	if (element.Flag() == ElementCardT::kMarkON)
	{
		if (ip ==0) {
				element.Dimension(0,fnstatev*NumIP());
		}
		Load(element,ip);
		
		const dSymMatrixT& strain = e();
		if (NumSD() == 2)
		{
			fStrain3D = 0;
			fStrain3D[0] = strain[0];
			fStrain3D[1] = strain[1];
			fStrain3D[5] = strain[2];
		}
		else 
			fStrain3D = strain;
		
		double I1 = fStrain3D[0]+fStrain3D[1]+fStrain3D[2];
		fStrain3D[0] -= third*I1;
		fStrain3D[1] -= third*I1;
		fStrain3D[2] -= third*I1;

		double muEQ = fMu[kEquilibrium];
		double kappaEQ = fKappa[kEquilibrium];

		double dt = fSSMatSupport->TimeStep();
		for (int k = 0; k < fnumprocess; k++)
		{
			fmeanQ[k] = 0.0;
			fmeanQ_n[k] = 0.0;
			fdevQ[k] = 0.0;
			fdevQ_n[k] = 0.0;

			double muNEQ = fMu[k+1];
			double kappaNEQ = fKappa[k+1];

			double taudtS = dt/ftauS[k];
			double taudtB = dt/ftauB[k];
		
			double alphaS = exp(-0.5*taudtS);
			double alphaB = exp(-0.5*taudtB);
			double betaS = exp(-taudtS);
			double betaB = exp(-taudtB);

			 dSymMatrixT& devQ_n = fdevQ_n[k];
			 dSymMatrixT& devSin_n = fdevSin_n[k];
			 dArrayT& meanQ_n = fmeanQ_n[k];
			 dArrayT& meanSin_n = fmeanSin_n[k];
		
			/*calculate out of plane component*/
			if (NumSD() == 2 && Constraint() == kPlaneStress)
			{
				double K1 = 3.0*(kappaEQ + kappaNEQ*alphaB) + 4.0*(muEQ + muNEQ*alphaS);
				double K2 = -3.0*(betaB*meanQ_n[0] + betaS*devQ_n[2] - alphaB*meanSin_n[0] - alphaS*devSin_n[2]);
				K2 += (-3.0*kappaNEQ*alphaB + 2.0*muNEQ*alphaS - 3.0*kappaEQ + 2.0*muEQ)*(strain[0]+strain[1]);
	
				fStrain3D[2] = K2/K1;
			}

			/*deviatoric part*/
			dSymMatrixT& devSin = fdevSin[k];
			dArrayT& meanSin = fmeanSin[k];
			
			devSin = fStrain3D;
			devSin *= 2.0*muNEQ;
			devSin_n = devSin;

			/*volumetric part*/
			meanSin[0] = kappaNEQ*I1;
			meanSin_n[0] = meanSin[0];
		}
		Store(element, ip);
	}
}
 
void SSViscoelasticityT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	int nip = NumIP();
	SS_Visc_Support::Update(element,nip);
	
}

void SSViscoelasticityT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	int nip = NumIP();
	SS_Visc_Support::Update(element,nip);
}

/* information about subordinate parameter lists */
void SSViscoelasticityT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("ss_eq_params", ParameterListT::Once);
	sub_list.AddSub("ss_neq_params", ParameterListT::OnePlus);

}

/* describe the parameters needed by the interface */
ParameterInterfaceT* SSViscoelasticityT::NewSub(const StringT& name) const
{

	/* common limit */
	LimitT positive(0.0, LimitT::Lower);

	if (name == "ss_eq_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT moduli("params_eq");
		moduli.AddParameter(ParameterT::Double, "mu_EQ");
		moduli.AddParameter(ParameterT::Double, "kappa_EQ");
		choice->AddSub(moduli, ParameterListT::Once, false);
		
		return choice;
	}
	else if (name == "ss_neq_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);

		ParameterContainerT moduli("params_neq");
		moduli.AddParameter(ParameterT::Double, "mu_NEQ");
		moduli.AddParameter(ParameterT::Double, "kappa_NEQ");
		moduli.AddParameter(ParameterT::Double, "tau_shear");
		moduli.AddParameter(ParameterT::Double, "tau_bulk");

		choice->AddSub(moduli, ParameterListT::Once, false);

		return choice;
	}
	else return (SSSolidMatT::NewSub(name));
}

/* accept parameter list */
void SSViscoelasticityT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SSViscoelasticityT::TakeParameterList";
	/* inherited */
	SSSolidMatT::TakeParameterList(list);

	fnumprocess = list.NumLists("ss_neq_params");

	if (NumSD() == 2 && Constraint() == SolidMaterialT::kPlaneStress && fnumprocess >1)
		ExceptionT::GeneralFail(caller, "Plane stress formulation only implemented for 1 neq process");
	
	fMu.Dimension(fnumprocess+1);
	fKappa.Dimension(fnumprocess+1);
	ftauS.Dimension(fnumprocess);
	ftauB.Dimension(fnumprocess);
	
	const ParameterListT& eq_params = list.GetListChoice(*this, "ss_eq_params");
	fMu[0] = eq_params.GetParameter("mu_EQ");
	fKappa[0] = eq_params.GetParameter("kappa_EQ");
	
	for (int i = 0; i < fnumprocess; i++)
	{
		const ParameterListT& neq_params = list.GetListChoice(*this, "ss_neq_params",i);
		fMu[i+1] = neq_params.GetParameter("mu_NEQ");
		fKappa[i+1] = neq_params.GetParameter("kappa_NEQ");
		ftauS[i] = neq_params.GetParameter("tau_shear");
		ftauB[i] = neq_params.GetParameter("tau_bulk");
	}
	
	/* dimension work space */
	fStrain3D.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
	fStress.Dimension(NumSD());

	SS_Visc_Support::SetSateVariables();
}

double SSViscoelasticityT::StrainEnergyDensity(void)
{
	/*get strains*/
	const dSymMatrixT& strain = e();
	fStress = s_ij();
	double energy = 0.5*strain.ScalarProduct(fStress);
	
	return(energy);
}

const dMatrixT& SSViscoelasticityT::C_IJKL(void) 
{ 
	return(c_ijkl());
}

const dSymMatrixT& SSViscoelasticityT::S_IJ(void)
{
	return(s_ij());
}

const dMatrixT& SSViscoelasticityT::c_ijkl(void)
{        
 	double dt = fSSMatSupport->TimeStep();
	if (NumSD()==2)
	{
		/*equilibrium*/
		SetModulus2D(fModulus, dt, Constraint(), -1);

		/*nonequilibrium*/
		for (int i = 0; i< fnumprocess; i++)
			SetModulus2D(fModulus,dt, Constraint(), i);
	}
	else
	{
		/*equilibrium*/
		SetModulus(fModulus,dt, -1);
	
		/*nonequilibrium*/
		for (int i = 0; i< fnumprocess; i++)
			SetModulus(fModulus,dt, i);
	}
	return(fModulus);
}

const dSymMatrixT& SSViscoelasticityT::s_ij(void)
{
	double dt = fSSMatSupport->TimeStep();
	const dSymMatrixT& strain = e();
	
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
	if (NumSD() == 2)
	{
		/*equilibrium components*/
		ComputeStress2D(strain, fStress, dt, Constraint(), -1);
	
		/*add nonequilibrium component*/
		if(fSSMatSupport->RunState() == GlobalT::kFormRHS)
		{
			for (int i = 0; i < fnumprocess; i++)
				ComputeStress2D(strain, fStress, dt, Constraint(), i);
			Store(element,CurrIP());
		}
		else
		{
			for (int i = 0; i < fnumprocess; i++)
			{
				const dSymMatrixT& devQ = fdevQ[i];
				const dArrayT& meanQ = fmeanQ[i];
		
				fStress[0] += devQ[0];
				fStress[1] += devQ[1];
				fStress[2] += devQ[5];

				fStress[0] += meanQ[0];
				fStress[1] += meanQ[0];
			}
		}
	}
	else
	{
		/*equilibrium components*/
		ComputeStress(strain, fStress, dt, -1);
	
		/*non-equilibrium components*/
		ElementCardT& element = CurrentElement();
		Load(element, CurrIP());

		/*nonequilibrium*/
		if(fSSMatSupport->RunState() == GlobalT::kFormRHS)
		{
			for (int i = 0; i< fnumprocess; i++)
				ComputeStress(strain, fStress, dt, i);

			Store(element,CurrIP());
		}
		else
		{
			for (int i = 0; i< fnumprocess; i++)
			{
				fStress += fdevQ[i];
				const dArrayT& meanQ = fmeanQ[i];
				fStress[0] += meanQ[0];
				fStress[1] += meanQ[0];
				fStress[2] += meanQ[0];
			}
		}
	}
	return(fStress);
}

int SSViscoelasticityT::NumOutputVariables() const {return kNumOutputVar;}

void SSViscoelasticityT::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}
	
void SSViscoelasticityT::ComputeOutput(dArrayT& output)
{
	/*non-equilibrium components*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());

	output[0] = 0.0;
	
	for (int i = 0; i< fnumprocess; i++)
	{
		double etaS = fMu[i+1]*ftauS[i];
		double etaB = fKappa[i+1]*ftauB[i];

		const dSymMatrixT& devQ = fdevQ[i];
		const dArrayT& meanQ = fmeanQ[i];
		 	
		output[0] += 0.5*(0.5/etaS*devQ.ScalarProduct() + 1.0/etaB*meanQ[0]*meanQ[0]); 
	}
}	

