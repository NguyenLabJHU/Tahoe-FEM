 /* $Id: SS_Visc_Support.cpp,v 1.1 2008-08-13 12:45:41 tdnguye Exp $ */
#include "SS_Visc_Support.h"
#include "ExceptionT.h"
#include "SolidMaterialT.h"

using namespace Tahoe;

const double third = 1.0/3.0;

SS_Visc_Support::SS_Visc_Support(void)
{

}	
void SS_Visc_Support::SetSateVariables(void)
{

	fdevQ.Dimension(fnumprocess);
	fdevSin.Dimension(fnumprocess);
	fmeanQ.Dimension(fnumprocess);
	fmeanSin.Dimension(fnumprocess);
	
	fdevQ_n.Dimension(fnumprocess);
	fdevSin_n.Dimension(fnumprocess);
	fmeanQ_n.Dimension(fnumprocess);
	fmeanSin_n.Dimension(fnumprocess);

	int ndof = 3;
	int numstress = (ndof*(ndof+1))/2;
		
	/*allocates storage for state variable array*/
	fnstatev = 0;
	fnstatev += numstress;           /*current deviatoric overstress*/
	fnstatev += numstress;           /*current deviatoric inelastic stress*/
	fnstatev ++;			 /*current mean over stress*/
	fnstatev ++; 			 /*current mean inelastic stress*/
	
	fnstatev += numstress;           /*preceding deviatoric overstress*/
	fnstatev += numstress;           /*preceding deviatoric inelastic stress*/
	fnstatev ++;			 /*preceding mean overstress*/
	fnstatev ++; 			 /*preceding mean inelastic stress*/
	
	fnstatev*=fnumprocess;
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	
	/* assign pointers to current and preceding blocks of state variable array */
	for (int i = 0; i < fnumprocess; i++)
	{
		fdevQ[i].Set(ndof, pstatev);        
		pstatev += numstress;
		fdevSin[i].Set(ndof, pstatev);	
		pstatev += numstress;
		fmeanQ[i].Set(1, pstatev);
		pstatev ++;
		fmeanSin[i].Set(1, pstatev);
		pstatev ++;
	
		fdevQ_n[i].Set(ndof, pstatev);        
		pstatev += numstress;
		fdevSin_n[i].Set(ndof, pstatev);	
		pstatev += numstress;
		fmeanQ_n[i].Set(1, pstatev);
		pstatev ++;
		fmeanSin_n[i].Set(1, pstatev);
		pstatev ++;
	}
	
	fMu.Dimension(fnumprocess+1);
	fKappa.Dimension(fnumprocess+1);
	
	ftauS.Dimension(fnumprocess);
	ftauB.Dimension(fnumprocess);
	
	fStress3D.Dimension(3);
	fStrain3D.Dimension(3);

}

void SS_Visc_Support::Update(ElementCardT& element,  int nip)
{
	/* current element */
	for (int ip = 0; ip < nip; ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		for (int i = 0; i < fnumprocess; i++)
		{
			fdevQ_n[i] = fdevQ[i];
			fdevSin_n[i] = fdevSin[i];
			fmeanQ_n[i] = fmeanQ[i];
			fmeanSin_n[i] = fmeanSin[i];
		}
		
		/* write to storage */
		Store(element, ip);
	}
}

void SS_Visc_Support::Reset(ElementCardT& element,  int nip)
{
	/* current element */
	for (int ip = 0; ip < nip; ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		for (int i = 0; i < fnumprocess; i++)
		{
			fdevQ[i] = fdevQ_n[i];
			fdevSin[i] = fdevSin_n[i];
			fmeanQ[i] = fmeanQ_n[i];
			fmeanSin[i] = fmeanSin_n[i];
		}
		/* write to storage */
		Store(element, ip);
	}
}

void SS_Visc_Support::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}

void SS_Visc_Support::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;
}

/*************** protected *****************/

void SS_Visc_Support::ComputeStress(const dSymMatrixT& strain, dSymMatrixT& stress, double dt, int pindex)
{

	fStrain3D = strain;
	double I1 = fStrain3D[0]+fStrain3D[1]+fStrain3D[2]; 

	fStrain3D[0] -= third*I1;
	fStrain3D[1] -= third*I1;
	fStrain3D[2] -= third*I1;
	
	if (pindex == -1)
	{
		/*equilibrium components*/
		double mu = fMu[pindex+1];
		double kappa = fKappa[pindex+1];

		/*deviatoric part*/
		stress = fStrain3D;
		stress *= 2.0*mu;

		/*volumetric part*/
		stress[0] += kappa*I1;
		stress[1] += kappa*I1;
		stress[2] += kappa*I1;
	}
	else
	{
		/*nonequilibrium components.  Stress is accumulated*/
		double mu = fMu[pindex+1];
		double kappa = fKappa[pindex+1];

		double taudtS = dt/ftauS[pindex];
		double taudtB = dt/ftauB[pindex];

		double alphaS = exp(-0.5*taudtS);
		double alphaB = exp(-0.5*taudtB);
		double betaS = exp(-taudtS);
		double betaB = exp(-taudtB);

		dSymMatrixT& devQ = fdevQ[pindex];
		dSymMatrixT& devQ_n = fdevQ_n[pindex];
		dSymMatrixT& devSin = fdevSin[pindex];
		dSymMatrixT& devSin_n = fdevSin_n[pindex];
		
		/*deviatoric part*/       
		devSin = fStrain3D;
		devSin *= 2.0*mu;
	
		devQ[0] = betaS*devQ_n[0] + alphaS*(devSin[0]-devSin_n[0]);
		devQ[1] = betaS*devQ_n[1] + alphaS*(devSin[1]-devSin_n[1]);
		devQ[2] = betaS*devQ_n[2] + alphaS*(devSin[2]-devSin_n[2]);
		devQ[3] = betaS*devQ_n[3] + alphaS*(devSin[3]-devSin_n[3]);
		devQ[4] = betaS*devQ_n[4] + alphaS*(devSin[4]-devSin_n[4]);
		devQ[5] = betaS*devQ_n[5] + alphaS*(devSin[5]-devSin_n[5]);
		
		dArrayT& meanSin = fmeanSin[pindex];
		dArrayT& meanSin_n = fmeanSin_n[pindex];
		dArrayT& meanQ = fmeanQ[pindex];
		dArrayT& meanQ_n = fmeanQ_n[pindex];
		
		/*volumetric part*/
		meanSin[0] = kappa*I1;
		meanQ[0] = betaB*meanQ_n[0] + alphaB * (meanSin[0]-meanSin_n[0]);
		
		stress += devQ;
		stress[0] += meanQ[0];
		stress[1] += meanQ[0];
		stress[2] += meanQ[0];
	}
	
}

void SS_Visc_Support::ComputeStress2D(const dSymMatrixT& strain, dSymMatrixT& stress, double dt, int constraint, int pindex)
{
	double muEQ = fMu[kEquilibrium];
	double kappaEQ = fKappa[kEquilibrium];

	fStrain3D = 0;
	fStrain3D[0] = strain[0];
	fStrain3D[1] = strain[1];
	fStrain3D[5] = strain[2];
	/*calculate out of plane component*/

	const dSymMatrixT& devQ_n = fdevQ_n[kNonEquilibrium-1];
	const dSymMatrixT& devSin_n = fdevSin_n[kNonEquilibrium-1];
	const dArrayT& meanSin_n = fmeanSin_n[kNonEquilibrium-1];
	const dArrayT& meanQ_n = fmeanQ_n[kNonEquilibrium-1];

	if (constraint == SolidMaterialT::kPlaneStress)
	{
		double muNEQ = fMu[kNonEquilibrium];
		double kappaNEQ = fKappa[kNonEquilibrium];

		double taudtS = dt/ftauS[kNonEquilibrium-1];
		double taudtB = dt/ftauB[kNonEquilibrium-1];
		double alphaS = exp(-0.5*taudtS);
		double alphaB = exp(-0.5*taudtB);
		double betaS = exp(-taudtS);
		double betaB = exp(-taudtB);
	
		double K1 = 3.0*(kappaEQ + kappaNEQ*alphaB) + 4.0*(muEQ + muNEQ*alphaS);
		double K2 = -3.0*(betaB*meanQ_n[0] + betaS*devQ_n[2] - alphaB*meanSin_n[0] - alphaS*devSin_n[2]);
		K2 += (-3.0*kappaNEQ*alphaB + 2.0*muNEQ*alphaS - 3.0*kappaEQ + 2.0*muEQ)*(strain[0]+strain[1]);

		fStrain3D[2] = K2/K1;
	}

	double I1 = fStrain3D[0]+fStrain3D[1]+fStrain3D[2];
	fStrain3D[0] -= third*I1;
	fStrain3D[1] -= third*I1;
	fStrain3D[2] -= third*I1;

	if (pindex == -1)
	{
		/*equilibrium components*/
		/*deviatoric part*/
		fStress3D = fStrain3D;
		fStress3D *= 2.0*muEQ;

		/*volumetric part*/
		fStress3D[0] += kappaEQ*I1;
		fStress3D[1] += kappaEQ*I1;
		fStress3D[2] += kappaEQ*I1;
		
		stress[0] = fStress3D[0];
		stress[1] = fStress3D[1];
		stress[2] = fStress3D[5];
	}
	else  /*nonequilibrium part.  stress is accumulated*/
	{
		double muNEQ = fMu[pindex+1];
		double kappaNEQ = fKappa[pindex+1];

		double taudtS = dt/ftauS[pindex];
		double taudtB = dt/ftauB[pindex];
		double alphaS = exp(-0.5*taudtS);
		double alphaB = exp(-0.5*taudtB);
		double betaB = exp(-taudtB);
		double betaS = exp(-taudtS);

		dSymMatrixT& devQ = fdevQ[pindex];
		dSymMatrixT& devSin = fdevSin[pindex];

		/*deviatoric part*/
		devSin = fStrain3D;
		devSin *= 2.0*muNEQ;
			
		devQ[0] = betaS*devQ_n[0] + alphaS*(devSin[0]-devSin_n[0]);
		devQ[1] = betaS*devQ_n[1] + alphaS*(devSin[1]-devSin_n[1]);
		devQ[2] = betaS*devQ_n[2] + alphaS*(devSin[2]-devSin_n[2]);
		devQ[3] = betaS*devQ_n[3] + alphaS*(devSin[3]-devSin_n[3]);
		devQ[4] = betaS*devQ_n[4] + alphaS*(devSin[4]-devSin_n[4]);
		devQ[5] = betaS*devQ_n[5] + alphaS*(devSin[5]-devSin_n[5]);

		stress[0] += devQ[0];
		stress[1] += devQ[1];
		stress[2] += devQ[5];
		
		dArrayT& meanSin = fmeanSin[pindex];
		dArrayT& meanQ = fmeanQ[pindex];

		/*volumetric part*/
		meanSin[0] = kappaNEQ*I1;
		meanQ[0] = betaB*meanQ_n[0] + alphaB*(meanSin[0]-meanSin_n[0]);
		
		stress[0] += meanQ[0];
		stress[1] += meanQ[0];
		
	}

}

void SS_Visc_Support::SetModulus(dMatrixT& modulus, double dt, int pindex)
{
	/*equilibrium component*/
	if (pindex == -1)
	{
		double mu = fMu[pindex+1];
		double kappa = fKappa[pindex+1];

        /*deviatoric part*/
		modulus = 0.0;
		modulus(0,0) = modulus(1,1) =  modulus(2,2) = 2.0*mu*(1.0 - third);
		modulus(3,3) = modulus(4,4) =  modulus(5,5) = mu;
		modulus(0,1) = modulus(0,2) =  modulus(1,2) = -2.0*mu*third;
		modulus(1,0) = modulus(2,0) =  modulus(2,1) = -2.0*mu*third;

		/*volumetric part*/
		modulus(0,0) += kappa; modulus(1,1) += kappa; modulus(2,2) += kappa;
		modulus(0,1) += kappa; modulus(0,2) += kappa; modulus(1,2) += kappa;
		modulus(1,0) += kappa; modulus(2,0) += kappa; modulus(2,1) += kappa;
	}
	else 		/*non-equilibrium component moduli are accumulated*/
	{
		double taudtS = dt/ftauS[pindex];
		double taudtB = dt/ftauB[pindex];
	
		double alphaS, alphaB;

		alphaS = exp(-0.5*taudtS);
		alphaB = exp(-0.5*taudtB);

		double mu = fMu[pindex+1];
		double kappa = fKappa[pindex+1];

		/*deviatoric part*/
		modulus(0,0) += 2.0*mu*alphaS*(1.0 - third);
		modulus(1,1) += 2.0*mu*alphaS*(1.0 - third);
		modulus(2,2) += 2.0*mu*alphaS*(1.0 - third);	
		modulus(3,3) += mu*alphaS;
		modulus(4,4) += mu*alphaS;
		modulus(5,5) += mu*alphaS;
	
		modulus(0,1) += -2.0*mu*alphaS*third;
		modulus(0,2) += -2.0*mu*alphaS*third;
		modulus(1,2) += -2.0*mu*alphaS*third;
		modulus(1,0) += -2.0*mu*alphaS*third;
		modulus(2,0) += -2.0*mu*alphaS*third;
		modulus(2,1) += -2.0*mu*alphaS*third;
	
	
		/*volumetric part*/
		modulus(0,0) += kappa*alphaB;
		modulus(1,1) += kappa*alphaB;
		modulus(2,2) += kappa*alphaB;

		modulus(0,1) += kappa*alphaB;
		modulus(0,2) += kappa*alphaB;
		modulus(1,2) += kappa*alphaB;
		modulus(1,0) += kappa*alphaB;
		modulus(2,0) += kappa*alphaB;
		modulus(2,1) += kappa*alphaB;
	}
}

void SS_Visc_Support::SetModulus2D(dMatrixT& modulus, double dt, int constraint, int pindex)
{

        double muEQ = fMu[kEquilibrium];
        double kappaEQ = fKappa[kEquilibrium];

        double r = 0.0;
        if (constraint == SolidMaterialT::kPlaneStress)
        {
			double taudtS = dt/ftauS[kNonEquilibrium-1];
			double taudtB = dt/ftauB[kNonEquilibrium-1];

			double alphaS, alphaB;
		
			alphaS = exp(-0.5*taudtS);
			alphaB = exp(-0.5*taudtB);
			
			double muNEQ = fMu[kNonEquilibrium];
			double kappaNEQ = fKappa[kNonEquilibrium];

			double K1 = 3.0*(kappaEQ + kappaNEQ*alphaB) + 4.0*(muEQ + muNEQ*alphaS);
			double K2 = 3.0*(kappaEQ + kappaNEQ*alphaB) - 2.0*(muEQ + muNEQ*alphaS);
			r = K2/K1;
        }

        /*equilibrium component*/
		if (pindex == -1)
		{
			/*deviatoric part*/
			modulus = 0.0;
			modulus(0,0) = modulus(1,1) = 2.0*muEQ*third*(2.0+r);
			modulus(2,2) = muEQ;
			modulus(0,1) = modulus(1,0) = -2.0*muEQ*third*(1.0-r);
	
			/*volumetric part*/
			modulus(0,0) += kappaEQ*(1.0-r); modulus(1,1) += kappaEQ*(1.0-r);
			modulus(0,1) += kappaEQ*(1.0-r); modulus(1,0) += kappaEQ*(1.0-r);
		}
		else
		{
			double taudtS = dt/ftauS[pindex];
			double taudtB = dt/ftauB[pindex];

			double alphaS, alphaB;
		
			alphaS = exp(-0.5*taudtS);
			alphaB = exp(-0.5*taudtB);
			double muNEQ = fMu[pindex+1];
			double kappaNEQ = fKappa[pindex+1];

			/*nonequilibrium component: moduli are accumulated*/
			/*deviatoric part*/
			modulus(0,0) +=  2.0*muNEQ*alphaS*third*(2.0+r);
			modulus(1,1) +=  2.0*muNEQ*alphaS*third*(2.0+r);
			modulus(2,2) +=	muNEQ*alphaS;
	
			modulus(0,1) += -2.0*muNEQ*alphaS*third*(1.0-r);
			modulus(1,0) += -2.0*muNEQ*alphaS*third*(1.0-r);
	
			/*volumetric part*/
			modulus(0,0) += kappaNEQ*alphaB*(1.0-r);
			modulus(1,1) += kappaNEQ*alphaB*(1.0-r);

			modulus(0,1) += kappaNEQ*alphaB*(1.0-r);
			modulus(1,0) += kappaNEQ*alphaB*(1.0-r);
	}
}
