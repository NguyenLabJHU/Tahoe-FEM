/* $Id: SSJ2LinHardBaseT.cpp,v 1.11 2004-07-15 08:28:34 paklein Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* H(a) = (1 - ftheta) fH_bar a                                           */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                        */

#include "SSJ2LinHardBaseT.h"

#include <iostream.h>
#include <math.h>

#include "ifstreamT.h"
#include "SSMatSupportT.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
SSJ2LinHardBaseT::SSJ2LinHardBaseT(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("SSJ2LinHardBaseT"),
//	SSSolidMatT(in, support),
	fthird(1.0/3.0),
	fsqrt23(sqrt(2.0/3.0)),
	fplastic(false)
{
ExceptionT::GeneralFail("SSJ2LinHardBaseT::SSJ2LinHardBaseT", "out of date");
#if 0
	/* read parameters */	
	in >> fMu;
	in >> fKappa;
	in >> fYield;	if (fYield <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fH_bar;	if (fH_bar < 0.0) throw ExceptionT::kBadInputValue;
	in >> ftheta;	if (ftheta < 0.0 || ftheta > 1.0) throw ExceptionT::kBadInputValue;
        IsotropicT::Set_mu_kappa(fMu,fKappa);

	/*dimension workspace*/
        int nsd = (NumSD()==2) ? 4 :3;
	int numstress = dSymMatrixT::NumValues(nsd);

	fStressCorr.Dimension(nsd);
	fUnitNorm.Dimension(nsd);
	fRelStress.Dimension(nsd);
	fDevStrain.Dimension(nsd);
	fModuliCorr.Dimension(numstress);
	fTensorTemp.Dimension(numstress);

	/*set internal dofs*/
	fInternalDOF.Dimension(3);
	fInternalDOF[0] = 1;
	fInternalDOF[1] = numstress;
	fInternalDOF[2] = numstress;

        fInternalStressVars.Dimension(2*numstress+1);
        fInternalStrainVars.Dimension(2*numstress+1);

	/*allocates storage for history variables*/
	fnstatev = 0;
    
	/* previous time step*/       
	fnstatev ++;               // alpha: equivalent plastic strain
	fnstatev += numstress;     // beta:  back stress
	fnstatev += numstress;     // plastic strain;
	
	/* current time step*/
	fnstatev ++;               // alpha: equivalent plastic strain
	fnstatev += numstress;     // beta:  back stress
	fnstatev += numstress;     // plastic strain;
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
    
	/*assign pointers to current and preceding blocks of state variable array*/	
	/*current time step*/
	falpha.Set(1,pstatev);
	pstatev ++;
	fBeta.Set(nsd,pstatev);
	pstatev += numstress;
	fPlasticStrain.Set(nsd,pstatev);
	pstatev += numstress;
	
	/*previous time step*/
	falpha_n.Set(1,pstatev);
	pstatev ++;
	fBeta_n.Set(nsd,pstatev);
	pstatev += numstress;
	fPlasticStrain_n.Set(nsd,pstatev);
#endif
}

/* output parameters to stream */
void SSJ2LinHardBaseT::Print(ostream& out) const
{
    out << " Shear modulus . . . . . . . . . . . . . . . . . = " << fMu << '\n'; 
    out << " Bulk modulus. . . . . . . . . . . . . . . . . . = " << fKappa << '\n'; 
	out << " Initial yield stress. . . . . . . . . . . . . . = " << fYield << '\n';
	out << " Hardening parameter . . . . . . . . . . . . . . = " << fH_bar << '\n';
	out << " Isotropic/kinematic mixity. . . . . . . . . . . = " << ftheta << '\n';
}

void SSJ2LinHardBaseT::PrintName(ostream& out) const
{
	/* inherited */
	out << "    J2 Isotropic/Kinematic\n";
	out << "    Hardening with Radial Return\n";
	out << "    Small Strain\n";
}

void SSJ2LinHardBaseT::PointInitialize(void)
{
	/* allocate element storage */
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
	/* initialize internal variables to 0.0*/
		element.DoubleData() = 0.0;
	}

	/* store results as last converged */
	if (CurrIP() == NumIP() - 1) 
		UpdateHistory();
}
 
void SSJ2LinHardBaseT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		falpha_n = falpha;
		fBeta_n = fBeta;
		fPlasticStrain_n = fPlasticStrain;
		
		/* write to storage */
		Store(element, ip);
	}
}

void SSJ2LinHardBaseT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		falpha = falpha_n;
		fBeta = fBeta_n;
		fPlasticStrain = fPlasticStrain_n;
		
		/* write to storage */
		Store(element, ip);
	}
}

void SSJ2LinHardBaseT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}

void SSJ2LinHardBaseT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;
}

/**************************protected************************************/
/* returns the value value of the yield function given the
* relative stress vector and state variables, where  alpha
* represents isotropic hardening.  NOTE: the relative stress
* should already contain the correction for any kinematic
* hardening. */
double SSJ2LinHardBaseT::YieldCondition(const dSymMatrixT& relstress, double alpha) const
{
	return sqrt(relstress.ScalarProduct()) - fsqrt23*K(alpha);
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface. Given the trialstrain, evaluates ftrial, fRelStress,
* fUnitNorm using fBeta_n and falpha_n*/
bool SSJ2LinHardBaseT::PlasticLoading(const dSymMatrixT& devtrialstress)
{	
  /*evaluate trial stress*/
  fRelStress = devtrialstress;
  fRelStress -= fBeta_n;
  
  ftrial = YieldCondition(fRelStress, falpha_n[0]);
  /* plastic */
  if (ftrial > kSmall)
  {		
    fdgamma = ftrial/(2.0*fMu + 2.0*fH_bar/3.0);
    fUnitNorm = fRelStress;
    double norm = sqrt(fRelStress.ScalarProduct());
    fUnitNorm /= norm;
    return true;
  }
  /* elastic */
  else return false;
}	

/* return the correction to stress vector computed by the mapping the
* stress back to the yield surface, if needed */
const dSymMatrixT& SSJ2LinHardBaseT::StressCorrection(const dSymMatrixT& devtrialstress)
{
	/* check consistency and initialize plastic element */
	/*note PlasticLoading evaluates ftrial, fUnitNorm*/
	
  fplastic = PlasticLoading(devtrialstress);
  if (fplastic)
  {
    if (fSSMatSupport->RunState() == GlobalT::kFormRHS)
    {

      /* plastic increment */
      falpha[0] = falpha_n[0] + fsqrt23*fdgamma;

      fBeta = fUnitNorm;
      fBeta *= fsqrt23*(H(falpha[0]) - H(falpha_n[0]));
      fBeta += fBeta_n;

      fPlasticStrain = fUnitNorm;
      fPlasticStrain *= fdgamma;
      fPlasticStrain += fPlasticStrain_n;
    }
    /* plastic increment stress correction */
    fStressCorr = fUnitNorm;
    fStressCorr *= -2.0*fMu*fdgamma;
  }
  else fStressCorr = 0.0;
  
  return fStressCorr;
}	

/* return the correction to moduli due to plasticity (if any)
*
* Note: Return mapping occurs during the call to StressCorrection.
*       The element passed in is already assumed to carry current
*       internal variable values */
const dMatrixT& SSJ2LinHardBaseT::ModuliCorrection(void)
{
	/* initialize */
        fModuliCorr = 0.0;

	if (fplastic)
	{
	  /* compute constants */
	  double norm = sqrt(fRelStress.ScalarProduct());
	  double thetahat = 2.0*fMu*fdgamma/norm;
	  double hardmod = (dK()+dH())/(3.0*fMu);
	  double thetabar = (1.0 / (1.0 + hardmod )) - thetahat;
	
	  int nsd = (NumSD()==2) ? 4 : 3;

	  /* moduli corrections */
	  fTensorTemp = 0.0;

	  if (nsd == 4){
	    fModuliCorr(0,0)=fModuliCorr(1,1)=fModuliCorr(3,3)=2.0*fthird;
	    fModuliCorr(0,1)=fModuliCorr(0,3)= -fthird;
	    fModuliCorr(1,0)=fModuliCorr(1,3)= -fthird;
	    fModuliCorr(3,0)=fModuliCorr(3,1)= -fthird;
	    fModuliCorr(2,2)=0.5;
	  }
	  else {
	    fModuliCorr(0,0)=fModuliCorr(1,1)=fModuliCorr(2,2)=2.0*fthird;
	    fModuliCorr(0,1)=fModuliCorr(0,2)= -fthird;
	    fModuliCorr(1,0)=fModuliCorr(1,2)= -fthird;
	    fModuliCorr(2,0)=fModuliCorr(2,1)= -fthird;
	    fModuliCorr(3,3)=fModuliCorr(4,4)=fModuliCorr(5,5)=0.5;
	  }
	  fModuliCorr *= 2.0*fMu*(-thetahat);
	  
	  fTensorTemp.Outer(fUnitNorm,fUnitNorm);
	  fTensorTemp *= 2.0*fMu*thetabar;
	  fModuliCorr -= fTensorTemp;
	}

	return fModuliCorr;
}	

