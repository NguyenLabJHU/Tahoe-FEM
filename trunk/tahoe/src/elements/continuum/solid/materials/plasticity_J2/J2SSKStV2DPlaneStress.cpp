/* $Id: J2SSKStV2DPlaneStress.cpp,v 1.1 2006-07-21 20:03:19 tdnguye Exp $ */
/* created: paklein (06/18/1997) */
#include "J2SSKStV2DPlaneStress.h"
#include "SSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

const int kNumOutput = 3;
static const char* Labels[kNumOutput] = {
	"alpha",  // equivalent plastic strain
	   "VM",  // Von Mises stress
	"press"}; // pressure

const double kYieldTol    = 1.0e-10;
const double third		  = 1.0/3.0;
const double sqrt23       = sqrt(2.0/3.0);

const int    kNumInternal = 4; // number of internal variables
const int kNSD = 2;

/* constructor */
J2SSKStV2DPlaneStress::J2SSKStV2DPlaneStress(void):
	ParameterInterfaceT("small_strain_StVenant_J2_PlaneStress"),
	HookeanMatT(2),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fPlasticStrain2D(2),
	fRelStress2D(2),
	fTrialRelStress2D(2),
	fBeta2D(2),
	fElasticStrain2D(2),
	fVec1(3),
	fVec2(3),
	fModuliCorr2D(dSymMatrixT::NumValues(2)),
	fP(dSymMatrixT::NumValues(2)),
	fModTemp(dSymMatrixT::NumValues(2)),
	fModTheta(dSymMatrixT::NumValues(2))
	{

	/* reset default value */
	fConstraint = kPlaneStress;

	/* acccount for thickness */
//	fDensity *= fThickness;
	
	/*Define Projection Matrix*/
	fP = 0.0;
	fP(0,0) = fP(1,1) = 2.0*third;
	fP(1,0) = fP(0,1) = -third;
	fP(2,2) = 2.0;

	fVec1 = 0.0;
	fVec2 = 0.0;
	fModTemp = 0.0;
	fModTheta = 0.0;
}

double J2SSKStV2DPlaneStress::YieldCondition(const dSymMatrixT& relstress,
	double alpha) const
{
	double meanstress = third*(relstress[0]+relstress[1]);
	double s0 = relstress[0]-meanstress;
	double s1 = relstress[1]-meanstress;
	double s2 = -meanstress;
	double s01 = relstress[2];
	
	return sqrt(s0*s0+s1*s1+s2*s2+2.0*s01*s01) - sqrt23*K(alpha);
}

/* returns elastic strain (3D) */
const dSymMatrixT& J2SSKStV2DPlaneStress::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{
	/* remove plastic strain */
	if (element.IsAllocated())
	{
		/* load internal variables */
		LoadData(element, ip);

		/* compute elastic strain */
		fElasticStrain2D.DiffOf(totalstrain, fPlasticStrain2D);
	
		return fElasticStrain2D;
	}	
	/* no plastic strain */
	else	
		return totalstrain;
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface */
int J2SSKStV2DPlaneStress::PlasticLoading(const dSymMatrixT& trialstrain,
	ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return( YieldCondition(RelativeStress(trialstrain, element), 0.0)
					> kYieldTol );
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
			
		/* load internal variables */
		LoadData(element, ip);
		
		fInternal[kftrial] = YieldCondition(RelativeStress(trialstrain,element),
						    fInternal[kalpha]);

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{		
		
			/* set flag */
			Flags[ip] = kIsPlastic;
	
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
			Flags[ip] = kIsElastic;
			
			return 0;
		}
	}
}	

/* return the correction to stress vector computed by the mapping the
* stress back to the yield surface, if needed */
const dSymMatrixT& J2SSKStV2DPlaneStress::StressCorrection(const dSymMatrixT& trialstrain, 
	ElementCardT& element, int ip)
{
	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, element, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element);
					
		/* initialize element data */
		PlasticLoading(trialstrain, element, ip);
	}

	/* initialize */
	fModTheta = HookeanMatT::Modulus();
	const dSymMatrixT& trialstress = RelativeStress(trialstrain, element);
	
	if (element.IsAllocated())
	{		
		LoadData(element,ip);

		/* fetch data */
		double  ftrial = fInternal[kftrial];
		double& dgamma = fInternal[kdgamma];
		double& fbar2  = fInternal[kstressnorm];
		double alpha   = fInternal[kalpha];

		/* return mapping (single step) */
		if (ftrial > kYieldTol)
		{
			double smean_tr = trialstress[0] + trialstress[1];
			double sdiff_tr = trialstress[0] - trialstress[1];
			double s12_tr = trialstress[2];
		
			double E = IsotropicT::Young();
			double nu = IsotropicT::Poisson();
			double mu = IsotropicT::Mu();
					
			dgamma = 0.0;
			
			double d1 = (1.0 +third*E/(1.0-nu)*dgamma);
			double d2 = (1.0 + 2.0*mu*dgamma);
			fbar2 = 0.5*third*smean_tr*smean_tr/(d1*d1) 
					+ (0.5*sdiff_tr*sdiff_tr+2.0*s12_tr*s12_tr)/(d2*d2);

			double alphan1 = alpha+sqrt23*sqrt(fbar2)*dgamma;
			double R2 = third*K(alphan1)*K(alphan1);					
				
			double f2 = 0.5*fbar2-R2;
			int count = 0;
			while (fabs(f2) > kYieldTol && ++count <= fmax_iteration)
			{
				/* stiffness */
				double dfbar2 = -third*third*smean_tr*smean_tr * E/((1.0-nu)*d1*d1*d1)
					-4.0*mu*(0.5*sdiff_tr*sdiff_tr+2.0*s12_tr*s12_tr)/(d2*d2*d2);
				double dR2 = 2.0*third*K(alphan1)*dK(alphan1)*sqrt23/sqrt(fbar2)*(fbar2+0.5*dgamma*dfbar2);
				double df2 = 0.5*dfbar2 - dR2;
				
				if (-df2 < kSmall)
				{
					cout << "\nJ2SSKStV2DPlaneStressPlaneStressT::StressCorrection: consistency function is nonconvex"<< endl;
					throw ExceptionT::kGeneralFail;
				}
				
					/* increment update */
					dgamma -= f2/df2;
				
					d1 = (1.0 +third*E/(1.0-nu)*dgamma);
					d2 = (1.0 + 2.0*mu*dgamma);
					fbar2 = 0.5*third*smean_tr*smean_tr/(d1*d1) 
						+ (0.5*sdiff_tr*sdiff_tr+2.0*s12_tr*s12_tr)/(d2*d2);

					alphan1 = alpha+sqrt23*sqrt(fbar2)*dgamma;
					R2 = third*K(alphan1)*K(alphan1);					

					f2 = 0.5*fbar2 - R2;
				}
			/* check for failure */
			if (count == fmax_iteration)
			{
				cout << "\n J2SSC0HardeningT::StressCorrection: local iteration failed after " 
					 << fmax_iteration << " iterations" << endl;
				throw ExceptionT::kGeneralFail;
			}
			fModTemp = fP;
			fModTemp *= dgamma;

			fModTheta.Inverse();
			fModTheta += fModTemp;
			fModTheta.Inverse();
			
			fModTemp = HookeanMatT::Modulus();
			fModTemp.Inverse();
			
			/*update stress*/	
//			s = Theta C^-1 strial
			fModTemp.MultTx(trialstress,fVec1);
			fModTheta.MultTx(fVec1,fRelStress2D);
		}
		else 
		{
			dgamma = 0.0;			
			fRelStress2D = trialstress;
		}
		return fRelStress2D;
	}
	else 
	{
		fTrialRelStress2D = trialstress;
		return fTrialRelStress2D;
	}
}	

/* return the correction to moduli due to plasticity (if any)
*
* Note: Return mapping occurs during the call to StressCorrection.
*       The element passed in is already assumed to carry current
*       internal variable values */
const dMatrixT& J2SSKStV2DPlaneStress::ModuliCorrection(const ElementCardT& element,
	int ip)
{
	/* initialize */
	fModTheta = HookeanMatT::Modulus();
	fModuliCorr2D = 0.0;
	if (element.IsAllocated() && (element.IntegerData())[ip] == kIsPlastic)
	{
		/* load internal variables */
		LoadData(element,ip);
			
		/* compute constants */
		double alpha    = fInternal[kalpha];
		double dgamma   = fInternal[kdgamma];
		const dSymMatrixT& rel = fRelStress2D;
		
		double theta1 = 1.0 + 2.0*third*H(alpha)*dgamma;
		double theta2 = 1.0 - 2.0*third*dK(alpha)*dgamma;

		/*norm of deviatoric part of the relative stress*/
		double norm = 2.0*third*(rel[0]*rel[0] -rel[0]*rel[1]+rel[1]*rel[1]+3.0*rel[2]*rel[2]);
		double thetabar = 2.0*third*theta1/theta2*(dK(alpha)*theta1 + H(alpha)*theta2)*norm;

		/*calculate Theta*/
		fModTemp = fP;
		fModTemp *= dgamma;
		fModTheta.Inverse();
		fModTheta += fModTemp;
		fModTheta.Inverse();

		/*P.rel*/
		fVec2[0] = third*(2.0*rel[0]-rel[1]);
		fVec2[1] = third*(2.0*rel[1]-rel[0]);
		fVec2[2] = 2.0*rel[2];	

		fModTheta.MultTx(fVec2,fVec1);

		double denom = (fVec2[0]*fVec1[0]+fVec2[1]*fVec1[1]+fVec2[2]*fVec1[2])+thetabar;
		
		fModuliCorr2D.Outer(fVec1,fVec1);
		fModuliCorr2D /= -denom;
	}
	return fModuliCorr2D;
}

/* return a pointer to a new plastic element object constructed with
* the data from element */
void J2SSKStV2DPlaneStress::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int ip = NumIP();
	int i_size = 0;
	i_size += ip; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*ip; //fPlasticStrain2D
	d_size += dSymMatrixT::NumValues(kNSD)*ip; //fRelStress2D
	d_size += dSymMatrixT::NumValues(kNSD)*ip; //fBeta2D
	d_size += kNumInternal*ip;          //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;
}

void J2SSKStV2DPlaneStress::LoadData(const ElementCardT& element, int ip)
{
	/* check */
	if (!element.IsAllocated()) { 
		cout << "\nElement is not allocated: ";	
		throw ExceptionT::kGeneralFail;
	}

	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*NumIP();
	int dex       = ip*stressdim;
	
	fPlasticStrain2D.Alias(         dim, &d_array[           dex]);
	     fRelStress2D.Alias(         dim, &d_array[  offset + dex]);
	         fBeta2D.Alias(         dim, &d_array[2*offset + dex]);
	     fInternal.Alias(kNumInternal, &d_array[3*offset + ip*kNumInternal]);     	
}

/* update internal variables */
void J2SSKStV2DPlaneStress::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void J2SSKStV2DPlaneStress::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element, NumIP());
}

/* element level data */
void J2SSKStV2DPlaneStress::Update(ElementCardT& element)
{
	/* get flags */
	iArrayT& flags = element.IntegerData();

	/* check if reset state (is same for all ip) */
	if (flags[0] == kReset)
	{
		flags = kIsElastic; //don't update again
		return;
	}

	/* update plastic variables */
	for (int ip = 0; ip < NumIP(); ip++)
	{
		if (flags[ip] ==  kIsPlastic) /* plastic update */
		{
			/* do not repeat if called again */
			flags[ip] = kReset;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* fetch element data */
			LoadData(element, ip);
		
			const double& dgamma = fInternal[kdgamma];
			const double& fbar2 = fInternal[kstressnorm];

			/* plastic increment */
			fInternal[kalpha] += sqrt23*sqrt(fbar2)*dgamma;

			fP.MultTx(fRelStress2D,fVec2);
			fVec2[2] *= 0.5;
			fPlasticStrain2D.AddScaled(dgamma, fVec2);	
		}
	}
}


/* moduli */
const dMatrixT& J2SSKStV2DPlaneStress::c_ijkl(void)
{
	fModulus2D =  ModuliCorrection(CurrentElement(), CurrIP());
	fModulus2D += fModTheta;
	return fModulus2D;
}

/* stress */
const dSymMatrixT& J2SSKStV2DPlaneStress::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);	
	fStress2D = StressCorrection(e_els, element, ip);
	fStress2D += fBeta2D;
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double J2SSKStV2DPlaneStress::StrainEnergyDensity(void)
{
        fElasticStrain2D = e();
	const dSymMatrixT& strain = e();
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated() && (element.IntegerData())[CurrIP()] == kIsPlastic) 
	{
		/* load internal variables */

		LoadData(element, CurrIP());
	        fElasticStrain2D -= fPlasticStrain2D;
                if (element.IntegerData()[CurrIP()] == kIsPlastic)
		{
		  const double& dgamma = fInternal[kdgamma];
		  fP.MultTx(fRelStress2D,fVec2);
		  fVec2[2] *= 0.5;

		  fElasticStrain2D[0] -= dgamma*fVec2[0];
		  fElasticStrain2D[1] -= dgamma*fVec2[1];
		  fElasticStrain2D[2] -= dgamma*fVec2[2];
		}
	}	

	double energy = HookeanEnergy(fElasticStrain2D);
	return energy;
}

/***********************************************************************
* Protected
***********************************************************************/


/* set (material) tangent modulus */
void J2SSKStV2DPlaneStress::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli2D(modulus, fConstraint);
//	modulus *= fThickness;
}


dSymMatrixT& J2SSKStV2DPlaneStress::RelativeStress(const dSymMatrixT& trialstrain, const ElementCardT& element)
{
	HookeanStress(trialstrain,fTrialRelStress2D);
	/* element is plastic - element data loaded */
	if (element.IsAllocated())
	{
		/* kinematic hardening */
		fTrialRelStress2D -= fBeta2D;
	}
	return fTrialRelStress2D;
}

int J2SSKStV2DPlaneStress::NumOutputVariables(void) const  { return kNumOutput; }
void J2SSKStV2DPlaneStress::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

const iArrayT& J2SSKStV2DPlaneStress::InternalDOF(void) const
{
  return(fInternalDOF);
}

const dArrayT& J2SSKStV2DPlaneStress::InternalStrainVars(void)
{
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated() )
	{
		/* load internal variables */
		LoadData(element, CurrIP());

		double* p = fInternalStrainVars.Pointer();

		if (element.IntegerData()[CurrIP()] == kIsPlastic)
		{
		  const double& dgamma = fInternal[kdgamma];
		  const double& fbar2 = fInternal[kstressnorm]; 
		  fP.MultTx(fRelStress2D,fVec2);
		  fVec2[2] *= 0.5;

		  *p++ = fInternal[kalpha]+sqrt23*sqrt(fbar2)*dgamma;

		  *p++ = fPlasticStrain2D[0]+dgamma*fVec2[0];
		  *p++ = fPlasticStrain2D[1]+dgamma*fVec2[1];
		  *p++ = -fPlasticStrain2D[0]-fPlasticStrain2D[1]-dgamma*fVec2[0]-dgamma*fVec2[1];
		  *p++ = 0.0;
		  *p++ = 0.0;
		  *p++ = fPlasticStrain2D[2]+dgamma*fVec2[2];

		  *p++ = fPlasticStrain2D[0]+dgamma*fVec2[0];
		  *p++ = fPlasticStrain2D[1]+dgamma*fVec2[1];
		  *p++ = -fPlasticStrain2D[0]-fPlasticStrain2D[1]-dgamma*fVec2[0]-dgamma*fVec2[1];
		  *p++ = 0.0;
		  *p++ = 0.0;
		  *p++ = fPlasticStrain2D[2]+dgamma*fVec2[2];
	        }
		else {
		  *p++ = fInternal[kalpha];
		  
		  *p++ = fPlasticStrain2D[0];
		  *p++ = fPlasticStrain2D[1];
		  *p++ = -fPlasticStrain2D[0]-fPlasticStrain2D[1];
		  *p++ = 0.0;
		  *p++ = 0.0;
		  *p++ = fPlasticStrain2D[2];

		  *p++ = fPlasticStrain2D[0];
		  *p++ = fPlasticStrain2D[1];
		  *p++ = -fPlasticStrain2D[0]-fPlasticStrain2D[1];
		  *p++ = 0.0;
		  *p++ = 0.0;
		  *p++ = fPlasticStrain2D[2];

		}
	}
	else fInternalStrainVars = 0.0;

	return(fInternalStrainVars);
}

const dArrayT& J2SSKStV2DPlaneStress::InternalStressVars(void)
{
	ElementCardT& element = CurrentElement();
	fElasticStrain2D = e();

	double* p = fInternalStressVars.Pointer();
	if (element.IsAllocated())
	{
		/* load internal variables */
		LoadData(element, CurrIP());

		const double& dgamma = fInternal[kdgamma];
		const double& fbar2 = fInternal[kstressnorm];
		
		if (element.IntegerData()[CurrIP()] == kIsPlastic) 
		  *p++ = -K(fInternal[kalpha]+sqrt23*sqrt(fbar2)*dgamma);
		else
		  *p++ = -K(fInternal[kalpha]);

		*p++ = -fBeta2D[0];
		*p++ = -fBeta2D[1];
		*p++ = 0.0;
		*p++ = 0.0;
		*p++ = 0.0;
		*p++ = -fBeta2D[2];

		fElasticStrain2D -= fPlasticStrain2D;
		if (element.IntegerData()[CurrIP()] == kIsPlastic) 
		{
		  const double& dgamma = fInternal[kdgamma];
		  const double& fbar2 = fInternal[kstressnorm]; 
		  fP.MultTx(fRelStress2D,fVec2);
		  fVec2[2] *= 0.5;

		  fElasticStrain2D[0] -= dgamma*fVec2[0];
		  fElasticStrain2D[1] -= dgamma*fVec2[1];
		  fElasticStrain2D[2] -= dgamma*fVec2[2];
		}
	}
	else {
		*p++ = 0.0;

		*p++ = 0.0;
		*p++ = 0.0;
		*p++ = 0.0;
		*p++ = 0.0;
		*p++ = 0.0;
		*p++ = 0.0;
	}
	
	HookeanStress(fElasticStrain2D,fStress2D);
 	*p++ = fStress2D[0];
	*p++ = fStress2D[1];
	*p++ = 0.0;
	*p++ = 0.0;
	*p++ = 0.0;
	*p++ = fStress2D[2];

	return(fInternalStressVars);
}

void J2SSKStV2DPlaneStress::ComputeOutput(dArrayT& output)
{
	/* pressure */
	double meanstress = fStress2D.Trace()/3.0;
	output[2] = meanstress;
	
	/* deviatoric Von Mises stress */
	fStress2D[0] -= meanstress;
	fStress2D[1] -= meanstress;
	double J2 = fStress2D.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

	output[0] = 0.0;
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		/* load internal variables */
		LoadData(element, CurrIP());

		/* plastic strain */
		output[0] = fInternal[kalpha];
		
		const double& dgamma = fInternal[kdgamma];
		const double& fbar2 = fInternal[kstressnorm];
		/* status flags */
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
		{
			output[0] += sqrt23*sqrt(fbar2)*dgamma;
		}
	}
}

/* describe the parameters needed by the interface */
void J2SSKStV2DPlaneStress::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SSSolidMatT::DefineParameters(list);
	IsotropicT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
	J2SSC0HardeningT::DefineParameters(list);

	/* max number of iterations */
	ParameterT max_its(fmax_iteration, "max_iteration");
	max_its.SetDefault(fmax_iteration);
	list.AddParameter(max_its);
}

/* information about subordinate parameter lists */
void J2SSKStV2DPlaneStress::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	J2SSC0HardeningT::DefineSubs(sub_list);
}

/* return the description of the given inline subordinate parameter list */
void J2SSKStV2DPlaneStress::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	SSSolidMatT::DefineInlineSub(name, order, sub_lists);
	IsotropicT::DefineInlineSub(name, order, sub_lists);
	HookeanMatT::DefineInlineSub(name, order, sub_lists);
	J2SSC0HardeningT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* J2SSKStV2DPlaneStress::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = NULL;

	/* try each base class */
	sub = SSSolidMatT::NewSub(name);
	if (sub) return sub;

	sub = IsotropicT::NewSub(name);
	if (sub) return sub;

	sub = HookeanMatT::NewSub(name);
	if (sub) return sub;
	
	return J2SSC0HardeningT::NewSub(name);
}

/* accept parameter list */
void J2SSKStV2DPlaneStress::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
	J2SSC0HardeningT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);

	/* number elastic iterations */
	fmax_iteration = list.GetParameter("max_iteration");
}
