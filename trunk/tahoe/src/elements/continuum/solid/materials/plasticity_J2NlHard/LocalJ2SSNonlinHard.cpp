/* $Id: LocalJ2SSNonlinHard.cpp,v 1.11 2004-07-15 08:29:01 paklein Exp $ */
#include "LocalJ2SSNonlinHard.h"

#include "ifstreamT.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

#include "SSMatSupportT.h"

using namespace Tahoe;

/* parameters */
const int    kNumInternal = 3;
const double sqrt23       = sqrt(2.0/3.0);
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;

/* element output data */
const int    kNumOutput = 6;
static const char* Labels[kNumOutput] = {
	"VMStrss",   // Von Mises stress
        "Prssure",   // pressure
        "EqPStrn",   // equivalent plastic strain
	"dEqPStrn",  // increment of equivalent plastic strain
	"PlsMult",   // plastic multiplier
        "IsoHard"};  // isotropic hardening

/* constructor */
LocalJ2SSNonlinHard::LocalJ2SSNonlinHard(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_J2_local"),
	HookeanMatT  (kNSD),
	fNumIP       (NumIP()),
	fmu          (Mu()),

	/* return values */
	fElasticStrain (kNSD),
	fStress        (kNSD),
	fModulus       (dSymMatrixT::NumValues(kNSD)),
	fModuliCorr    (dSymMatrixT::NumValues(kNSD)),

	/* general workspaces */
	fRelStress (kNSD),
	fsymmatx1  (kNSD),
	fmatx1     (kNSD,kNSD),
	fmatx2     (kNSD,kNSD),
	fmatx3     (kNSD,kNSD),
	ftnsr1     (dSymMatrixT::NumValues(kNSD))
{
        /* obtain hardening coefficients */
        in >> yield >> k1 >> k2 >> k3 >> k4;

	if (yield < 0 )
	{
                cout << "\n GradJ2SSNonlinHard: yield <0" << endl;
		throw ExceptionT::kBadInputValue;
	}
	if (k1 < 0 || k2 < 0 || k3 < 0 || k4 <0)
	{
	        cout << "\n GradJ2SSNonlinHard: bad hardening parameter k1, k2, k3, or k4" << endl;
/*		throw ExceptionT::kBadInputValue;  */
	}
}

/* initialization */
void LocalJ2SSNonlinHard::Initialize(void)
{
	/* void */
	HookeanMatT::Initialize();

	// allocate space for all elements
	AllocateAllElements();
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT LocalJ2SSNonlinHard::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* update internal variables */
void LocalJ2SSNonlinHard::UpdateHistory(void)
{
	ElementCardT& element = CurrentElement();

        /* update plastic variables */
        for (int ip = 0; ip < fNumIP; ip++)
        {
	        LoadData(element, ip);
		
	        /* update state */
	        fStress_n   = fStress;
	        fPlstStrn_n = fPlstStrn;
	        fUnitNorm_n = fUnitNorm;
	        fKineHard_n = fKineHard;
	        fInternal_n = fInternal;
        }
}

/* reset internal variables to last converged solution */
void LocalJ2SSNonlinHard::ResetHistory(void)
{
	ElementCardT& element = CurrentElement();

	/* status flags */
	iArrayT& flags = element.IntegerData();

        for (int ip = 0; ip < fNumIP; ip++)
	{
	        LoadData(element, ip);

	        /* reset state */
	        fStress   = fStress_n;
	        fPlstStrn = fPlstStrn_n;
	        fUnitNorm = fUnitNorm_n;
	        fKineHard = fKineHard_n;
	        fInternal = fInternal_n;
		flags[ip]   = kIsElastic;
        }
}

#if 0
        out << " Hardening coefficients:\n";
	out << "     k1 = " << k1 << "  (Kinematic Hardening          )" << endl;
	out << "     k2 = " << k2 << "  (Isotropic Hardening          )" << endl;
	out << "     k3 = " << k3 << "  (Nonlinear Kinematic Hardening)" << endl;
	out << "     k4 = " << k4 << "  (Nonlinear Isotropic Hardening)" << endl;
#endif

/* modulus */
const dMatrixT& LocalJ2SSNonlinHard::c_ijkl(void)
{
        /* gather element/ip information */
	int fCurrIP = CurrIP();
	ElementCardT& element = CurrentElement();

	LoadData(element, fCurrIP);

	return fModulus;
}

/* stress */
const dSymMatrixT& LocalJ2SSNonlinHard::s_ij(void)
{
        /* gather element/ip information */
	int fCurrIP = CurrIP();
	ElementCardT& element = CurrentElement();

	int iteration = fSSMatSupport->IterationNumber();

	if (fSSMatSupport->RunState() == GlobalT::kFormRHS)
	{
	        if (iteration > -1)
	                /* solve state at current integration point */
		        SolveState(element, fCurrIP);
		else
		{
		        /* load internal variables */
	                LoadData(element, fCurrIP);

		        /* compute elastic stress */
			const dSymMatrixT& e_tot = e();
			const dSymMatrixT& e_els = ElasticStrain(e_tot, element, fCurrIP);
		        HookeanStress(e_els, fStress);

			/* compute elastic moduli */
			fModulus = HookeanMatT::Modulus();
		}
	}

	LoadData(element, fCurrIP);

	return fStress;	
}

/* returns the strain energy density for the specified strain */
double LocalJ2SSNonlinHard::StrainEnergyDensity(void)
{
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, CurrentElement(), CurrIP());
	return HookeanEnergy(e_els);		
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int LocalJ2SSNonlinHard::NumOutputVariables(void) const  { return kNumOutput; }
void LocalJ2SSNonlinHard::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void LocalJ2SSNonlinHard::ComputeOutput(dArrayT& output)
{
        /* gather element/integ point information */
        ElementCardT& element = CurrentElement();
        int ip = CurrIP();

        /* load element data */
        LoadData(element, ip);

	/* deviatoric Von Mises stress */
	fsymmatx1.Deviatoric(fStress);
	double J2 = fsymmatx1.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[0] = sqrt(3.0*J2);

	/* pressure */
	output[1] = fStress.Trace()/3.0;

	/* equivalent plastic strain */
        output[2] = sqrt23*sqrt(fPlstStrn.ScalarProduct());

	/* plastic multiplier */
	output[3] = fInternal[kdelLmbda];

	/* evolution equivalent plastic strain */
	fsymmatx1.SetToScaled(fInternal[kdelLmbda], fUnitNorm_n);
	output[4] = sqrt23*sqrt(fsymmatx1.ScalarProduct());

	/* isotropic hardening */
	output[5] = fInternal[kIsotHard];
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void LocalJ2SSNonlinHard::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

/* returns elastic strain */
const dSymMatrixT& LocalJ2SSNonlinHard::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{	
	/* load internal variables */
	LoadData(element, ip);

	/* compute elastic strain */
	fElasticStrain.DiffOf(totalstrain, fPlstStrn);

	return fElasticStrain;
}	

/* solve for the state at ip */
void LocalJ2SSNonlinHard::SolveState(ElementCardT& element, int fCurrIP)
{
        /* load internal variables */
        LoadData(element, fCurrIP);

	/* status flags */
	iArrayT& flags = element.IntegerData();

	/* step 1. set initial values of plastic strain & internal variables to 
  	           converged values at end of previous time step */
	if (fCurrIP == 0) ResetHistory();

        /* load internal variables */
        LoadData(element, fCurrIP);

	/* step 2. evaluate elastic trial stresses */
        const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, fCurrIP);
        HookeanStress(e_els, fStress);

	/* step 3. initialize unit normal and yield criteria */
	UpdateState();

	if (fInternal[kYieldCrt] > kYieldTol)
	        flags[fCurrIP] = kIsPlastic;
	else 
	        flags[fCurrIP] = kIsElastic;

	/* check for inelastic processes */
	if (flags[fCurrIP] == kIsPlastic)
	{
		/* local Newton iteration */
	       	int max_iteration = 15;
	       	int count = 0;

		/* step 4. zero the increment in plasticity parameter */
		fInternal[kdelLmbda] = 0.;

	       	while (fInternal[kYieldCrt] > kYieldTol && ++count <= max_iteration)
	       	{
		        double varLambda;

		        /* step 5. increment plasticity parameter */
		        IncrementPlasticParameter(varLambda);

			/* step 6. increment stress and state variables */
			IncrementState(varLambda);

			/* step 7. update unit normal and yield criteria */
			UpdateState();
	       	}

	       	/* check for failure */
	       	if (count == max_iteration)
	       	{
	       		cout << "\n LocalJ2SSNonlinHard::SolveState: local iteration failed after " 
	       		     << max_iteration << " iterations" << endl;
	       		throw ExceptionT::kGeneralFail;
	       	}

		/* step 8. compute consistent tangent moduli */
		TangentModuli();
	}
	/* elastic process */
	else if (flags[fCurrIP] == kIsElastic)
	{
	        /* step 4. compute elastic modulus */
	        fModulus = HookeanMatT::Modulus();
	}
	else
	{
	        cout << "\n LocalJ2SSNonlinHard::SolveState: bad flag value " ;
       		throw ExceptionT::kGeneralFail;
       	}
}

/* return a pointer to a new element object constructed with
* the data from element */
void LocalJ2SSNonlinHard::AllocateAllElements(void)
{
	/* determine storage */
	int d_size = 0;
	int dim = dSymMatrixT::NumValues(kNSD);

	d_size += dim;          //fStress
	d_size += dim;          //fStress_n
	d_size += dim;          //fPlstStrn
	d_size += dim;          //fPlstStrn_n
	d_size += dim;          //fUnitNorm
	d_size += dim;          //fUnitNorm_n
	d_size += dim;          //fKineHard
	d_size += dim;          //fKineHard_n
	d_size += kNumInternal; //fInternal
	d_size += kNumInternal; //fInternal_n
	d_size += dim*dim;      //fModulus

	d_size *= fNumIP;

	/* allocate space for all elements */
	for (int el = 0; el < NumElements(); el++)
	{
	        /* get pointer to element el */
		ElementCardT& element = ElementCard(el);

	        /* construct new element */
		element.Dimension(fNumIP, d_size);
	
		/* initialize values */
		element.IntegerData() = kIsElastic;
		element.DoubleData()  = 0.0;
	}
}

double LocalJ2SSNonlinHard::YieldCondition(const dSymMatrixT& relstress,
	double isotropic) const
{
	return sqrt(relstress.ScalarProduct()) - sqrt23*(yield+isotropic);
}

/***********************************************************************
* Private
***********************************************************************/

/* load element data for the specified integration point */
void LocalJ2SSNonlinHard::LoadData(const ElementCardT& element, int fCurrIP)
{
	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT sdim = dSymMatrixT::int2DimensionT(kNSD);
	int dim   = dSymMatrixT::NumValues(kNSD);
	int block = 8*dim + 2*kNumInternal + dim*dim;
	int dex   = fCurrIP*block;

	fStress.Alias     (sdim,         &d_array[dex                ]);
	fStress_n.Alias   (sdim,         &d_array[dex += dim         ]);
	fPlstStrn.Alias   (sdim,         &d_array[dex += dim         ]);
	fPlstStrn_n.Alias (sdim,         &d_array[dex += dim         ]);
	fUnitNorm.Alias   (sdim,         &d_array[dex += dim         ]);
	fUnitNorm_n.Alias (sdim,         &d_array[dex += dim         ]);
	fKineHard.Alias   (sdim,         &d_array[dex += dim         ]);
	fKineHard_n.Alias (sdim,         &d_array[dex += dim         ]);
	fInternal.Alias   (kNumInternal, &d_array[dex += dim         ]);
	fInternal_n.Alias (kNumInternal, &d_array[dex += kNumInternal]);
	fModulus.Alias    (dim,dim,      &d_array[dex += kNumInternal]);
}

/* computes the increment in the plasticity parameter */
void LocalJ2SSNonlinHard::IncrementPlasticParameter(double& varLambda)
{
        /* operations to compute dot product for varLambda */
        fUnitNorm.ToMatrix(fmatx1);
        fUnitNorm_n.ToMatrix(fmatx2);
        fKineHard_n.ToMatrix(fmatx3);
	double cnn = dMatrixT::Dot(fmatx1, fmatx2);
	double cnx = dMatrixT::Dot(fmatx1, fmatx3);

	/* stiffness */
	double dYieldCrt = 2*fmu*cnn + (k1*cnn - k3*cnx)
			     + sqrt23*(sqrt23*k2-k4*fInternal_n[kIsotHard]);

	if (dYieldCrt < kSmall)
	{
		cout << "\n LocalJ2SSNonlinHardT::StressCorrection: consistency function is nonconvex" << endl;
		throw ExceptionT::kGeneralFail;
	}
		
	/* variation of plasticity multiplier */
	varLambda = fInternal[kYieldCrt]/dYieldCrt;

	/* increment of plasticity */
	fInternal[kdelLmbda] += varLambda;
}

/* computes the increments in the stress and internal variables */
void LocalJ2SSNonlinHard::IncrementState(const double& varLambda)
{
	/* increment stress */
	fsymmatx1.SetToScaled(-2.0*fmu*varLambda, fUnitNorm_n);
	fStress += fsymmatx1;

	/* increment kinematic hardening */
	fsymmatx1.SetToScaled(k1*varLambda, fUnitNorm_n);
	fsymmatx1.AddScaled(-1.0*k3*varLambda, fKineHard_n);
	fKineHard += fsymmatx1;

	/* increment isotropic hardening */
	fInternal[kIsotHard] += varLambda*(k2*sqrt23-k4*fInternal_n[kIsotHard]);

	/* increment plastic strain */
	fsymmatx1.SetToScaled(varLambda, fUnitNorm_n);
	fPlstStrn += fsymmatx1;
}

/* computes the unit normal and the yield condition */
void LocalJ2SSNonlinHard::UpdateState()
{
        /* compute relative stress */
	fRelStress.Deviatoric(fStress);
	fRelStress.AddScaled(-1.0, fKineHard);

	/* compute unit normal to yield surface */
	fUnitNorm.SetToScaled(1.0/ sqrt(fRelStress.ScalarProduct()), fRelStress);

	/* compute yield criteria */ 
	fInternal[kYieldCrt] = YieldCondition(fRelStress,fInternal[kIsotHard]);
}

/* computes the consistent tangent moduli */
void LocalJ2SSNonlinHard::TangentModuli()
{
	/* initialize moduli correction */
	fModuliCorr = 0.0;

	/* compute corrections to elastic moduli */
	fUnitNorm.ToMatrix(fmatx1);
	fUnitNorm_n.ToMatrix(fmatx2);
        fKineHard_n.ToMatrix(fmatx3);
	double cnn = dMatrixT::Dot(fmatx1, fmatx2);
	double cnx = dMatrixT::Dot(fmatx1, fmatx3);

	ftnsr1.Outer(fUnitNorm_n,fUnitNorm);
	double h = 2.0*fmu*cnn + (k1*cnn - k3*cnx) + sqrt23*(k2*sqrt23 - k4*fInternal_n[kIsotHard]);
	fModuliCorr.AddScaled(-4*fmu*fmu/h,ftnsr1);

	/* make corrections to elastic moduli */
	fModulus.SumOf(HookeanMatT::Modulus(), fModuliCorr);
}
