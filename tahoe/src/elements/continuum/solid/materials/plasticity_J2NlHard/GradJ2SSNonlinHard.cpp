#include "GradJ2SSNonlinHard.h"

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "ifstreamT.h"
#include "ContinuumElementT.h"

using namespace Tahoe;

/* parameters */
const int    kNumInternal = 4;
const double sqrt23       = sqrt(2.0/3.0);
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;

/* element output data */
const int    kNumOutput = 5;
static const char* Labels[kNumOutput] = {
        "EqPStrn",  // equivalent plastic strain
	"VMStrss",  // Von Mises stress
        "Prssure",  // pressure
        "IsoHard",  // isotropic hardening
        "NlIsoHard"}; // nonlocal isotropic hardening

/* constructor */
GradJ2SSNonlinHard::GradJ2SSNonlinHard(ifstreamT& in, const SmallStrainT& element):
	SSStructMatT (in, element),
	IsotropicT   (in),
	HookeanMatT  (kNSD),
	fStatus      (ContinuumElement().RunState()),
        fNumIP       (NumIP()),
        fmu          (Mu()),
	fNumNodes    (ContinuumElement().InitialCoordinates().NumberOfNodes()),

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
        in >> yield >> k1 >> k2 >> k3 >> k4  >> c1 >> c2;

	if (yield < 0 )
	{
                cout << "\n GradJ2SSNonlinHard: yield <0" << endl;
		throw eBadInputValue;
	}
	if (k1 < 0 || k2 < 0 || k3 < 0 || k4 <0)
	{
	        cout << "\n GradJ2SSNonlinHard: bad hardening parameter k1, k2, k3, or k4" << endl;
		throw eBadInputValue;
	}
}

/* initialization */
void GradJ2SSNonlinHard::Initialize(void)
{
	/* void */
	HookeanMatT::Initialize();

	// allocate space for all elements
	AllocateAllElements();
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT GradJ2SSNonlinHard::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* update internal variables */
void GradJ2SSNonlinHard::UpdateHistory(void)
{
	ElementCardT& element = CurrentElement();

        /* update plastic variables */
        for (int ip = 0; ip < fNumIP; ip++)
        {
	        LoadData(element, ip);
		
	        /* update state */
	        fStress_n     = fStress;
	        fPlstStrn_n   = fPlstStrn;
	        fUnitNorm_n   = fUnitNorm;
	        fKineHard_n   = fKineHard;
	        fNlKineHard_n = fNlKineHard;
	        fInternal_n   = fInternal;
        }
}

/* reset internal variables to last converged solution */
void GradJ2SSNonlinHard::ResetHistory(void)
{
	ElementCardT& element = CurrentElement();

	/* status flags */
	iArrayT& flags = element.IntegerData();

        for (int ip = 0; ip < fNumIP; ip++)
	{
	        LoadData(element, ip);

	        /* reset state */
	        fStress     = fStress_n;
	        fPlstStrn   = fPlstStrn_n;
	        fUnitNorm   = fUnitNorm_n;
	        fKineHard   = fKineHard_n;
		fNlKineHard = fNlKineHard_n;
	        fInternal   = fInternal_n;
		flags[ip]   = kIsElastic;
        }
}

/* print parameters */
void GradJ2SSNonlinHard::Print(ostream& out) const
{
	/* inherited */
	SSStructMatT::Print(out);
	IsotropicT::Print(out);

        /* hardening coefficients */
        out << " Hardening coefficients:\n";
	out << "     k1 = " << k1 << "  (Kinematic Hardening          )" << endl;
	out << "     k2 = " << k2 << "  (Isotropic Hardening          )" << endl;
	out << "     k3 = " << k3 << "  (Nonlinear Kinematic Hardening)" << endl;
	out << "     k4 = " << k4 << "  (Nonlinear Isotropic Hardening)" << endl;

	/* nonlocal coefficients */
	out << " Nonlocal coefficients:\n";
	out << "     c1 = " << c1 << "  (Nonlocal Kinematic Hardening )" << endl;
	out << "     c2 = " << c2 << "  (Nonlocal Isotropic Hardening )" << endl;
}

/* print name */
void GradJ2SSNonlinHard::PrintName(ostream& out) const
{
	/* inherited */
	SSStructMatT::PrintName(out);
	out << "    Nonlocal small strain J2 plasticity\n";
	out << "      with nonlinear gradient dependent\n";
        out << "      isotropic/kinematic hardening\n";
        out << "    Semi-implicit backward Euler\n";
	out << "      integrating scheme\n";
}

/* modulus */
const dMatrixT& GradJ2SSNonlinHard::c_ijkl(void)
{
        /* gather element/ip information */
	int fCurrIP = CurrIP();
	ElementCardT& element = CurrentElement();

	LoadData(element, fCurrIP);

	return fModulus;
}

/* stress */
const dSymMatrixT& GradJ2SSNonlinHard::s_ij(void)
{
        /* gather element/ip information */
	int fCurrIP = CurrIP();
	ElementCardT& element = CurrentElement();

	int iteration = ContinuumElement().IterationNumber();

	if (fStatus == GlobalT::kFormRHS && fCurrIP == 0)
	{
	        if (iteration > -1)
		        /* solve state at each integration point (all at once) */
		        SolveState(element);
		else
		{
		        for (int ip = 0; ip < fNumIP; ip++)
			{
			        /* load internal variables */
			        LoadData(element, ip);

				/* compute elastic stress */
				const dSymMatrixT& e_tot = e(ip);
				const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);
				HookeanStress(e_els, fStress);

				/* compute elastic moduli */
				fModulus = HookeanMatT::Modulus();
			}
		}
	}

	LoadData(element, fCurrIP);

	return fStress;	
}

/* returns the strain energy density for the specified strain */
double GradJ2SSNonlinHard::StrainEnergyDensity(void)
{
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, CurrentElement(), CurrIP());
	return HookeanEnergy(e_els);		
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int GradJ2SSNonlinHard::NumOutputVariables(void) const  { return kNumOutput; }
void GradJ2SSNonlinHard::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Allocate(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void GradJ2SSNonlinHard::ComputeOutput(dArrayT& output)
{
        /* gather element/integ point information */
        ElementCardT& element = CurrentElement();
        int ip = CurrIP();

        /* load element data */
        LoadData(element, ip);

	/* pressure */
	output[2] = fStress.Trace()/3.0;

	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

	/* equivalent plastic strain */
        output[0] = sqrt23*sqrt(fPlstStrn.ScalarProduct());

	/* isotropic hardening */
	output[3] = fInternal[kIsotHard];

	output[4] = fInternal[kNlIsotHard];
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void GradJ2SSNonlinHard::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

/* returns elastic strain */
const dSymMatrixT& GradJ2SSNonlinHard::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{	
	/* load internal variables */
	LoadData(element, ip);

	/* compute elastic strain */
	fElasticStrain.DiffOf(totalstrain, fPlstStrn);

	return fElasticStrain;
}	

/* solve for the state at ip */
void GradJ2SSNonlinHard::SolveState(ElementCardT& element)
{
	/* status flags */
	iArrayT& flags = element.IntegerData();

	/* step 1. set initial values of plastic strain & internal variables to 
  	           converged values at end of previous time step */
	ResetHistory();

        for (int ip = 0; ip < fNumIP; ip ++)
	{
	        /* load internal variables */
	        LoadData(element, ip);

		/* step 2. evaluate elastic trial stresses */
		const dSymMatrixT& e_tot = e(ip);
		const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);
		HookeanStress(e_els, fStress);

		/* step 3. initialize unit normal and yield criteria */
		UpdateState();

		if (fInternal[kYieldCrt] > kYieldTol)
		        flags[ip] = kIsPlastic;
	}

        /* check for inelastic processes at any ip in element */
        bool Converged = CheckElementState(element);

	if (!Converged)
	{
	        /* local Newton iteration */
	        int max_iteration = 30;
	        int count = 0;

	        for (int ip = 0; ip < fNumIP; ip ++)
		{
		        /* load internal variables */
		        LoadData(element, ip);

		        /* step 4. zero the increment in plasticity parameter */
		        fInternal[kdelLmbda] = 0.;
		}

		while (!Converged && ++count <= max_iteration)
		{
			/*array of IsotHard at all ip in current element used to compute Laplacian*/
			dArrayT  ip_IsotHard(fNumIP);
			ip_IsotHard = 0.;

		        for (int ip = 0; ip < fNumIP; ip ++)
			{
			        /* load internal variables */
			        LoadData(element, ip);

				/* check for inelastic processes */
				if (flags[ip] == kIsPlastic && fInternal[kYieldCrt] > kYieldTol)
				{
				        double varLambda;

					/* step 5. increment plasticity parameter */
					IncrementPlasticParameter(varLambda);

					/* step 6. increment stress and state variables */
					IncrementState(varLambda);

					ip_IsotHard[ip] = fInternal[kIsotHard];
				}
			}

			dArrayT ip_LapIsotHard(fNumIP);
			ip_LapIsotHard = Laplacian(ip_IsotHard, 1);

		        for (int ip = 0; ip < fNumIP; ip ++)
			{
			        /* load internal variables */
			        LoadData(element, ip);

			        /* update nonlocal variables */
			        fInternal[kNlIsotHard] = fInternal[kIsotHard] + c2 * ip_LapIsotHard[ip];
				fNlKineHard = fKineHard;

				/* step 7. update unit normal and yield criteria */
				UpdateState();

				if (fInternal[kYieldCrt] > kYieldTol)
				        flags[ip] = kIsPlastic;
			}

			/* check if stress state for all ip is elastic or returned to yield surface */
			Converged = CheckElementState(element);
		}

		/* check for failure */
		if (count == max_iteration)
		{
		        cout << "\n GradJ2SSNonlinHard::SolveState: local iteration failed after " 
			     << max_iteration << " iterations" << endl;
			throw eGeneralFail;
		}
	}

        for (int ip = 0; ip < fNumIP; ip ++)
	{
	        /* load internal variables */
	        LoadData(element, ip);

		/* check for inelastic processes */
		if (flags[ip] == kIsPlastic)
		        /* step 8. compute consistent tangent moduli */
		        TangentModuli();

		else if (flags[ip] == kIsElastic)
		        /* step 4. compute elastic modulus */
		        fModulus = HookeanMatT::Modulus();
		else
	        {
		        cout << "\n GradJ2SSNonlinHard::SolveState: bad flag value " ;
		        throw eGeneralFail;
		}
	}
}

/* return a pointer to a new element object constructed with
* the data from element */
void GradJ2SSNonlinHard::AllocateAllElements(void)
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
	d_size += dim;          //fNlKineHard
	d_size += dim;          //fNlKineHard_n
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
		element.Allocate(fNumIP, d_size);
	
		/* initialize values */
		element.IntegerData() = kIsElastic;
		element.DoubleData()  = 0.0;
	}
}

double GradJ2SSNonlinHard::YieldCondition(const dSymMatrixT& relstress,
	double isotropic) const
{
	return sqrt(relstress.ScalarProduct()) - sqrt23*(yield+isotropic);
}

/***********************************************************************
* Private
***********************************************************************/

/* load element data for the specified integration point */
void GradJ2SSNonlinHard::LoadData(const ElementCardT& element, int fCurrIP)
{
	/* fetch arrays */
	dArrayT& d_array = element.DoubleData();
	
	/* decode */
	int dim   = dSymMatrixT::NumValues(kNSD);
	int block = 10*dim + 2*kNumInternal + dim*dim;
	int dex   = fCurrIP*block;

        fStress.Set      (kNSD,         &d_array[dex                ]);
        fStress_n.Set    (kNSD,         &d_array[dex += dim         ]);
        fPlstStrn.Set    (kNSD,         &d_array[dex += dim         ]);
        fPlstStrn_n.Set  (kNSD,         &d_array[dex += dim         ]);
        fUnitNorm.Set    (kNSD,         &d_array[dex += dim         ]);
        fUnitNorm_n.Set  (kNSD,         &d_array[dex += dim         ]);
        fKineHard.Set    (kNSD,         &d_array[dex += dim         ]);
        fKineHard_n.Set  (kNSD,         &d_array[dex += dim         ]);
        fNlKineHard.Set  (kNSD,         &d_array[dex += dim         ]);
        fNlKineHard_n.Set(kNSD,         &d_array[dex += dim         ]);
        fInternal.Set    (kNumInternal, &d_array[dex += dim         ]);
        fInternal_n.Set  (kNumInternal, &d_array[dex += kNumInternal]);
        fModulus.Set     (dim,dim,      &d_array[dex += kNumInternal]);
}

/* computes the increment in the plasticity parameter */
void GradJ2SSNonlinHard::IncrementPlasticParameter(double& varLambda)
{
        /* operations to compute dot product for varLambda */
        fUnitNorm.ToMatrix(fmatx1);
        fUnitNorm_n.ToMatrix(fmatx2);
        fNlKineHard_n.ToMatrix(fmatx3);
	double cnn = dMatrixT::Dot(fmatx1, fmatx2);
	double cnx = dMatrixT::Dot(fmatx1, fmatx3);

	/* stiffness */
	double dYieldCrt = (2*fmu+k1)*cnn - k1*k3*cnx
			     + k2*(sqrt23-k4*fInternal_n[kNlIsotHard]);

	if (dYieldCrt < kSmall)
	{
		cout << "\n GradJ2SSNonlinHardT::StressCorrection: consistency function is nonconvex" << endl;
		throw eGeneralFail;
	}
		
	/* variation of plasticity multiplier */
	varLambda = fInternal[kYieldCrt]/dYieldCrt;

	/* increment of plasticity */
	fInternal[kdelLmbda] += varLambda;
}

/* computes the increments in the stress and internal variables */
void GradJ2SSNonlinHard::IncrementState(const double& varLambda)
{
	/* increment stress */
	fsymmatx1.SetToScaled(-2.0*fmu*varLambda, fUnitNorm_n);
	fStress += fsymmatx1;

	/* increment kinematic hardening */
	fsymmatx1.SetToScaled(k1*varLambda, fUnitNorm_n);
	fsymmatx1.AddScaled(-1.0*k1*k3*varLambda, fNlKineHard_n);
	fKineHard += fsymmatx1;

	/* increment isotropic hardening */
	fInternal[kIsotHard] += k2*varLambda*(sqrt23-k4*fInternal_n[kNlIsotHard]);

	/* increment plastic strain */
	fPlstStrn = fPlstStrn_n;
	fPlstStrn.AddScaled(fInternal[kdelLmbda], fUnitNorm_n);
}

/* computes the unit normal and the yield condition */
void GradJ2SSNonlinHard::UpdateState()
{
        /* compute relative stress */
	fRelStress.Deviatoric(fStress);
	fRelStress.AddScaled(-1.0, fNlKineHard);

	/* compute unit normal to yield surface */
	fUnitNorm.SetToScaled(1.0/ sqrt(fRelStress.ScalarProduct()), fRelStress);

	/* compute yield criteria */ 
	fInternal[kYieldCrt] = YieldCondition(fRelStress,fInternal[kNlIsotHard]);
}

/* computes the consistent tangent moduli */
void GradJ2SSNonlinHard::TangentModuli()
{
	/* initialize moduli correction */
	fModuliCorr = 0.0;

	/* compute corrections to elastic moduli */
	fUnitNorm.ToMatrix(fmatx1);
	fUnitNorm_n.ToMatrix(fmatx2);
        fNlKineHard_n.ToMatrix(fmatx3);
	double cnn = dMatrixT::Dot(fmatx1, fmatx2);
	double cnx = dMatrixT::Dot(fmatx1, fmatx3);

	ftnsr1.Outer(fUnitNorm_n,fUnitNorm);
	double h = 2.0*fmu*cnn + k1*(cnn - k3*cnx) + k2*(sqrt23 - k4*fInternal_n[kNlIsotHard]);
	fModuliCorr.AddScaled(-4*fmu*fmu/h,ftnsr1);

	/* make corrections to elastic moduli */
	fModulus.SumOf(HookeanMatT::Modulus(), fModuliCorr);
}

bool GradJ2SSNonlinHard::CheckElementState(const ElementCardT& element)
{
        int ip = 0;
	bool test = true;
	while (ip < fNumIP && test)
	{
	        LoadData(element, ip);
	        test = (fInternal[kYieldCrt] < kYieldTol);
		ip++;
	}
	return test;
}

dArrayT GradJ2SSNonlinHard::Laplacian(const dArrayT& ip_field, int field_length)
{
	dArrayT     dA_ip_lap_field(fNumIP);
	int fNumSD = NumSD();

	/* initialize laplacian */
	dA_ip_lap_field = 0.;

        if (field_length == 1)
	{
                LocalArrayT LA_nd_field(LocalArrayT::kUnspecified,fNumNodes,field_length);
                LocalArrayT LA_nd_grad_field(LocalArrayT::kUnspecified,fNumNodes,1);

		dArrayT     dA_nd_field(fNumNodes);
		dArrayT     dA_nd_grad_field(fNumNodes);

		dMatrixT    dM_ip_grad_field(field_length,fNumSD);
		dMatrixT    dM_ip_secgrad_field(field_length,fNumSD);

		ArrayT<dArrayT> A_dA_ip_grad_field(fNumSD);
		for (int sd = 0; sd < fNumSD; sd++)
		        A_dA_ip_grad_field[sd].Allocate(fNumIP);
	
		/* extrapolate values of field from ip to nodes */
		ContinuumElement().IP_ExtrapolateAll(ip_field,dA_nd_field);

		/* move nodal data from dArrayT to LocalArrayT */
		LA_nd_field.Copy(fNumNodes, 1, dA_nd_field);

		for (int ip = 0; ip < fNumIP; ip ++)
		{
		        /* compute gradient of nodal field at ip */
		        ContinuumElement().IP_ComputeGradient(LA_nd_field,dM_ip_grad_field,ip);

			for (int sd = 0; sd < fNumSD; sd ++)
			{
			        /* store the field of derivatives wrt sd at ip in a dArray */
			        A_dA_ip_grad_field[sd][ip] = dM_ip_grad_field[sd];
			}
		}

		for (int sd = 0; sd < fNumSD; sd ++)
		{ 
		        /* extrapolate the field of derivatives wrt sd from ips to nodes */
		        ContinuumElement().IP_ExtrapolateAll(A_dA_ip_grad_field[sd],dA_nd_grad_field);

			/* move nodal data from dArrayT to LocalArrayT */
			LA_nd_grad_field.Copy(fNumNodes, 1, dA_nd_grad_field);

			for (int ip = 0; ip < fNumIP; ip ++)
			{
			        /* compute gradient of nodal field at ip */
			        ContinuumElement().IP_ComputeGradient(LA_nd_grad_field,dM_ip_secgrad_field,ip);

				/* add the second derivative wrt sd at ip to the laplacian at ip */
				dA_ip_lap_field[ip] += dM_ip_secgrad_field[sd];
			}
		}

			
        }
	else
	{
	        cout << "\n GradJ2SSNonlinHardT::Laplacian: laplacian of multi-dimensional array not yet implemented" << endl;
		throw eGeneralFail;
	}

	return dA_ip_lap_field;
}
