/* $Id: GradJ2SSNonlinHard.cpp,v 1.15 2004-07-22 21:10:18 paklein Exp $ */
#include "GradJ2SSNonlinHard.h"

#include "ifstreamT.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "SSMatSupportT.h"

#include "ContinuumElementT.h" //needed for global information about nodes

using namespace Tahoe;

/* parameters */
const int    kNumInternal = 7;
const double sqrt23       = sqrt(2.0/3.0);
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;
const double fStateTol    = 1.0e-1;

/* element output data */
const int    kNumOutput = 10;
static const char* Labels[kNumOutput] = {
	"VMStrss",   // Von Mises stress
        "Prssure",   // pressure
        "EqPStrn",   // equivalent plastic strain
	"dEqPStrn",  // increment of equivalent plastic strain
	"PlsMult",   // plastic multiplier
        "IsoHardCF", // isotropic hardening conjugate force
        "LapIsoCF", // laplacian of isotropic hardening conjugate force
        "IsoHard",   // isotropic hardening
        "LapIso",   // isotropic hardening
        "NLIsoHard"};  // nonlocal isotropic hardening

/* constructor */
GradJ2SSNonlinHard::GradJ2SSNonlinHard(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_J2_nonlocal"),
	HookeanMatT  (kNSD),
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
        in >> yield >> k1 >> k2 >> k3 >> k4 >> c1 >> c2;

	if (yield < 0 )
	{
                cout << "\n GradJ2SSNonlinHard: yield <0" << endl;
		throw ExceptionT::kBadInputValue;
	}
	if (k1 < 0 || k2 < 0 || k3 < 0 || k4 <0)
	        cout << "\n GradJ2SSNonlinHard: Warning, nonnegative hardening parameter k1, k2, k3, or k4" << endl;
	if (k2 == 0)
	{
                cout << "\n GradJ2SSNonlinHard: k2 = 0" << endl;
		throw ExceptionT::kBadInputValue;
	}
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
	        fKinHardCF_n   = fKinHardCF;
		fNLKinHardCF_n = fNLKinHardCF;
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
	        fKinHardCF   = fKinHardCF_n;
		fNLKinHardCF = fNLKinHardCF_n;
	        fInternal   = fInternal_n;
		flags[ip]   = kIsElastic;
        }
}

#if 0
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
#endif

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

	int iteration = fSSMatSupport->IterationNumber();

	if (fSSMatSupport->RunState() == GlobalT::kFormRHS && fCurrIP == 0)
	{
	        if (iteration > -1)
		        /* solve state at each integration point (all at once) */
		        SolveState(element);
		else
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
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void GradJ2SSNonlinHard::ComputeOutput(dArrayT& output)
{
        /* gather element/integ point information */
        ElementCardT& element = CurrentElement();

        int elem = CurrElementNumber();
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

	/* isotropic hardening conjugate force*/
	output[5] = fInternal[kIsoHardCF];

	/* laplacian of isotropic hardening conjugate force*/
	output[6] = fInternal[kLapIsoCF];

	/* isotropic hardening */
	output[7] = fInternal[kIsoHardCF] / k2;

	/* laplacian of isotropic hardening */
	output[8] = fInternal[kLapIsoCF] / k2;

	/* nonlocal isotropic hardening */
	output[9] = fInternal[kNLIsoHardCF] / k2;
}

/* information about subordinate parameter lists */
void GradJ2SSNonlinHard::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	 SSIsotropicMatT::DefineSubs(sub_list);
	 HookeanMatT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradJ2SSNonlinHard::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = SSIsotropicMatT::NewSub(name);
	if (sub)
		return sub;
	else
		return HookeanMatT::NewSub(name);
}

/* accept parameter list */
void GradJ2SSNonlinHard::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	 SSIsotropicMatT::TakeParameterList(list);
	 HookeanMatT::TakeParameterList(list);

	/* allocate space for all elements */
	AllocateAllElements();
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

/* solve for the state at ip: only solved for fCurrIP = 0 */
void GradJ2SSNonlinHard::SolveState(ElementCardT& element)
{
        int elem = CurrElementNumber();

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
		HookeanStress(e_els, fStress_n);
		fStress = fStress_n;

		/* step 3. initialize unit normal and yield criteria */
		UpdateState();

		if (fInternal[kYieldCrt] > kYieldTol)
		        flags[ip] = kIsPlastic;
		else
		        flags[ip] = kIsElastic;

	        /* step 4. zero the increment in plasticity parameter */
	        fInternal[kdelLmbda] = 0.;
	}

        /* check for inelastic processes at any ip in element */
        bool Converged = CheckElementState(element);

        /* local Newton iteration */
        int max_iteration = 60;
        int count = 0;

	while (!Converged)
	{
		/*array of IsotHard at all ip in current element used to compute Laplacian*/
		dArrayT  ip_IsoHardCF(fNumIP);
		ip_IsoHardCF = 0.;

		/*array of the Laplacian of IsoHard at all ip in current element*/
		dArrayT ip_LapIsoCF(fNumIP);
		ip_LapIsoCF = 0.;

	        for (int ip = 0; ip < fNumIP; ip ++)
		{
		        /* load internal variables */
		        LoadData(element, ip);

		        if (fInternal[kYieldCrt] > kYieldTol)
			{
			        double varLambda;

				/* step 5. increment plasticity parameter */
				IncrementPlasticParameter(varLambda);

				/* step 6. increment stress and state variables */
				IncrementState(varLambda);
			}

			/* store isotropic hardening at ip in array */
			ip_IsoHardCF[ip] = fInternal[kIsoHardCF];
		}

		Laplacian(ip_LapIsoCF,ip_IsoHardCF, 1);

	        for (int ip = 0; ip < fNumIP; ip ++)
		{
		        /* load internal variables */
		        LoadData(element, ip);

			fInternal[kLapIsoCF] = ip_LapIsoCF[ip];

			/* update nonlocal variables */
			fInternal[kNLIsoHardCF] = fInternal[kIsoHardCF] + c2 * fInternal[kLapIsoCF];
			fNLKinHardCF = fKinHardCF;

			/* step 7. update unit normal and yield criteria */
			UpdateState();

			if (flags[ip] == kIsElastic && fInternal[kYieldCrt] > kYieldTol)
			       flags[ip] = kIsPlastic;
		}

		/* check if stress state for all ip is elastic or returned to yield surface */
		Converged = CheckElementState(element);

		/* check for failure */
		if (++count == max_iteration)
		{
		        cout << "\n GradJ2SSNonlinHard::SolveState: local iteration failed after " 
			     << max_iteration << " iterations" << endl;
			throw ExceptionT::kGeneralFail;
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
		        /* step 5. compute elastic modulus */
		        fModulus = HookeanMatT::Modulus();
		else
	        {
		        cout << "\n GradJ2SSNonlinHard::SolveState: bad flag value " ;
		        throw ExceptionT::kGeneralFail;
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
	d_size += dim;          //fKinHardCF
	d_size += dim;          //fKinHardCF_n
	d_size += dim;          //fNLKinHardCF
	d_size += dim;          //fNLKinHardCF_n
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
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT sdim = dSymMatrixT::int2DimensionT(kNSD);
	int dim   = dSymMatrixT::NumValues(kNSD);
	int block = 10*dim + 2*kNumInternal + dim*dim;
	int dex   = fCurrIP*block;

	fStress.Alias       (sdim,         &d_array[dex                ]);
	fStress_n.Alias     (sdim,         &d_array[dex += dim         ]);
	fPlstStrn.Alias     (sdim,         &d_array[dex += dim         ]);
	fPlstStrn_n.Alias   (sdim,         &d_array[dex += dim         ]);
	fUnitNorm.Alias     (sdim,         &d_array[dex += dim         ]);
	fUnitNorm_n.Alias   (sdim,         &d_array[dex += dim         ]);
	fKinHardCF.Alias    (sdim,         &d_array[dex += dim         ]);
	fKinHardCF_n.Alias  (sdim,         &d_array[dex += dim         ]);
	fNLKinHardCF.Alias  (sdim,         &d_array[dex += dim         ]);
	fNLKinHardCF_n.Alias(sdim,         &d_array[dex += dim         ]);
	fInternal.Alias     (kNumInternal, &d_array[dex += dim         ]);
	fInternal_n.Alias   (kNumInternal, &d_array[dex += kNumInternal]);
	fModulus.Alias      (dim,dim,      &d_array[dex += kNumInternal]);
}

/* computes the increment in the plasticity parameter */
void GradJ2SSNonlinHard::IncrementPlasticParameter(double& varLambda)
{
        /* operations to compute dot product for variation of plasticity multiplier */
        fUnitNorm.ToMatrix(fmatx1);
        fUnitNorm_n.ToMatrix(fmatx2);
        fNLKinHardCF_n.ToMatrix(fmatx3);
	double cnn = dMatrixT::Dot(fmatx1, fmatx2);
	double cnx = dMatrixT::Dot(fmatx1, fmatx3);

	/* stiffness */
	double dYieldCrt = 2*fmu*cnn + (k1*cnn - k3*cnx)
			     + sqrt23*(sqrt23*k2 - k4*fInternal_n[kNLIsoHardCF]);

	//	if (dYieldCrt < kSmall)
	//	{
	//		cout << "\n GradJ2SSNonlinHardT::StressCorrection: consistency function is nonconvex" << endl;
	//		throw ExceptionT::kGeneralFail;
	//	}
		
	varLambda = fInternal[kYieldCrt]/dYieldCrt;

	/* increment of plasticity */
	fInternal[kdelLmbda] += varLambda;

	if (fInternal[kdelLmbda] < 0.0)
	        fInternal[kdelLmbda] = 0.0;
}

/* computes the increments in the stress and internal variables */
void GradJ2SSNonlinHard::IncrementState(const double& varLambda)
{
#pragma unused(varLambda)
        //	/* increment stress */
        //        fsymmatx1.SetToScaled(-2.0*fmu*varLambda, fUnitNorm_n);
        //	fStress += fsymmatx1;
	//
	//	/* increment kinematic hardening */
	//	fsymmatx1.SetToScaled(k1*varLambda, fUnitNorm_n);
	//	fsymmatx1.AddScaled(-1.0*k3*varLambda, fNLKinHardCF_n);
	//	fKinHardCF += fsymmatx1;
	//
	//	/* increment isotropic hardening */
	//	fInternal[kIsoHardCF] += varLambda*(k2*sqrt23-k4*fInternal_n[kNLIsoHardCF]);
	//
	//	/* increment plastic strain */
	//	fsymmatx1.SetToScaled(varLambda, fUnitNorm_n);
	//	fPlstStrn += fsymmatx1;
	//

	/* increment stress */
        fsymmatx1.SetToScaled(-2.0*fmu*fInternal[kdelLmbda], fUnitNorm_n);
	fStress = fStress_n;
	fStress += fsymmatx1;

	/* increment kinematic hardening */
	fsymmatx1.SetToScaled(k1*fInternal[kdelLmbda], fUnitNorm_n);
	fsymmatx1.AddScaled(-1.0*k3*fInternal[kdelLmbda], fNLKinHardCF_n);
	fKinHardCF = fKinHardCF_n;
	fKinHardCF += fsymmatx1;

	/* increment isotropic hardening */
	fInternal[kIsoHardCF] = fInternal_n[kIsoHardCF];
	fInternal[kIsoHardCF] += fInternal[kdelLmbda]*(k2*sqrt23-k4*fInternal_n[kNLIsoHardCF]);

	/* increment plastic strain */
	fsymmatx1.SetToScaled(fInternal[kdelLmbda], fUnitNorm_n);
	fPlstStrn = fPlstStrn_n;
	fPlstStrn += fsymmatx1;
}

/* computes the unit normal and the yield condition */
void GradJ2SSNonlinHard::UpdateState()
{
        /* compute relative stress */
	fRelStress.Deviatoric(fStress);
	fRelStress.AddScaled(-1.0, fNLKinHardCF);

	/* compute unit normal to yield surface */
	fUnitNorm.SetToScaled(1.0/ sqrt(fRelStress.ScalarProduct()), fRelStress);

	/* compute yield criteria */ 
	fInternal[kYieldCrt] = YieldCondition(fRelStress,fInternal[kNLIsoHardCF]);
}

/* computes the consistent tangent moduli */
void GradJ2SSNonlinHard::TangentModuli()
{
	/* initialize moduli correction */
	fModuliCorr = 0.0;

	/* compute corrections to elastic moduli */
	fUnitNorm.ToMatrix(fmatx1);
	fUnitNorm_n.ToMatrix(fmatx2);
        fNLKinHardCF_n.ToMatrix(fmatx3);
	double cnn = dMatrixT::Dot(fmatx1, fmatx2);
	double cnx = dMatrixT::Dot(fmatx1, fmatx3);

	ftnsr1.Outer(fUnitNorm_n,fUnitNorm);
	double h = 2.0*fmu*cnn + (k1*cnn - k3*cnx) + sqrt23*(k2*sqrt23 - k4*fInternal_n[kNLIsoHardCF]);
	fModuliCorr.AddScaled(-4*fmu*fmu/h,ftnsr1);

	/* make corrections to elastic moduli */
	fModulus.SumOf(HookeanMatT::Modulus(), fModuliCorr);
}

bool GradJ2SSNonlinHard::CheckElementState(const ElementCardT& element)
{
	bool test = true;

        for (int ip = 0; ip < fNumIP; ip ++)
	{
	        LoadData(element, ip);
		if (fInternal[kYieldCrt] > kYieldTol)
		        test = false;
	}

	return test;
}

void GradJ2SSNonlinHard::Laplacian(dArrayT& ip_laplacian_field, const dArrayT& ip_field, int field_length)
{
        int fNumIP = NumIP();
	int fNumSD = NumSD();
	int fNumNodes = ContinuumElement().InitialCoordinates().NumberOfNodes();

        int elem = CurrElementNumber();

	ip_laplacian_field = 0.;

	//	dArrayT ip_coord(fNumSD);
	//	dArrayT ip_coord_x(fNumIP);
	//	dArrayT ip_coord_y(fNumIP);
	//	dArrayT nd_coord_x(fNumNodes);
	//	dArrayT nd_coord_y(fNumNodes);
	//	ip_coord = 0.;
	//
	//	for (int ip = 0; ip < fNumIP; ip ++)
	//	{
	//	        ContinuumElement().IP_Coords(ip_coord, ip); 
	//		ip_coord_x[ip] = ip_coord[0];
	//		ip_coord_y[ip] = ip_coord[1];
	//	}

	//	/* extrapolate values of field from ip to nodes */
	//	ContinuumElement().IP_ExtrapolateAll(ip_coord_x,nd_coord_x);
	//	ContinuumElement().IP_ExtrapolateAll(ip_coord_y,nd_coord_y);

        if (field_length == 1)
	{
                LocalArrayT LA_nd_field(LocalArrayT::kUnspecified,fNumNodes,field_length);
                LocalArrayT LA_nd_grad_field(LocalArrayT::kUnspecified,fNumNodes,field_length);

		dArrayT     dA_nd_field(fNumNodes);
		dArrayT     dA_nd_grad_field(fNumNodes);
		dArrayT     dA_nd_lap_field(fNumNodes);

		dMatrixT    dM_ip_grad_field(field_length,fNumSD);
		dMatrixT    dM_ip_secgrad_field(field_length,fNumSD);

		ArrayT<dArrayT> A_dA_ip_grad_field(fNumSD);
		for (int sd = 0; sd < fNumSD; sd++)
		        A_dA_ip_grad_field[sd].Dimension(fNumIP);
	
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
				ip_laplacian_field[ip] += dM_ip_secgrad_field[sd];
			}
		}

		ContinuumElement().IP_ExtrapolateAll(ip_laplacian_field,dA_nd_lap_field);
        }
	else
	{
	        cout << "\n GradJ2SSNonlinHardT::Laplacian: laplacian of multi-dimensional array not yet implemented" << endl;
		throw ExceptionT::kGeneralFail;
	}
}
