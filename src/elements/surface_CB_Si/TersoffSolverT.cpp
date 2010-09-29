/* $Id: TersoffSolverT.cpp,v 1.7 2010-09-29 14:50:13 hspark Exp $ */
#include "TersoffSolverT.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"
#include "Tersoff_inc.h"
#include "PDM_inc.h"

using namespace Tahoe;

const int kNSD       = 3;
const int kNumDOF    = 3;
const int kNumAngles = 12;
const int kStressDim =dSymMatrixT::NumValues(kNSD);

/* set pair numbers */
static int pairdata[kNumAngles*2] =
{2,  3,
					1,  3,
					1,  2,
					0,  3,
					0,  2,
					0,  1,
					5,  6,
					4,  6,
					4,  5,
					0,  6,
					0,  5,
					0,  4};

/* Constructor */
TersoffSolverT::TersoffSolverT(const ThermalDilatationT* thermal):
	ParameterInterfaceT("Tersoff_CB_solver"),
	fEquilibrate(true),
	fThermal(thermal),
	fPairs(kNumAngles, 2, pairdata),
//	fGeometry(NULL),
	f_A(0.0),
	f_B(0.0),
	f_lambda(0.0),
	f_mu(0.0),
	f_beta(0.0),
	f_n(0.0),
	f_c(0.0),
	f_d(0.0),
	f_h(0.0),
	f_chi(0.0),
	f_R(0.0),
	f_S(0.0),
	f_omega0(0.0),
	f_ex(0.0),
	f_ey(0.0),
	f_ez(0.0),
	f_econv(0.0),
	f_alphatot(0.0),
	f_alpha1(0.0)
{

}

/* Destructor */
TersoffSolverT::~TersoffSolverT(void)
{
//	delete fGeometry;
}

/* moduli - assume Xsi already determined */
void TersoffSolverT::SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	/* compute second derivatives wrt {C,C} and {C,Xsi} for Tersoff/mechanical part */
	ddC_driver(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dCdC_hat.Pointer(), dCdXsi_hat.Pointer());

	/* compute second derivatives wrt {C,C} and {C,Xsi} for PDM/electrostatic part */
	get_ddC_pdm(fParams_pdm.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dCdC_hat_pdm.Pointer(), dCdXsi_hat_pdm.Pointer());	

	/* Compute moduli */
//	moduli = dCdC_hat;
	moduli = dCdC_hat;
	moduli+=dCdC_hat_pdm;
	dCdXsi_hat_tot = dCdXsi_hat;
	dCdXsi_hat_tot+=dCdXsi_hat_pdm;

	if (fEquilibrate)
	{
		dXsidXsi_tot.Inverse();
		fTempMixed.MultAB(dCdXsi_hat_tot, dXsidXsi_tot);
		fTempRank4.MultABT(fTempMixed,dCdXsi_hat_tot);
		moduli -= fTempRank4;
	}
	moduli *= 4.0;
	moduli *= f_omega0;
}

//for now return symmetric matrix
void TersoffSolverT::SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress)
{
	/* set internal equilibrium */
	if (fEquilibrate) Equilibrate(CIJ, Xsi);

	/* call C function for bulk mechanical stress */
	get_dUdC(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(),
		stress.Pointer()); 
		
	/* call C function for bulk electrostatic/PDM stress */
	get_dUdC_pdm(fParams_pdm.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(),
		stress_pdm.Pointer()); 
		
	stress+=stress_pdm;
	stress *= 2.0;
	stress *= f_omega0;
	
#if 0
	cout << "stress: " << stress.no_wrap() << endl;
#endif
}

/* strain energy density */
double TersoffSolverT::StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);
	
	/* compute bulk mechanical energy */
	double energy = get_energy(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer());

	/* compute bulk electrostatic/PDM energy */
	double energy_elec = get_energy_pdm(fParams_pdm.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer());
	
	energy += energy_elec;
	energy *= f_omega0;		// scale by atomic volume
	return energy;
}

/* describe the parameters needed by the interface */
void TersoffSolverT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT equilibrate(ParameterT::Boolean, "equilibrate");
	equilibrate.SetDefault(true);
	list.AddParameter(equilibrate);

	/* lattice parameter */
	ParameterT a0(f_a0, "a0");
	a0.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(a0);

	/* Revised from TersoffPairT.cpp */
	ParameterT mass(fMass, "mass");
	mass.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mass);

	ParameterT A(f_A, "rep_energy_scale_Aij");
	A.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(A);
	
	ParameterT B(f_B, "attr_energy_scale_Bij");
	B.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(B);
	
	ParameterT lambda(f_lambda, "rep_energy_exponent_lambdaij");
	lambda.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(lambda);
	
	ParameterT mu(f_mu, "attr_energy_exponent_muij");
	mu.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mu);
	
	ParameterT beta(f_beta, "bond_order_coeff1_betai");
	beta.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(beta);
	
	ParameterT n(f_n, "bond_order_exponent_ni");
	n.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(n);
	
	ParameterT c(f_c, "bond_order_coeff2_ci");
	c.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c);
	
	ParameterT d(f_d, "bond_order_coeff3_di");
	d.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(d);
	
	ParameterT h(f_h, "bond_order_coeff4_hi");
	list.AddParameter(h);
	
	ParameterT chi(f_chi, "bond_order_scaling_chi_ij");
	chi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi);
	
	ParameterT R(f_R, "cutoff_func_length_1_Rij");
	R.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(R);
	
	ParameterT S(f_S, "cutoff_func_length_2_Sij");
	S.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(S);
	
	ParameterT ex(f_ex, "x_direction_efield");
	ex.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(ex);
	
	ParameterT ey(f_ey, "y_direction_efield");
	ey.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(ey);
	
	ParameterT ez(f_ez, "z_direction_efield");
	ez.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(ez);
	
	ParameterT econv(f_econv, "energy_conversion_parameter");
	econv.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(econv);
	
	ParameterT alphatot(f_alphatot, "thole_s");
	alphatot.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alphatot);
	
	ParameterT alpha1(f_alpha1, "thole_alpha1");
	alpha1.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha1);
}

/* accept parameter list */
void TersoffSolverT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* dimension work space */
	dXsi.Dimension(kNumDOF);
	dXsi_pdm.Dimension(kNumDOF);
	dXsi_tot.Dimension(kNumDOF);
	dXsidXsi.Dimension(kNumDOF);
	dXsidXsi_pdm.Dimension(kNumDOF);
	dXsidXsi_tot.Dimension(kNumDOF);
	dCdC_hat.Dimension(kStressDim);
	dCdC_hat_pdm.Dimension(kStressDim);
	dCdC_hat_tot.Dimension(kStressDim);
	dCdXsi_hat.Dimension(kStressDim,kNumDOF);
	dCdXsi_hat_pdm.Dimension(kStressDim,kNumDOF);
	dCdXsi_hat_tot.Dimension(kStressDim,kNumDOF);
	fTempRank4.Dimension(kStressDim);
	fTempMixed.Dimension(kStressDim, kNumDOF);
	stress_pdm.Dimension(kNumDOF);

#if 0
	fMatrices.Dimension(kNumDOF);
	fMat2.Dimension(kNumDOF);
	fGradl_i.Dimension(3,kNumDOF); 
	fSymMat1.Dimension(kNSD);
	fGradl_C.Dimension(3,kStressDim);
#endif
	fMat1.Dimension(kNumDOF); 
	fVec.Dimension(kNumDOF);

	/* unit cell coordinates */
	fUnitCellCoords.Dimension(5, 3); /* [5 atoms] x [3 dim]: first atom is 'center' */
	fUnitCellCoords(0,0) = 0.25;
	fUnitCellCoords(1,0) = 0.00;
	fUnitCellCoords(2,0) = 0.50;
	fUnitCellCoords(3,0) = 0.50;
	fUnitCellCoords(4,0) = 0.00;

	fUnitCellCoords(0,1) = 0.25;
	fUnitCellCoords(1,1) = 0.00;
	fUnitCellCoords(2,1) = 0.50;
	fUnitCellCoords(3,1) = 0.00;
	fUnitCellCoords(4,1) = 0.50;

	fUnitCellCoords(0,2) = 0.25;
	fUnitCellCoords(1,2) = 0.00;
	fUnitCellCoords(2,2) = 0.00;
	fUnitCellCoords(3,2) = 0.50;
	fUnitCellCoords(4,2) = 0.50;

	/* flag */
	fEquilibrate = list.GetParameter("equilibrate");

	/* All parameters required for Tersoff Si */
	f_a0 = list.GetParameter("a0");
	fMass = list.GetParameter("mass");
	f_A = list.GetParameter("rep_energy_scale_Aij");
	f_B = list.GetParameter("attr_energy_scale_Bij");
	f_lambda = list.GetParameter("rep_energy_exponent_lambdaij");
	f_mu = list.GetParameter("attr_energy_exponent_muij");
	f_beta = list.GetParameter("bond_order_coeff1_betai");
	f_n = list.GetParameter("bond_order_exponent_ni");
	f_c = list.GetParameter("bond_order_coeff2_ci");
	f_d = list.GetParameter("bond_order_coeff3_di");
	f_h = list.GetParameter("bond_order_coeff4_hi");
	f_chi = list.GetParameter("bond_order_scaling_chi_ij");
	f_R = list.GetParameter("cutoff_func_length_1_Rij");
	f_S = list.GetParameter("cutoff_func_length_2_Sij");
	
	/* Parameters required for Thole-based PDM */
	f_ex = list.GetParameter("x_direction_efield");
	f_ey = list.GetParameter("y_direction_efield");
	f_ez = list.GetParameter("z_direction_efield");
	f_econv = list.GetParameter("energy_conversion_parameter");
	f_alphatot = list.GetParameter("thole_s");
	f_alpha1 = list.GetParameter("thole_alpha1");
	
	/* scale unit cell coordinates */
	fUnitCellCoords *= f_a0;

	/* Calculate atomic volume = a^{3}/8 */
	double asdf = f_a0*f_a0*f_a0/8.0;
	f_omega0 = 1.0/asdf;

	/* write into vector to pass to C code for Tersoff Si */
	fParams.Dimension(13);
	fParams[ 0] = f_A;
	fParams[ 1] = f_B;
	fParams[ 2] = fMass;
	fParams[ 3] = f_lambda;
	fParams[ 4] = f_mu;
	fParams[ 5] = f_beta;
	fParams[ 6] = f_n;
	fParams[ 7] = f_c;
	fParams[ 8] = f_d;
	fParams[ 9] = f_h;
	fParams[10] = f_chi;
	fParams[11] = f_R;
	fParams[12] = f_S;	
	
	/* write into vector to pass to C code for electrostatic PDM */
	fParams_pdm.Dimension(6);
	fParams_pdm[ 0] = f_ex;
	fParams_pdm[ 1] = f_ey;
	fParams_pdm[ 2] = f_ez;
	fParams_pdm[ 3] = f_econv;
	fParams_pdm[ 4] = f_alphatot;
	fParams_pdm[ 5] = f_alpha1;
}

double TersoffSolverT::Density(void) const
{
	return 8.0*fMass*1.0365e-4/(f_a0*f_a0*f_a0);
}

/**********************************************************************
 * Private
 **********************************************************************/

/* Minimize the energy wrt Xsi using the initial value passed */
void TersoffSolverT::Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi)
{
	bool debug = false;

	/* check initial value */
	SetdXsi(CIJ, Xsi);

	int count = 0;
	if (debug) cout << dXsi_tot.Magnitude() << '\n';
	while (count++ < 15 && dXsi_tot.Magnitude() > 1.0e-12)
	{
		fMat1.Inverse(dXsidXsi_tot);
		fMat1.Multx(dXsi_tot, fVec);
		
		Xsi -= fVec;
		
		/* recompute */
		SetdXsi(CIJ, Xsi);
		if (debug) cout << dXsi.Magnitude() << '\n';
	}

	/* assume not converged */
	if (count == 15) ExceptionT::GeneralFail("TersoffSolverT::Equilibrate", "failed");
}

/* set free dof - triggers recomputation */
void TersoffSolverT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* call C function for Tersoff */
	get_dXsi(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dXsi.Pointer(), dXsidXsi.Pointer());

	/* call C function for electrostatic PDM */
	get_dXsi_pdm(fParams_pdm.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dXsi_pdm.Pointer(), dXsidXsi_pdm.Pointer());

	/* Total dXsi and total dXsidXsi */
	dXsi_tot = dXsi;
	dXsi_tot+=dXsi_pdm;
	dXsidXsi_tot = dXsidXsi;
	dXsidXsi_tot+=dXsidXsi_pdm;
	
//debugging
#if 0
cout << dXsi.no_wrap() << ":" << dXsidXsi.no_wrap() << endl;
#endif
}
