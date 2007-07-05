/* $Id: TersoffSolverT_surf.cpp,v 1.4 2007-07-05 14:35:15 hspark Exp $ */
#include "TersoffSolverT_surf.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"
#include "Tersoff_inc_surf.h"

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
TersoffSolverT_surf::TersoffSolverT_surf(const ThermalDilatationT* thermal, int normal):
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
	fSurfaceThickness(0.0),
	fNormalCode(normal),
	f_area0(0.0)
{
//	SetName("Tersoff_CB_solver");
}

/* Destructor */
TersoffSolverT_surf::~TersoffSolverT_surf(void)
{
//	delete fGeometry;
}

/* moduli - assume Xsi already determined */
void TersoffSolverT_surf::SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	/* compute second derivatives wrt {C,C} and {C,Xsi} */
	get_ddC_surf(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dCdC_hat.Pointer(), dCdXsi_hat.Pointer());

	/* Compute moduli */
	moduli = dCdC_hat;

	if (fEquilibrate)
	{
		dXsidXsi.Inverse();
		fTempMixed.MultAB(dCdXsi_hat, dXsidXsi);
		fTempRank4.MultABT(fTempMixed,dCdXsi_hat);
		moduli -= fTempRank4;
	}
	moduli *= 4.0;
	moduli *= f_omega0;
}

//for now return symmetric matrix
void TersoffSolverT_surf::SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress)
{
	/* set internal equilibrium */
	if (fEquilibrate) Equilibrate(CIJ, Xsi);

	/* call C function */
	get_dUdC_surf(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(),
		stress.Pointer()); 
	stress *= 2.0;
	stress *= f_omega0;
	
#if 0
	cout << "stress: " << stress.no_wrap() << endl;
#endif
}

/* strain energy density */
double TersoffSolverT_surf::StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi)
{
#if 0
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

// 	return( (f2Body->Phi()).Sum() + (f3Body->Phi()).Sum() );
#endif

//not implemented
return 0.0;
}

/* describe the parameters needed by the interface */
void TersoffSolverT_surf::DefineParameters(ParameterListT& list) const
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
}

/* information about subordinate parameter lists */
//void TersoffSolverT_surf::DefineSubs(SubListT& sub_list) const
//{
	/* inherited */
//	ParameterInterfaceT::DefineSubs(sub_list);

	/* crystal orientation */
//	sub_list.AddSub("FCC_lattice_orientation", ParameterListT::Once, true);

	/* choice of potentials */
//	sub_list.AddSub("DC_potential_choice", ParameterListT::Once, true);
//}

/* a pointer to the ParameterInterfaceT of the given subordinate */
//ParameterInterfaceT* TersoffSolverT_surf::NewSub(const StringT& name) const
//{
// 	if (name == "DC_potential_choice")
// 	{
// 		ParameterContainerT* choice = new ParameterContainerT(name);
// 		choice->SetSubSource(this);
// 		choice->SetListOrder(ParameterListT::Choice);
// 	
// 		choice->AddSub("Stillinger-Weber");
// 
// 		ParameterContainerT PTHT("PTHT");
// 		PTHT.AddParameter(ParameterT::Double, "A");
// 		PTHT.AddParameter(ParameterT::Double, "A1");
// 		PTHT.AddParameter(ParameterT::Double, "A2");
// 		
// 		PTHT.AddParameter(ParameterT::Double, "B");
// 		PTHT.AddParameter(ParameterT::Double, "Z");
// 		choice->AddSub(PTHT);
// 
// 		//choice->AddSub(ParameterContainerT("Tersoff"));
// 
// 		return choice;
// 	}
// 	else if (name == "FCC_lattice_orientation")
// 	{
// 		FCCLatticeT lattice(0);
// 		return lattice.NewSub(name);
// 	}
// 	else if (name == "Stillinger-Weber")
// 		return new SWDataT;
// 	else /* inherited */
// 		return ParameterInterfaceT::NewSub(name);
//}

/* accept parameter list */
void TersoffSolverT_surf::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* dimension work space */
	dXsi.Dimension(kNumDOF);
	dXsidXsi.Dimension(kNumDOF);
	dCdC_hat.Dimension(kStressDim);
	dCdXsi_hat.Dimension(kStressDim,kNumDOF);
	fTempRank4.Dimension(kStressDim);
	fTempMixed.Dimension(kStressDim, kNumDOF);

#if 0
	fMatrices.Dimension(kNumDOF);
	fMat2.Dimension(kNumDOF);
	fGradl_i.Dimension(3,kNumDOF); 
	fSymMat1.Dimension(kNSD);
	fGradl_C.Dimension(3,kStressDim);
#endif
	fMat1.Dimension(kNumDOF); 
	fVec.Dimension(kNumDOF);

	/* unit cell coordinates - FOR IDEAL, UNRECONSTRUCTED {100} SURFACE */
	/* NORMAL IS IN [001] DIRECTION */
	fUnitCellCoords.Dimension(9, 3); /* [9 atoms] x [3 dim]: first atom is 'center' */
	fUnitCellCoords(0,0) = 0.00;
	fUnitCellCoords(1,0) = 0.25;
	fUnitCellCoords(2,0) = 0.00;
	fUnitCellCoords(3,0) = 0.50;
	fUnitCellCoords(4,0) = 0.50;
	fUnitCellCoords(5,0) = -0.25;
	fUnitCellCoords(6,0) = -0.50;
	fUnitCellCoords(7,0) = -0.50;
	fUnitCellCoords(8,0) = 0.00;

	fUnitCellCoords(0,1) = 0.00;
	fUnitCellCoords(1,1) = -0.25;
	fUnitCellCoords(2,1) = -0.50;
	fUnitCellCoords(3,1) = 0.00;
	fUnitCellCoords(4,1) = -0.50;
	fUnitCellCoords(5,1) = 0.25;
	fUnitCellCoords(6,1) = 0.50;
	fUnitCellCoords(7,1) = 0.00;
	fUnitCellCoords(8,1) = 0.50;

	fUnitCellCoords(0,2) = 0.00;
	fUnitCellCoords(1,2) = -0.25;
	fUnitCellCoords(2,2) = -0.50;
	fUnitCellCoords(3,2) = -0.50;
	fUnitCellCoords(4,2) = 0.00;
	fUnitCellCoords(5,2) = -0.25;
	fUnitCellCoords(6,2) = 0.00;
	fUnitCellCoords(7,2) = -0.50;
	fUnitCellCoords(8,2) = -0.50;

	/* flag */
	fEquilibrate = list.GetParameter("equilibrate");

	/* resolve orientation */
// 	FCCLatticeT lattice(0);
// 	const ParameterListT& orientation = list.GetListChoice(lattice, "FCC_lattice_orientation");
// 	dMatrixT Q;
// 	FCCLatticeT::SetQ(orientation, Q);
// 	
// 	/* construct bond lattice */
// 	fGeometry = new LengthsAndAnglesT(Q,fPairs);

	/* set potentials */
//	const ParameterListT& potential = list.GetListChoice(*this, "DC_potential_choice");

	/* All parameters required */
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
	
	/* scale unit cell coordinates */
	fUnitCellCoords *= f_a0;

	/* Calculate atomic volume = a^{3}/8 */
	double asdf = f_a0*f_a0*f_a0/8.0;
	f_omega0 = 1.0/asdf;

	/* CALCULATE ATOMIC AREA!! */
	f_area0 = 1.0;

	/* Calculate surface thickness - TEMPORARY VALUE! */
	fSurfaceThickness = 0.375*f_a0;

	/* write into vector to pass to C code */
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
}

/**********************************************************************
 * Private
 **********************************************************************/

/* Minimize the energy wrt Xsi using the initial value passed */
void TersoffSolverT_surf::Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi)
{
	bool debug = false;

	/* check initial value */
	SetdXsi(CIJ, Xsi);

	int count = 0;
	if (debug) cout << dXsi.Magnitude() << '\n';
	while (count++ < 15 && dXsi.Magnitude() > 1.0e-12)
	{
		fMat1.Inverse(dXsidXsi);
		fMat1.Multx(dXsi, fVec);
		
		Xsi -= fVec;
		
		/* recompute */
		SetdXsi(CIJ, Xsi);
		if (debug) cout << dXsi.Magnitude() << '\n';
	}

	/* assume not converged */
	if (count == 15) ExceptionT::GeneralFail("TersoffSolverT_surf::Equilibrate", "failed");
}

/* set free dof - triggers recomputation */
void TersoffSolverT_surf::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* call C function */
	get_dXsi_surf(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dXsi.Pointer(), dXsidXsi.Pointer());

//debugging
#if 0
cout << dXsi.no_wrap() << ":" << dXsidXsi.no_wrap() << endl;
#endif
}
