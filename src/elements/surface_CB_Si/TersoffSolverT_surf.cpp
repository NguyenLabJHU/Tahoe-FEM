/* $Id: TersoffSolverT_surf.cpp,v 1.4 2008-02-15 15:34:03 hspark Exp $ */
#include "TersoffSolverT_surf.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"
#include "Tersoff_inc_surf.h"

using namespace Tahoe;

const int kNSD       = 3;
const int kNumDOF    = 3;
const int kNumAngles = 12;
const int kStressDim =dSymMatrixT::NumValues(kNSD);
const double piby2 = 4.0 * atan(1.0) / 2.0;

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
	moduli *= f_area0;
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
	stress *= f_area0;
	
#if 0
	cout << "stress: " << stress.no_wrap() << endl;
#endif
}

/* strain energy density */
double TersoffSolverT_surf::StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	double surf_energy, surf_energy2;
	surf_energy = get_energy_surf(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer());
	surf_energy2 = get_energy_surf_bulk(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer());
	surf_energy *= f_area0;
	surf_energy2 *= f_area0;	// scale this further by 1/2 or something
	return surf_energy;
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
	/* DOUBLE DIMENSIONS OF XSI-RELATED MATRICES DUE TO 2 SETS OF
	XDOFS FOR SURFACES */
//	dXsi.Dimension(2*kNumDOF);
//	dXsidXsi.Dimension(2*kNumDOF);
	dXsi.Dimension(kNumDOF);
	dXsidXsi.Dimension(kNumDOF);
	dCdC_hat.Dimension(kStressDim);
//	dCdXsi_hat.Dimension(kStressDim,kNumDOF*2);
	dCdXsi_hat.Dimension(kStressDim,kNumDOF);
	fTempRank4.Dimension(kStressDim);
//	fTempMixed.Dimension(kStressDim, kNumDOF*2);
	fTempMixed.Dimension(kStressDim, kNumDOF);

#if 0
	fMatrices.Dimension(kNumDOF);
	fMat2.Dimension(kNumDOF);
	fGradl_i.Dimension(3,kNumDOF); 
	fSymMat1.Dimension(kNSD);
	fGradl_C.Dimension(3,kStressDim);
#endif
	/* DOUBLE DIMENSIONS FOR 2 SETS OF XDOFS */
//	fMat1.Dimension(kNumDOF*2); 
//	fVec.Dimension(kNumDOF*2);
	fMat1.Dimension(kNumDOF); 
	fVec.Dimension(kNumDOF);

	/* unit cell coordinates - FOR IDEAL, UNRECONSTRUCTED {100} SURFACE */
	/* NORMAL IS IN [001] DIRECTION */
	dMatrixT tempUnitCellCoords, temp_bonds2;
	temp_bonds2.Dimension(9, 3);
	tempUnitCellCoords.Dimension(9, 3);	// temporary until bond rotation
	fUnitCellCoords.Dimension(9, 3); /* [9 atoms] x [3 dim]: first atom is 'center' */
	tempUnitCellCoords(0,0) = 0.00;
	tempUnitCellCoords(1,0) = 0.25;
	tempUnitCellCoords(2,0) = 0.00;
	tempUnitCellCoords(3,0) = 0.50;
	tempUnitCellCoords(4,0) = 0.50;
	tempUnitCellCoords(5,0) = -0.25;
	tempUnitCellCoords(6,0) = -0.50;
	tempUnitCellCoords(7,0) = -0.50;
	tempUnitCellCoords(8,0) = 0.00;

	tempUnitCellCoords(0,1) = 0.00;
	tempUnitCellCoords(1,1) = -0.25;
	tempUnitCellCoords(2,1) = -0.50;
	tempUnitCellCoords(3,1) = 0.00;
	tempUnitCellCoords(4,1) = -0.50;
	tempUnitCellCoords(5,1) = 0.25;
	tempUnitCellCoords(6,1) = 0.50;
	tempUnitCellCoords(7,1) = 0.00;
	tempUnitCellCoords(8,1) = 0.50;

	tempUnitCellCoords(0,2) = 0.00;
	tempUnitCellCoords(1,2) = -0.25;
	tempUnitCellCoords(2,2) = -0.50;
	tempUnitCellCoords(3,2) = -0.50;
	tempUnitCellCoords(4,2) = 0.00;
	tempUnitCellCoords(5,2) = -0.25;
	tempUnitCellCoords(6,2) = 0.00;
	tempUnitCellCoords(7,2) = -0.50;
	tempUnitCellCoords(8,2) = -0.50;
	
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

	/* Calculate atomic volume = a^{3}/8:  8 atoms per unit volume */
	double asdf = f_a0*f_a0*f_a0/8.0;
	f_omega0 = 1.0/asdf;

	/* Area calculated from RamstadPRB1995, Figure 1 */
	double asdf2 = f_a0*f_a0/2.0;
	f_area0 = 1.0/asdf2;
	
	/* Calculate surface thickness */
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
	
	/* ROTATE UNIT CELL COORDS DEPENDING ON UNIT NORMAL */
	dMatrixT blah1(3);
	dArrayT temp(9);
	
	if (fNormalCode == 0)	// rotate [0,0,1] to [1,0,0]
	{
		/* CHECKED HSP 7/12/07 - FINE */
		fUnitCellCoords(0,0) = 0.00;
		fUnitCellCoords(1,0) = -0.25;
		fUnitCellCoords(2,0) = -0.50;
		fUnitCellCoords(3,0) = -0.50;
		fUnitCellCoords(4,0) = 0.00;
		fUnitCellCoords(5,0) = -0.25;
		fUnitCellCoords(6,0) = 0.00;
		fUnitCellCoords(7,0) = -0.50;
		fUnitCellCoords(8,0) = -0.50;

		fUnitCellCoords(0,1) = 0.00;
		fUnitCellCoords(1,1) = -0.25;
		fUnitCellCoords(2,1) = 0.00;
		fUnitCellCoords(3,1) = -0.50;
		fUnitCellCoords(4,1) = -0.50;
		fUnitCellCoords(5,1) = 0.25;
		fUnitCellCoords(6,1) = 0.50;
		fUnitCellCoords(7,1) = 0.50;
		fUnitCellCoords(8,1) = 0.00;

		fUnitCellCoords(0,2) = 0.00;
		fUnitCellCoords(1,2) = 0.25;
		fUnitCellCoords(2,2) = 0.50;
		fUnitCellCoords(3,2) = 0.00;
		fUnitCellCoords(4,2) = 0.50;
		fUnitCellCoords(5,2) = -0.25;
		fUnitCellCoords(6,2) = -0.50;
		fUnitCellCoords(7,2) = 0.00;
		fUnitCellCoords(8,2) = -0.50;
	}
	else if (fNormalCode == 1) // rotate [0,0,1] to [-1,0,0]
	{
		/* CHECKED HSP 7/12/07 - FINE */
		fUnitCellCoords(0,0) = 0.00;
		fUnitCellCoords(1,0) = 0.25;
		fUnitCellCoords(2,0) = 0.50;
		fUnitCellCoords(3,0) = 0.50;
		fUnitCellCoords(4,0) = 0.00;
		fUnitCellCoords(5,0) = 0.25;
		fUnitCellCoords(6,0) = 0.00;
		fUnitCellCoords(7,0) = 0.50;
		fUnitCellCoords(8,0) = 0.50;

		fUnitCellCoords(0,1) = 0.00;
		fUnitCellCoords(1,1) = 0.25;
		fUnitCellCoords(2,1) = 0.00;
		fUnitCellCoords(3,1) = 0.50;
		fUnitCellCoords(4,1) = 0.50;
		fUnitCellCoords(5,1) = -0.25;
		fUnitCellCoords(6,1) = -0.50;
		fUnitCellCoords(7,1) = -0.50;
		fUnitCellCoords(8,1) = 0.00;

		fUnitCellCoords(0,2) = 0.00;
		fUnitCellCoords(1,2) = 0.25;
		fUnitCellCoords(2,2) = 0.50;
		fUnitCellCoords(3,2) = 0.00;
		fUnitCellCoords(4,2) = 0.50;
		fUnitCellCoords(5,2) = -0.25;
		fUnitCellCoords(6,2) = -0.50;
		fUnitCellCoords(7,2) = 0.00;
		fUnitCellCoords(8,2) = -0.50;
	}
	else if (fNormalCode == 2)	// rotate [0,0,1] to [0,1,0]
	{
		/* CHECKED HSP 7/12/07 - FINE */
		temp_bonds2 = tempUnitCellCoords;
		blah1 = RotationMatrixB(-piby2);
		fUnitCellCoords.MultAB(temp_bonds2,blah1);
		fUnitCellCoords.CopyColumn(2,temp);
		temp *= -1.0;
		fUnitCellCoords.SetCol(2,temp);
	}
	else if (fNormalCode == 3)	// rotate [0,0,1] to [0,-1,0]
	{
		/* CHECKED HSP 7/12/07 - FINE */
		temp_bonds2 = tempUnitCellCoords;
		blah1 = RotationMatrixB(piby2);
		fUnitCellCoords.MultAB(temp_bonds2,blah1);
		fUnitCellCoords.CopyColumn(0,temp);
		temp *= -1.0;
		fUnitCellCoords.SetCol(0,temp);
	}
	else if (fNormalCode == 4)	// this is the default orientation
	{
		fUnitCellCoords = tempUnitCellCoords;
	}	
	else if (fNormalCode == 5)	// rotate [0,0,1] to [0,0,-1]
	{
		/* HSP - KEEP Y CONSTANT */
		temp_bonds2 = tempUnitCellCoords;
		fUnitCellCoords = temp_bonds2;
		fUnitCellCoords *= -1.0;
		fUnitCellCoords.CopyColumn(1,temp);
		temp *= -1.0;
		fUnitCellCoords.SetCol(1,temp);
	}	

	/* scale to correct lattice parameter */				     		
	fUnitCellCoords *= f_a0;
}

double TersoffSolverT_surf::Density(void) const
{
	return 8.0*fMass*1.0365e-4/(f_a0*f_a0*f_a0);
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

 /* Rotate bonds with [0,0,1] normal to bonds with [1,0,0]-type normals */
 /* Use positive piby2 to go to [-1,0,0], -piby2 to go to [1,0,0] */
dMatrixT TersoffSolverT_surf::RotationMatrixA(const double angle)
 {
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = cos(angle);
	rmatrix(0,2) = sin(angle);
	rmatrix(1,1) = 1.0;
	rmatrix(2,0) = -sin(angle);
	rmatrix(2,2) = cos(angle);

	return rmatrix;
 }
 
/* Rotate bonds with [0,0,1] normal to bonds with [0,1,0]-type normals */
/* Use positive piby2 to go to [0,-1,0], -piby2 to go to [0,1,0] */
dMatrixT TersoffSolverT_surf::RotationMatrixB(const double angle)
{
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = 1.0;
	rmatrix(1,1) = cos(angle);
	rmatrix(1,2) = sin(angle);
	rmatrix(2,1) = -sin(angle);
	rmatrix(2,2) = cos(angle);
	
	return rmatrix;
}

