/* $Id: TersoffSolverT.cpp,v 1.5 2007-06-24 18:11:17 paklein Exp $ */
#include "TersoffSolverT.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"
#include "Tersoff_inc.h"

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
	f_S(0.0)
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
#if 0
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	/* Compute all needed derivatives */
	SetAll(CIJ);

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
#endif
	ExceptionT::GeneralFail("TersoffSolverT::SetModuli", "not implemented");
}

//for now return symmetric matrix
void TersoffSolverT::SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress)
{
	/* set internal equilibrium */
	if (fEquilibrate) Equilibrate(CIJ, Xsi);

	/* call C function */
	get_dUdC(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(),
		stress.Pointer()); 
	stress *= 2.0;
	
#if 0
	cout << "stress: " << stress.no_wrap() << endl;
#endif
}

/* strain energy density */
double TersoffSolverT::StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi)
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
}

/* information about subordinate parameter lists */
//void TersoffSolverT::DefineSubs(SubListT& sub_list) const
//{
	/* inherited */
//	ParameterInterfaceT::DefineSubs(sub_list);

	/* crystal orientation */
//	sub_list.AddSub("FCC_lattice_orientation", ParameterListT::Once, true);

	/* choice of potentials */
//	sub_list.AddSub("DC_potential_choice", ParameterListT::Once, true);
//}

/* a pointer to the ParameterInterfaceT of the given subordinate */
//ParameterInterfaceT* TersoffSolverT::NewSub(const StringT& name) const
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
void TersoffSolverT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* dimension work space */
	dXsi.Dimension(kNumDOF);
	dXsidXsi.Dimension(kNumDOF);
#if 0
	dCdC_hat.Dimension(kStressDim);
	dCdXsi_hat.Dimension(kStressDim,kNumDOF);
	fMatrices.Dimension(kNumDOF);
	fMat2.Dimension(kNumDOF);
	fGradl_i.Dimension(3,kNumDOF); 
	fSymMat1.Dimension(kNSD);
	fTempRank4.Dimension(kStressDim);
	fTempMixed.Dimension(kStressDim, kNumDOF);
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
void TersoffSolverT::Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi)
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
	if (count == 15) ExceptionT::GeneralFail("TersoffSolverT::Equilibrate", "failed");
}

/* set free dof - triggers recomputation */
void TersoffSolverT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* call C function */
	get_dXsi(fParams.Pointer(), Xsi.Pointer(), 
		fUnitCellCoords(0), fUnitCellCoords(1), fUnitCellCoords(2), 
		CIJ.Pointer(), 
		dXsi.Pointer(), dXsidXsi.Pointer());

//debugging
#if 0
cout << dXsi.no_wrap() << ":" << dXsidXsi.no_wrap() << endl;
#endif
}

/* set free dof - triggers recomputation */
void TersoffSolverT::SetAll(const dMatrixT& CIJ)
{
	/* set geometry */
//	fGeometry->SetAll(CIJ);

#if 0	
	/* Initialize */
	dCdC_hat   = 0.0;
	dCdXsi_hat = 0.0;

	/* scalar derivatives */
	const dArray2DT& dl_dXsi      = fGeometry->dl_dXsi();
	const dArray2DT& d2l_dXsidXsi = fGeometry->d2l_dXsidXsi();

	const dArray2DT& dl_dC        = fGeometry->dl_hat_dC();
	const dArray2DT& d2l_dCdC     = fGeometry->d2l_hat_dCdC();
	const dArray2DT& d2l_dCdXsi   = fGeometry->d2l_hat_dCdXsi();

	const dArray2DT& dc_dXsi      = fGeometry->dCos_dXsi();
	const dArray2DT& d2c_dXsidXsi = fGeometry->d2Cos_dXsidXsi();

	const dArray2DT& dc_dC        = fGeometry->dCos_hat_dC();
	const dArray2DT& d2c_dCdC     = fGeometry->d2Cos_hat_dCdC();
	const dArray2DT& d2c_dCdXsi   = fGeometry->d2Cos_hat_dCdXsi();
		
	/* shallow work temps */
	dMatrixT	d2ldCdC, dldC;
	dArrayT		dldXsi;
	dMatrixT	d2ldCdXsi;
#endif

	/* 2-body derivatives */
// 	const dArrayT& dPhi_2  = f2Body->dPhi();
// 	const dArrayT& ddPhi_2 = f2Body->ddPhi();	
// 	for (int i = 0 ; i < dPhi_2.Length(); i++)
// 	{
// 		/* d2/dCdC */
// 		dldC.Alias(kNSD, kNSD, dl_dC(i));
// 		fSymMat1.FromMatrix(dldC);
// 		fTempRank4.Outer(fSymMat1,fSymMat1);
// 
// 		d2ldCdC.Alias(kStressDim, kStressDim, d2l_dCdC(i));		
// 
// 		dCdC_hat.AddCombination(dPhi_2[i], d2ldCdC, ddPhi_2[i], fTempRank4);
// 	
// 		/* d2/dCdXsi */
// 		dldXsi.Alias(kNumDOF, dl_dXsi(i));
// 		fTempMixed.Outer(fSymMat1,dldXsi);
// 		
// 		d2ldCdXsi.Alias(kStressDim, kNumDOF, d2l_dCdXsi(i));
// 			
// 		dCdXsi_hat.AddCombination(dPhi_2[i], d2ldCdXsi, ddPhi_2[i], fTempMixed);
// 	}
// 
// 	/* for the linear combinations */
// 	dArrayT  coeffs;
// 	dMatrixT ddl1, ddl2, ddc12;
// 	fMatrices[0] = &ddl1;
// 	fMatrices[1] = &ddl2;
// 	fMatrices[2] = &ddc12;
// 	
// 	/* shallow temps */
// 	dMatrixT ddPhi3;
// 	dArrayT	 dl1dXsi, dl2dXsi, dCosdXsi;
// 	dMatrixT dl1dC  , dl2dC  , dCosdC;
// 
// 	/* 3-body derivatives */
// 	const dArray2DT& dPhi_3  = f3Body->dPhi();
// 	const dArray2DT& ddPhi_3 = f3Body->ddPhi();
// 	for (int j = 0 ; j < dPhi_3.MajorDim(); j++)
// 	{
// 		int n1 = fPairs(j,0);
// 		int n2 = fPairs(j,1);
// 
// 		coeffs.Alias(kNumDOF, dPhi_3(j));
// 	
// 		/* d2/dCdC */
// 		ddPhi3.Alias(kNumDOF, kNumDOF, ddPhi_3(j));
// 		
// 		dl1dC.Alias(kNSD, kNSD, dl_dC(n1));		
// 		dl2dC.Alias(kNSD, kNSD, dl_dC(n2));
// 		dCosdC.Alias(kNSD, kNSD, dc_dC(j));
// 	
// 		fSymMat1.FromMatrix(dl1dC);
// 		fGradl_C.SetRow(0, fSymMat1);
// 		fSymMat1.FromMatrix(dl2dC);
// 		fGradl_C.SetRow(1, fSymMat1 );
// 		fSymMat1.FromMatrix(dCosdC);
// 		fGradl_C.SetRow(2, fSymMat1);
// 
// 		fTempMixed.MultATB(fGradl_C,ddPhi3);
// 		fTempRank4.MultAB(fTempMixed,fGradl_C);
// 		
// 		//testing
// 		//dCdC_hat += fTempRank4;
// 		
// 		ddl1.Alias(kStressDim, kStressDim, d2l_dCdC(n1));
// 		ddl2.Alias(kStressDim, kStressDim, d2l_dCdC(n2));
// 		ddc12.Alias(kStressDim, kStressDim, d2c_dCdC(j));
// 		
// 		//testing
// 		//dCdC_hat.AddCombination(coeffs,fMatrices);
// 		
// 		dCdC_hat.AddCombination(1.0,fTempRank4,coeffs[0],ddl1);
// 		dCdC_hat.AddCombination(coeffs[1],ddl2,
// 		                        coeffs[2],ddc12);		
// 				
// 		/* d2/dCdXsi */
// 		dl1dXsi.Alias(kNumDOF, dl_dXsi(n1));
// 		dl2dXsi.Alias(kNumDOF, dl_dXsi(n2));
// 		dCosdXsi.Alias(kNumDOF, dc_dXsi(j));
// 		
// 		fGradl_i.SetRow(0, dl1dXsi);
// 		fGradl_i.SetRow(1, dl2dXsi );
// 		fGradl_i.SetRow(2, dCosdXsi);
// 
// 		fMat1.MultAB(ddPhi3,fGradl_i);
// 		fTempMixed.MultATB(fGradl_C,fMat1);
// 	
// 		//testing
// 		//dCdXsi_hat += fTempMixed;
// 		
// 		ddl1.Alias(kStressDim , kNumDOF, d2l_dCdXsi(n1));
// 		ddl2.Alias(kStressDim , kNumDOF, d2l_dCdXsi(n2));
// 		ddc12.Alias(kStressDim, kNumDOF, d2c_dCdXsi(j));
// 		
// 		//testing
// 		//dCdC_hat.AddCombination(coeffs,fMatrices);
// 		
// 		dCdXsi_hat.AddCombination(1.0, fTempMixed, coeffs[0],ddl1);
// 		dCdXsi_hat.AddCombination(coeffs[1],ddl2,
// 		                          coeffs[2],ddc12);
// 	}
}
