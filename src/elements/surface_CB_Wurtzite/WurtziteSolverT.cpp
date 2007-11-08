/* $Id: WurtziteSolverT.cpp,v 1.1 2007-11-08 19:37:46 hspark Exp $ */
#include "WurtziteSolverT.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"
#include "Wurtzite_inc.h"

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
WurtziteSolverT::WurtziteSolverT(const ThermalDilatationT* thermal):
	ParameterInterfaceT("Wurtzite_CB_solver"),
	fEquilibrate(true),
	fThermal(thermal),
	fPairs(kNumAngles, 2, pairdata),
//	fGeometry(NULL),
	f_Caxis(0.0),
	f_Aaxis(0,0),
	f_D0(0.0),
	f_S0(0.0),
	f_r0(0.0),
	f_beta(0.0),
	f_gamma(0.0),
	f_c(0.0),
	f_d(0.0),
	f_h(0.0),
	f_R(0.0),
	f_D(0.0),
	f_omega0(0.0)
{

}

/* Destructor */
WurtziteSolverT::~WurtziteSolverT(void)
{
//	delete fGeometry;
}

/* moduli - assume Xsi already determined */
void WurtziteSolverT::SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	/* compute second derivatives wrt {C,C} and {C,Xsi} */
	get_ddC(fParams.Pointer(), Xsi.Pointer(), 
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
void WurtziteSolverT::SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress)
{
	/* set internal equilibrium */
	if (fEquilibrate) Equilibrate(CIJ, Xsi);

	/* call C function */
	get_dUdC(fParams.Pointer(), Xsi.Pointer(), 
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
double WurtziteSolverT::StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi)
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
void WurtziteSolverT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT equilibrate(ParameterT::Boolean, "equilibrate");
	equilibrate.SetDefault(true);
	list.AddParameter(equilibrate);

	/* lattice parameters */
	ParameterT Caxis(f_Caxis, "Caxis");
	Caxis.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(Caxis);

	ParameterT Aaxis(f_Aaxis, "Aaxis");
	Aaxis.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(Aaxis);

	/* Revised from WurtzitePairT.cpp */
	ParameterT mass(fMass, "mass");
	mass.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mass);

	ParameterT D0(f_D0, "energy_scale_d0");
	D0.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(D0);
	
	ParameterT S0(f_S0, "energy_scale_s0");
	S0.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(S0);
	
	ParameterT r0(f_r0, "dimer_length_r0");
	r0.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(r0);
	
	ParameterT beta(f_beta, "bond_order_coeff1_betai");
	beta.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(beta);
	
	ParameterT gamma(f_gamma, "bond_order_coeff0_gamma");
	gamma.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(gamma);
	
	ParameterT c(f_c, "bond_order_coeff2_ci");
	c.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c);
	
	ParameterT d(f_d, "bond_order_coeff3_di");
	d.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(d);
	
	ParameterT h(f_h, "bond_order_coeff4_hi");
	list.AddParameter(h);
	
	ParameterT R(f_R, "cutoff_func_length_1_Rij");
	R.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(R);
	
	ParameterT D(f_D, "cutoff_func_length_2_Dij");
	D.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(D);
}

/* information about subordinate parameter lists */
//void WurtziteSolverT::DefineSubs(SubListT& sub_list) const
//{
	/* inherited */
//	ParameterInterfaceT::DefineSubs(sub_list);

	/* crystal orientation */
//	sub_list.AddSub("FCC_lattice_orientation", ParameterListT::Once, true);

	/* choice of potentials */
//	sub_list.AddSub("DC_potential_choice", ParameterListT::Once, true);
//}

/* a pointer to the ParameterInterfaceT of the given subordinate */
//ParameterInterfaceT* WurtziteSolverT::NewSub(const StringT& name) const
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
// 		//choice->AddSub(ParameterContainerT("Wurtzite"));
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
void WurtziteSolverT::TakeParameterList(const ParameterListT& list)
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
	f_Caxis = list.GetParameter("Caxis");
	f_Aaxis = list.GetParameter("Aaxis");
	fMass = list.GetParameter("mass");
	f_D0 = list.GetParameter("energy_scale_d0");
	f_S0 = list.GetParameter("energy_scale_s0");
	f_r0 = list.GetParameter("dimer_length_r0");
	f_beta = list.GetParameter("bond_order_coeff1_betai");
	f_gamma = list.GetParameter("bond_order_coeff0_gamma");
	f_c = list.GetParameter("bond_order_coeff2_ci");
	f_d = list.GetParameter("bond_order_coeff3_di");
	f_h = list.GetParameter("bond_order_coeff4_hi");
	f_R = list.GetParameter("cutoff_func_length_1_Rij");
	f_D = list.GetParameter("cutoff_func_length_2_Dij");
	
	/* scale unit cell coordinates */
	fUnitCellCoords *= f_a0;

	/* Calculate atomic volume - FIX FOR WURTZITE! */
	double asdf = f_a0*f_a0*f_a0/8.0;
	f_omega0 = 1.0/asdf;

	/* write into vector to pass to C code */
	fParams.Dimension(13);
	fParams[ 0] = f_Caxis;
	fParams[ 1] = f_Aaxis;
	fParams[ 2] = fMass;
	fParams[ 3] = f_D0;
	fParams[ 4] = f_S0;
	fParams[ 5] = f_r0;
	fParams[ 6] = f_beta;
	fParams[ 7] = f_gamma;
	fParams[ 8] = f_c;
	fParams[ 9] = f_d;
	fParams[10] = f_h;
	fParams[11] = f_R;
	fParams[12] = f_D;	
}

/**********************************************************************
 * Private
 **********************************************************************/

/* Minimize the energy wrt Xsi using the initial value passed */
void WurtziteSolverT::Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi)
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
	if (count == 15) ExceptionT::GeneralFail("WurtziteSolverT::Equilibrate", "failed");
}

/* set free dof - triggers recomputation */
void WurtziteSolverT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
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
