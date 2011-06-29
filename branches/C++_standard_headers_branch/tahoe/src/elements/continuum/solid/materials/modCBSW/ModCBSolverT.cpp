/* $Id: ModCBSolverT.cpp,v 1.6 2004-07-15 08:28:36 paklein Exp $ */
/* created: paklein (05/27/1997) */
#include "ModCBSolverT.h"

#include "dSymMatrixT.h"
#include "SW2BodyT.h"
#include "SW3BodyT.h"
#include "PTHT2BodyT.h"
#include "PTHT3BodyT.h"
#include "ParameterContainerT.h"
#include "FCCLatticeT.h" /* needed for lattice orientation */

using namespace Tahoe;

const int kNSD       = 3;
const int kNumDOF    = 3;
const int kNumAngles = 12;
const int kStressDim =dSymMatrixT::NumValues(kNSD);

/* set pair numbers */
int pairdata[kNumAngles*2] =
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
ModCBSolverT::ModCBSolverT(const ThermalDilatationT* thermal):
	ParameterInterfaceT("mod_Cauchy-Born_solver"),
	fEquilibrate(true),
	fThermal(thermal),
	fPairs(kNumAngles, 2, pairdata),
	fGeometry(NULL),
	f2Body(NULL),
	f3Body(NULL)
{

}

/* Destructor */
ModCBSolverT::~ModCBSolverT(void)
{
	delete f2Body;
	delete f3Body;
	delete fGeometry;
}

/* moduli - assume Xsi already determined */
void ModCBSolverT::SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli)
{
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
}

//for now return symmetric matrix
void ModCBSolverT::SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	/* Compute all needed derivatives */
	fGeometry->SetdC(CIJ);

	/* Initialize stress */
	stress = 0.0;

	/* scalar derivatives */
	const dArray2DT& dlh_dC   = fGeometry->dl_hat_dC();
	const dArray2DT& dCosh_dC = fGeometry->dCos_hat_dC();

	/* shallow work temps */
	dMatrixT dl1hdC, dl2hdC, dCoshdC;

	/* 2-body derivatives */
	const dArrayT& dPhi_2 = f2Body->dPhi();
	for (int i = 0 ; i < dPhi_2.Length(); i++)
	{
		/* stress */
		dl1hdC.Alias(kNSD, kNSD, dlh_dC(i));

		stress.AddScaled(dPhi_2[i], dl1hdC);	
	}

	/* for the linear combinations */
	dArrayT  coeffs;

	fMatrices[0] = &dl1hdC;
	fMatrices[1] = &dl2hdC;
	fMatrices[2] = &dCoshdC;

	/* 3-body derivatives */
	const dArray2DT& dPhi_3  = f3Body->dPhi();
	for (int j = 0 ; j < dPhi_3.MajorDim(); j++)
	{
		int n1 = fPairs(j,0);
		int n2 = fPairs(j,1);

		coeffs.Alias(kNumDOF, dPhi_3(j));
	
		/* stress */
		dl1hdC.Alias(kNSD, kNSD, dlh_dC(n1));
		dl2hdC.Alias(kNSD, kNSD, dlh_dC(n2));
		dCoshdC.Alias(kNSD, kNSD, dCosh_dC(j));
	
		stress.AddScaled(coeffs[0],dl1hdC);
		stress.AddCombination(coeffs[1],dl2hdC,
		                      coeffs[2],dCoshdC);
	}
	
	/* factor of 2 to get to PK2 */
	stress *= 2.0;
}

/* strain energy density */
double ModCBSolverT::StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi)
{
	/* set internal equilibrium */
	if (fEquilibrate)
		Equilibrate(CIJ, Xsi);
	else
		SetdXsi(CIJ, Xsi);

	return( (f2Body->Phi()).Sum() + (f3Body->Phi()).Sum() );
}

/* describe the parameters needed by the interface */
void ModCBSolverT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT equilibrate(ParameterT::Boolean, "equilibrate");
	equilibrate.SetDefault(true);
	list.AddParameter(equilibrate);
}

/* information about subordinate parameter lists */
void ModCBSolverT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* crystal orientation */
	sub_list.AddSub("FCC_lattice_orientation", ParameterListT::Once, true);

	/* choice of potentials */
	sub_list.AddSub("DC_potential_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ModCBSolverT::NewSub(const StringT& name) const
{
	if (name == "DC_potential_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetSubSource(this);
		choice->SetListOrder(ParameterListT::Choice);
	
		choice->AddSub("Stillinger-Weber");

		ParameterContainerT PTHT("PTHT");
		PTHT.AddParameter(ParameterT::Double, "A");
		PTHT.AddParameter(ParameterT::Double, "A1");
		PTHT.AddParameter(ParameterT::Double, "A2");
		
		PTHT.AddParameter(ParameterT::Double, "B");
		PTHT.AddParameter(ParameterT::Double, "Z");
		choice->AddSub(PTHT);

		//choice->AddSub(ParameterContainerT("Tersoff"));

		return choice;
	}
	else if (name == "FCC_lattice_orientation")
	{
		FCCLatticeT lattice(0);
		return lattice.NewSub(name);
	}
	else if (name == "Stillinger-Weber")
		return new SWDataT;
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void ModCBSolverT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* dimension work space */
	dXsi.Dimension(kNumDOF);
	dXsidXsi.Dimension(kNumDOF);
	dCdC_hat.Dimension(kStressDim);
	dCdXsi_hat.Dimension(kStressDim,kNumDOF);
	fMatrices.Dimension(kNumDOF);
	fMat1.Dimension(kNumDOF); 
	fMat2.Dimension(kNumDOF);
	fGradl_i.Dimension(3,kNumDOF); 
	fVec.Dimension(kNumDOF);
	fSymMat1.Dimension(kNSD);
	fTempRank4.Dimension(kStressDim);
	fTempMixed.Dimension(kStressDim, kNumDOF);
	fGradl_C.Dimension(3,kStressDim);

	/* flag */
	fEquilibrate = list.GetParameter("equilibrate");

	/* resolve orientation */
	FCCLatticeT lattice(0);
	const ParameterListT& orientation = list.GetListChoice(lattice, "FCC_lattice_orientation");
	dMatrixT Q;
	FCCLatticeT::SetQ(orientation, Q);
	
	/* construct bond lattice */
	fGeometry = new LengthsAndAnglesT(Q,fPairs);

	/* set potentials */
	const ParameterListT& potential = list.GetListChoice(*this, "DC_potential_choice");
	if (potential.Name() == "Stillinger-Weber")
	{
		/* extract parameters */
		fSW.TakeParameterList(potential);	
		
		/* construct potentials */
		f2Body = new SW2BodyT(fGeometry->Lengths(), fThermal, fSW);
		f3Body = new SW3BodyT(fGeometry->Lengths(), fGeometry->Cosines(), fPairs, fThermal, fSW);
	}
	else if (potential.Name() == "Stillinger-Weber")
	{
		/* extract parameters */
		double A = list.GetParameter("A");
		double A1 = list.GetParameter("A1");
		double A2 = list.GetParameter("A2");
		double B = list.GetParameter("B");
		double Z = list.GetParameter("Z");	

		/* construct potentials */
		f2Body = new PTHT2BodyT(fGeometry->Lengths(), fThermal, A, A1, A2);
		f3Body = new PTHT3BodyT(fGeometry->Lengths(), fGeometry->Cosines(), fPairs, fThermal, B, Z);
	}
	else
		ExceptionT::BadInputValue("ModCBSolverT::TakeParameterList",
			"unknown potential \"%s\"", potential.Name().Pointer());
}

/**********************************************************************
 * Private
 **********************************************************************/

/* Minimize the energy wrt Xsi using the initial value passed */
void ModCBSolverT::Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi)
{
	/* check initial value */
	SetdXsi(CIJ, Xsi);

	int count = 0;	
	while (count++ < 15 && dXsi.Magnitude() > 1.0e-12)
	{
		fMat1.Inverse(dXsidXsi);
		fMat1.Multx(dXsi, fVec);
		
		Xsi -= fVec;
		
		/* recompute */
		SetdXsi(CIJ, Xsi);
	}

	/* assume not converged */
	if (count == 15) ExceptionT::GeneralFail("ModCBSolverT::Equilibrate", "failed");
}

/* set free dof - triggers recomputation */
void ModCBSolverT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* set geometry */
	fGeometry->SetdXsi(CIJ,Xsi);

	/* potentials and derivatives */
	f2Body->Set();
	f3Body->Set();
	
	/* initialize */
	dXsi     = 0.0;
	dXsidXsi = 0.0;
		
	/* scalar derivatives */
	const dArray2DT& dl_dXsi      = fGeometry->dl_dXsi();
	const dArray2DT& d2l_dXsidXsi = fGeometry->d2l_dXsidXsi();

	const dArray2DT& dc_dXsi      = fGeometry->dCos_dXsi();
	const dArray2DT& d2c_dXsidXsi = fGeometry->d2Cos_dXsidXsi();
		
	/* shallow work temps */
	dArrayT dl1dXsi, dl2dXsi, dCosdXsi;
	dMatrixT d2ldXsidXsi;

	/* 2-body derivatives */
	const dArrayT& dPhi_2  = f2Body->dPhi();
	const dArrayT& ddPhi_2 = f2Body->ddPhi();	
	for (int i = 0 ; i < dPhi_2.Length(); i++)
	{
		/* gradient */
		dl1dXsi.Alias(kNumDOF, dl_dXsi(i));

		dXsi.AddScaled(dPhi_2[i], dl1dXsi);
	
		/* hessian */
		d2ldXsidXsi.Alias(kNumDOF,kNumDOF,d2l_dXsidXsi(i));
		fMat1.Outer(dl1dXsi,dl1dXsi);
	
		dXsidXsi.AddCombination(ddPhi_2[i], fMat1, dPhi_2[i], d2ldXsidXsi);
	}

	/* for the linear combinations */
	dArrayT  coeffs;
	dMatrixT ddl1, ddl2, ddc12;
	fMatrices[0] = &ddl1;
	fMatrices[1] = &ddl2;
	fMatrices[2] = &ddc12;
	
	/* shallow temps */
	dMatrixT ddPhi3;

	/* 3-body derivatives */
	const dArray2DT& dPhi_3  = f3Body->dPhi();
	const dArray2DT& ddPhi_3 = f3Body->ddPhi();
	for (int j = 0 ; j < dPhi_3.MajorDim(); j++)
	{
		int n1 = fPairs(j,0);
		int n2 = fPairs(j,1);

		coeffs.Alias(kNumDOF, dPhi_3(j));
	
		/* gradient */
		dl1dXsi.Alias(kNumDOF, dl_dXsi(n1));
		dl2dXsi.Alias(kNumDOF, dl_dXsi(n2));
		dCosdXsi.Alias(kNumDOF, dc_dXsi(j));
	
		fGradl_i.SetRow(0, dl1dXsi );
		fGradl_i.SetRow(1, dl2dXsi );
		fGradl_i.SetRow(2, dCosdXsi);

		fGradl_i.MultTx(coeffs, fVec);
		
		dXsi += fVec;
		
		/* hessian */
		ddPhi3.Alias(kNumDOF,kNumDOF,ddPhi_3(j));
		fMat1.MultATB(fGradl_i,ddPhi3);
		fMat2.MultAB(fMat1,fGradl_i);
	
		//testing
		//dXsidXsi += fMat2;
		
		ddl1.Alias(kNumDOF , kNumDOF, d2l_dXsidXsi(n1));
		ddl2.Alias(kNumDOF , kNumDOF, d2l_dXsidXsi(n2));
		ddc12.Alias(kNumDOF, kNumDOF, d2c_dXsidXsi(j));
		
		//testing
		//dXsidXsi.AddCombination(coeffs, fMatrices);
		
		dXsidXsi.AddCombination(1.0,fMat2, coeffs[0], ddl1);
		dXsidXsi.AddCombination(coeffs[1], ddl2,
		                        coeffs[2], ddc12);
	}
}

/* set free dof - triggers recomputation */
void ModCBSolverT::SetAll(const dMatrixT& CIJ)
{
	/* set geometry */
	fGeometry->SetAll(CIJ);
	
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

	/* 2-body derivatives */
	const dArrayT& dPhi_2  = f2Body->dPhi();
	const dArrayT& ddPhi_2 = f2Body->ddPhi();	
	for (int i = 0 ; i < dPhi_2.Length(); i++)
	{
		/* d2/dCdC */
		dldC.Alias(kNSD, kNSD, dl_dC(i));
		fSymMat1.FromMatrix(dldC);
		fTempRank4.Outer(fSymMat1,fSymMat1);

		d2ldCdC.Alias(kStressDim, kStressDim, d2l_dCdC(i));		

		dCdC_hat.AddCombination(dPhi_2[i], d2ldCdC, ddPhi_2[i], fTempRank4);
	
		/* d2/dCdXsi */
		dldXsi.Alias(kNumDOF, dl_dXsi(i));
		fTempMixed.Outer(fSymMat1,dldXsi);
		
		d2ldCdXsi.Alias(kStressDim, kNumDOF, d2l_dCdXsi(i));
			
		dCdXsi_hat.AddCombination(dPhi_2[i], d2ldCdXsi, ddPhi_2[i], fTempMixed);
	}

	/* for the linear combinations */
	dArrayT  coeffs;
	dMatrixT ddl1, ddl2, ddc12;
	fMatrices[0] = &ddl1;
	fMatrices[1] = &ddl2;
	fMatrices[2] = &ddc12;
	
	/* shallow temps */
	dMatrixT ddPhi3;
	dArrayT	 dl1dXsi, dl2dXsi, dCosdXsi;
	dMatrixT dl1dC  , dl2dC  , dCosdC;

	/* 3-body derivatives */
	const dArray2DT& dPhi_3  = f3Body->dPhi();
	const dArray2DT& ddPhi_3 = f3Body->ddPhi();
	for (int j = 0 ; j < dPhi_3.MajorDim(); j++)
	{
		int n1 = fPairs(j,0);
		int n2 = fPairs(j,1);

		coeffs.Alias(kNumDOF, dPhi_3(j));
	
		/* d2/dCdC */
		ddPhi3.Alias(kNumDOF, kNumDOF, ddPhi_3(j));
		
		dl1dC.Alias(kNSD, kNSD, dl_dC(n1));		
		dl2dC.Alias(kNSD, kNSD, dl_dC(n2));
		dCosdC.Alias(kNSD, kNSD, dc_dC(j));
	
		fSymMat1.FromMatrix(dl1dC);
		fGradl_C.SetRow(0, fSymMat1);
		fSymMat1.FromMatrix(dl2dC);
		fGradl_C.SetRow(1, fSymMat1 );
		fSymMat1.FromMatrix(dCosdC);
		fGradl_C.SetRow(2, fSymMat1);

		fTempMixed.MultATB(fGradl_C,ddPhi3);
		fTempRank4.MultAB(fTempMixed,fGradl_C);
		
		//testing
		//dCdC_hat += fTempRank4;
		
		ddl1.Alias(kStressDim, kStressDim, d2l_dCdC(n1));
		ddl2.Alias(kStressDim, kStressDim, d2l_dCdC(n2));
		ddc12.Alias(kStressDim, kStressDim, d2c_dCdC(j));
		
		//testing
		//dCdC_hat.AddCombination(coeffs,fMatrices);
		
		dCdC_hat.AddCombination(1.0,fTempRank4,coeffs[0],ddl1);
		dCdC_hat.AddCombination(coeffs[1],ddl2,
		                        coeffs[2],ddc12);		
				
		/* d2/dCdXsi */
		dl1dXsi.Alias(kNumDOF, dl_dXsi(n1));
		dl2dXsi.Alias(kNumDOF, dl_dXsi(n2));
		dCosdXsi.Alias(kNumDOF, dc_dXsi(j));
		
		fGradl_i.SetRow(0, dl1dXsi);
		fGradl_i.SetRow(1, dl2dXsi );
		fGradl_i.SetRow(2, dCosdXsi);

		fMat1.MultAB(ddPhi3,fGradl_i);
		fTempMixed.MultATB(fGradl_C,fMat1);
	
		//testing
		//dCdXsi_hat += fTempMixed;
		
		ddl1.Alias(kStressDim , kNumDOF, d2l_dCdXsi(n1));
		ddl2.Alias(kStressDim , kNumDOF, d2l_dCdXsi(n2));
		ddc12.Alias(kStressDim, kNumDOF, d2c_dCdXsi(j));
		
		//testing
		//dCdC_hat.AddCombination(coeffs,fMatrices);
		
		dCdXsi_hat.AddCombination(1.0, fTempMixed, coeffs[0],ddl1);
		dCdXsi_hat.AddCombination(coeffs[1],ddl2,
		                          coeffs[2],ddc12);
	}
}
