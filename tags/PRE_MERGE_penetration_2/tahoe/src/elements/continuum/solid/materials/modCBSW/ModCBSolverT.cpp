/* $Id: ModCBSolverT.cpp,v 1.5 2004-06-17 07:41:03 paklein Exp $ */
/* created: paklein (05/27/1997)                                          */
/* Q defines the orientation of the crystals' natural coordinates         */
/* and the global coordinate frame. Q is defined as:                      */
/* 			Q = d x_natural / d x_global                                        */
/* So that the vectors are transformed by:                                */
/* 			r_global = Transpose[Q].r_natural                                   */

#include "ModCBSolverT.h"

#include <iostream.h>

#include "ExceptionT.h"

#include "ifstreamT.h"
#include "dSymMatrixT.h"
#include "SW2BodyT.h"
#include "SW3BodyT.h"
#include "PTHT2BodyT.h"
#include "PTHT3BodyT.h"


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
ModCBSolverT::ModCBSolverT(const dMatrixT& Q,
	const ThermalDilatationT* thermal, ifstreamT& in, bool equilibrate):
	fEquilibrate(equilibrate),
	fPairs(kNumAngles,2,pairdata),
	fGeometry(Q,fPairs),
	dXsi(kNumDOF),
	dXsidXsi(kNumDOF),
	dCdC_hat(kStressDim),
	dCdXsi_hat(kStressDim,kNumDOF),
	fMatrices(kNumDOF),
	fMat1(kNumDOF), fMat2(kNumDOF),
	fGradl_i(3,kNumDOF), fVec(kNumDOF),
	fSymMat1(kNSD),
	fTempRank4(kStressDim),
	fTempMixed(kStressDim, kNumDOF),
	fGradl_C(3,kStressDim)
{
	/* check */
	if (fEquilibrate != 0 && fEquilibrate != 1)
		throw ExceptionT::kBadInputValue;

	/* set potentials */
	in >> fPotential;
	switch (fPotential)
	{
		case kSW:
		
			/* read SW data */
			fSW.Read(in);
		
			f2Body = new SW2BodyT(fGeometry.Lengths(), thermal, fSW);
			f3Body = new SW3BodyT(fGeometry.Lengths(), fGeometry.Cosines(),
							fPairs, thermal, fSW);
			break;
			
		case kPTHT:
		
			f2Body = new PTHT2BodyT(fGeometry.Lengths(), thermal, in);
			f3Body = new PTHT3BodyT(fGeometry.Lengths(), fGeometry.Cosines(),
							fPairs, thermal, in);
			break;
			
		case kTersoff:
		
			throw ExceptionT::kGeneralFail; //not yet implemented
			break;
			
		default:

			cout << "\nModCBSolverT::ModCBSolverT: unknown potential code:";
			cout << fPotential << endl;
			throw ExceptionT::kBadInputValue;
	}

	if (!f2Body || !f3Body) throw ExceptionT::kOutOfMemory;
}

/* Destructor */
ModCBSolverT::~ModCBSolverT(void)
{
	delete f2Body;
	delete f3Body;
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
	fGeometry.SetdC(CIJ);

	/* Initialize stress */
	stress = 0.0;

	/* scalar derivatives */
	const dArray2DT& dlh_dC      = fGeometry.dl_hat_dC();
	const dArray2DT& dCosh_dC    = fGeometry.dCos_hat_dC();

	/* shallow work temps */
	dMatrixT dl1hdC, dl2hdC, dCoshdC;

	/* 2-body derivatives */
	const dArrayT& dPhi_2  = f2Body->dPhi();
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

/*
* Print parameters.
*/
void ModCBSolverT::Print(ostream& out) const
{
	/* print potential data */
	if (fPotential == kSW) fSW.Write(out);

	//printing not implemented for other potentials

	out << " Number of internal DOF. . . . . . . . . . . . . = ";
	out << ((fEquilibrate == 1) ? dXsi.Length() : 0) << '\n';
}

void ModCBSolverT::PrintName(ostream& out) const
{
	const char* potentials[] = {"Stillinger-Weber",
	                            "PTHT",
	                            "Tersoff"};
	
	out << "    " << potentials[fPotential] << '\n';
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
	if (count == 15)
	{
		cout << "\n ModCBSolverT::Equilibrate: could not find internal equilibrium" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* set free dof - triggers recomputation */
void ModCBSolverT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* set geometry */
	fGeometry.SetdXsi(CIJ,Xsi);

	/* potentials and derivatives */
	f2Body->Set();
	f3Body->Set();
	
	/* initialize */
	dXsi     = 0.0;
	dXsidXsi = 0.0;
		
	/* scalar derivatives */
	const dArray2DT& dl_dXsi      = fGeometry.dl_dXsi();
	const dArray2DT& d2l_dXsidXsi = fGeometry.d2l_dXsidXsi();

	const dArray2DT& dc_dXsi      = fGeometry.dCos_dXsi();
	const dArray2DT& d2c_dXsidXsi = fGeometry.d2Cos_dXsidXsi();
		
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
	fGeometry.SetAll(CIJ);
	
	/* Initialize */
	dCdC_hat   = 0.0;
	dCdXsi_hat = 0.0;
	
		/* scalar derivatives */
	const dArray2DT& dl_dXsi      = fGeometry.dl_dXsi();
	const dArray2DT& d2l_dXsidXsi = fGeometry.d2l_dXsidXsi();

	const dArray2DT& dl_dC        = fGeometry.dl_hat_dC();
	const dArray2DT& d2l_dCdC     = fGeometry.d2l_hat_dCdC();
	const dArray2DT& d2l_dCdXsi   = fGeometry.d2l_hat_dCdXsi();

	const dArray2DT& dc_dXsi      = fGeometry.dCos_dXsi();
	const dArray2DT& d2c_dXsidXsi = fGeometry.d2Cos_dXsidXsi();

	const dArray2DT& dc_dC        = fGeometry.dCos_hat_dC();
	const dArray2DT& d2c_dCdC     = fGeometry.d2Cos_hat_dCdC();
	const dArray2DT& d2c_dCdXsi   = fGeometry.d2Cos_hat_dCdXsi();
		
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
