/* $Id: TotalLagrangianT.cpp,v 1.3 2001-07-10 07:29:54 paklein Exp $ */
/* created: paklein (09/07/1998)                                          */

#include "TotalLagrangianT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "StructuralMaterialT.h"
#include "MaterialListT.h"
#include "ShapeFunctionT.h"

/* constructor */
TotalLagrangianT::TotalLagrangianT(FEManagerT& fe_manager):
	FiniteStrainT(fe_manager),
	fStressMat(fNumSD),
	fTempMat1(fNumSD),
	fTempMat2(fNumSD)
{
	/* disable any strain-displacement options */
	if (fStrainDispOpt != 0)
	{
		cout << "\n TotalLagrangianT::TotalLagrangianT: no strain-displacement options\n" << endl;
		fStrainDispOpt = 0;
	}

}

/* data initialization */
void TotalLagrangianT::Initialize(void)
{
	/* inherited */
	FiniteStrainT::Initialize();

//TEMP
	if (fMaterialList->HasThermalStrains())
	{
		cout << "\n TotalLagrangianT::Initialize: element group ";
		cout << fFEManager.ElementGroupNumber(this) + 1 << ": total Lagrangian\n";
		cout << "     formulation has a bug associated with the multiplicative treatment\n";
		cout << "     of thermal strains. The solutions appear correct, but the convergence\n";
		cout << "     is not quadratic. Use updated Lagrangian formulation." << endl;
		throw eGeneralFail;
	}	

	/* dimension */
	fGradNa.Allocate(fNumSD, fNumElemNodes);
	fStressStiff.Allocate(fNumElemNodes);
	fTemp2.Allocate(fNumElemNodes*fNumDOF);
}

/***********************************************************************
* Protected
***********************************************************************/

/* form the element stiffness matrix */
void TotalLagrangianT::FormStiffness(double constK)
{		
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
	/* S T R E S S   S T I F F N E S S */
				
		/* PK2 stress (and set deformation gradient) */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
	
		/* change rule shape function derivatives */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();
		fTempMat1.Inverse();
		fShapes->SetChainRule(fTempMat1, fDNa_x);

		/* get shape function gradients matrix */
		fShapes->GradNa(fDNa_x, fGradNa);

		/* scale factor */
		double scale = constK*(*Det++)*(*Weight++)*J;

		/* integration constants */		
		fStressMat *= scale;
	
		/* using the stress symmetry */
		fStressStiff.MultQTBQ(fGradNa, fStressMat, format,
			dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
	
		/* strain displacement matrix */
		fShapes->B(fDNa_x, fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, fNumDOF);
}

//DEV - Rayleigh damping should be added to the constitutive level
#if 0
/*
* Compute the effective acceleration and velocities based
* on the algorithmic flags formXx and the given constants
* constXx.
*
*      acc_eff  = constMa acc  + constCv a vel
*      vel_eff  = constCv b vel;
*      disp_eff = constKd disp
*
* where a and b are the Rayleigh damping coefficients.
*
*        ***The effective displacement does not include
*           velocity since the internal force is a nonlinear
*           function of the displacements
*/
void TotalLagrangianT::ComputeEffectiveDVA(int formBody,
	int formMa, double constMa, int formCv, double constCv,
	int formKd, double constKd)
{
	/* acceleration */
	if (formMa || formBody)
	{
		if (formMa)
			SetLocalU(fLocAcc);
		else
			fLocAcc = 0.0;
		
		if (formBody) AddBodyForce(fLocAcc);

		fLocAcc *= constMa;	
	}
	else
		fLocAcc = 0.0;
	
	/* displacement */
	if (formKd)
	{
		SetLocalU(fLocDisp);
		fLocDisp *= constKd;	
	}
	else
		fLocDisp = 0.0;
	
	/* Rayleigh damping */
	if (formCv)
	{
		SetLocalU(fLocVel);
		fLocVel *= constCv;
		
		/* effective a */
		fLocAcc.AddScaled(fCurrMaterial->MassDamping(), fLocVel);
		
		/* effective v */
		fLocVel *= fCurrMaterial->StiffnessDamping();
	}
	else
		fLocVel = 0.0;
}	
#endif

/* calculate the internal force contribution ("-k*d") */
void TotalLagrangianT::FormKd(double constK)
{
	/* matrix alias to fTemp */
	dMatrixT fWP(fNumSD, fStressStiff.Rows(), fNEEvec.Pointer());

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
//TEMP: compute 1st PK stress based on Cauchy stress since most materials
//      have not been optimized to compute PK2 directly.	
	
		/* get Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F^(-1) */
		fTempMat2 = DeformationGradient();
		double J = fTempMat2.Det();
		if (J <= 0.0)
		{
			cout << "\n TotalLagrangianT::FormKd: negative jacobian determinant" << endl;
			throw eBadJacobianDet;
		}
		else
			fTempMat2.Inverse();

		/* compute PK1 */
		fStressMat.MultABT(fTempMat1, fTempMat2);

		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* Wi,J PiJ */
		fWP.MultAB(fStressMat, fGradNa);

		/* accumulate */
		fRHS.AddScaled(J*constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}
