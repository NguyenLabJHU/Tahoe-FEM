/* $Id: SimoFiniteStrainT.cpp,v 1.1 2001-07-11 01:02:15 paklein Exp $ */
#include "SimoFiniteStrainT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "StructuralMaterialT.h"
#include "MaterialListT.h" //NOTE - only needed for check in Initialize?
#include "SimoShapeFunctionT.h"

/* constructor */
SimoFiniteStrainT::SimoFiniteStrainT(FEManagerT& fe_manager):
	FiniteStrainT(fe_manager),
	fEnhancedShapes(NULL),
	fStressMat(fNumSD),
	fTempMat1(fNumSD),
	fTempMat2(fNumSD)
{
	/* disable any strain-displacement options */
	if (fStrainDispOpt != 0)
	{
		cout << "\n SimoFiniteStrainT::SimoFiniteStrainT: no strain-displacement options\n" << endl;
		fStrainDispOpt = 0;
	}

	/* read parameters */
	ifstreamT& in = fe_manager.Input();
	int inc_mode = -1;
	in >> inc_mode;
	in >> fLocalIterationMax;

	/* resolve */
	if (inc_mode != 0 && inc_mode != 1) throw eBadInputValue;
	if (fLocalIterationMax < 1) throw eGeneralFail;
	fIncompressibleMode = (inc_mode == 1);
	if (NumSD() == 2) fIncompressibleMode = false;
	
	/* set number of mode shapes */
	if (NumSD() == 2)
		fNumModeShapes = 2;
	else if (NumSD() == 3)
		fNumModeShapes = (fIncompressibleMode) ? 4 : 3;
	else throw eGeneralFail;
}

/* data initialization */
void SimoFiniteStrainT::Initialize(void)
{
	/* inherited */
	FiniteStrainT::Initialize();

//TEMP
	if (fMaterialList->HasThermalStrains())
	{
		cout << "\n SimoFiniteStrainT::Initialize: element group ";
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

/* write element parameter to out */
void SimoFiniteStrainT::PrintControlData(ostream& out) const
{
	/* inherited */
	FiniteStrainT::PrintControlData(out);
	
	/* parameters */
	out << " Include incompressible mode . . . . . . . . . . = " << fIncompressibleMode << '\n';
	out << "    Number of enhanced mode shapes = " << fNumModeShapes << '\n';
	out << " Maximum number of local sub-iterations. . . . . = " << fLocalIterationMax << '\n';
}

/* construct shape function */
void SimoFiniteStrainT::SetShape(void)
{
	/* allocate memory for incompressible modes */
	fElementModes.Allocate(NumElements(), fNumModeShapes*NumSD());
	fElementModes = 0;
	
	/* dimension before sending to the shape functions */
	fCurrElementModes.Set(fNumModeShapes, NumSD(), fElementModes.Pointer());

	/* construct shape functions */
	fEnhancedShapes = new SimoShapeFunctionT(fGeometryCode, fNumIP,
		fLocInitCoords, fCurrElementModes);
	if (!fEnhancedShapes) throw eOutOfMemory;

	/* initialize */
	fEnhancedShapes->Initialize();
	
	/* set base class pointer */
	fShapes = fEnhancedShapes;
}

/* form shape functions and derivatives */
void SimoFiniteStrainT::SetGlobalShape(void)
{
//NOTE: there's a loop in here somewhere to solve for the
//      for the enhanced modes

	/* 3D with incompressible mode */
	if (fIncompressibleMode)
		ModifiedEnhancedDeformation();
	else /* modification is strictly additive */
	{
		/* inherited */
		FiniteStrainT::SetGlobalShape();
	
	
	
	}
}

/* form the element stiffness matrix */
void SimoFiniteStrainT::FormStiffness(double constK)
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

/* calculate the internal force contribution ("-k*d") */
void SimoFiniteStrainT::FormKd(double constK)
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
			cout << "\n SimoFiniteStrainT::FormKd: negative jacobian determinant" << endl;
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

/***********************************************************************
* Private
***********************************************************************/

/* compute modified, enhanced deformation gradient */
void SimoFiniteStrainT::ModifiedEnhancedDeformation(void)
{
	/* check */
	if (NumSD() != 3)
	{
		cout << "\n SimoFiniteStrainT::ModifiedEnhancedDeformation: for 3D only" << endl;
		throw eGeneralFail;
	}
	
	/* skip base class implementation because the deformation gradient
	 * modification are not simply additive */
	ElasticT::SetGlobalShape();
	
	cout << "\n SimoFiniteStrainT::ModifiedEnhancedDeformation: not done" << endl;
	throw;
}
