/* $Id: UpdatedLagrangianT.cpp,v 1.3.4.3 2002-05-17 01:27:22 paklein Exp $ */
/* created: paklein (07/03/1996) */

#include "UpdatedLagrangianT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "StructuralMaterialT.h"
#include "ShapeFunctionT.h"

/* constructor */
UpdatedLagrangianT::UpdatedLagrangianT(const ElementSupportT& support, const FieldT& field):
	FiniteStrainT(support, field),
	fCurrShapes(NULL),
	fCauchyStress(NumSD()),
	fLocCurrCoords(LocalArrayT::kCurrCoords)
{
	/* disable any strain-displacement options */
	if (fStrainDispOpt != 0)
	{
		cout << "\nUpLag_FDElasticT::UpdatedLagrangianT: no strain-displacement options\n" << endl;
		fStrainDispOpt = 0;
	}

	/* consistency check */
	if (ElementSupport().Analysis() == GlobalT::kLinStatic ||
	    ElementSupport().Analysis() == GlobalT::kLinDynamic)
	{
		cout << "\nUpLag_FDElasticT::UpdatedLagrangianT: no current coordinates required\n" << endl;
		fLocCurrCoords.SetType(LocalArrayT::kInitCoords);
	}	
}

/* destructors */
UpdatedLagrangianT::~UpdatedLagrangianT(void)
{
	delete fCurrShapes;
	fCurrShapes = NULL;
}

/* data initialization */
void UpdatedLagrangianT::Initialize(void)
{
	/* inherited */
	FiniteStrainT::Initialize();

	/* dimension */
	fGradNa.Allocate(NumSD(), NumElementNodes());
	fStressStiff.Allocate(NumElementNodes());
	fTemp2.Allocate(NumElementNodes()*NumDOF());
}

/***********************************************************************
* Protected
***********************************************************************/

/* initialize local arrays */
void UpdatedLagrangianT::SetLocalArrays(void)
{
	/* inherited */
	FiniteStrainT::SetLocalArrays();

	/* allocate and set source */
	fLocCurrCoords.Allocate(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}

/* initialization functions */
void UpdatedLagrangianT::SetShape(void)
{
	/* inherited */
	FiniteStrainT::SetShape();

	/* linked shape functions */
	fCurrShapes = new ShapeFunctionT(*fShapes, fLocCurrCoords);
	if (!fCurrShapes) throw eOutOfMemory ;

	fCurrShapes->Initialize();
}

/* form shape functions and derivatives */
void UpdatedLagrangianT::SetGlobalShape(void)
{
	/* inherited */
	FiniteStrainT::SetGlobalShape();

	/* shape function wrt current config */
	SetLocalX(fLocCurrCoords);
	fCurrShapes->SetDerivatives();
}

/* form the element stiffness matrix */
void UpdatedLagrangianT::FormStiffness(double constK)
{		
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* double scale factor */
		double scale = constK*(*Det++)*(*Weight++);
	
	/* S T R E S S   S T I F F N E S S */			
		/* compute Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fCauchyStress);
	
		/* integration constants */		
		fCauchyStress *= scale;

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
	
		/* using the stress symmetry */
		fStressStiff.MultQTBQ(fGradNa, fCauchyStress,
			format, dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
		/* strain displacement matrix */
		fCurrShapes->B(fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF());
}

/* calculate the internal force contribution ("-k*d") */
void UpdatedLagrangianT::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		/* get strain-displacement matrix */
		fCurrShapes->B(fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}
