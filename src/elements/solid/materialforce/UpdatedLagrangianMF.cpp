#include "UpdatedLagrangianMF.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "toolboxConstants.h"
#include "SolidMaterialT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
UpdatedLagrangianMF::UpdatedLagrangianMF(const ElementSupportT& support, const FieldT& field):
	FiniteStrainMF(support, field),
	fCauchyStress(NumSD()),
	fLocCurrCoords(LocalArrayT::kCurrCoords)
{
	/* consistency check */
	if (ElementSupport().Analysis() == GlobalT::kLinStatic ||
	    ElementSupport().Analysis() == GlobalT::kLinDynamic)
	{
		cout << "\nUpLag_FDElasticT::UpdatedLagrangianT: no current coordinates required\n" << endl;
		fLocCurrCoords.SetType(LocalArrayT::kInitCoords);
	}	
}

/* destructors */
UpdatedLagrangianMF::~UpdatedLagrangianMF(void)
{
	delete fCurrShapes;
	fCurrShapes = NULL;
}

/* data initialization */
void UpdatedLagrangianMF::Initialize(void)
{
	/* inherited */
	FiniteStrainT::Initialize();

	/* dimension */
	fGradNa.Dimension(NumSD(), NumElementNodes());
	fStressStiff.Dimension(NumElementNodes());
}

/***********************************************************************
* Protected
***********************************************************************/

/* initialize local arrays */
void UpdatedLagrangianMF::SetLocalArrays(void)
{
	/* inherited */
	FiniteStrainT::SetLocalArrays();

	/* allocate and set source */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}

/* initialization functions */
void UpdatedLagrangianMF::SetShape(void)
{
	/* inherited */
	FiniteStrainT::SetShape();

	/* linked shape functions */
	fCurrShapes = new ShapeFunctionT(*fShapes, fLocCurrCoords);
	if (!fCurrShapes) throw ExceptionT::kOutOfMemory ;

	fCurrShapes->Initialize();
}

/* form shape functions and derivatives */
void UpdatedLagrangianMF::SetGlobalShape(void)
{
	/* inherited */
	FiniteStrainT::SetGlobalShape();

	/* shape function wrt current config */
	SetLocalX(fLocCurrCoords);
	fCurrShapes->SetDerivatives();
}

/* form the element stiffness matrix */
void UpdatedLagrangianMF::FormStiffness(double constK)
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
		Set_B(fCurrShapes->Derivatives_U(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
}

/* calculate the internal force contribution ("-k*d") */
void UpdatedLagrangianMF::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		/* strain displacement matrix */
		Set_B(fCurrShapes->Derivatives_U(), fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}
