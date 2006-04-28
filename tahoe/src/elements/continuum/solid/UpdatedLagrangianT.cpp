/* $Id: UpdatedLagrangianT.cpp,v 1.16 2006-04-28 16:37:13 tdnguye Exp $ */
/* created: paklein (07/03/1996) */
#include "UpdatedLagrangianT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "SolidMaterialT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
UpdatedLagrangianT::UpdatedLagrangianT(const ElementSupportT& support):
	FiniteStrainT(support),
	fLocCurrCoords(LocalArrayT::kCurrCoords)
{
	SetName("updated_lagrangian");
}

/* destructors */
UpdatedLagrangianT::~UpdatedLagrangianT(void)
{
	delete fCurrShapes;
	fCurrShapes = NULL;
}

/* describe the parameters needed by the interface */
void UpdatedLagrangianT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FiniteStrainT::DefineParameters(list);

	/* remove option to store shape functions */
	list.RemoveParameter("store_shapefunctions");
}

/* accept parameter list */
void UpdatedLagrangianT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FiniteStrainT::TakeParameterList(list);
	
	/* allocate workspace */
	int nsd = NumSD();
	int nen = NumElementNodes();
	fCauchyStress.Dimension(nsd);
	fGradNa.Dimension(nsd, nen);
	fStressStiff.Dimension(nen);
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
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}

/* initialization functions */
void UpdatedLagrangianT::SetShape(void)
{
	/* inherited */
	FiniteStrainT::SetShape();

	/* linked shape functions */
	fCurrShapes = new ShapeFunctionT(*fShapes, fLocCurrCoords);
	if (!fCurrShapes) throw ExceptionT::kOutOfMemory ;

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
void UpdatedLagrangianT::FormKd(double constK)
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

	    const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
//		cout << "\n sij tot: "<< cauchy<< endl;
		/* B^T * Cauchy stress */
		fB.MultTx(cauchy, fNEEvec);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}
