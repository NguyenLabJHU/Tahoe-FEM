/* $Id: UpdatedLagrangianAxiT.cpp,v 1.2 2004-06-28 22:41:14 hspark Exp $ */
#include "UpdatedLagrangianAxiT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "toolboxConstants.h"
#include "SolidMaterialT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

const double Pi2 = 2.0*acos(-1.0);
const int kRadialDirection = 0; /* x <-> r */

/* constructor */
UpdatedLagrangianAxiT::UpdatedLagrangianAxiT(const ElementSupportT& support, const FieldT& field):
	FiniteStrainAxiT(support, field),
	fStress2D_axi(dSymMatrixT::k3D_plane),
	fStressMat(3)
{
	/* consistency check */
	if (ElementSupport().Analysis() == GlobalT::kLinStatic ||
	    ElementSupport().Analysis() == GlobalT::kLinDynamic)
	{
		cout << "\nUpLag_FDElasticT::UpdatedLagrangianAxiT: no current coordinates required\n" << endl;
		fLocCurrCoords.SetType(LocalArrayT::kInitCoords);
	}	
}

/* destructors */
UpdatedLagrangianAxiT::~UpdatedLagrangianAxiT(void)
{
	delete fCurrShapes;
	fCurrShapes = NULL;
}

/* data initialization */
void UpdatedLagrangianAxiT::Initialize(void)
{
	/* inherited */
	FiniteStrainAxiT::Initialize();

	/* dimension */
	fGradNa.Dimension(NumSD(), NumElementNodes());
	fStressStiff.Dimension(NumElementNodes());
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialization functions */
void UpdatedLagrangianAxiT::SetShape(void)
{
	/* inherited */
	FiniteStrainAxiT::SetShape();

	/* linked shape functions */
	fCurrShapes = new ShapeFunctionT(*fShapes, fLocCurrCoords);
	if (!fCurrShapes) throw ExceptionT::kOutOfMemory ;

	fCurrShapes->Initialize();
}

/* form shape functions and derivatives */
void UpdatedLagrangianAxiT::SetGlobalShape(void)
{
	/* inherited */
	FiniteStrainAxiT::SetGlobalShape();

	/* shape function wrt current config */
	fCurrShapes->SetDerivatives();
}

/* form the element stiffness matrix */
void UpdatedLagrangianAxiT::FormStiffness(double constK)
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
	fNEEvec = 0.0;

	int ndof = NumDOF();
	int nen  = NumElementNodes();
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		int ip = fShapes->CurrIP();
		double r = fRadius_x[ip];

		/* double scale factor */
		double scale = Pi2*r*constK*(*Det++)*(*Weight++);

		/* collect array of nodal shape functions */
		const double* Na_u = fCurrShapes->IPShapeU();
		fIPShape.Alias(nen, Na_u);
		double* u_r = fNEEvec.Pointer(kRadialDirection);
		for (int a = 0; a < nen; a++) {
			*u_r = *Na_u++;
			u_r += ndof;
		}
	
	/* S T R E S S   S T I F F N E S S */			
		/* compute Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
	
		/* integration constants */		
		fStressMat *= scale;

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
	
		/* using the stress symmetry */
		fMat2D.Rank2ReduceFrom3D(fStressMat);
		fStressStiff.MultQTBQ(fGradNa, fMat2D, format, dMatrixT::kAccumulate);

		/* contribution from out-of-plane stress */
		fLHS.Outer(fNEEvec, fNEEvec, fStressMat(2,2)/(r*r), dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
		/* strain displacement matrix */
		Set_B_axi(fIPShape, fCurrShapes->Derivatives_U(), r, fB);

		/* get D matrix */
		fD.Rank4ReduceFrom3D(fCurrMaterial->c_ijkl());
		fD *= scale;
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, ndof, dMatrixT::kAccumulate);
}

/* calculate the internal force contribution ("-k*d") */
void UpdatedLagrangianAxiT::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	int nen = NumElementNodes();
	fCurrShapes->TopIP();
	while (fCurrShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		double r = fRadius_x[ip];
	
		/* collect array of nodal shape functions */
		fIPShape.Alias(nen, fShapes->IPShapeU());
	
		/* strain displacement matrix */
		Set_B_axi(fIPShape, fCurrShapes->Derivatives_U(), r, fB);

		/* translate to axisymmetric */
		fStress2D_axi.ReduceFrom3D(fCurrMaterial->s_ij());

		/* B^T * Cauchy stress */
		fB.MultTx(fStress2D_axi, fNEEvec);

		/* accumulate */
		fRHS.AddScaled(Pi2*r*constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}
