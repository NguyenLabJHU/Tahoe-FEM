/* $Id: TotalLagrangianT.cpp,v 1.14 2004-07-15 08:26:27 paklein Exp $ */
/* created: paklein (09/07/1998) */
#include "TotalLagrangianT.h"

#include "toolboxConstants.h"
#include "SolidMaterialT.h"
#include "MaterialListT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
TotalLagrangianT::TotalLagrangianT(const ElementSupportT& support):
	FiniteStrainT(support)
{
	SetName("total_lagrangian");
}

/* accept parameter list */
void TotalLagrangianT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FiniteStrainT::TakeParameterList(list);
	
	/* dimension workspace */
	fStressMat.Dimension(NumSD());
	fTempMat1.Dimension(NumSD());
	fTempMat2.Dimension(NumSD());

	fGradNa.Dimension(NumSD(), NumElementNodes());
	fStressStiff.Dimension(NumElementNodes());
	fTemp2.Dimension(NumElementNodes()*NumDOF());
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form the element stiffness matrix */
void TotalLagrangianT::FormStiffness(double constK)
{		
//NOTE: Because most materials have been optimized to calculate the
//      Cauchy stress s_ij and the material tangent modulus c_ijkl,
//      the derivatives with respect to the reference coordinates X 
//      are transformed to derivatives with respect to the current
//      coordinates x using the inverse of the deformation gradient.
//      Alternately, the stress and modulus could have been transferred
//      to their material representations, i.e., pulled back.

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
				
		/* Cauchy stress (and set deformation gradient) */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
	
		/* chain rule shape function derivatives */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();
		fTempMat1.Inverse();
		fShapes->TransformDerivatives(fTempMat1, fDNa_x);

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
		Set_B(fDNa_x, fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
}

/* calculate the internal force contribution ("-k*d") */
void TotalLagrangianT::FormKd(double constK)
{
//NOTE: compute 1st P-K stress based on Cauchy stress since most materials
//      have not been optimized to compute PK2 directly.	

	/* matrix alias to fTemp */
	dMatrixT fWP(NumSD(), fStressStiff.Rows(), fNEEvec.Pointer());

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* get Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F^(-1) */
		fTempMat2 = DeformationGradient();
		double J = fTempMat2.Det();
		if (J <= 0.0)
			ExceptionT::BadJacobianDet("TotalLagrangianT::FormKd");
		else
			fTempMat2.Inverse();

		/* compute PK1/J */
		fStressMat.MultABT(fTempMat1, fTempMat2);

		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* Wi,J PiJ */
		fWP.MultAB(fStressMat, fGradNa);

		/* accumulate */
		fRHS.AddScaled(J*constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}
