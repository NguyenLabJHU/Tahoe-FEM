/* $Id: TotalLagrangianAxiT.cpp,v 1.2 2004-02-03 08:24:57 paklein Exp $ */
#include "TotalLagrangianAxiT.h"

#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"

const double Pi2 = 2.0*acos(-1.0);
const int kRadialDirection = 0; /* x <-> r */

using namespace Tahoe;

/* constructor */
TotalLagrangianAxiT::TotalLagrangianAxiT(const ElementSupportT& support, const FieldT& field):
	FiniteStrainAxiT(support, field),
	fStressMat(3),
	fTempMat1(3),
	fTempMat2(3)
{

}

/* data initialization */
void TotalLagrangianAxiT::Initialize(void)
{
	/* inherited */
	FiniteStrainAxiT::Initialize();

	/* dimension */
	fGradNa.Dimension(NumSD(), NumElementNodes());
	fStressStiff.Dimension(NumElementNodes());
	fTemp2.Dimension(NumElementNodes()*NumDOF());
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form the element stiffness matrix */
void TotalLagrangianAxiT::FormStiffness(double constK)
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
	fNEEvec = 0.0;
	
	int ndof = NumSD();
	int nun = fLocDisp.NumberOfNodes();
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		int ip = fShapes->CurrIP();
		double R = fRadius_X[ip];
		double r = fRadius_x[ip];

		/* collect array of nodal shape functions */
		const double* Na_u = fShapes->IPShapeU();
		fIPShape.Alias(nun, Na_u);
		double* u_r = fNEEvec.Pointer(kRadialDirection);
		for (int a = 0; a < nun; a++) {
			*u_r = *Na_u++;
			u_r += ndof;
		}
			
	/* S T R E S S   S T I F F N E S S */
				
		/* Cauchy stress (and set deformation gradient) */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);

		/* F */	
		fTempMat1 = DeformationGradient();
		double F_33 = fTempMat1(2,2);

		/* compute F^-1 in 2D */
		fMat2D.Rank2ReduceFrom3D(fTempMat1);
		double J = fMat2D.Det()*F_33;
		fMat2D.Inverse();
		fTempMat1.Rank2ExpandFrom2D(fMat2D);
		 fTempMat1(2,2) = 1.0/F_33;

		/* chain rule shape function derivatives */
		fShapes->TransformDerivatives(fMat2D, fDNa_x);

		/* get shape function gradients matrix */
		fShapes->GradNa(fDNa_x, fGradNa);

		/* scale factor */
		double scale = Pi2*R*constK*(*Det++)*(*Weight++)*J;

		/* integration constants */		
		fStressMat *= scale;
	
		/* using the stress symmetry */
		fMat2D.Rank2ReduceFrom3D(fStressMat);
		fStressStiff.MultQTBQ(fGradNa, fMat2D, format,
			dMatrixT::kAccumulate);

		/* contribution from out-of-plane stress */
		fLHS.Outer(fNEEvec, fNEEvec, fStressMat(2,2)/(r*r), dMatrixT::kAccumulate);
			
	/* M A T E R I A L   S T I F F N E S S */

		/* strain displacement matrix */
		Set_B_axi(fIPShape, fDNa_x, r, fB);

		/* get D matrix */
		fD.Rank4ReduceFrom3D(fCurrMaterial->c_ijkl());
		fD *= scale;
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
}

/* calculate the internal force contribution ("-k*d") */
void TotalLagrangianAxiT::FormKd(double constK)
{
	/* matrix alias to fTemp */
	dMatrixT fWP(NumSD(), fStressStiff.Rows(), fNEEvec.Pointer());

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	int ndof = NumSD();
	int nun  = fLocDisp.NumberOfNodes();
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		double R = fRadius_X[ip];
	
		/* get Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F */
		fTempMat2 = DeformationGradient();
		double F_33 = fTempMat2(2,2);

		/* compute F^-1 in 2D */
		fMat2D.Rank2ReduceFrom3D(fTempMat2);
		double J = fMat2D.Det()*F_33;
		if (J <= 0.0) ExceptionT::BadJacobianDet("TotalLagrangianAxiT::FormKd");
		fMat2D.Inverse();
		fTempMat2.Rank2ExpandFrom2D(fMat2D);
		fTempMat2(2,2) = 1.0/F_33;

		/* compute PK1/J */
		fStressMat.MultABT(fTempMat1, fTempMat2);
		fMat2D.Rank2ReduceFrom3D(fStressMat);

		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* Wi,J PiJ */
		fWP.MultAB(fMat2D, fGradNa);

		/* accumulate */
		double scale = Pi2*R*J*constK*(*Weight++)*(*Det++);
		fRHS.AddScaled(scale, fNEEvec);
		
		/* contribution from out-of-plane component: x <-> r */
		scale *= fStressMat(2,2)/R;
		const double* NaU = fShapes->IPShapeU();
		double* pRHS = fRHS.Pointer(kRadialDirection);
		for (int a = 0; a < nun; a++) {
			*pRHS += scale*(*NaU++);
			pRHS += ndof;
		}
	}	
}
