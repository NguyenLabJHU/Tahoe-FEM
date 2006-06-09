/* $Id: TotalLagrangianFlatT.cpp,v 1.3 2006-06-09 22:58:50 jzimmer Exp $ */
#include "TotalLagrangianFlatT.h"

#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "SolidMaterialT.h"
#include "MaterialListT.h"

using namespace Tahoe;

/* constructor */
TotalLagrangianFlatT::TotalLagrangianFlatT(const ElementSupportT& support):
	TotalLagrangianT(support)
{
	SetName("total_lagrangian_flat");
}

void TotalLagrangianFlatT::RHSDriver(void)
{
	/* inherited from SolidElementT */
	ContinuumElementT::RHSDriver();

	/* heat source if needed */
	const FieldT* temperature = ElementSupport().Field("temperature");

	/* initialize sources */
	if (temperature && fIncrementalHeat.Length() == 0) {
	
		/* allocate the element heat */
		fElementHeat.Dimension(fShapes->NumIP());
			
		/* initialize heat source arrays */
		fIncrementalHeat.Dimension(fBlockData.Length());
		for (int i = 0; i < fIncrementalHeat.Length(); i++)
		{
			/* dimension */
			fIncrementalHeat[i].Dimension(fBlockData[i].Dimension(), NumIP());

			/* register */
			temperature->RegisterSource(fBlockData[i].ID(), fIncrementalHeat[i]);
		}
	}
	
	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodySchedule && fBody.Magnitude() > kSmall))
	{	
		formBody = 1;
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;

	/* workspace */
	dMatrixT fWP;
	fWP.Alias(NumSD(), fStressStiff.Rows(), fNEEvec.Pointer());

	/* materials information */
	const MaterialListT& materials = MaterialsList();

	/* loop over elements */
	int block_count = 0, block_dex = 0;
	fElementCards.Top();
	while (fElementCards.Next())
	{
		/* advance to block [skip empty blocks] */		
		while (block_count == fBlockData[block_dex].Dimension()) {
			block_count = 0;
			block_dex++;
		}

		/* current element information */
		const ElementCardT& element = CurrentElement();
		const iArrayT& nodesU = element.NodesU(); /* nodes defining field over the element */
		const iArrayT& nodesX = element.NodesX(); /* nodes defining geometry over the element */

		/* set pointer to the material for the current element */
		ContinuumMaterialT* pcont_mat = materials[element.MaterialNumber()]; /* pull pointer out of proxy return value */
		fCurrMaterial = (SolidMaterialT*) pcont_mat; /* cast is safe since class contructs materials list */
		
		/* fetch (initial) coordinates for the element */
		fLocInitCoords.SetLocal(nodesX);
	
		/* compute shape functions and derivatives */
		fShapes->SetDerivatives();

		/* collect field values over the element */
		fLocDisp.SetLocal(nodesU);
		fLocLastDisp.SetLocal(nodesU);

		/* get nodal temperatures if available */
		if (fLocTemp) fLocTemp->SetLocal(nodesU);
		if (fLocTemp_last) fLocTemp_last->SetLocal(nodesU);

		/* initialize */
		fRHS = 0.0;
		fElementHeat = 0.0;

		/* integrate internal force */
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			int currIP = fShapes->CurrIP();

			/* deformation gradient */
			if (fF_List.Length() > 0) {
				dMatrixT& F_mat = fF_List[currIP];
				fShapes->GradU(fLocDisp, F_mat, currIP);
				F_mat.PlusIdentity();
			}

			/* "last" deformation gradient */
			if (fF_last_List.Length() > 0) {
				dMatrixT& F_mat_last = fF_last_List[currIP];
				fShapes->GradU(fLocLastDisp, F_mat_last, currIP);
				F_mat_last.PlusIdentity();
			}

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
			fRHS.AddScaled(-J*constKd*(*Weight++)*(*Det++), fNEEvec);
		}

		/* inertia forces */
		if (formMa || formBody)
		{
			/* nodal accelerations */
			if (formMa)
				SetLocalU(fLocAcc);
			else 
				fLocAcc = 0.0;
			
			/* body force contribution */
			if (formBody) AddBodyForce(fLocAcc);
		
			FormMa(fMassType, -constMa*fCurrMaterial->Density(), false, &fLocAcc, NULL, NULL);
		}
		
		/* store incremental heat */
		if (temperature)
			fIncrementalHeat[block_dex].SetRow(block_count, fElementHeat);

		/* assemble element force */
		ElementSupport().AssembleRHS(Group(), fRHS, element.Equations());

		/* next in block */
		block_count++;
	}
}
