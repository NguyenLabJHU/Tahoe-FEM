/* $Id: SmallStrainT.cpp,v 1.7 2002-09-23 06:58:25 paklein Exp $ */

#include "SmallStrainT.h"
#include "ShapeFunctionT.h"
#include "SSStructMatT.h"
#include "MaterialListT.h"

/* constructor */

using namespace Tahoe;

SmallStrainT::SmallStrainT(const ElementSupportT& support, const FieldT& field):
	ElasticT(support, field),
	fNeedsOffset(-1),
	fGradU(NumSD())
{

}

/* called immediately after constructor */
void SmallStrainT::Initialize(void)
{
	/* inherited */
	ElasticT::Initialize();

	/* what's needed */
	bool need_strain = false;
	bool need_strain_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_strain = need_strain || needs[fNeedsOffset + kstrain];
		need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
	}

	/* allocate deformation gradient list */
	if (need_strain)
	{
		fStrain_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_List[i].Allocate(NumSD());
	}
	
	/* allocate "last" deformation gradient list */
	if (need_strain_last)
	{
		fStrain_last_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_last_List[i].Allocate(NumSD());
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct list of materials from the input stream */
void SmallStrainT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	ElasticT::ReadMaterialData(in);

	/* offset to class needs flags */
	fNeedsOffset = fMaterialNeeds[0].Length();
	
	/* set material needs */
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* needs array */
		ArrayT<bool>& needs = fMaterialNeeds[i];

		/* resize array */
		needs.Resize(needs.Length() + 2, true);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SSStructMatT* mat = (SSStructMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kstrain     ] = mat->Need_Strain();
		needs[fNeedsOffset + kstrain_last] = mat->Need_Strain_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kstrain];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kstrain_last];
	}
}

/* initialize local field arrays. Allocate B-bar workspace if needed. */
void SmallStrainT::SetLocalArrays(void)
{
	/* inherited */
	ElasticT::SetLocalArrays();

	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar) {
		fLocDispTranspose.Dimension(fLocDisp.Length());
		fMeanGradient.Dimension(NumSD(), NumElementNodes());
	}
}

/* calculate the internal force contribution ("-k*d") */
void SmallStrainT::FormKd(double constK)
{
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();
	
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);
		
		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
		
		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}

/* form the element stiffness matrix */
void SmallStrainT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

/* compute the measures of strain/deformation over the element */
void SmallStrainT::SetGlobalShape(void)
{
	/* inherited */
	ElasticT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		/* compute mean of shape function gradients */
		SetMeanGradient(fMeanGradient);

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* set B-bar */
			int ip = fShapes->CurrIP();
			Set_B_bar(fShapes->Derivatives_U(ip), fMeanGradient, fB);
	
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* transpose displacement array */
				fLocDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* transpose displacement array */
				fLocLastDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_last_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}
		}		
	}
	else
	{
		/* loop over integration points */
		for (int i = 0; i < NumIP(); i++)
		{
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* displacement gradient */
				fShapes->GradU(fLocDisp, fGradU, i);

				/* symmetric part */
				 fStrain_List[i].Symmetrize(fGradU);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* displacement gradient */
				fShapes->GradU(fLocLastDisp, fGradU, i);

				/* symmetric part */
				 fStrain_last_List[i].Symmetrize(fGradU);
			}
		}
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean shape function gradient, Hughes (4.5.23) */
void SmallStrainT::SetMeanGradient(dArray2DT& mean_gradient) const
{
	int nip = NumIP();
	const double* det = fShapes->IPDets();
	const double*   w = fShapes->IPWeights();

	/* volume */
	double vol = 0.0;
	for (int i = 0; i < nip; i++)
		vol += w[i]*det[i];

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/vol, fShapes->Derivatives_U(i));
}