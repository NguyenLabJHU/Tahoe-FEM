/* $Id: SmallStrainT.cpp,v 1.5 2002-07-02 19:55:24 cjkimme Exp $ */

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

/* compute field gradients */
void SmallStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
#pragma unused(u)
#pragma unused(grad_u)
	cout << "\n SmallStrainT::ComputeGradient: not implemented" << endl;
	throw;
}

/* compute field gradients */
void SmallStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	cout << "\n SmallStrainT::ComputeGradient: not implemented" << endl;
	throw;
}

/* compute field gradients from the end of the previous time step */
void SmallStrainT::ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u) const
{
#pragma unused(u)
#pragma unused(grad_u)
	cout << "\n SmallStrainT::ComputeGradient_last: not implemented" << endl;
	throw;
}

/* compute field gradients from the end of the previous time step */
void SmallStrainT::ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	cout << "\n SmallStrainT::ComputeGradient_last: not implemented" << endl;
	throw;
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
	if (fStrainDispOpt == ShapeFunctionT::kMeanDilBbar)
		fLocDispTranspose.Allocate(fLocDisp.Length());
}

/* increment current element */
void SmallStrainT::SetGlobalShape(void)
{
	/* inherited */
	ElasticT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == ShapeFunctionT::kMeanDilBbar)
	{
		/* compute mean dilatation */
		fShapes->SetMeanDilatation();

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* get B-bar */
			fShapes->B(fB);
		
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* transpose displacement array */
				fLocDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_List[fShapes->CurrIP()];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* transpose displacement array */
				fLocLastDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_last_List[fShapes->CurrIP()];
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
