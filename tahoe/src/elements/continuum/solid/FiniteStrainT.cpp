/* $Id: FiniteStrainT.cpp,v 1.2 2001-07-03 01:34:50 paklein Exp $ */

#include "FiniteStrainT.h"
#include "ShapeFunctionT.h"
#include "FDStructMatT.h"
#include "MaterialListT.h"

/* constructor */
FiniteStrainT::FiniteStrainT(FEManagerT& fe_manager):
	ElasticT(fe_manager),
	fNeedsOffset(-1)
{

}

/* called immediately after constructor */
void FiniteStrainT::Initialize(void)
{
	/* inherited */
	ElasticT::Initialize();

	/* what's needed */
	bool need_F      = false;
	bool need_F_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_F = need_F || needs[fNeedsOffset + kF];
		need_F_last = need_F_last || needs[fNeedsOffset + kF_last];
	}

	/* allocate deformation gradient list */
	if (need_F)
	{
		fF_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fF_List[i].Allocate(NumSD());
	}
	
	/* allocate "last" deformation gradient list */
	if (need_F_last)
	{
		fF_last_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fF_last_List[i].Allocate(NumSD());
	}
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	/* field gradient */
	fShapes->GradU(u, grad_u);
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
	/* field gradient */
	fShapes->GradU(u, grad_u, ip);
}

/* compute field gradients with respect to reference coordinates */
void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, 
	dMatrixT& grad_u) const
{
#pragma unused(u)
#pragma unused(grad_u)
	cout << "\n FiniteStrainT::ComputeGradient_reference: not implemented" << endl;
	throw;
}

void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	cout << "\n FiniteStrainT::ComputeGradient_reference: not implemented" << endl;
	throw;
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct list of materials from the input stream */
void FiniteStrainT::ReadMaterialData(ifstreamT& in)
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
		FDStructMatT* mat = (FDStructMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kF     ] = mat->Need_F();
		needs[fNeedsOffset + kF_last] = mat->Need_F_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kF];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kF_last];
	}
}

/* increment current element */
void FiniteStrainT::SetGlobalShape(void)
{
	/* inherited */
	ElasticT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs[fNeedsOffset + kF])
		{
			dMatrixT& mat = fF_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}

		/* "last" deformation gradient */
		if (needs[fNeedsOffset + kF_last])
		{
			dMatrixT& mat = fF_last_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocLastDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}
	}
}
