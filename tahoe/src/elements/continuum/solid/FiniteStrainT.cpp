/* $Id: FiniteStrainT.cpp,v 1.19.26.5 2004-07-13 16:42:32 paklein Exp $ */
#include "FiniteStrainT.h"

#include "ShapeFunctionT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"
#include "ParameterContainerT.h"

/* materials lists */
#include "FSSolidMatList1DT.h"
#include "FSSolidMatList2DT.h"
#include "FSSolidMatList3DT.h"

using namespace Tahoe;

/* constructor */
FiniteStrainT::FiniteStrainT(const ElementSupportT& support):
	SolidElementT(support),
	fNeedsOffset(-1),
	fCurrShapes(NULL),
	fFSMatSupport(NULL)
{
	SetName("large_strain");
}

/* destructor */
FiniteStrainT::~FiniteStrainT(void) {
	delete fFSMatSupport;
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	if (fCurrShapes)
	{
		/* field gradient */
		fCurrShapes->GradU(u, grad_u);
	}
	else
		ExceptionT::GeneralFail("FiniteStrainT::ComputeGradient", "shape functions wrt current coords not defined");
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
	if (fCurrShapes)
	{
		/* field gradient */
		fCurrShapes->GradU(u, grad_u, ip);
	}
	else
		ExceptionT::GeneralFail("FiniteStrainT::ComputeGradient", "shape functions wrt current coords not defined");
}

/* information about subordinate parameter lists */
void FiniteStrainT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	

	/* element block/material specification */
	sub_list.AddSub("large_strain_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list. */
void FiniteStrainT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "large_strain_material_choice")
	{
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("large_strain_material_1D");
		sub_lists.AddSub("large_strain_material_2D");
		sub_lists.AddSub("large_strain_material_3D");
	}
	else /* inherited */
		SolidElementT::DefineInlineSub(name, order, sub_lists);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* FiniteStrainT::NewSub(const StringT& name) const
{
	if (name == "large_strain_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("large_strain_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);

}

/* accept parameter list */
void FiniteStrainT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidElementT::TakeParameterList(list);

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
		FSSolidMatT* mat = (FSSolidMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kF     ] = mat->Need_F();
		needs[fNeedsOffset + kF_last] = mat->Need_F_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kF];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kF_last];
	}

	/* what's needed */
	bool need_F = false;
	bool need_F_last = false;
	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		need_F = need_F || Needs_F(i);		
		need_F_last = need_F_last || Needs_F_last(i);
	}	

	/* allocate deformation gradient list */
	if (need_F)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_all.Dimension(nip*nsd*nsd);
		fF_List.Dimension(nip);
		for (int i = 0; i < nip; i++)
			fF_List[i].Set(nsd, nsd, fF_all.Pointer(i*nsd*nsd));
	}
	
	/* allocate "last" deformation gradient list */
	if (need_F_last)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_last_all.Dimension(nip*nsd*nsd);
		fF_last_List.Dimension(nip);
		for (int i = 0; i < nip; i++)
			fF_last_List[i].Set(nsd, nsd, fF_last_all.Pointer(i*nsd*nsd));
	}
}

/* extract the list of material parameters */
void FiniteStrainT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "SmallStrainT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("large_strain_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("large_strain_element_block", i);
		
		/* resolve material list name */
		if (i == 0) {
			const ParameterListT& mat_list_params = block.GetListChoice(*this, "large_strain_material_choice");
			mat_params.SetName(mat_list_params.	Name());
		}
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* FiniteStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new FSMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	SolidElementT::NewMaterialSupport(p);
	
	/* set FiniteStrainT fields */
	FSMatSupportT* ps = TB_DYNAMIC_CAST(FSMatSupportT*, p);
	if (ps) {
		ps->SetDeformationGradient(&fF_List);
		ps->SetDeformationGradient_last(&fF_last_List);
	}

	return p;
}

/* construct materials manager and read data */
MaterialListT* FiniteStrainT::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	int nsd = -1;
	if (name == "large_strain_material_1D")
		nsd = 1;
	else if (name == "large_strain_material_2D")
		nsd = 2;
	else if (name == "large_strain_material_3D")
		nsd = 3;
	
	/* no match */
	if (nsd == -1) return NULL;

	if (size > 0)
	{
		/* material support */
		if (!fFSMatSupport) {
			fFSMatSupport = TB_DYNAMIC_CAST(FSMatSupportT*, NewMaterialSupport());
			if (!fFSMatSupport) ExceptionT::GeneralFail("FiniteStrainT::NewMaterialList");
		}

		if (nsd == 1)
			return new FSSolidMatList1DT(size, *fFSMatSupport);
		else if (nsd == 2)
			return new FSSolidMatList2DT(size, *fFSMatSupport);
		else if (nsd == 3)
			return new FSSolidMatList3DT(size, *fFSMatSupport);
	}
	else
	 {
	 	if (nsd == 1)
	 		return new FSSolidMatList1DT;
		else if (nsd == 2)
			return new FSSolidMatList2DT;
		else if (nsd == 3)
			return new FSSolidMatList3DT;
	}
	
	/* no match */
	return NULL;
}

/* form shape functions and derivatives */
void FiniteStrainT::SetGlobalShape(void)
{
	/* inherited */
	SolidElementT::SetGlobalShape();

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			dMatrixT& mat = fF_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			dMatrixT& mat = fF_last_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocLastDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}
	}
}

/* write all current element information to the stream */
void FiniteStrainT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	SolidElementT::CurrElementInfo(out);
	
	/* write deformation gradients */
	out << "\n i.p. deformation gradients:\n";
	for (int i = 0; i < fF_List.Length(); i++)
		out << " ip: " << i+1 << '\n'
		    << fF_List[i] << '\n';
	out << '\n';
}
