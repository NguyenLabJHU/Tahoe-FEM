/* $Id: SmallStrainT.cpp,v 1.13.2.1 2004-01-21 19:09:58 paklein Exp $ */
#include "SmallStrainT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
SmallStrainT::SmallStrainT(const ElementSupportT& support, const FieldT& field):
	SolidElementT(support, field),
	fNeedsOffset(-1),
	fGradU(NumSD()),
	fSSMatSupport(NULL)
{
	SetName("small_strain");
}

SmallStrainT::SmallStrainT(const ElementSupportT& support):
	SolidElementT(support),
	fNeedsOffset(-1),
	fSSMatSupport(NULL)
{
	SetName("small_strain");
}

/* destructor */
SmallStrainT::~SmallStrainT(void)
{
	delete fSSMatSupport;
}

/* called immediately after constructor */
void SmallStrainT::Initialize(void)
{
	/* inherited */
	SolidElementT::Initialize();

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
		fStrain_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_List[i].Dimension(NumSD());
	}
	
	/* allocate "last" deformation gradient list */
	if (need_strain_last)
	{
		fStrain_last_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_last_List[i].Dimension(NumSD());
	}
}

/* TEMPORARY */
void SmallStrainT::InitialCondition(void)
{
	/* inherited */
	SolidElementT::InitialCondition();
	
	/* set the source for the iteration number */
	fSSMatSupport->SetIterationNumber(ElementSupport().IterationNumber(Group()));
}

/* implementation of the ParameterInterfaceT interface */
void SmallStrainT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SolidElementT::DefineParameters(list);

	/* strain-displacement relation */
	ParameterT strain_displacement(ParameterT::Enumeration, "strain_displacement");
	strain_displacement.AddEnumeration("standard", kStandardB);
    strain_displacement.AddEnumeration("B_bar", kMeanDilBbar);
    strain_displacement.SetDefault(kStandardB);
	list.AddParameter(strain_displacement);
}

/* information about subordinate parameter lists */
void SmallStrainT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	

	/* element block/material specification */
	sub_list.AddSub("small_strain_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* SmallStrainT::NewSub(const StringT& list_name) const
{
	/* create non-const this */
	SmallStrainT* non_const_this = const_cast<SmallStrainT*>(this);

	if (list_name == "small_strain_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT("small_strain_element_block");
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("small_strain_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else if (list_name == "small_strain_material_1D")
		return non_const_this->NewMaterialList(1,0);
	else if (list_name == "small_strain_material_2D")
		return non_const_this->NewMaterialList(2,0);
	else if (list_name == "small_strain_material_3D")
		return non_const_this->NewMaterialList(3,0);
	else /* inherited */
		return ContinuumElementT::NewSub(list_name);
}

/* return the description of the given inline subordinate parameter list. */
void SmallStrainT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "small_strain_material_choice")
	{
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_sub_list.AddSub("small_strain_material_1D");
		sub_sub_list.AddSub("small_strain_material_2D");
		sub_sub_list.AddSub("small_strain_material_3D");
	}
	else /* inherited */
		return SolidElementT::DefineInlineSub(sub, order, sub_sub_list);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SmallStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new SSMatSupportT(NumSD(), NumDOF(), NumIP());

	/* inherited initializations */
	SolidElementT::NewMaterialSupport(p);
	
	/* set SolidMatSupportT fields */
	SSMatSupportT* ps = TB_DYNAMIC_CAST(SSMatSupportT*, p);
	if (ps) {
		ps->SetLinearStrain(&fStrain_List);
		ps->SetLinearStrain_last(&fStrain_last_List);
	}

	return p;
}

/* return a pointer to a new material list */
MaterialListT* SmallStrainT::NewMaterialList(int nsd, int size)
{
	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSMatSupport) {
			fSSMatSupport = TB_DYNAMIC_CAST(SSMatSupportT*, NewMaterialSupport());
			if (!fSSMatSupport)
				ExceptionT::GeneralFail("SmallStrainT::NewMaterialList");
		}

		if (nsd == 1)
			return new SSSolidMatList1DT(size, *fSSMatSupport);
		else if (nsd == 2)
			return new SSSolidMatList2DT(size, *fSSMatSupport);
		else if (nsd == 3)
			return new SSSolidMatList3DT(size, *fSSMatSupport);
		else
			return NULL;
	}
	else
	{
		if (nsd == 1)
			return new SSSolidMatList1DT;
		else if (nsd == 2)
			return new SSSolidMatList2DT;
		else if (nsd == 3)
			return new SSSolidMatList3DT;
		else
			return NULL;
	}	
}

/* construct list of materials from the input stream */
void SmallStrainT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	SolidElementT::ReadMaterialData(in);

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
		SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;

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
	SolidElementT::SetLocalArrays();

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

	/********DEBUG*******/
	bool print = false; 
	int pos = fElementCards.Position(); 
	if (pos == 1&&0)  
	  print = true; 
	/*******************/
	
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
		if (print) cout << "\nmodulus: "<<fCurrMaterial->c_ijkl();
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

/* compute the measures of strain/deformation over the element */
void SmallStrainT::SetGlobalShape(void)
{
	/* inherited */
	SolidElementT::SetGlobalShape();

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
