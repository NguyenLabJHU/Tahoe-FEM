/* $Id: SmallStrainAxiT.cpp,v 1.3 2004-07-15 08:26:27 paklein Exp $ */
#include "SmallStrainAxiT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"

#include <math.h>

const double Pi = acos(-1.0);
const int kRadialDirection = 0; /* the x direction is radial */
const int kNSD = 2;

using namespace Tahoe;

/* constructor */
SmallStrainAxiT::SmallStrainAxiT(const ElementSupportT& support):
	SmallStrainT(support),
	fIPInterp(kNSD),
	fStrain2D(kNSD),
	fStress2D_axi(dSymMatrixT::k3D_plane)	
{
	SetName("small_strain_axi");
}

/* information about subordinate parameter lists */
void SmallStrainAxiT::DefineSubs(SubListT& sub_list) const
{
	const char caller[] = "SmallStrainAxiT::DefineSubs";

	/* inherited */
	SmallStrainT::DefineSubs(sub_list);

	/* remove previous block definition */
	const char old_block_def[] = "small_strain_element_block";
	if (!sub_list.RemoveSub(old_block_def))
		ExceptionT::GeneralFail(caller, "did not find \"%s\"", old_block_def);

	/* element block/material specification */
	sub_list.AddSub("small_strain_axi_element_block", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SmallStrainAxiT::NewSub(const StringT& name) const
{
	if (name == "small_strain_axi_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("small_strain_material_3D");
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return SmallStrainT::NewSub(name);
}

/* accept parameter list */
void SmallStrainAxiT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SmallStrainAxiT::TakeParameterList";

	/* inherited */
	SmallStrainT::TakeParameterList(list);

	//TEMP - B-bar not implemented
	if (fStrainDispOpt == kMeanDilBbar)
		ExceptionT::GeneralFail(caller, "no B-bar");

	/* redimension with out-of-plane component */
	int nstrs = dSymMatrixT::NumValues(kNSD) + 1;
	fD.Dimension(nstrs);
	fB.Dimension(nstrs, NumSD()*NumElementNodes());

	/* allocate strain list */
	for (int i = 0; i < fStrain_List.Length(); i++)
		fStrain_List[i].Dimension(dSymMatrixT::k3D);
	
	/* allocate "last" strain list */
	for (int i = 0; i < fStrain_last_List.Length(); i++)
		fStrain_last_List[i].Dimension(dSymMatrixT::k3D);	
}

/* extract the list of material parameters */
void SmallStrainAxiT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "SmallStrainT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	mat_params.SetName("small_strain_material_3D");
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("small_strain_axi_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("small_strain_axi_element_block", i);
		
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
MaterialSupportT* SmallStrainAxiT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* construct 3D support */
	if (!p) {
		p = new SSMatSupportT(NumDOF(), NumIP());
		p->SetNumSD(3);
	}

	/* inherited initializations */
	SmallStrainT::NewMaterialSupport(p);

	return p;
}

/* initialize local field arrays. Allocate B-bar workspace if needed. */
void SmallStrainAxiT::SetLocalArrays(void)
{
	/* inherited */
	SmallStrainT::SetLocalArrays();

	/* using B-bar - need average of out of plane strain*/
	if (fStrainDispOpt == kMeanDilBbar)
		fMeanGradient.Dimension(3, NumElementNodes());
}

/* calculate the internal force contribution ("-k*d") */
void SmallStrainAxiT::FormKd(double constK)
{
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	double Pi2 = Pi*2.0;
	int nen = NumElementNodes();
	
	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* coordinates of the current integration point */
		fShapes->IPCoords(fIPInterp);
		double r = fIPInterp[kRadialDirection];

		/* collect array of nodal shape functions */
		fIPShape.Alias(nen, fShapes->IPShapeX());

		/* strain displacement matrix */
		//if (fStrainDispOpt == kMeanDilBbar)
		//	Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		//else
		Set_B_axi(fIPShape, fShapes->Derivatives_U(), r, fB);

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

/* form the element stiffness matrix */
void SmallStrainAxiT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	double Pi2 = Pi*2.0;
	int nen = NumElementNodes();

	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* coordinates of the current integration point */
		fShapes->IPCoords(fIPInterp);
		double r = fIPInterp[kRadialDirection];

		/* collect array of nodal shape functions */
		fIPShape.Alias(nen, fShapes->IPShapeX());

		double scale = Pi2*r*constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
//		if (fStrainDispOpt == kMeanDilBbar)
//			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
//		else
		Set_B_axi(fIPShape, fShapes->Derivatives_U(), r, fB);

		/* get D matrix */
		fD.Rank4ReduceFrom3D(fCurrMaterial->c_ijkl());
		fD *= scale;
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

/* compute the measures of strain/deformation over the element */
void SmallStrainAxiT::SetGlobalShape(void)
{
	/* skip call to SmallStrainT::SetGlobalShape */
	SolidElementT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		ExceptionT::GeneralFail("SmallStrainAxiT::SetGlobalShape", "no B-bar");
	}
	else
	{
		/* loop over integration points */
		for (int i = 0; i < NumIP(); i++)
		{
			/* integration point coordinates */
			fShapes->IPCoords(fIPInterp, i);
			double r = fIPInterp[kRadialDirection];
		
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* displacement gradient */
				fShapes->GradU(fLocDisp, fGradU, i);

				/* symmetric part */
				fStrain2D.Symmetrize(fGradU);

				/* integration point displacement */
				fShapes->InterpolateU(fLocDisp, fIPInterp, i);			

				/* make axisymmetric */
				dSymMatrixT& strain_axi = fStrain_List[i];
				strain_axi.ExpandFrom2D(fStrain2D);
				strain_axi(2,2) = fIPInterp[kRadialDirection]/r;
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* displacement gradient */
				fShapes->GradU(fLocLastDisp, fGradU, i);

				/* symmetric part */
				fStrain2D.Symmetrize(fGradU);

				/* integration point displacement */
				fShapes->InterpolateU(fLocLastDisp, fIPInterp, i);			

				/* make axisymmetric */
				dSymMatrixT& strain_axi = fStrain_last_List[i];
				strain_axi.ExpandFrom2D(fStrain2D);
				strain_axi(2,2) = fIPInterp[kRadialDirection]/r;
			}
		}
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean shape function gradient, Hughes (4.5.23) */
void SmallStrainAxiT::SetMeanGradient(dArray2DT& mean_gradient) const
{
#pragma message("correct integration volume")

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
