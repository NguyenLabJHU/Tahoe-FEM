/* $Id: SmallStrainAxiT.cpp,v 1.2 2004-02-02 23:48:04 paklein Exp $ */
#include "SmallStrainAxiT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"
#include "SSMatSupportT.h"

/* materials lists */
#include "SolidMatList3DT.h" //TEMP - not needed with ParameterListT input

#include <math.h>

const double Pi = acos(-1.0);
const int kRadialDirection = 0; /* the x direction is radial */
const int kNSD = 2;

using namespace Tahoe;

/* constructor */
SmallStrainAxiT::SmallStrainAxiT(const ElementSupportT& support, const FieldT& field):
	SmallStrainT(support, field),
	fIPInterp(kNSD),
	fStrain2D(kNSD),
	fStress2D_axi(dSymMatrixT::k3D_plane)
{
	SetName("small_strain_axi");
}

SmallStrainAxiT::SmallStrainAxiT(const ElementSupportT& support):
	SmallStrainT(support),
	fIPInterp(kNSD),
	fStrain2D(kNSD),
	fStress2D_axi(dSymMatrixT::k3D_plane)	
{
	SetName("small_strain_axi");
}

/* called immediately after constructor */
void SmallStrainAxiT::Initialize(void)
{
	/* inherited */
	SmallStrainT::Initialize();

	//TEMP - B-bar not implemented
	if (fStrainDispOpt == kMeanDilBbar)
		ExceptionT::GeneralFail("SmallStrainAxiT::Initialize", "no B-bar");

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

/* information about subordinate parameter lists */
void SmallStrainAxiT::DefineSubs(SubListT& sub_list) const
{	
	/* inherited */
	SmallStrainT::DefineSubs(sub_list);

	/* remove choice of materials lists */
	sub_list.RemoveSub("solid_materials");

	/* 3D material list only */
	sub_list.AddSub("solid_materials_3D", ParameterListT::Once);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SmallStrainAxiT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate 3D support */
	if (!p) p = new SSMatSupportT(3, NumDOF(), NumIP());

	/* inherited initializations */
	SmallStrainT::NewMaterialSupport(p);

	return p;
}

/* return a pointer to a new material list */
MaterialListT* SmallStrainAxiT::NewMaterialList(int nsd, int size)
{
#pragma unused(nsd)
#pragma message("not needed with ParameterListT input")

	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSMatSupport) {
			fSSMatSupport = TB_DYNAMIC_CAST(SSMatSupportT*, NewMaterialSupport());
			if (!fSSMatSupport)
				ExceptionT::GeneralFail("SmallStrainAxiT::NewMaterialList");
		}

		return new SolidMatList3DT(size, *fSSMatSupport);
	}
	else
		return new SolidMatList3DT;
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
