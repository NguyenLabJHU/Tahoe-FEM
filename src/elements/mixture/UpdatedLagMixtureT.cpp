/* $Id: UpdatedLagMixtureT.cpp,v 1.7 2005-02-16 21:37:47 paklein Exp $ */
#include "UpdatedLagMixtureT.h"
#include "ShapeFunctionT.h"
#include "FSSolidMixtureT.h"
#include "SolidMatListT.h"
#include "ScheduleT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

/* constructor */
UpdatedLagMixtureT::UpdatedLagMixtureT(const ElementSupportT& support):
	UpdatedLagrangianT(support)
{
	SetName("updated_lagrangian_mixture");
}

/* resolve the species name into the index */
int UpdatedLagMixtureT::SpeciesIndex(const StringT& field_name) const
{
	/* get materials */
	const ContinuumMaterialT* pcont_mat = (*fMaterialList)[0]; /* just use first material */
	const FSSolidMixtureT* mixture = TB_DYNAMIC_CAST(const FSSolidMixtureT*, pcont_mat);
	if (!mixture) ExceptionT::GeneralFail("UpdatedLagMixtureT::SpeciesIndex", "material is not a mixture");

	/* resolve index */
	return mixture->SpeciesIndex(field_name);
}

/* project the given partial first Piola-Kirchoff stress to the nodes */
void UpdatedLagMixtureT::ProjectPartialStress(int i)
{
	const char caller[] = "UpdatedLagMixtureT::ProjectPartialStress";

	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);

	/* work space */
	dMatrixT P(nsd);
	dArrayT P_1D;
	P_1D.Alias(P);

	/* loop over elements */
	dArray2DT nodal_P(nen, nsd*nsd);
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get materials */
			FSSolidMixtureT* mixture = TB_DYNAMIC_CAST(FSSolidMixtureT*, fCurrMaterial);
			if (!mixture) ExceptionT::GeneralFail(caller, "material is not a mixture");
		
			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture->UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_P = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& cauchy = mixture->s_ij(i);
				
				/* Cauchy -> 1st PK stress */
				cauchy.ToMatrix(fStress);
				const dMatrixT& F = DeformationGradient();
				fF_inv.Inverse(F);
				P.MultABT(fStress, fF_inv);
				P *= F.Det();

				/* extrapolate to the nodes */
				fShapes->Extrapolate(P_1D, nodal_P);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_P);
		}
}

/* project the variation with concentration of the given partial first
 * Piola-Kirchoff stress to the nodes */
void UpdatedLagMixtureT::ProjectDPartialStress(int i)
{
	const char caller[] = "UpdatedLagMixtureT::ProjectDPartialStress";

	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);

	/* work space */
	dMatrixT P(nsd), F_inv(nsd), s(nsd);
	dArrayT P_1D;
	P_1D.Alias(P);

	/* loop over elements */
	dArray2DT nodal_P(nen, nsd*nsd);
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get materials */
			FSSolidMixtureT* mixture = TB_DYNAMIC_CAST(FSSolidMixtureT*, fCurrMaterial);
			if (!mixture) ExceptionT::GeneralFail(caller, "material is not a mixture");
		
			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture->UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_P = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture->ds_ij_dc(i);
				
				/* Cauchy -> 1st PK stress */
				dcauchy.ToMatrix(fStress);
				const dMatrixT& F = DeformationGradient();
				fF_inv.Inverse(F);
				P.MultABT(fStress, fF_inv);
				P *= F.Det();

				/* extrapolate to the nodes */
				fShapes->Extrapolate(P_1D, nodal_P);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_P);
		}
}

void UpdatedLagMixtureT::IP_PartialStress(int i, ArrayT<dMatrixT>* ip_stress, 
	ArrayT<dMatrixT>* ip_dstress)
{
	/* nothing wanted */
	if (!ip_stress && !ip_dstress)
		return;
	/* element is active */
	else if (CurrentElement().Flag() != ElementCardT::kOFF)
	{
		/* get materials */
		FSSolidMixtureT* mixture = TB_DYNAMIC_CAST(FSSolidMixtureT*, fCurrMaterial);
		if (!mixture) 
			ExceptionT::GeneralFail("UpdatedLagMixtureT::IP_PartialStress", 
				"material is not a mixture");
	
		/* global shape function values */
		SetGlobalShape();
		
		/* collect concentration */
		mixture->UpdateConcentrations(i);

		/* collect integration point element stresses */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* destination */
			int ip = fShapes->CurrIP();
		
			/* deformation gradient */
			const dMatrixT& F = DeformationGradient();
			fF_inv.Inverse(F);

			/* stress */
			if (ip_stress)
			{
				/* Cauchy stress */
				const dSymMatrixT& cauchy = mixture->s_ij(i);
				cauchy.ToMatrix(fStress);
			
				/* Cauchy -> 1st PK stress */
				dMatrixT& P = (*ip_stress)[ip];
				P.MultABT(fStress, fF_inv);
				P *= F.Det();
			}

			/* stress variation */
			if (ip_dstress)
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture->ds_ij_dc(i);
				dcauchy.ToMatrix(fStress);
			
				/* Cauchy -> 1st PK stress */
				dMatrixT& dP = (*ip_dstress)[ip];
				dP.MultABT(fStress, fF_inv);
				dP *= F.Det();
			}
		}
	}
	else /* zero them out */
	{
		if (ip_stress)
			for (int i = 0; i < ip_stress->Length(); i++)
				(*ip_stress)[i] = 0.0;

		if (ip_dstress)
			for (int i = 0; i < ip_dstress->Length(); i++)
				(*ip_dstress)[i] = 0.0;
	}
}

/* return the nodal accelerations over the current element */
void UpdatedLagMixtureT::Acceleration(LocalArrayT& acc)
{
	if (fIntegrator->Order() == 2) {
		SetLocalU(fLocAcc);
		acc = fLocAcc;
	}
	else
		acc = 0.0;
}

/* return the body force vector */
void UpdatedLagMixtureT::BodyForce(dArrayT& body_force) const
{
	if (fBodySchedule)
		body_force.SetToScaled(fBodySchedule->Value(), fBody);
	else
		body_force = 0.0;
}

/* accept parameter list */
void UpdatedLagMixtureT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	UpdatedLagrangianT::TakeParameterList(list);

	/* dimension work space */
	int nsd = NumSD();
	fF_inv.Dimension(nsd);
	fStress.Dimension(nsd);
}
