/* $Id: UpdatedLagMixtureT.cpp,v 1.4 2005-01-04 00:52:16 paklein Exp $ */
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

/* project the Cauchy stress for the given species to the nodes */
void UpdatedLagMixtureT::ProjectPartialStress(int i)
{
	const char caller[] = "UpdatedLagMixtureT::ProjectPartialStress";

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
		if (CurrentElement().Flag() != kOFF)
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
				cauchy.ToMatrix(s);
				const dMatrixT& F = DeformationGradient();
				F_inv.Inverse(F);
				P.MultABT(s, F_inv);
				P *= F.Det();

				/* extrapolate to the nodes */
				fShapes->Extrapolate(P_1D, nodal_P);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_P);
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
