/* $Id: UpdatedLagMixtureT.cpp,v 1.2 2004-11-07 17:09:26 paklein Exp $ */
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
	int nst = dSymMatrixT::NumValues(nsd);

	/* reset averaging workspace */
	ElementSupport().ResetAverage(nst);

	/* loop over elements */
	dArray2DT nodalstress(nen, nst);
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
			nodalstress = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* get Cauchy stress */
				const dSymMatrixT& cauchy = mixture->s_ij(i);

				/* extrapolate to the nodes */
				fShapes->Extrapolate(cauchy, nodalstress);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodalstress);
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
	body_force.SetToScaled(fBodySchedule->Value(), fBody);
}
