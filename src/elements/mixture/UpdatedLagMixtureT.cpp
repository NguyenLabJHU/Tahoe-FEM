/* $Id: UpdatedLagMixtureT.cpp,v 1.16 2006-01-04 17:40:39 thao Exp $ */
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
int UpdatedLagMixtureT::SpeciesIndex(const StringT& field_name) const {
	return FSSolidMixture().SpeciesIndex(field_name);
}

/* density of the given species */
double UpdatedLagMixtureT::Density(int i) {
	return FSSolidMixture().Density(i);
}

/* set concentration flag */
void UpdatedLagMixtureT::SetConcentration(int i, ConcentrationT conc)
{
	const char caller[] = "UpdatedLagMixtureT::SetConcentration";

	/* get material */
	FSSolidMixtureT& mixture = FSSolidMixture();

	/* set flag */
	if (conc == kReference)
		mixture.SetConcentration(i, FSSolidMixtureT::kReference);
	else if (conc == kCurrent)
		mixture.SetConcentration(i, FSSolidMixtureT::kCurrent);	
	else
		ExceptionT::GeneralFail(caller, "unrecognized flag %d", conc);
	
	/* (re-)set form of element stiffness matrix */
	GlobalT::SystemTypeT type = TangentType();
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);
}

/* project the given partial first Piola-Kirchoff stress to the nodes */
void UpdatedLagMixtureT::ProjectPartialStress(int i)
{
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
			/* get material */
			FSSolidMixtureT& mixture = FSSolidMixture();
		
			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_P = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& cauchy = mixture.s_ij(i);
				
				/* Cauchy -> 1st PK stress */
				cauchy.ToMatrix(fStress);
				const dMatrixT& F = DeformationGradient();
				fF_inv.Inverse(F);
				P.MultABT(fStress, fF_inv);
				P *= F.Det();
//                cout << "\n Elem: "<<CurrElementNumber()
//                     << "\t IP: "<<CurrIP()
//                     << "\n P: "<<P;
                     
				/* extrapolate to the nodes */
				fShapes->Extrapolate(P_1D, nodal_P);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_P);
		}
}

/* project the given partial first Piola-Kirchoff stress to the nodes */
void UpdatedLagMixtureT::ProjectPartialTau(int i)
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();
	
	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);
	
	/* work space */
	dMatrixT tau(nsd);
	dArrayT tau_1D;
	tau_1D.Alias(tau);

	/* loop over elements */
	dArray2DT nodal_tau(nen, nsd*nsd);

	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get material */
			FSSolidMixtureT& mixture = FSSolidMixture();
			
			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);
			
			/* extrapolate element stresses */
			nodal_tau = 0.0;
			fCurrShapes->TopIP();
			while (fCurrShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& stress = mixture.specific_tau_ij(i);
                stress.ToMatrix(tau);
                
                const dArrayT& conc = mixture.Get_IPConcentration();	
                dMatrixT P=tau;
                P *= conc[i];
//                cout << "\n Elem: "<<CurrElementNumber()
//                     << "\t IP: "<<CurrIP()
//                     << "\n tau: "<<tau
//                     << "\n P: "<<P;

				/* extrapolate to the nodes */
				fCurrShapes->Extrapolate(tau_1D, nodal_tau);
			}
			
			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_tau);
		}
}

/* project the variation with concentration of the given partial first
 * Piola-Kirchoff stress to the nodes */
void UpdatedLagMixtureT::ProjectDPartialStress(int i)
{
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
			FSSolidMixtureT& mixture = FSSolidMixture();

			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_P = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture.ds_ij_dc(i);
				
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
		FSSolidMixtureT& mixture = FSSolidMixture();

		/* collect concentration */
		mixture.UpdateConcentrations(i);

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
				const dSymMatrixT& cauchy = mixture.s_ij(i);
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
				const dSymMatrixT& dcauchy = mixture.ds_ij_dc(i);
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

void UpdatedLagMixtureT::IP_PartialTau(int i, ArrayT<dMatrixT>* ip_stress)
{
	/* nothing wanted */
	if (!ip_stress)
		return;
	/* element is active */
	else if (CurrentElement().Flag() != ElementCardT::kOFF)
	{
		/* get materials */
		FSSolidMixtureT& mixture = FSSolidMixture();
		
		/* collect concentration */
		mixture.UpdateConcentrations(i);
		
		/* collect integration point element stresses */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* destination */
			int ip = fShapes->CurrIP();			

        
            /* stress */
			if (ip_stress)
			{
				/* Cauchy stress */
				const dSymMatrixT& tau = mixture.specific_tau_ij(i);
 				dMatrixT& mat_stress = (*ip_stress)[ip];
				tau.ToMatrix(mat_stress);
			}
		}
	}
	else /* zero them out */
	{
		if (ip_stress)
			for (int i = 0; i < ip_stress->Length(); i++)
				(*ip_stress)[i] = 0.0;
	}
}

const dMatrixT& UpdatedLagMixtureT::IP_PartialTau(int i,int ip)
{
	/* element is active */
	if (CurrentElement().Flag() != ElementCardT::kOFF)
	{
		/* get materials */
		FSSolidMixtureT& mixture = FSSolidMixture();
		
		/* collect concentration */
		mixture.UpdateConcentrations(i);
		
        fShapes->SetIP(ip);

        /* Cauchy stress */
		mixture.specific_tau_ij(i).ToMatrix(fStress);
	}
	else fStress = 0.0;
	
	return (fStress);
}

/* return the nodal accelerations over the current element */
void UpdatedLagMixtureT::Acceleration(LocalArrayT& acc) const
{
	if (fIntegrator->Order() == 2) {
		acc.SetGlobal(fLocAcc.Global());
		SetLocalU(acc);
	}
	else
		acc = 0.0;
}

/* return the nodal velocities over the current element */
void UpdatedLagMixtureT::Velocity(LocalArrayT& vel) const
{
	/* should have been set during SetGlobalShape because FSSolidMixtureT
	 * "needs" velocity information */
	vel = fLocVel;
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

/***********************************************************************
 * Protected
 ***********************************************************************/
const FSSolidMixtureT& UpdatedLagMixtureT::FSSolidMixture(void) const 
{
	const char caller[] = "UpdatedLagMixtureT::FSSolidMixture";

#if __option(extended_errorcheck)
	if (fMaterialList->Length() > 1)
		ExceptionT::GeneralFail(caller, "expecting only 1 material %d",
			fMaterialList->Length());
#endif

	const ContinuumMaterialT* pcont_mat = (*fMaterialList)[0]; /* just use first material */
	if (!pcont_mat) 
		ExceptionT::GeneralFail(caller, "material 0 is NULL");

	const FSSolidMixtureT* mixture = TB_DYNAMIC_CAST(const FSSolidMixtureT*, pcont_mat);
	if (!mixture)
		ExceptionT::GeneralFail(caller, "material \"%s\" is not a mixture",
			pcont_mat->Name().Pointer());

	return *mixture;
}
