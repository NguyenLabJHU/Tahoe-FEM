/* $Id: MixtureSpeciesT.cpp,v 1.4 2004-11-10 02:00:13 paklein Exp $ */
#include "MixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "ShapeFunctionT.h"
#include "NLDiffusionMaterialT.h"

using namespace Tahoe;

/* constructor */
MixtureSpeciesT::MixtureSpeciesT(const ElementSupportT& support):
	NLDiffusionElementT(support),
	fUpdatedLagMixture(NULL),
	fBackgroundSpecies(NULL),
	fIndex(-1)
{
	SetName("mixture_species");
}

/* describe the parameters needed by the interface */
void MixtureSpeciesT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLDiffusionElementT::DefineParameters(list);

	/* associated solid element group */
	list.AddParameter(ParameterT::Integer, "solid_element_group");

#if 0
	/* velocity of species is calculated wrt this reference frame */
	ParameterT frame(ParameterT::Word, "reference_frame");
	frame.SetDefault("global");
	species->AddParameter(frame);
#endif
}

/* accept parameter list */
void MixtureSpeciesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MixtureSpeciesT::TakeParameterList";

	/* inherited */
	NLDiffusionElementT::TakeParameterList(list);

	/* resolve background solid element group */
	int solid_element_group = list.GetParameter("solid_element_group");
	solid_element_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(solid_element_group);	
	fUpdatedLagMixture = TB_DYNAMIC_CAST(UpdatedLagMixtureT*, &element);
	if (!fUpdatedLagMixture)
		ExceptionT::GeneralFail(caller, "group %d is not a mixture", solid_element_group+1);
	
	/* checks */
	if (fUpdatedLagMixture->NumElements() != NumElements() ||
		fUpdatedLagMixture->NumElementNodes() != NumElementNodes() ||
		fUpdatedLagMixture->NumIP() != NumIP())
		ExceptionT::SizeMismatch(caller);

	/* resolve species index */
	const ArrayT<StringT>& labels = Field().Labels();
	if (labels.Length() != 1) 
		ExceptionT::GeneralFail(caller, "field dim 1 not %d", labels.Length());
	fIndex = fUpdatedLagMixture->SpeciesIndex(labels[0]);

	/* dimension */
	fFluxVelocity.Dimension(NumElements(), NumIP()*NumSD());
	fMassFlux.Dimension(NumElements(), NumIP()*NumSD());	
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form the residual force vector */
void MixtureSpeciesT::RHSDriver(void)
{
	/* compute the flux velocities */
	ComputeMassFlux();

	/* inherited */
	NLDiffusionElementT::RHSDriver();
}

/* calculate the internal force contribution ("-k*d") */
void MixtureSpeciesT::FormKd(double constK)
{
/* include effects from:
 * (1) source/sink terms
 * (2) divergence of flux
 */

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	int nsd = NumSD();
	dMatrixT grad;
	dArrayT field;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Set(1, fField_list.Pointer(CurrIP()));
		IP_Interpolate(fLocDisp, field);

		/* get strain-displacement matrix */
		B(fShapes->CurrIP(), fB);

		/* compute heat flow */
//		fB.MultTx(fCurrMaterial->q_i(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}

/* compute the flux velocities */
void MixtureSpeciesT::ComputeMassFlux(void)
{
	const char caller[] = "MixtureSpeciesT::ComputeFluxVelocities";

	/* project partial stresses to the nodes */
	fUpdatedLagMixture->ProjectPartialStress(fIndex);

	/* get the array of nodal stresses (PK1) */
	const dArray2DT& P_avg = ElementSupport().OutputAverage();
	
	/* work space */
	int nsd = NumSD();
	int nip = NumIP();
	dArrayT force(nsd);
	dMatrixT F_inv(nsd);
	
	/* get the body force */
	dArrayT body_force(nsd);
	fUpdatedLagMixture->BodyForce(body_force);

	/* element values */
	LocalArrayT acc(LocalArrayT::kAcc, NumElementNodes(), nsd);
	LocalArrayT P(LocalArrayT::kUnspecified, NumElementNodes(), nsd*nsd);
	P.SetGlobal(P_avg);

	/* integration point values */
	dArrayT ip_conc(1);
	dArrayT ip_acc(nsd);
	dMatrixT ip_Grad_P(nsd*nsd, nsd), ip_Grad_P_j;
	
	dArray2DT V_e, M_e;
	dArrayT V, M;

	Top();
	fUpdatedLagMixture->Top();
	while (NextElement()) 
	{
		int e = CurrElementNumber();
	
		/* set solid element */
		fUpdatedLagMixture->NextElement();
		fUpdatedLagMixture->SetGlobalShape();
	
		/* collect nodal accelerations */
		fUpdatedLagMixture->Acceleration(acc);

		/* collect nodal concentrations */
		SetLocalU(fLocDisp);

		/* collect nodal stresses */
		SetLocalU(P);
		
		/* mass flux and velocity */
		V_e.Alias(nip, nsd, fFluxVelocity(e));
		M_e.Alias(nip, nsd, fMassFlux(e));

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			int ip = fShapes->CurrIP();
		
			/* ip values */
			IP_Interpolate(fLocDisp, ip_conc);		
			IP_Interpolate(acc, ip_acc);		
		
			/* inertial forces */
			force.SetToCombination(-1.0, ip_acc, 1.0, body_force);
						
			/* stress divergence */
			IP_ComputeGradient(P, ip_Grad_P);
			for (int j = 0; j < nsd; j++) {
				ip_Grad_P_j.Alias(nsd, nsd, ip_Grad_P(j));
				for (int i = 0; i < nsd; i++)
					force[i] += ip_Grad_P_j(i,j)/ip_conc[0];
			}
			
			
			/* compute (scaled) relative flux velocity */
			const dMatrixT& D = fCurrMaterial->k_ij();
			V_e.RowAlias(ip, V);
			D.Multx(force, V); /* c*V */

			/* compute (scaled) flux velocity */
//			V.AddScaled(ip_conc[0], [velocity of background]);
// add motion of background

			/* compute mass flux */
			const dMatrixT& F = fUpdatedLagMixture->DeformationGradient();
			F_inv.Inverse(F);
			M_e.RowAlias(ip, M);
			F_inv.Multx(V, M);
			
			/* compute velocity */
			V /= ip_conc[0];
		}
	}	
}
