/* $Id: MixtureSpeciesT.cpp,v 1.9 2005-01-14 00:20:30 paklein Exp $ */
#include "MixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "ShapeFunctionT.h"
#include "NLDiffusionMaterialT.h"

//DEBUG
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
MixtureSpeciesT::MixtureSpeciesT(const ElementSupportT& support):
	NLDiffusionElementT(support),
	fGradientOption(kGlobalProjection),
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

	/* gradient option */
	ParameterT grad_opt(ParameterT::Enumeration, "stress_gradient_option");
	grad_opt.AddEnumeration("global_projection", kGlobalProjection);
	grad_opt.AddEnumeration("element_projection", kElementProjection);
	grad_opt.SetDefault(fGradientOption);
	list.AddParameter(grad_opt);
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

	/* method used to compute stress gradient */
	int grad_opt = list.GetParameter("stress_gradient_option");
	if (grad_opt == kGlobalProjection)
		fGradientOption = kGlobalProjection;
	else if (grad_opt == kElementProjection)
		fGradientOption = kElementProjection;
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"stress_gradient_option\" %d", grad_opt);

	/* allocate work space for element-by-element stress projection */
	if (fGradientOption == kElementProjection) 
	{
		/* parent domain information */
		const ParentDomainT& parent_domain = fShapes->ParentDomain();
	
		/* dimensions */
		int nip = NumIP();
		int nsd = NumSD();
		
		/* allocate/initialize */
		fP_ip.Dimension(nip);
		fdP_ip.Dimension(nip);
		fip_gradient.Dimension(nip);
		for (int i = 0; i < nip; i++) {
			fP_ip[i].Dimension(nsd);
			fdP_ip[i].Dimension(nsd);
			fip_gradient[i].Dimension(nsd,nip);
			parent_domain.IPGradientTransform(i, fip_gradient[i]);
		}
	}

	/* resolve species index */
	fIndex = fUpdatedLagMixture->SpeciesIndex(Field().FieldName());
	if (fIndex < 0)
		ExceptionT::GeneralFail(caller, "could not resolve index of field \"%s\"",
			Field().FieldName().Pointer());

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
	ComputeMassFlux(false);

	/* inherited */
	NLDiffusionElementT::RHSDriver();
}

/* form group contribution to the stiffness matrix */
void MixtureSpeciesT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	/* compute the variation in flux velocities */
	ComputeMassFlux(true);

	/* inherited */
	NLDiffusionElementT::LHSDriver(sys_type);
}

/* calculate the internal force contribution ("-k*d") */
void MixtureSpeciesT::FormKd(double constK)
{
/* include effects from:
 * (1) source/sink terms
 * (2) divergence of flux
 */

	/* dimensions */
	int nsd = NumSD();
	int nip = NumIP();

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* mass flux */
	dArray2DT M_e(nip, nsd, fMassFlux(CurrElementNumber()));
	dArrayT M;
	
	dMatrixT grad;
	dArrayT field;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
	
		/* retrieve the mass flux */
		M_e.RowAlias(ip, M);
	
		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Set(1, fField_list.Pointer(CurrIP()));
		IP_Interpolate(fLocDisp, field);

		/* get strain-displacement matrix */
		B(fShapes->CurrIP(), fB);

		/* (div) flux contribution */
		fB.MultTx(M, fNEEvec);

		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
}

/* form the element stiffness matrix */
void MixtureSpeciesT::FormStiffness(double constK)
{
	/* must be nonsymmetric */
	if (fLHS.Format() != ElementMatrixT::kNonSymmetric)
		ExceptionT::GeneralFail("MixtureSpeciesT::FormStiffness",
			"LHS matrix must be nonsymmetric");

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* dimensions */
	int nsd = NumSD();
	int nip = NumIP();
	int nen = NumElementNodes();	

	/* (linearization) mass flux */
	dArray2DT DM_e(nip, nsd, fDMassFlux(CurrElementNumber()));
	dArrayT DM;
	
	/* integrate element stiffness */
	dMatrixT grad;
	dArrayT field;
	dArrayT Na;
	fShapes->TopIP();
	dArrayT dfield(1);
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		
		/* retrieve the mass flux derivative */
		DM_e.RowAlias(ip, DM);
	
		double scale = constK*(*Det++)*(*Weight++);

		/* set field gradient */
		grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Set(1, fField_list.Pointer(CurrIP()));
		IP_Interpolate(fLocDisp, field);

		/* shape function array */
		Na.Set(nen, (double*) fShapes->IPShapeU());

		/* strain displacement matrix */
		B(fShapes->CurrIP(), fB);
	
		/* (divergence) mass flux contribution */
		fB.MultTx(DM, fNEEvec);
		fLHS.Outer(fNEEvec, Na, scale, dMatrixT::kAccumulate);
	}
}

/* compute the flux velocities */
void MixtureSpeciesT::ComputeMassFlux(bool compute_dmass_flux)
{
	const char caller[] = "MixtureSpeciesT::ComputeMassFlux";
	
	/* work space */
	int nsd = NumSD();
	int nip = NumIP();
	dArrayT force(nsd);
	dMatrixT F_inv(nsd);

	/* project partial stresses to the nodes */
	LocalArrayT P(LocalArrayT::kUnspecified), dP(LocalArrayT::kUnspecified);
	dArray2DT dP_avg;
	dMatrixT ip_Grad_P, ip_Grad_P_j;
	dMatrixT ip_grad_X;
	if (fGradientOption == kGlobalProjection)
	{	
		/* get nodal stresses (PK1) */
		fUpdatedLagMixture->ProjectPartialStress(fIndex);
		fP_avg = ElementSupport().OutputAverage();
		P.Dimension(NumElementNodes(), nsd*nsd);		
		P.SetGlobal(fP_avg);

		/* project variation in partial stresses to the nodes */
		if (compute_dmass_flux) {
			fUpdatedLagMixture->ProjectDPartialStress(fIndex);
			dP_avg.Alias(ElementSupport().OutputAverage());
			dP.Dimension(NumElementNodes(), nsd*nsd);		
			dP.SetGlobal(dP_avg);
		}

		ip_Grad_P.Dimension(nsd*nsd, nsd);
	}
	else {
		ip_grad_X.Dimension(nsd, nip);
	}
	
	/* get the body force */
	dArrayT body_force(nsd), divP(nsd), vec(nsd);
	fUpdatedLagMixture->BodyForce(body_force);

	/* dimension work space */
	if (compute_dmass_flux)
		fDMassFlux.Dimension(NumElements(), nip*nsd);	

	/* element values */
	LocalArrayT acc(LocalArrayT::kAcc, NumElementNodes(), nsd);

	/* integration point values */
	dArrayT ip_conc(1);
	dArrayT ip_acc(nsd);
	
	dArray2DT V_e, M_e, dM_e;
	dArrayT V, M, dM;

	Top();
	fUpdatedLagMixture->Top();
	while (NextElement()) 
	{
		int e = CurrElementNumber();

		/* global shape function values */
		SetGlobalShape();
	
		/* set solid element */
		fUpdatedLagMixture->NextElement();
		if (fGradientOption == kGlobalProjection)
		{
			/* compute shape functions */
			fUpdatedLagMixture->SetGlobalShape();
			
			/* collect nodal stresses */
			SetLocalU(P);
			if (compute_dmass_flux) SetLocalU(dP);
		}

		/* collect integration point stresses - sets shapes functions over the element */
		if (fGradientOption == kElementProjection)
			fUpdatedLagMixture->IP_PartialStress(fIndex, &fP_ip, (compute_dmass_flux) ? &fdP_ip : NULL);

		/* collect nodal accelerations */
		fUpdatedLagMixture->Acceleration(acc);

		/* collect nodal concentrations */
		SetLocalU(fLocDisp);
		
		/* mass flux and velocity */
		V_e.Alias(nip, nsd, fFluxVelocity(e));
		M_e.Alias(nip, nsd, fMassFlux(e));
		if (compute_dmass_flux) dM_e.Alias(nip, nsd, fDMassFlux(e));

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
						
			/* compute stress divergence */
			if (fGradientOption == kGlobalProjection)
			{
				divP = 0.0;
				IP_ComputeGradient(P, ip_Grad_P);
				for (int j = 0; j < nsd; j++) {
					ip_Grad_P_j.Alias(nsd, nsd, ip_Grad_P(j));
					for (int i = 0; i < nsd; i++)
						divP[i] += ip_Grad_P_j(i,j);
				}
			}
			else /* element-by-element gradient calculation */
			{
				/* transform gradient matrix to element (ref) coordinates */
				fShapes->ParentDomain().DomainJacobian(fLocInitCoords, ip, F_inv);
				F_inv.Inverse();
				ip_grad_X.MultATB(F_inv, fip_gradient[ip]);

				/* compute divergence of P */
				ComputeDivergence(ip_grad_X, fP_ip, divP);
			}
			
			/* add to force */
			force.AddScaled(1.0/ip_conc[0], divP);
			
			/* compute (scaled) relative flux velocity */
			const dMatrixT& D = fCurrMaterial->k_ij();
			V_e.RowAlias(ip, V);
			D.Multx(force, V); /* c*V */

			/* compute (scaled) flux velocity */
//			V.AddScaled(ip_conc[0], [velocity of background]);
// add motion of background

			/* compute mass flux */
			const dMatrixT& F = fUpdatedLagMixture->DeformationGradient(ip);
			F_inv.Inverse(F);
			M_e.RowAlias(ip, M);
			F_inv.Multx(V, M);
			
			/* compute velocity */
			V /= ip_conc[0];

			/* mass flux variation */
			if (compute_dmass_flux)
			{
				dM_e.RowAlias(ip, dM);
			
				/* contribution from changing diffusivity */
				const dMatrixT& dD = fCurrMaterial->dk_ij();
				dD.Multx(force, vec);
				F_inv.Multx(vec, dM);

				/* stress variation divergence */
				if (fGradientOption == kGlobalProjection)
				{
					divP /= -ip_conc[0]*ip_conc[0];
					IP_ComputeGradient(dP, ip_Grad_P);
					for (int j = 0; j < nsd; j++) {
						ip_Grad_P_j.Alias(nsd, nsd, ip_Grad_P(j));
						for (int i = 0; i < nsd; i++)
							divP[i] += ip_Grad_P_j(i,j)/ip_conc[0];
					}
					D.Multx(divP, vec);
				}
				else /* element-by-element gradient calculation */
				{
					/* div P contribution */
					D.Multx(divP, vec, -1.0/(ip_conc[0]*ip_conc[0]));

					/* compute divergence of P variation */
					ComputeDivergence(ip_grad_X, fdP_ip, divP);
				
					/* div dP contribution */
					D.Multx(divP, vec, 1.0/ip_conc[0], dMatrixT::kAccumulate);				
				}

				/* accumulate contritbution from stress variation */
				F_inv.Multx(vec, dM, 1.0, dMatrixT::kAccumulate);
			}
		}
	}	
}

#if 0
/* compute the flux velocities */
void MixtureSpeciesT::ComputeDMassFlux(void)
{
	const char caller[] = "MixtureSpeciesT::ComputeDMassFlux";

	/* global projections */
	dArray2DT dP_avg;
	if (fGradientOption == kGlobalProjection) 
	{
		/* project partial stresses to the nodes */
		fUpdatedLagMixture->ProjectPartialStress(fIndex);
		fP_avg = ElementSupport().OutputAverage();

		/* project variation in partial stresses to the nodes */
		fUpdatedLagMixture->ProjectDPartialStress(fIndex);
		dP_avg.Alias(ElementSupport().OutputAverage());
	}

	/* work space */
	int nsd = NumSD();
	int nip = NumIP();
	int nen = NumElementNodes();
	dArrayT force(nsd), vec(nsd), divP(nsd);
	dMatrixT F_inv(nsd);

	/* dimension work space */
	fDMassFlux.Dimension(NumElements(), nip*nsd);	
	
	/* get the body force */
	dArrayT body_force(nsd);
	fUpdatedLagMixture->BodyForce(body_force);

	/* element values */
	LocalArrayT acc(LocalArrayT::kAcc, nen, nsd);
	LocalArrayT P(LocalArrayT::kUnspecified, nen, nsd*nsd);
	P.SetGlobal(fP_avg);
	LocalArrayT dP(LocalArrayT::kUnspecified, nen, nsd*nsd);
	dP.SetGlobal(dP_avg);

	/* integration point values */
	dArrayT ip_conc(1);
	dArrayT ip_acc(nsd);
	dMatrixT ip_Grad_P(nsd*nsd, nsd), ip_Grad_P_j;
	
	dArray2DT V_e, M_e, dM_e;
	dArrayT V, M, dM;

	Top();
	fUpdatedLagMixture->Top();
	while (NextElement()) 
	{
		int e = CurrElementNumber();

		/* global shape function values */
		SetGlobalShape();
	
		/* set solid element */
		fUpdatedLagMixture->NextElement();
		fUpdatedLagMixture->SetGlobalShape();
	
		/* collect nodal accelerations */
		fUpdatedLagMixture->Acceleration(acc);

		/* collect nodal concentrations */
		SetLocalU(fLocDisp);

		/* collect nodal stress and variation */
		SetLocalU(P);
		SetLocalU(dP);
		
		/* mass flux and velocity */
		V_e.Alias(nip, nsd, fFluxVelocity(e));
		M_e.Alias(nip, nsd, fMassFlux(e));
		dM_e.Alias(nip, nsd, fDMassFlux(e));

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
			divP = 0.0;
			IP_ComputeGradient(P, ip_Grad_P);
			for (int j = 0; j < nsd; j++) {
				ip_Grad_P_j.Alias(nsd, nsd, ip_Grad_P(j));
				for (int i = 0; i < nsd; i++)
					divP[i] += ip_Grad_P_j(i,j);
			}
			
			/* stress divergence contribution */
			force.AddScaled(1.0/ip_conc[0], divP);

			/* compute (scaled) relative flux velocity */
			const dMatrixT& D = fCurrMaterial->k_ij();
			V_e.RowAlias(ip, V);
			D.Multx(force, V); /* c*V */

			/* compute (scaled) flux velocity */
//			V.AddScaled(ip_conc[0], [velocity of background]);
// add motion of background

			/* compute mass flux */
			const dMatrixT& F = fUpdatedLagMixture->DeformationGradient(ip);
			F_inv.Inverse(F);
			M_e.RowAlias(ip, M);
			F_inv.Multx(V, M);

			/* compute velocity */
			V /= ip_conc[0];			
			
			/* mass flux variation */
			dM_e.RowAlias(ip, dM);
			
			/* contribution from changing diffusivity */
			const dMatrixT& dD = fCurrMaterial->dk_ij();
			dD.Multx(force, vec);
			F_inv.Multx(vec, dM);

			/* stress variation divergence */
			divP /= -ip_conc[0]*ip_conc[0];
			IP_ComputeGradient(dP, ip_Grad_P);
			for (int j = 0; j < nsd; j++) {
				ip_Grad_P_j.Alias(nsd, nsd, ip_Grad_P(j));
				for (int i = 0; i < nsd; i++)
					divP[i] += ip_Grad_P_j(i,j)/ip_conc[0];
			}
			D.Multx(divP, vec);
			F_inv.Multx(vec, dM, 1.0, dMatrixT::kAccumulate);
		}
	}	
}
#endif

/* compute the divergence tensor field given the values at the integration points */
void MixtureSpeciesT::ComputeDivergence(const dMatrixT& ip_grad_transform, 
	const ArrayT<dMatrixT>& tensor_ip, dArrayT& div) const
{
	/* dimensions */
	int nsd = ip_grad_transform.Rows();
	int nip = ip_grad_transform.Cols();

	div = 0.0;
	for (int k = 0; k < nip; k++) {
		const dMatrixT& A_k = tensor_ip[k];
		for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++)
				div[i] += ip_grad_transform(j,k)*A_k(i,j);			
	}
}
