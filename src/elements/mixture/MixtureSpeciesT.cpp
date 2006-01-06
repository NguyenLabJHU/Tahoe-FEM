/* $Id: MixtureSpeciesT.cpp,v 1.23 2006-01-06 02:55:57 thao Exp $ */
#include "MixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "Q1P0MixtureT.h"
#include "ShapeFunctionT.h"
#include "NLDiffusionMaterialT.h"
#include "MaterialListT.h"

//DEBUG
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
MixtureSpeciesT::MixtureSpeciesT(const ElementSupportT& support):
	NLDiffusionElementT(support),
	fGradientOption(kGlobalProjection),
	fConcentration(kReference),
	fOutputMass(false),
	fUpdatedLagMixture(NULL),
	fQ1P0Mixture(NULL),
	fBackgroundSpecies(NULL),
	fIndex(-1),
	fLocCurrCoords(LocalArrayT::kCurrCoords)
{
	SetName("mixture_species");
}

/* write element output */
void MixtureSpeciesT::WriteOutput(void)
{
	/* inherited */
	NLDiffusionElementT::WriteOutput();

	if (fOutputMass)
	{
		/* compute total species mass */
		double mass = 0.0;
		dArrayT ip_conc(1);
		Top();
		while (NextElement()) 
		{
			/* global shape function values */
			SetGlobalShape();
	
			/* collect nodal concentrations */
			SetLocalU(fLocDisp);
		
			/* loop over integration points */
			const double* j = fShapes->IPDets();
			const double* w = fShapes->IPWeights();
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* ip values */
				IP_Interpolate(fLocDisp, ip_conc);		

				/* accumulate */
				mass += (*j++)*(*w++)*ip_conc[0];
			}
		}
	
		/* output */
		ofstreamT& out = ElementSupport().Output();
		int d_width = OutputWidth(out, &mass);
		out << '\n'
		    << setw(d_width) << ElementSupport().Time()
		    << setw(d_width) << mass
		    << ": time, mass of \"" << Field().FieldName() << "\"" << '\n';
	}
}

/* describe the parameters needed by the interface */
void MixtureSpeciesT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLDiffusionElementT::DefineParameters(list);

	/* associated solid element group */
	list.AddParameter(ParameterT::Integer, "solid_element_group");

	/* velocity of species is calculated wrt this reference frame */
	ParameterT frame(ParameterT::Word, "background_species");
	list.AddParameter(frame, ParameterListT::ZeroOrOnce);

	/* gradient option */
	ParameterT grad_opt(ParameterT::Enumeration, "stress_gradient_option");
	grad_opt.AddEnumeration("global_projection", kGlobalProjection);
	grad_opt.AddEnumeration("element_projection", kElementProjection);
	grad_opt.SetDefault(fGradientOption);
	list.AddParameter(grad_opt);

	/* concentration type */
	ParameterT conc_opt(ParameterT::Enumeration, "concentration");
	conc_opt.AddEnumeration("reference_concentration", kReference);
	conc_opt.AddEnumeration("current_concentration", kCurrent);
	conc_opt.SetDefault(fConcentration);
	list.AddParameter(conc_opt);

	/* output total species mass */
	ParameterT output_mass(fOutputMass, "output_mass");
	output_mass.SetDefault(fOutputMass);
	list.AddParameter(output_mass);
}

/* accept parameter list */
void MixtureSpeciesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MixtureSpeciesT::TakeParameterList";

	/* concentration */
	int conc_opt = list.GetParameter("concentration");
	if (conc_opt == kReference)
		fConcentration = kReference;
	else if (conc_opt == kCurrent)
		fConcentration = kCurrent;
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"concentration\" %d", conc_opt);

	/* inherited */
	NLDiffusionElementT::TakeParameterList(list);

	/* resolve background solid element group */
	int solid_element_group = list.GetParameter("solid_element_group");
	solid_element_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(solid_element_group);	
	fUpdatedLagMixture = TB_DYNAMIC_CAST(UpdatedLagMixtureT*, &element);
	fQ1P0Mixture       = TB_DYNAMIC_CAST(Q1P0MixtureT*, &element);
	if (!fUpdatedLagMixture && !fQ1P0Mixture)
		ExceptionT::GeneralFail(caller, "group %d \"%s\" is not a mixture", 
			solid_element_group+1, element.Name().Pointer());
	
	/* checks */
	if (fUpdatedLagMixture) {
		if (fUpdatedLagMixture->NumElements() != NumElements() ||
			fUpdatedLagMixture->NumElementNodes() != NumElementNodes() ||
			fUpdatedLagMixture->NumIP() != NumIP())
			ExceptionT::SizeMismatch(caller);
	} else {
		if (fQ1P0Mixture->NumElements() != NumElements() ||
			fQ1P0Mixture->NumElementNodes() != NumElementNodes() ||
			fQ1P0Mixture->NumIP() != NumIP())
			ExceptionT::SizeMismatch(caller);	
	}

	/* resolve background species */
	const ParameterT* bg_species = list.Parameter("background_species");
	if (bg_species)
	{
		/* resolve bg field name */
		StringT bg_species_name = *bg_species;
		if (bg_species_name == Field().FieldName())
			ExceptionT::GeneralFail(caller, "background_species must differ from this species \"%s\"",
				Field().FieldName().Pointer());
			
		/* find name */
		int num_groups = ElementSupport().NumElementGroups();
		for (int i = 0; !fBackgroundSpecies && i < num_groups; i++) {
			ElementBaseT& element = ElementSupport().ElementGroup(i);
			if (element.Field().FieldName() == bg_species_name) {
				fBackgroundSpecies = TB_DYNAMIC_CAST(MixtureSpeciesT*, &element);
			}
		}
		if (!fBackgroundSpecies)
			ExceptionT::GeneralFail(caller, "could not resolve background_species \"%s\"",
				bg_species_name.Pointer());
	}
	
	/* method used to compute stress gradient */
	int grad_opt = list.GetParameter("stress_gradient_option");
	if (grad_opt == kGlobalProjection)
		fGradientOption = kGlobalProjection;
	else if (grad_opt == kElementProjection)
		fGradientOption = kElementProjection;
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"stress_gradient_option\" %d", grad_opt);

	/* output the total system mass */
	fOutputMass = list.GetParameter("output_mass");

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
	fIndex = (fUpdatedLagMixture) ? 
		fUpdatedLagMixture->SpeciesIndex(Field().FieldName()) :
		fQ1P0Mixture->SpeciesIndex(Field().FieldName());
	if (fIndex < 0)
		ExceptionT::GeneralFail(caller, "could not resolve index of species \"%s\" in \"%s\"",
			Field().FieldName().Pointer(),
			((fUpdatedLagMixture) ? 
				fUpdatedLagMixture->Name().Pointer() : 
				fQ1P0Mixture->Name().Pointer()));

#if 0
	/* consistency check */
	double solid_density = fUpdatedLagMixture->Density(fIndex);
	for (int i = 0; i < fMaterialList->Length(); i++) {
		const ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		const DiffusionMaterialT* pdiff_mat = TB_DYNAMIC_CAST(const DiffusionMaterialT*, pcont_mat);
		if (!pdiff_mat) ExceptionT::GeneralFail(caller, "error resolving density of material %d", i+1);
		if (fabs(pdiff_mat->Density() - solid_density) > kSmall)
			ExceptionT::GeneralFail(caller, 
				"density %g of the species at index %d of the solid mixture differs from %g",
					solid_density, fIndex+1, pdiff_mat->Density());
	}
#endif

	/* set concentration type */
	if (fUpdatedLagMixture) {
		if (fConcentration == kReference)
			fUpdatedLagMixture->SetConcentration(fIndex, UpdatedLagMixtureT::kReference);
		else
			fUpdatedLagMixture->SetConcentration(fIndex, UpdatedLagMixtureT::kCurrent);
	} else {
		if (fConcentration == kReference)
			fQ1P0Mixture->SetConcentration(fIndex, Q1P0MixtureT::kReference);
		else
			fQ1P0Mixture->SetConcentration(fIndex, Q1P0MixtureT::kCurrent);	
	}

	/* dimension */
	fFluxVelocity.Dimension(NumElements(), NumIP()*NumSD());
	fMassFlux.Dimension(NumElements(), NumIP()*NumSD());
	fNEEmat.Dimension(NumElementNodes());
	fNSDmat1.Dimension(NumSD());
	fNSDmat2.Dimension(NumSD());
	fNSDmat3.Dimension(NumSD());
	
	/* initialize */
	fFluxVelocity = 0.0;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* allocate and initialize shape function objects */
void MixtureSpeciesT::SetShape(void)
{
	/* clear existing */
	delete fShapes;

	/* configuration */
	const LocalArrayT& coords = (fConcentration == kCurrent) ?
		fLocCurrCoords : fLocInitCoords;

	/* construct shape functions */
	fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), coords);
	fShapes->Initialize();
}

/* allocate and initialize local arrays */
void MixtureSpeciesT::SetLocalArrays(void)
{
	/* inherited */
	NLDiffusionElementT::SetLocalArrays();
	
	/* set up array of current nodal coordinates */
	if (fConcentration == kCurrent) {
		fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
		ElementSupport().RegisterCoordinates(fLocCurrCoords);
	}
}

/* compute shape functions and derivatives */
void MixtureSpeciesT::SetGlobalShape(void)
{
	/* collect coordinates */
	if (fConcentration == kCurrent) SetLocalX(fLocCurrCoords);
	
	/* inherited */
	NLDiffusionElementT::SetGlobalShape();
	
	/* will need deformation gradient */
	if (fConcentration == kCurrent || fGradientOption == kGlobalProjection) {
		if (fUpdatedLagMixture)	
			fUpdatedLagMixture->SetGlobalShape();
		else
			fQ1P0Mixture->SetGlobalShape();
	}
}

/* reset loop */
void MixtureSpeciesT::Top(void)
{
	/* inherited */
	NLDiffusionElementT::Top();

	/* synchronize solid element group */
	if (fUpdatedLagMixture)
		fUpdatedLagMixture->Top();
	else
		fQ1P0Mixture->Top();
}
	
/* advance to next element */ 
bool MixtureSpeciesT::NextElement(void)
{
	/* inherited */
	bool next = NLDiffusionElementT::NextElement();

	/* synchronize solid element group */
	if (fUpdatedLagMixture)
		return fUpdatedLagMixture->NextElement() && next;
	else
		return fQ1P0Mixture->NextElement() && next;
}

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
		grad.Alias(1, nsd, fGradient_list[ip].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
		field.Alias(1, fField_list.Pointer(ip));
		IP_Interpolate(fLocDisp, field);

		/* get strain-displacement matrix */
		B(ip, fB);

		/* (div) flux contribution */
		fB.MultTx(M, fNEEvec);

		/* c div(v) contribution */
		if (fConcentration == kCurrent) 
		{
			/* deformation gradients */
			const dMatrixT& F = (fUpdatedLagMixture) ? 
				fUpdatedLagMixture->DeformationGradient(ip) : 
				fQ1P0Mixture->DeformationGradient(ip);
			const dMatrixT& F_last = (fUpdatedLagMixture) ?
				fUpdatedLagMixture->DeformationGradient_last(ip) :
				fQ1P0Mixture->DeformationGradient_last(ip);
			
			/* compute h and h^T h */
			fNSDmat1.DiffOf(F, F_last);           /* GRAD U (8.1.9) */
			fNSDmat2.Inverse(F);                  /* F^-1 */
			fNSDmat3.MultAB(fNSDmat1, fNSDmat2);  /* h (8.1.7) */
			fNSDmat1.MultATB(fNSDmat3, fNSDmat3); /* h^T h */

			/* compute velocity gradient (Simo: 8.1.22) and (Simo: 8.3.13) */
			double dt = ElementSupport().TimeStep();
			double by_dt = (dt > kSmall) ? 1.0/dt : 0.0;
            /*rate of deformation tensor, symmetric part of velocity gradient*/
            fNSDmat2.SetToCombination(by_dt, fNSDmat3, -0.5*by_dt, fNSDmat1); 
			/* c div(v) */
			double c_div_v = field[0]*fNSDmat2.Trace();
			fNEEvec.AddScaled(-c_div_v, fShapes->IPShapeU());
		}

		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
    if (CurrElementNumber() == 0)
        cout << "\nfRHS: "<<fRHS;
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
	dArray2DT DM_e(nip, nsd*nen, fDMassFlux(CurrElementNumber()));
	dMatrixT DM;
	
	/* integrate element stiffness */
	dMatrixT grad;
	dArrayT field;
	dArrayT Na;
	fShapes->TopIP();
	dArrayT dfield(1);
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		double scale = constK*(*Det++)*(*Weight++);
		
		/* retrieve the mass flux derivative */
		DM.Alias(nsd, nen, DM_e(ip));

		/* set field gradient */
//		grad.Alias(1, nsd, fGradient_list[CurrIP()].Pointer());
//		IP_ComputeGradient(fLocDisp, grad);
		
		/* interpolate field */
//		field.Alias(1, fField_list.Pointer(CurrIP()));
//		IP_Interpolate(fLocDisp, field);

		/* strain displacement matrix */
		B(ip, fB);
	
		/* (divergence) mass flux contribution */
		fNEEmat.MultATB(fB, DM);
		fLHS.AddScaled(scale, fNEEmat);
	}
    if (CurrElementNumber() == 0)
        cout << "\nfLHS: "<<fLHS;
}

/* compute the flux velocities */
void MixtureSpeciesT::ComputeMassFlux(bool compute_dmass_flux)
{
	const char caller[] = "MixtureSpeciesT::ComputeMassFlux";
	
	/* work space */
	int nen = NumElementNodes();
	int nsd = NumSD();
	int nip = NumIP();
	dArrayT force(nsd);
	dMatrixT F_inv(nsd);
	dMatrixT jacobian(nsd);

	/* project partial stresses to the nodes */
	LocalArrayT P(LocalArrayT::kUnspecified), dP(LocalArrayT::kUnspecified);
	dArray2DT dP_avg;
	dMatrixT ip_Grad_P, ip_Grad_P_j;
	dMatrixT ip_grad_X;
	if (fGradientOption == kGlobalProjection)
	{	
		/* get nodal stresses (PK1) */
		if (fUpdatedLagMixture) fUpdatedLagMixture->ProjectPartialStress(fIndex);
		else fQ1P0Mixture->ProjectPartialStress(fIndex);
		
		fP_avg = ElementSupport().OutputAverage();
		P.Dimension(NumElementNodes(), nsd*nsd);		
		P.SetGlobal(fP_avg);
//        cout<<"\nnodal P: "<<P;

		/* project variation in partial stresses to the nodes */
		if (compute_dmass_flux) {
			if (fUpdatedLagMixture) fUpdatedLagMixture->ProjectDPartialStress(fIndex);
			else fQ1P0Mixture->ProjectDPartialStress(fIndex);
			dP_avg.Alias(ElementSupport().OutputAverage());
			dP.Dimension(NumElementNodes(), nsd*nsd);		
			dP.SetGlobal(dP_avg);
            
//            cout << "\nnodal dP: "<<dP;
		}

		ip_Grad_P.Dimension(nsd*nsd, nsd);
	}
	else {
		ip_grad_X.Dimension(nsd, nip);
	}
	    
	/* get the body force */
	dArrayT body_force(nsd), divP(nsd), vec(nsd);
	dArrayT vec1(nsd), vec2(nsd);
	if (fUpdatedLagMixture) fUpdatedLagMixture->BodyForce(body_force);
	else fQ1P0Mixture->BodyForce(body_force);
	dMatrixT d_divP(nsd,nen), mat(nsd,nen), matnsd(nsd);

	/* initialize mass flux variation */
	if (compute_dmass_flux)
		fDMassFlux.Dimension(NumElements(), nip*nsd*nen);
	fDMassFlux = 0.0;

	/* element values */
	LocalArrayT acc(LocalArrayT::kAcc, NumElementNodes(), nsd);

	/* integration point values */
	dArrayT ip_conc(1);
	dArrayT ip_acc(nsd);
	
	dArray2DT V_e, M_e, dM_e, V_e_bg;
	dArrayT V, M, Na, V_bg;
	dMatrixT dM;

	Top();
	while (NextElement()) 
	{
		int e = CurrElementNumber();

		/* global shape function values */
		SetGlobalShape();
	
		/* collect nodal stresses */
		if (fGradientOption == kGlobalProjection) {
			SetLocalU(P);
			if (compute_dmass_flux) SetLocalU(dP);
		}

		/* collect integration point stresses - sets shapes functions over the element */
		if (fGradientOption == kElementProjection) {
			if (fUpdatedLagMixture)
				fUpdatedLagMixture->IP_PartialStress(fIndex, &fP_ip, (compute_dmass_flux) ? &fdP_ip : NULL);
			else
				fQ1P0Mixture->IP_PartialStress(fIndex, &fP_ip, (compute_dmass_flux) ? &fdP_ip : NULL);
		}

		/* collect nodal accelerations */
		if (fUpdatedLagMixture)
			fUpdatedLagMixture->Acceleration(acc);
		else
			fQ1P0Mixture->Acceleration(acc);

		/* collect nodal concentrations */
		SetLocalU(fLocDisp);
		
		/* mass flux and velocity */
		V_e.Alias(nip, nsd, fFluxVelocity(e));
		M_e.Alias(nip, nsd, fMassFlux(e));
		if (compute_dmass_flux) dM_e.Alias(nip, nen*nsd, fDMassFlux(e));

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			int ip = fShapes->CurrIP();

			/* deformation gradient */
			const dMatrixT& F = (fUpdatedLagMixture) ?
				fUpdatedLagMixture->DeformationGradient(ip) :
				fQ1P0Mixture->DeformationGradient(ip);
			double detF = F.Det();
			F_inv.Inverse(F);

			/* ip values */
			IP_Interpolate(fLocDisp, ip_conc);
			IP_Interpolate(acc, ip_acc);		

			/* convert to reference concentration */
			if (fConcentration == kCurrent)	ip_conc *= detF;
		
			/* inertial forces */
			force.SetToCombination(-1.0, ip_acc, 1.0, body_force);
                        
			/* compute stress divergence */
			if (fGradientOption == kGlobalProjection)
			{
				divP = 0.0;
//                cout << "\nP: "<<P;
//                cout << "\nnodal conc: "<<fLocDisp;
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
				fShapes->ParentDomain().DomainJacobian(fLocInitCoords, ip, jacobian);
				jacobian.Inverse();
				ip_grad_X.MultATB(jacobian, fip_gradient[ip]);

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
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				V_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
				V_e_bg.RowAlias(ip, V_bg);

				/* add (scaled) velocity of background */
				V.AddScaled(ip_conc[0], V_bg);
			}

			/* compute mass flux */
			M_e.RowAlias(ip, M);
			F_inv.Multx(V, M);
			
			/* compute velocity */
			V /= ip_conc[0];
//            cout << "\nvelocity: "<<V;
//            cout << "\nmass flux: "<<M;

/*            if (CurrElementNumber() == 0)
            {
                cout << "\n ip: "<< ip
                     << "\n ip_conc: "<<ip_conc[0]
                     << "\n force/conc: "<<force;
            }
*/
			/* mass flux variation */
			if (compute_dmass_flux)
			{
				dM.Alias(nsd, nen, dM_e(ip));

				/* shape function array */
				Na.Alias(nen, fShapes->IPShapeU());			
			
				/* contribution from changing diffusivity */
				const dMatrixT& dD = fCurrMaterial->dk_ij();
				dD.Multx(force, vec);

				/* div P contribution */
				D.Multx(divP, vec, -1.0/(ip_conc[0]*ip_conc[0]), dMatrixT::kAccumulate);
					
				/* accumulate */
				F_inv.Multx(vec, vec1);
				dM.Outer(vec1, Na, -1.0, dMatrixT::kAccumulate);

				/* stress variation divergence */
				if (fGradientOption == kGlobalProjection)
				{
					/* compute divergence of P variation */
					ComputeDDivergence(dP, d_divP, matnsd);
#if 0
					divP /= -ip_conc[0]*ip_conc[0];
					IP_ComputeGradient(dP, ip_Grad_P);
					for (int j = 0; j < nsd; j++) {
						ip_Grad_P_j.Alias(nsd, nsd, ip_Grad_P(j));
						for (int i = 0; i < nsd; i++)
							divP[i] += ip_Grad_P_j(i,j)/ip_conc[0];
					}
					D.Multx(divP, vec);
#endif
				}
				else /* element-by-element gradient calculation */
				{
					/* compute divergence of P variation */
					ComputeDDivergence(ip_grad_X, fdP_ip, d_divP);
				}

				/* accumulate */
				mat.MultABC(F_inv, D, d_divP);
				dM.AddScaled(-1.0/ip_conc[0], mat);
				
				/* transform to current configuration */
				if (fConcentration == kCurrent) {
					mat.SetToScaled(1.0/detF, dM);
					dM.MultAB(F, mat);
				}
			}
			
			/* transform to current configuration */
			if (fConcentration == kCurrent) {
				vec1.SetToScaled(1.0/detF, M);
				F.Multx(vec1, M);
			}
		}
	}	
}

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
			for (int j = 0; j < nsd; j++) /* div */
				div[i] += ip_grad_transform(j,k)*A_k(i,j);			
	}
}

/* compute variation of divergence with respect to the nodal values */
void MixtureSpeciesT::ComputeDDivergence(const dMatrixT& ip_grad_transform, 
	const ArrayT<dMatrixT>& tensor_ip, dMatrixT& d_div) const
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = ip_grad_transform.Rows();
	int nip = ip_grad_transform.Cols();	
	
	dArrayT Na;
	d_div = 0.0;
    dArrayT divdP(nsd);
    divdP = 0.0;
	for (int k = 0; k < nip; k++) 
	{
		/* shape function array */
		Na.Alias(nen, fShapes->IPShapeU(k));			
	
		const dMatrixT& A_k = tensor_ip[k];
		for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++) /* div */
				for (int n = 0; n < nen; n++)
					d_div(i,n) += ip_grad_transform(j,k)*A_k(i,j)*Na[n];

//        if (CurrElementNumber() == 0)
//            cout << "\n div dP ip: "<<A_k;
        for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++) /* div */
					divdP[i] += ip_grad_transform(j,k)*A_k(i,j);
	}
        if (CurrElementNumber() == 0) {
            cout << "\ndivdP: "<<divdP;
            cout << "\nd_div: "<<d_div;
        }
}

void MixtureSpeciesT::ComputeDDivergence(const LocalArrayT& nodal_dP, dMatrixT& d_div,
	dMatrixT& dP_ip) const
{
	/* dimensions */
	int nen = NumElementNodes();
	int nip = fShapes->NumIP();
	int nsd = NumSD();
	
	/* extrapolation matrix */
	const dMatrixT& extrap = fShapes->Extrapolation();
//    cout << "\nextrap: "<<extrap;
	/* shape function derivatives (at the current integration point) */
	const dArray2DT& DNa = fShapes->Derivatives_U();

	dArrayT Na;
	d_div = 0.0;
	for (int k = 0; k < nip; k++) 
	{
		/* shape function array */
		Na.Alias(nen, fShapes->IPShapeU(k));			

		/* get interpolation point stresses */
		fShapes->InterpolateU(nodal_dP, dP_ip, k);

		for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++) /* div */
				for (int n = 0; n < nen; n++)
					for (int m = 0; m < nen; m++)
 						d_div(i,n) += DNa(j,m)*extrap(m,k)*dP_ip(i,j)*Na[n];
 	}
}
