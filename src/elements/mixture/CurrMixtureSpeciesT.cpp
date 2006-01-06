/* $Id: CurrMixtureSpeciesT.cpp,v 1.6 2006-01-06 02:55:57 thao Exp $ */
#include "CurrMixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "Q1P0MixtureT.h"
#include "ShapeFunctionT.h"
#include "NLDiffusionMaterialT.h"
#include "MaterialListT.h"

//DEBUG
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
CurrMixtureSpeciesT::CurrMixtureSpeciesT(const ElementSupportT& support):
	NLDiffusionElementT(support),
	fGradientOption(kGlobalProjection),
	fConcentration(kCurrent),
	fOutputMass(false),
	fUpdatedLagMixture(NULL),
	fQ1P0Mixture(NULL),
	fBackgroundSpecies(NULL),
	fIndex(-1),
	fLocCurrCoords(LocalArrayT::kCurrCoords)
{
	SetName("current_mixture_species");
}

/* write element output */
void CurrMixtureSpeciesT::WriteOutput(void)
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
				fShapes->InterpolateU(fLocDisp, ip_conc);		

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
void CurrMixtureSpeciesT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLDiffusionElementT::DefineParameters(list);

	/* associated solid element group */
	list.AddParameter(ParameterT::Integer, "solid_element_group");

	/* species type */
	ParameterT species_opt(ParameterT::Enumeration, "species_type");
	species_opt.AddEnumeration("solid", kSolid);
	species_opt.AddEnumeration("fluid", kFluid);
	species_opt.AddEnumeration("solute", kSolute);
	species_opt.SetDefault(fSpecies);
	list.AddParameter(species_opt);

	/* velocity of species is calculated wrt this reference frame */
	ParameterT frame(ParameterT::Word, "background_species");
	list.AddParameter(frame, ParameterListT::ZeroOrOnce);

	/* gradient option */
	ParameterT grad_opt(ParameterT::Enumeration, "stress_gradient_option");
	grad_opt.AddEnumeration("global_projection", kGlobalProjection);
	grad_opt.AddEnumeration("element_projection", kElementProjection);
	grad_opt.SetDefault(fGradientOption);
	list.AddParameter(grad_opt);

	/* output total species mass */
	ParameterT output_mass(fOutputMass, "output_mass");
	output_mass.SetDefault(fOutputMass);
	list.AddParameter(output_mass);
}

/* accept parameter list */
void CurrMixtureSpeciesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CurrMixtureSpeciesT::TakeParameterList";
 
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

	/* species type */
	int species_opt = list.GetParameter("species_type");
	if (species_opt == kSolid)
		fSpecies = kSolid;
	else if (species_opt == kFluid)
		fSpecies = kFluid;
	else if (species_opt == kSolute)
		fSpecies = kSolute;
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"species type\" %d", species_opt);
	
	/* resolve background species */
	const ParameterT* bg_species = list.Parameter("background_species");
	if (bg_species)
	{
		StringT bg_species_name = *bg_species;
		if (bg_species_name == Field().FieldName())
			ExceptionT::GeneralFail(caller, "background_species must differ from this species \"%s\"",
				Field().FieldName().Pointer());
			
		int num_groups = ElementSupport().NumElementGroups();
		for (int i = 0; !fBackgroundSpecies && i < num_groups; i++) {
			ElementBaseT& element = ElementSupport().ElementGroup(i);
			if (element.Field().FieldName() == bg_species_name) {
				fBackgroundSpecies = TB_DYNAMIC_CAST(CurrMixtureSpeciesT*, &element);
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

    /* dimensions */
    int nip = NumIP();
    int nsd = NumSD();

    /* allocate/initialize */
    fcauchy_ip.Dimension(nip);
    fdcauchy_ip.Dimension(nip);
    for (int i = 0; i < nip; i++) {
        fcauchy_ip[i].Dimension(nsd);
        fdcauchy_ip[i].Dimension(nsd);
    }
    
	/* allocate work space for element-by-element stress projection */
	if (fGradientOption == kElementProjection) 
	{
		/* parent domain information */
		const ParentDomainT& parent_domain = fShapes->ParentDomain();
	
		/* allocate/initialize */
		fip_gradient.Dimension(nip);
		for (int i = 0; i < nip; i++) {
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
	if (fUpdatedLagMixture)
		fUpdatedLagMixture->SetConcentration(fIndex, UpdatedLagMixtureT::kCurrent);
	else
		fQ1P0Mixture->SetConcentration(fIndex, Q1P0MixtureT::kCurrent);	

	/* dimension */
	fFluxVelocity.Dimension(NumElements(), NumIP()*NumSD());
	fMassFlux.Dimension(NumElements(), NumIP()*NumSD());
	fDMassFlux.Dimension(NumElements(), NumIP()*NumSD());
	fDivBackgroundVel.Dimension(NumElements(), NumIP());
	fNEEmat.Dimension(NumElementNodes());
	fNSDmat1.Dimension(NumSD());
	fNSDmat2.Dimension(NumSD());
	fNSDmat3.Dimension(NumSD());
	
	/* initialize */
	fFluxVelocity = 0.0;
    fMassFlux = 0.0;
    fDMassFlux = 0.0;
	fDivBackgroundVel = 0.0;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* allocate and initialize shape function objects */
void CurrMixtureSpeciesT::SetShape(void)
{
	/* clear existing */
	delete fShapes;

	/*set current configuration */
	const LocalArrayT& coords = fLocCurrCoords;

	/* construct shape functions in current coordinates*/
	fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), coords);
	fShapes->Initialize();
}

/* allocate and initialize local arrays */
void CurrMixtureSpeciesT::SetLocalArrays(void)
{
	/* inherited */
	NLDiffusionElementT::SetLocalArrays();
	
	/* set up array of current nodal coordinates */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}

/* compute shape functions and derivatives */
void CurrMixtureSpeciesT::SetGlobalShape(void)
{
	/* collect current coordinates */
	SetLocalX(fLocCurrCoords);
	
	/* inherited */
	NLDiffusionElementT::SetGlobalShape();
	
	/* will need deformation gradient */
	if (fUpdatedLagMixture)	
		fUpdatedLagMixture->SetGlobalShape();
	else
		fQ1P0Mixture->SetGlobalShape();
}

/* reset loop */
void CurrMixtureSpeciesT::Top(void)
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
bool CurrMixtureSpeciesT::NextElement(void)
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
void CurrMixtureSpeciesT::RHSDriver(void)
{
    /* compute the flux velocities */
	ComputeMassFlux(false);
    
	/* inherited */
	NLDiffusionElementT::RHSDriver();
}

/* form group contribution to the stiffness matrix */
void CurrMixtureSpeciesT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	/* compute the variation in flux velocities */
	ComputeMassFlux(true);

	/* inherited */
	NLDiffusionElementT::LHSDriver(sys_type);
}

/* calculate the internal force contribution ("-k*d") */
void CurrMixtureSpeciesT::FormKd(double constK)
{
	/* dimensions */
	int nsd = NumSD();
	int nip = NumIP();

	int e = CurrElementNumber();
	/* integration parameters fShapes is in current coordinates*/
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* mass flux */
	dArray2DT m_e(nip, nsd, fMassFlux(e));
	dArrayT m;
	
	dMatrixT grad_conc;  /*gradient of concetration at integration point in current coords*/
	dArrayT conc;  /*concetration at integration point*/
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
	
		/* retrieve the mass flux */
		m_e.RowAlias(ip, m);

		/* set field gradient */
		grad_conc.Alias(1, nsd, fGradient_list[ip].Pointer());
		fShapes->GradU(fLocDisp, grad_conc);
		
		/* interpolate field */
		conc.Alias(1, fField_list.Pointer(ip));
		fShapes->InterpolateU(fLocDisp, conc);

		/* get strain-displacement matrix */
		B(ip, fB);

		/* (div) flux contribution */
		fB.MultTx(m, fNEEvec);
		
		/* c div(v) contribution */
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
		fNSDmat2.SetToCombination(by_dt, fNSDmat3, -0.5*by_dt, fNSDmat1);

		/* c div(v) */
		double c_div_v = conc[0]*fNSDmat2.Trace();
		fNEEvec.AddScaled(-c_div_v, fShapes->IPShapeU());

		if (fSpecies == kSolute) 
		{
			/*c div(v_f)*/
			dArray2DT v_e_bg;
			dArrayT v_bg(nsd);
			fNEEvec.AddScaled(-conc[0]*fDivBackgroundVel(e,ip), fShapes->IPShapeU());
			
			/*grad_x c dot v_f*/
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				v_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
			}
			v_e_bg.RowAlias(ip,v_bg);
			dArrayT vec(1);
			grad_conc.Multx(v_bg,vec);
			fNEEvec.AddScaled(-vec[0],fShapes->IPShapeU());
		}
		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}	
    if (CurrElementNumber() == 0)
        cout << "\nfRHS: "<<fRHS;
}

/* form the element stiffness matrix */
void CurrMixtureSpeciesT::FormStiffness(double constK)
{
	/* must be nonsymmetric */
	if (fLHS.Format() != ElementMatrixT::kNonSymmetric)
		ExceptionT::GeneralFail("CurrMixtureSpeciesT::FormStiffness",
			"LHS matrix must be nonsymmetric");

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* dimensions */
	int nsd = NumSD();
	int nip = NumIP();
	int nen = NumElementNodes();	

	/* flux velocity */
	dArray2DT dm_e(nip, nsd, fFluxVelocity(CurrElementNumber()));
	dArrayT dm(nsd);

	/* integrate element stiffness */
	dMatrixT grad_conc;
	dArrayT conc;
	
	/*work space*/
	dArrayT Na;
    	
 	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();

		double scale = constK*(*Det++)*(*Weight++);

		/* set field gradient */
		grad_conc.Alias(1, nsd, fGradient_list[ip].Pointer());
		fShapes->GradU(fLocDisp, grad_conc);
		
		/* interpolate field */
		conc.Alias(1, fField_list.Pointer(ip));
		fShapes->InterpolateU(fLocDisp, conc);
        
		/* strain displacement matrix */
		B(ip, fB);

		/* shape function array */
		Na.Alias(nen, fShapes->IPShapeU());			

		/*div(v)*/
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
		fNSDmat2.SetToCombination(by_dt, fNSDmat3, -0.5*by_dt, fNSDmat1);
		
		double div_v = fNSDmat2.Trace();
		fLHS.Outer(Na, Na, scale*div_v, dMatrixT::kAccumulate);
//        cout << "\nLHS1: "<<fLHS;

		if(fSpecies == kFluid)  
		{
			/*calculate Ba^T dm/dc Nb part of the stiffness matrix*/
            /*get derivative of mass flux wrt to conc*/
            dm_e.RowAlias(ip,dm);
		
			fB.MultTx(dm, fNEEvec);
            cout << "\ndm: "<<dm;
            cout << "\nNa: "<<Na;
			fLHS.Outer(fNEEvec, Na, scale, dMatrixT::kAccumulate);
        }
		else if (fSpecies == kSolute)
		{
			/*c div(v_f)*/
			dArray2DT v_e_bg;
			dArrayT v_bg(nsd);
			fLHS.Outer(Na, Na, scale*fDivBackgroundVel(CurrElementNumber(),ip), dMatrixT::kAccumulate);
			
			/*grad_x c dot v_f*/
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				v_e_bg.Alias(nip, nsd, bg_flux_velocity(CurrElementNumber()));
			}
			v_e_bg.RowAlias(ip,v_bg);
			fB.MultTx(v_bg,fNEEvec);
			fLHS.Outer(Na, fNEEvec, scale, dMatrixT::kAccumulate);
			
			/*Diffusion*/
			fD.SetToScaled(scale, fCurrMaterial->k_ij());
			fLHS.MultQTBQ(fB, fD, dMatrixT::kWhole, dMatrixT::kAccumulate);
            
            /*get derivative of mass flux wrt to conc (dD*grad_c) */
            dm_e.RowAlias(ip,dm);
            fB.MultTx(dm,fNEEvec);
			fLHS.Outer(fNEEvec, Na, scale, dMatrixT::kAccumulate);            
		}
	}
    if (CurrElementNumber() == 0)
        cout << "\nfLHS: "<<fLHS;
}

/* compute the relative mass flux and velocities. Implemented now only for fluid species  for momentum driving forces only and for solute species for diffusion only */
void CurrMixtureSpeciesT::ComputeMassFlux(bool compute_dmass_flux)
{
	const char caller[] = "CurrMixtureSpeciesT::ComputeMassFlux";
	
	/* work space */
	int nen = NumElementNodes();
	int nsd = NumSD();
	int nip = NumIP();

	dMatrixT ip_grad_x;                   /*ip gradient operator*/
	ip_grad_x.Dimension(nsd,nip);

	dArrayT ip_conc(1);
	dMatrixT grad_conc(1,nsd); /*current coords*/

	dArray2DT v_e, m_e,dm_e;
	dArrayT v, m, dm;

	if(fSpecies == kFluid) 
	{
		/*work space for projecting partial stresses to the nodes*/
        /*element local array of nodal stresses and variation of nodal stresses with concentration*/
		LocalArrayT cauchy(LocalArrayT::kUnspecified), dcauchy(LocalArrayT::kUnspecified); 
        dArray2DT dcauchy_avg;
        
        /*workspaces for calculating divergence*/
		dMatrixT ip_Grad_val, ip_Grad_val_j;
        dArrayT div_val(nsd);	
		if (fGradientOption == kGlobalProjection)
		{	
			/* get concentration specific nodal Kirchhoff stresses (tau) */
			if (fUpdatedLagMixture) fUpdatedLagMixture->ProjectPartialCauchy(fIndex);
			else ExceptionT::GeneralFail(caller, "Not implemented for Q1P0");

			fcauchy_avg = ElementSupport().OutputAverage();
			cauchy.Dimension(NumElementNodes(), nsd*nsd);		
			cauchy.SetGlobal(fcauchy_avg);
//            cout << "\nnodal cauchy: "<<cauchy;

            /* project variation in partial stresses to the nodes */
            if (compute_dmass_flux) {
                if (fUpdatedLagMixture) fUpdatedLagMixture->ProjectDPartialCauchy(fIndex);
                else ExceptionT::GeneralFail(caller, "Not implemented for Q1P0");
                dcauchy_avg.Alias(ElementSupport().OutputAverage());
                dcauchy.Dimension(NumElementNodes(), nsd*nsd);		
                dcauchy.SetGlobal(dcauchy_avg);
//                cout << "\nnodal dcauchy: "<<dcauchy;
            }
			ip_Grad_val.Dimension(nsd*nsd, nsd);
		}
        else {
            ip_grad_x.Dimension(nsd, nip);
        }
	
		/*calculate momentum driving force*/
		/* get the body force */
		dArrayT body_force(nsd);
		if (fUpdatedLagMixture) fUpdatedLagMixture->BodyForce(body_force);
		else fQ1P0Mixture->BodyForce(body_force);
        
		LocalArrayT acc(LocalArrayT::kAcc, NumElementNodes(), nsd);
		dArrayT ip_acc(nsd);
		
		dArrayT force(nsd); /*v=D*force, force = -phi*/;

		Top();
		while (NextElement()) 
		{
			int e = CurrElementNumber();
			/* global shape function values */
			SetGlobalShape();
	
			/* collect nodal stresses */
			if (fGradientOption == kGlobalProjection) {
                SetLocalU(cauchy);
                if (compute_dmass_flux) SetLocalU(dcauchy);
            }

            /* collect integration point stresses - sets shapes functions over the element */
            if (fGradientOption == kElementProjection) {
                if (fUpdatedLagMixture)
                    fUpdatedLagMixture->IP_PartialStress(fIndex, &fcauchy_ip, (compute_dmass_flux) ? &fdcauchy_ip : NULL);
                else
				ExceptionT::GeneralFail(caller, "Not implemented for Q1P0");
            }

            /* collect nodal accelerations */
            if (fUpdatedLagMixture)
                fUpdatedLagMixture->Acceleration(acc);
            else
                fQ1P0Mixture->Acceleration(acc);

			/* collect nodal concentrations */
			SetLocalU(fLocDisp);
		
			/* mass flux and velocity */
			v_e.Alias(nip, nsd, fFluxVelocity(e));
			m_e.Alias(nip, nsd, fMassFlux(e));

			/* loop over integration points */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				int ip = fShapes->CurrIP();

//                cout << "\nelement: "<< e
//                     << "\tip: "<< ip;
                     
				/* ip values of current concentration and acceleration*/
				fShapes->InterpolateU(fLocDisp, ip_conc);
				fShapes->InterpolateU(acc, ip_acc);		

				/* inertial forces */
				force.SetToCombination(-ip_conc[0], ip_acc, ip_conc[0], body_force);
 
                /* compute stress divergence */
				if (fGradientOption == kGlobalProjection)
				{
					div_val = 0.0;
					if (fUpdatedLagMixture)
						fShapes->GradU(cauchy, ip_Grad_val);   /*gradient wrt to current configuration*/
					else ExceptionT::GeneralFail(caller, "Not implemented for Q1P0");
					
					for (int j = 0; j < nsd; j++) {
						ip_Grad_val_j.Alias(nsd, nsd, ip_Grad_val(j));
						for (int i = 0; i < nsd; i++)
							div_val[i] += ip_Grad_val_j(i,j);
					}
				}
				else /* element-by-element gradient calculation */
				{
					/* transform gradient matrix to element (curr) coordinates */
					dMatrixT jacobian(nsd);
					fShapes->ParentDomain().DomainJacobian(fLocCurrCoords, ip, jacobian);
					jacobian.Inverse();
					ip_grad_x.MultATB(jacobian, fip_gradient[ip]);

					/* compute divergence of tau */
                    /* dimensions */
                    int nsd = ip_grad_x.Rows();
                    int nip = ip_grad_x.Cols();

                    div_val = 0.0;
                    for (int k = 0; k < nip; k++) {
                        const dMatrixT& ip_stress = fcauchy_ip[k];
                        for (int i = 0 ; i < nsd; i++)
                            for (int j = 0; j < nsd; j++) /* div */
                                div_val[i] += ip_grad_x(j,k)*ip_stress(i,j);			
                    }

				}
            
				/* add divergence of stress term to force */
				force += div_val;

				v_e.RowAlias(ip, v);
            
				const dMatrixT& D = fCurrMaterial->k_ij();
				D.Multx(force, v); /* c*V */

				/* compute relative mass flux */
				m_e.RowAlias(ip, m);
				m.SetToScaled(ip_conc[0],v);
  
/*                if (CurrElementNumber() == 0)
                {
                    cout << "\n ip: "<< ip
                        << "\n ip_conc: "<<ip_conc[0]
                        << "\n force: "<<force;
                }
*/
                /* mass flux variation */
                if (compute_dmass_flux)
                {
                    dm_e.Alias(nip, nsd, fDMassFlux(e));
                    dm_e.RowAlias(ip, dm);
    
                    /* contribution from changing diffusivity -cf dB/dc phi*/
                    const dMatrixT& dD = fCurrMaterial->dk_ij();
                    dD.Multx(force, dm, ip_conc[0]);

                    /*contribution form -B(2 cf dv/dt - 2 cf g - div sigma)*/
                    force.AddScaled(-ip_conc[0], ip_acc);
                    force.AddScaled(ip_conc[0], body_force);
                    const dMatrixT& D = fCurrMaterial->k_ij();
                    D.Multx(force, dm, 1.0, dMatrixT::kAccumulate);
                
                
                    /* compute divergence of dstress/dconc */
                    if (fGradientOption == kGlobalProjection)
                    {
                        div_val = 0.0;
                        if (fUpdatedLagMixture)
                            fShapes->GradU(dcauchy, ip_Grad_val);   /*gradient wrt to current configuration*/
                        else ExceptionT::GeneralFail(caller, "Not implemented for Q1P0");
					
                        for (int j = 0; j < nsd; j++) {
                            ip_Grad_val_j.Alias(nsd, nsd, ip_Grad_val(j));
                            for (int i = 0; i < nsd; i++)
                                div_val[i] += ip_Grad_val_j(i,j);
                        }
                    }
                    else /* element-by-element gradient calculation */
                    {
                        /* transform gradient matrix to element (curr) coordinates */
                        dMatrixT jacobian(nsd);
                        fShapes->ParentDomain().DomainJacobian(fLocCurrCoords, ip, jacobian);
                        jacobian.Inverse();
                        ip_grad_x.MultATB(jacobian, fip_gradient[ip]);

                        /* compute divergence of tau */
                        /* dimensions */
                        int nsd = ip_grad_x.Rows();
                        int nip = ip_grad_x.Cols();

                        div_val = 0.0;
                        for (int k = 0; k < nip; k++) {
                            const dMatrixT& ip_val= fdcauchy_ip[k];

/*                    if (CurrElementNumber() == 0)
                    {
                        cout << "\n div dcauchy ip: "<<ip_val;
                    }
                    */

                            for (int i = 0 ; i < nsd; i++)
                                for (int j = 0; j < nsd; j++) /* div */
                                    div_val[i] += ip_grad_x(j,k)*ip_val(i,j);			
                        }
                    }
                    D.Multx(div_val, dm, ip_conc[0], dMatrixT::kAccumulate); 

/*                    if (CurrElementNumber() == 0)
                    {
                        cout << "\n div dcauchy: "<<div_val;
                    } */
                
                }
            }
        }
    }
	else if (fSpecies == kSolute) 
	{		
		if (fGradientOption == kGlobalProjection){
			ProjectV();
			fv_bg_avg = ElementSupport().OutputAverage();
		}

		Top();
		while (NextElement()) 
		{
			int e = CurrElementNumber();

			/* global shape function values */
			SetGlobalShape();

			/* collect nodal concentrations */
			SetLocalU(fLocDisp);
			
			m_e.Alias(nip, nsd, fMassFlux(e));
			v_e.Alias(nip, nsd, fFluxVelocity(e));

			/* loop over integration points */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				int ip = fShapes->CurrIP();

				/* ip values of current concentration*/  //confirm that fLocDisp is current concentration*/
				fShapes->InterpolateU(fLocDisp, ip_conc);
				
				if (fGradientOption == kGlobalProjection){
					LocalArrayT v_bg(LocalArrayT::kUnspecified);
					dMatrixT ip_grad_v(nsd);

					v_bg.Dimension(NumElementNodes(), nsd);		
					v_bg.SetGlobal(fv_bg_avg);
					fShapes->GradU(v_bg, ip_grad_v);
					for (int i = 0; i < nsd; i++)
						fDivBackgroundVel(e,ip) += ip_grad_v(i,i);
				}
				else{
					dArray2DT v_e_bg;
					dArrayT v_bg(nsd);
					/* retrieve background fluid velocity */
					if (fBackgroundSpecies)
					{
						/* background flux velocity */
						const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
						v_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
					}
					else ExceptionT::GeneralFail(caller, "Background velocity not set");	/* work space */
												 
					dMatrixT jacobian(nsd);
					fShapes->ParentDomain().DomainJacobian(fLocCurrCoords, ip, jacobian);
					jacobian.Inverse();
					ip_grad_x.MultATB(jacobian, fip_gradient[ip]);
					
					for (int k = 0; k < nip; k++) {
						v_e_bg.RowAlias(k,v_bg);
						for (int i = 0 ; i < nsd; i++)
							fDivBackgroundVel(e,ip) += ip_grad_x(i,k)*v_bg[i];
					}
				}

				/*ip values spatial gradient of current concentration*/
				fShapes->GradU(fLocDisp, grad_conc);   /*gradient wrt current coordinates*/
				
				/*calculate diffusive flux and velocities*/
				const dMatrixT& D = fCurrMaterial->k_ij();
				m_e.RowAlias(ip,m);
				D.Multx(grad_conc,m, -1.0);
				
				v_e.RowAlias(ip,v);
				v.SetToScaled(1.0/ip_conc[0],m);
                
                /* mass flux variation with respect to conc: only dD grad_c implemented for now */
                if (compute_dmass_flux)
                {
                    dm_e.Alias(nip, nsd, fDMassFlux(e));
                    dm_e.RowAlias(ip, dm);
    
                    /* contribution from changing diffusivity -cf dB/dc phi*/
                    const dMatrixT& dD = fCurrMaterial->dk_ij();
                    dD.Multx(grad_conc, dm);
                }
			}
		}
	}
}


/* project the given partial first Piola-Kirchoff stress to the nodes */
void CurrMixtureSpeciesT::ProjectV(void)
{
	const char caller[] = "CurrMixtureSpeciesT::ProjectV";

	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();
	int nip = NumIP();
	
	dArray2DT nodal_v_bg(nen, nsd);
	dArray2DT v_e_bg;
	dArrayT v_bg;
	
	
	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd);
	
	Top();
	while (NextElement()) 
	{
		int e = CurrElementNumber();
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* global shape function values */
			SetGlobalShape();
			
			/* retrieve background fluid velocity */
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				v_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
			}
			else ExceptionT::GeneralFail(caller, "Background velocity not set");	/* work space */
										 
			/* extrapolate element background velocities */
			nodal_v_bg = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				int ip = fShapes->CurrIP();
				v_e_bg.RowAlias(ip, v_bg);
				
				/* extrapolate to the nodes */
				fShapes->Extrapolate(v_bg, nodal_v_bg);
			}
			
			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_v_bg);
		}
	}
}
