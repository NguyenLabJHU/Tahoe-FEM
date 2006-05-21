/* $Id: InelasticDuctile_RP2DT.cpp,v 1.8 2006-05-21 17:47:59 paklein Exp $  */
#include "InelasticDuctile_RP2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "dArrayT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* class parameters */
const int      knumDOF = 2;
const int knumInternal = 2;

/* state variable information */
const int kNumState = 
    knumDOF /* \Delta_n */
  + knumDOF /* T_n */
        + 2 /* dq = {dphi, dphi_s} */
        + 1 /* \kappa */
        + 1 /* \phi_n */
        + 1 /* \phi^*_n */
        + 1 /* \int \mathbf{T} \cdot d\boldsymbol{\Delta} */
		+ 1 /* \mathbf{T} \cdot d\boldsymbol{\Delta} */
		+ 1 /* status flag */
		+ 1 /* tied node flag (last slot) */
		+ (knumDOF + 1); /* 0.0/1.0 for active local equations (the traction and phi) */
        
/* indicies in state variable array */
const int         k_dex_T_n = 0;
const int         k_dex_phi = 2;
const int       k_dex_phi_s = 3;
const int       k_dex_kappa = 4;
const int k_dex_dissipation = 5;
const int   k_dex_incr_diss = 6;
const int     k_dex_Delta_n = 7;
const int          k_dex_dq = 9;
const int     k_status_flag = 11;
const int       k_tied_flag = 12;
const int    k_active_flags = 13; 

/* values in the BCJ output array - defined in BCJHypoIsoDamageYC3D.cpp */
//const int kBCJ_kappa_dex = 0;
//const int   kBCJ_phi_dex = 3;

/* values in the BCJ output array - defined in ABAQUS_BCJ_ISO.cpp */
const int kBCJ_kappa_dex = 0;
const int   kBCJ_phi_dex = 2;

/* code for material output */
const int kMaterialOutputCode = 6;

/* local iteration tolerances */
const int max_iter = 25;

//TEMP - enable local debugging code
#undef DEBUG

/* constructor */
InelasticDuctile_RP2DT::InelasticDuctile_RP2DT(void): 
	SurfacePotentialT(knumDOF),
	fTimeStep(NULL),
	fArea(NULL),
	fIteration(NULL),
	
	//TEMP
	//fOut(out),
	
	fw_0(-1.0),
	feps_0(-1.0),
	fphi_init(-1.0),
	fFixedIterationCount(0),
	fdD(2),
	fReversible(false),
	fUpdateIterations(1),
	
	/* work space */
	fR(knumDOF + 1), /* the traction and phi */
	fK(knumDOF + 1),

	fR_man(0, fR_temp),
	fK_man(0, fK_temp),

	/* state variable data */
	fState(kNumState),
	fDelta(2, fState.Pointer(k_dex_Delta_n)),
	fkappa(fState[k_dex_kappa]),
	fphi(fState[k_dex_phi]),
	fphi_s(fState[k_dex_phi_s]),
	fdissipation(fState[k_dex_dissipation]),
	fdq(2, fState.Pointer(k_dex_dq)),
	feq_active(knumDOF + 1, fState.Pointer(k_active_flags))
{
	const char caller[] = "InelasticDuctile_RP2DT::InelasticDuctile_RP2DT";

	/* reset memory in traction return value */
	fTraction.Alias(2, fState.Pointer(k_dex_T_n));
	
	SetName("rigid-inelastic_BCJ_2D");

#if 0
	/* traction potential parameters */
	int reverse = -1;
	in >> fw_0; if (fw_0 < 0.0) ExceptionT::BadInputValue(caller);
	in >> fphi_init; if (fphi_init < 0.0) ExceptionT::BadInputValue(caller);
	in >> fphi_max; if (fphi_init > 1.0) ExceptionT::BadInputValue(caller);
	in >> fdamage_exp;
	in >> fabs_tol;
	in >> frel_tol;

//	in >> fkappa_scale; if (fkappa_scale < 0.0) ExceptionT::BadInputValue(caller);
	in >> reverse; if (reverse != 0 && reverse != 1) ExceptionT::BadInputValue(caller);
	fReversible = bool(reverse);
	in >> fConstraintStiffness; if (fConstraintStiffness < 0.0) ExceptionT::BadInputValue(caller);
	in >> fUpdateIterations;

	/* BCJ model kinetic parameters */
	in >> fTemperature
	   >> fC1
	   >> fC2
	   >> fC3
	   >> fC4
	   >> fC5
	   >> fC6
	   >> fC19
	   >> fC20
	   >> fC21
	   >> fphi_0;

	/* BCJ model derived parameters */
	fV = fC1*exp(-fC2/fTemperature);
	if (fabs(fV) < kSmall) 
		fV = 1.0;
	
	fY = fC3/(fC21 + exp(-fC4/fTemperature));
	if (fabs(fC19) > kSmall) /* modified yield function */
		fY *= 0.5*(1.0 + tanh(fC19*(fC20 - fTemperature)));

	feps_0 = fC5*exp(-fC6/fTemperature);

	/* bulk group information for TiedPotentialBaseT */
	int num_bulk_groups = -99;
	in >> num_bulk_groups;
	if (num_bulk_groups < 0)
		ExceptionT::BadInputValue(caller, "bad number of bulk groups %d", num_bulk_groups);
	iBulkGroups.Dimension(num_bulk_groups);
	for (int i = 0; i < iBulkGroups.Length(); i++)
		in >> iBulkGroups[i];
	iBulkGroups--;
#endif
}

/* create an alias to the traction vector within the state variable array */
void InelasticDuctile_RP2DT::GetTraction(const ArrayT<double>& state, dArrayT& traction)
{
	traction.Alias(knumDOF, state.Pointer(k_dex_T_n));
}

/* make alias to the active flags */
void InelasticDuctile_RP2DT::GetActiveFlags(const ArrayT<double>& state, dArrayT& active)
{
	active.Alias(knumDOF, state.Pointer(k_active_flags));
}

/* update state variables with the given bulk data */
void InelasticDuctile_RP2DT::UpdateState(const dArrayT& bulk_nodal_data, ArrayT<double>& state) const
{
	/* don't assume new damage is OK */
	double new_damage = bulk_nodal_data[kBCJ_phi_dex];
	new_damage = (new_damage > 0.0) ? new_damage : 0.0;

	state[k_dex_phi] = new_damage;
	state[k_dex_kappa] = bulk_nodal_data[kBCJ_kappa_dex];

	if (fReversible)
		state[k_dex_phi_s] = new_damage;
	else
		state[k_dex_phi_s] = (state[k_dex_phi_s] > new_damage) ? 
			state[k_dex_phi_s] : new_damage;
}

void InelasticDuctile_RP2DT::RigidQ(const ArrayT<double>& state, ArrayT<bool>& rigid) const
{
	/* state variables */
	double kappa = state[k_dex_kappa];
	double   phi = state[k_dex_phi];
	double phi_s = state[k_dex_phi_s];

	double k1 = 1.0 - phi;
	double k2 = 1.0 - phi_s;

	/* tractions */
	double T_t = state[k_dex_T_n];
	double T_n = state[k_dex_T_n + 1];
	double T_eff = T_n + fabs(T_t);

	/* sinh flow arguments */
	double   arg_dq = (fabs(T_eff)/k2 - (kappa + fY))/fV;
	double arg_dD_t = (fabs(T_t)/k2 - (kappa + fY))/fV;

	/* evaluate constraints */
	rigid[0] = (arg_dD_t < 0.0);
	rigid[1] = (arg_dq < 0.0);
}

/* return the number of state variables needed by the model */
int InelasticDuctile_RP2DT::NumStateVariables(void) const { return kNumState; }

/* location in state variable array of the state flag */
int InelasticDuctile_RP2DT::TiedStatusPosition(void) const { return k_tied_flag; }

/* describe the parameters needed by the interface */
void InelasticDuctile_RP2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	/* model parameters */
	ParameterT w_0(fw_0, "zone_width");
	w_0.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(w_0);

	ParameterT phi_init(fphi_init, "phi_init");
	phi_init.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_init);

	ParameterT phi_max(fphi_max, "phi_max");
	phi_max.AddLimit(1.0, LimitT::UpperInclusive);
	list.AddParameter(phi_max);

	ParameterT phi_exp(fdamage_exp, "phi_exp");
	list.AddParameter(phi_exp);

	ParameterT abs_tol(fabs_tol, "abs_tol");
	list.AddParameter(abs_tol);
	ParameterT rel_tol(frel_tol, "rel_tol");
	list.AddParameter(rel_tol);
	
	ParameterT reverse_damage(fReversible, "reversible_damage");
	reverse_damage.SetDefault(fReversible);
	list.AddParameter(reverse_damage);

	ParameterT update_iter(fUpdateIterations, "update_iterations");
	update_iter.SetDefault(fUpdateIterations);
	list.AddParameter(update_iter);

	/* BCJ model kinetic parameters */
	list.AddParameter(fTemperature, "temperature");
	list.AddParameter(fC1, "C_1");
	list.AddParameter(fC2, "C_2");
	list.AddParameter(fC3, "C_3");
	list.AddParameter(fC4, "C_4");
	list.AddParameter(fC5, "C_5");
	list.AddParameter(fC6, "C_6");
	list.AddParameter(fC19, "C_19");
	list.AddParameter(fC20, "C_20");
	list.AddParameter(fC21, "C_21");
	list.AddParameter(fphi_0, "phi_0");
}

/* information about subordinate parameter lists */
void InelasticDuctile_RP2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);

	/* bulk group information for TiedPotentialBaseT */
	sub_list.AddSub("bulk_element_groups");
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* InelasticDuctile_RP2DT::NewSub(const StringT& name) const
{
	if (name == "bulk_element_groups")
		return new IntegerListT(name);
	else /* inherited */
		return SurfacePotentialT::NewSub(name);
}

/* accept parameter list */
void InelasticDuctile_RP2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);

	/* extract parameters */
	fw_0 = list.GetParameter("zone_width");
	fphi_init = list.GetParameter("phi_init");
	fphi_max = list.GetParameter("phi_max");
	fdamage_exp = list.GetParameter("phi_exp");
	fabs_tol = list.GetParameter("abs_tol");
	frel_tol = list.GetParameter("rel_tol");
	fReversible = list.GetParameter("reversible_damage");
	fUpdateIterations = list.GetParameter("update_iterations");

	/* BCJ model kinetic parameters */
	fTemperature = list.GetParameter("temperature");
	fC1 = list.GetParameter("C_1");
	fC2 = list.GetParameter("C_2");
	fC3 = list.GetParameter("C_3");
	fC4 = list.GetParameter("C_4");
	fC5 = list.GetParameter("C_5");
	fC6 = list.GetParameter("C_6");
	fC19 = list.GetParameter("C_19");
	fC20 = list.GetParameter("C_20");
	fC21 = list.GetParameter("C_21");
	fphi_0 = list.GetParameter("phi_0");

	/* BCJ model derived parameters */
	fV = fC1*exp(-fC2/fTemperature);
	if (fabs(fV) < kSmall) 
		fV = 1.0;
	
	fY = fC3/(fC21 + exp(-fC4/fTemperature));
	if (fabs(fC19) > kSmall) /* modified yield function */
		fY *= 0.5*(1.0 + tanh(fC19*(fC20 - fTemperature)));

	feps_0 = fC5*exp(-fC6/fTemperature);

	/* bulk group information for TiedPotentialBaseT */
	const ParameterListT& bulk_groups = list.GetList("bulk_element_groups");
	iBulkGroups.Dimension(bulk_groups.NumLists("Integer"));
	for (int i = 0; i < iBulkGroups.Length(); i++) {
		iBulkGroups[i] = bulk_groups.GetList("Integer", i).GetParameter("value");
	}
	iBulkGroups--;
}

/* initialize the state variable array. By default, initialization
 * involves only setting the array to zero. */
void InelasticDuctile_RP2DT::InitStateVariables(ArrayT<double>& state)
{
	state = 0.0;

	state[k_dex_phi]   = fphi_0;
	state[k_dex_phi_s] = fphi_0;
	state[k_dex_kappa] = 0.0;
	
	/* collecting initiation data */
	if (iBulkGroups.Length() > 0)
		state[k_tied_flag] = kTiedNode;
	else
		state[k_tied_flag] = kFreeNode;
}

/* incremental heat */
double InelasticDuctile_RP2DT::IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state)
{
#pragma unused(jump)
	return state[k_dex_incr_diss];
}

/* surface potential */ 
double InelasticDuctile_RP2DT::FractureEnergy(const ArrayT<double>& state) 
{
	return state[k_dex_dissipation]; 
}

double InelasticDuctile_RP2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
	return 0.0;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& InelasticDuctile_RP2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
	const char caller[] = "InelasticDuctile_RP2DT::Traction";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
	if ((*fTimeStep) < 0.0) ExceptionT::BadInputValue(caller, "expecting non-zero time increment: %g", (*fTimeStep));
#endif

	/* check pointers */
	if (!fArea) ExceptionT::GeneralFail(caller, "area pointer not set");
	if (!(*fTimeStep)) ExceptionT::GeneralFail(caller, "time step pointer not set");

// the material state should be taken from the surrounding bulk up until the point
// at which the zone initiates - handle this with the state flag ????????

	if (state[k_dex_phi_s] > fphi_max) /* failed */
		fTraction = 0.0;
	else /* calculate tractions */
	{
		/* copy state variables */
		fState = state;
		
		/* traction from the previous time increment */
		dArrayT traction_last(2, state.Pointer(k_dex_T_n));
	
		bool update_traction = true;
		if (fUpdateIterations > 1)
		{
			if (!fIteration) ExceptionT::GeneralFail(caller, "iteration pointer required for %d subiterations", fUpdateIterations);
			int mod = *fIteration % fUpdateIterations;
			update_traction = (mod == 0);
		}
	
		/* integrate traction and state variables */
		if (qIntegrate && update_traction)
		{
			dArrayT state_temp;
			state_temp.Alias(fState);
//TEMP
#ifdef DEBUG
fOut << " state_in =\n" << state_temp.wrap(5) << '\n';
fOut << "     dD = " << jump_u.no_wrap() << '\n';
fOut << "dD_last = " << fDelta.no_wrap() << '\n';
fOut << "D_t_rate = " << (jump_u[0] - fDelta[0])/(*fTimeStep) << '\n';
fOut << "D_n_rate = " << (jump_u[1] - fDelta[1])/(*fTimeStep) << '\n';
#endif

			double kappa = state[k_dex_kappa];
			double phi_s = state[k_dex_phi_s];
			double& phi_n = state[k_dex_phi];

			/* compute rates */
			Rates(fState, jump_u, fTraction, fdD, fdq, feq_active);

			bool t_free = feq_active[0] > 0.5; /* tangential opening */
			bool n_free = feq_active[1] > 0.5; /* normal opening */

			/* first two set externally, but phi evolution is not */
			feq_active[2] = feq_active[1]; 

			/* compute residual */
			fR[0] = -fdD[0]*feq_active[0];
			fR[1] = -fdD[1]*feq_active[1];
			fR[2] = -fdq[0]*feq_active[2];
			if (fabs((*fTimeStep)) > kSmall) {
				fR[0] += feq_active[0]*(jump_u[0] - fDelta[0])/(*fTimeStep);
				fR[1] += feq_active[1]*(jump_u[1] - fDelta[1])/(*fTimeStep);
				fR[2] += feq_active[2]*(fphi - phi_n)/(*fTimeStep);
			}

			double error, error0;
			error = error0 = fR.Magnitude();
			
			int count = 0;
			double phi_s_in = fphi_s;
			while (count++ < max_iter && error > fabs_tol && error/error0 > frel_tol) 
			{
//TEMP
#ifdef DEBUG
fOut << "error/error0 = " << error/error0 << " t = " << fTraction.no_wrap() << " phi = " << fphi << endl;
#endif
			
				/* compute Jacobian */
				Jacobian(fState, jump_u, fTraction, fdq, fK);		

				/* collect active parts */
				Assemble(fR, fK, feq_active);

				/* solve update */
				fK_temp.LinearSolve(fR_temp);

				/* limit T_t evolution */
				if (t_free) {
					int T_t_dex = 0;

					/* traction changing too fast - 5x flow stress */
					double max_T_t_jump = 1.0*(kappa + fY);
					if (fabs(fR_temp[T_t_dex]) > max_T_t_jump) {
						double scale = max_T_t_jump/fabs(fR_temp[T_t_dex]);
						fR_temp *= scale;
					}
				}

				/* limit T_n and phi evolution */
				if (n_free) {
					int T_n_dex = (t_free) ? 1 : 0;
					int phi_dex = T_n_dex + 1;

					/* bigger than 5% change */
					if (fabs(fR_temp[phi_dex]) > 0.05) {
						double scale = 0.05/fabs(fR_temp[phi_dex]);
						fR_temp *= scale;
					}

					/* exceeds failure max damage - go half way to max */
					double new_phi = fphi + fR_temp[phi_dex];
					if (new_phi > 1.0) {
						double scale = 0.5*(1.0 - fphi)/fR_temp[phi_dex];
						fR_temp *= scale;
					}
#if 0
					/* exceeds max damage */
					double new_phi = fphi + fR_temp[phi_dex];
					if (new_phi > fphi_max) {
						double scale = (fphi_max - fphi)/fR_temp[phi_dex];
						fR_temp *= scale;
					}
#endif

					/* less than min damage */
					else if (new_phi < 0.0) {
						double scale = 0.5*(0.0 - fphi)/fR_temp[phi_dex];
						fR_temp *= scale;
					}
#if 0
					/* less than initiation damage */
					else if (new_phi < fphi_init) {
						double scale = (fphi_init - fphi)/fR_temp[phi_dex];
						fR_temp *= scale;
					}
#endif

#if 1
					/* traction changing too fast - 1x flow stress */
					double max_T_n_jump = 1.0*(kappa + fY);
					if (fabs(fR_temp[T_n_dex]) > max_T_n_jump) {
						double scale = max_T_n_jump/fabs(fR_temp[T_n_dex]);
						fR_temp *= scale;
					}
#endif
				}
//TEMP
#if 0
int old_prec = fOut.precision();
fOut.precision(12);
fOut << "sc upd = " << fR_temp.no_wrap() << endl;
fOut.precision(old_prec);
#endif
				/* updates */
				int R_dex = 0;
				if (t_free) fTraction[0] += fR_temp[R_dex++];
				if (n_free) {
					fTraction[1] += fR_temp[R_dex++];
					fphi += fR_temp[R_dex++];
				}
				
				/* check for bad phi */
				if (fphi < 0)
					ExceptionT::GeneralFail(caller, "invalid phi: %g", fphi);

				/* irreversible damage */
				if (fReversible || fphi > phi_s_in) 
					fphi_s = fphi;
				else
					fphi_s = phi_s_in;
		
				/* new rates and error */
				Rates(fState, jump_u, fTraction, fdD, fdq, feq_active);

				/* compute residual */
				fR[0] = -fdD[0]*feq_active[0];
				fR[1] = -fdD[1]*feq_active[1];
				fR[2] = -fdq[0]*feq_active[2];
				if (fabs((*fTimeStep)) > kSmall) {
					fR[0] += feq_active[0]*(jump_u[0] - fDelta[0])/(*fTimeStep);
					fR[1] += feq_active[1]*(jump_u[1] - fDelta[1])/(*fTimeStep);
					fR[2] += feq_active[2]*(fphi - phi_n)/(*fTimeStep);
				}

				error = fR.Magnitude();
			}

			/* not converged */
			if (count >= max_iter)
				ExceptionT::BadJacobianDet(caller, "not converged after %d its: e_rel = %g > %g", max_iter, error/error0, frel_tol);

			/* incremental and integrated dissipation */
			double* T_n = state.Pointer(k_dex_T_n);
			fState[k_dex_incr_diss] = 0.5*(T_n[0] + fTraction[0])*(jump_u[0] - fDelta[0]) + 
			                          0.5*(T_n[1] + fTraction[1])*(jump_u[1] - fDelta[1]);
			fState[k_dex_dissipation] += fState[k_dex_incr_diss];

			/* update state variables */
			fDelta = jump_u;
			state = fState;

//TEMP
#ifdef DEBUG
fOut << "t_free = " << t_free << '\n';
fOut << "n_free = " << n_free << '\n';
fOut << "traction = " << fTraction.no_wrap() << '\n';
fOut << "state_out =\n" << state_temp.wrap(5) << '\n';
#endif
		}
	}	

	return fTraction;
}

/* traction vector given displacement jump vector */	
const dArrayT& InelasticDuctile_RP2DT::Traction_penalty(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
	const char caller[] = "InelasticDuctile_RP2DT::Traction_penalty";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
	if ((*fTimeStep) < 0.0) ExceptionT::BadInputValue(caller, "expecting non-zero time increment: %g", (*fTimeStep));
#endif

	/* check pointers */
	if (!fArea) ExceptionT::GeneralFail(caller, "area pointer not set");
	if (!fTimeStep) ExceptionT::GeneralFail(caller, "time step pointer not set");

// the material state should be taken from the surrounding bulk up until the point
// at which the zone initiates - handle this with the state flag ????????

	/* not free */
	if (fabs(state[k_tied_flag] - kFreeNode) > kSmall)
	{
		/* no traction */
		fTraction = 0.0;
	
		/* have nodal parameters */
		if (sigma.Length() > 0) {
			state[k_dex_phi]   = sigma[kBCJ_phi_dex];
			state[k_dex_phi_s] = sigma[kBCJ_phi_dex];
			state[k_dex_kappa] = sigma[kBCJ_kappa_dex];
		}
	}
	else /* calculate tractions */
	{
		/* copy state variables */
		fState = state;
		
		/* traction from the previous time increment */
		dArrayT traction_last(2, state.Pointer(k_dex_T_n));
	
		bool update_traction = true;
		if (fUpdateIterations > 1)
		{
			if (!fIteration) ExceptionT::GeneralFail(caller, "iteration pointer required for %d subiterations", fUpdateIterations);
			int mod = *fIteration % fUpdateIterations;
			update_traction = (mod == 0);
		}
	
		/* integrate traction and state variables */
		if (qIntegrate && update_traction)
		{
			dArrayT state_temp;
			state_temp.Alias(fState);
//TEMP
#ifdef DEBUG
fOut << " state_in = " << state_temp.wrap(5) << '\n';
fOut << "D_t_rate = " << (jump_u[0] - fDelta[0])/(*fTimeStep) << '\n';
fOut << "D_n_rate = " << (jump_u[1] - fDelta[1])/(*fTimeStep) << '\n';
#endif

			double kappa = state[k_dex_kappa];
			double phi_s = state[k_dex_phi_s];
			double& phi_n = state[k_dex_phi];

			/* supply initial guess (for T_t) */
//			if (fabs(fTraction[0]) < kSmall)
//				fTraction[0] = 1.0001*(kappa + fY)*(1.0 - phi_s); // just a little more than the resistance

			/* supply initial guess (for T_n) */
//			if (fabs(fTraction[1]) < kSmall)
//				fTraction[1] = 1.01*(kappa + fY)*(1.0 - phi_s); // just a little more than the resistance
		
			/* compute rates */
			Rates(fState, jump_u, fTraction, fdD, fdq, feq_active);

			/* check for constraints */
			bool t_free = feq_active[0] > 0.5; /* tangential opening */
			bool n_free = feq_active[1] > 0.5; /* normal opening */
			if (!t_free) fTraction[0] = traction_last[0] + fConstraintStiffness*(jump_u[0] - fDelta[0])/(*fArea);
			if (!n_free) fTraction[1] = traction_last[1] + fConstraintStiffness*(jump_u[1] - fDelta[1])/(*fArea);

			/* compute residual */
			fR[0] = -fdD[0]*feq_active[0];
			fR[1] = -fdD[1]*feq_active[1];
			fR[2] = -fdq[0]*feq_active[2];
			if (fabs((*fTimeStep)) > kSmall) {
				fR[0] += feq_active[0]*(jump_u[0] - fDelta[0])/(*fTimeStep);
				fR[1] += feq_active[1]*(jump_u[1] - fDelta[1])/(*fTimeStep);
				fR[2] += feq_active[2]*(fphi - phi_n)/(*fTimeStep);
			}

			double error, error0;
			error = error0 = fR.Magnitude();

			int count = 0;
			double phi_s_in = fphi_s;
			while ((t_free || n_free) && count++ < max_iter && error > fabs_tol && error/error0 > frel_tol) 
			{
//TEMP
#ifdef DEBUG
fOut << "error/error0 = " << error/error0 << endl;
#endif
			
				/* compute Jacobian */
				Jacobian(fState, jump_u, fTraction, fdq, fK);		

				/* collect active parts */
				Assemble(fR, fK, feq_active);

				/* solve update */
				fK_temp.LinearSolve(fR_temp);

				/* limit phi evolution */
				if (n_free) {
					int phi_dex = (t_free) ? 2 : 1;
					double new_phi = fphi + fR_temp[phi_dex];
					if (new_phi > fphi_max) {
						double scale = 0.9*(fphi_max - fphi)/fR_temp[phi_dex];
						fR_temp *= scale;
					}
					else if (new_phi < fphi_init) {
						double scale = (fphi_init - fphi)/fR_temp[phi_dex];
						fR_temp *= scale;
					}
				}

				/* updates */
				int R_dex = 0;
				if (t_free) fTraction[0] += fR_temp[R_dex++];
				if (n_free) {
					fTraction[1] += fR_temp[R_dex++];
					fphi += fR_temp[R_dex++];
				}

				/* irreversible damage */
				if (fReversible || fphi > phi_s_in) 
					fphi_s = fphi;
				else
					fphi_s = phi_s_in;
		
				/* max damage */
				if (fphi_s > fphi_max)
				{
					fphi_s = fphi_max;

					/* no traction */
					fTraction[0] = 0.0;
					fTraction[1] = 0.0;
		
					error = 0.0;
				}
				/* min damage */
				else if (fabs(fphi_init - fphi) < kSmall)
				{
					t_free = n_free = false;
					fTraction[0] = traction_last[0] + fConstraintStiffness*(jump_u[0] - fDelta[0])/(*fArea);
					fTraction[1] = traction_last[1] + fConstraintStiffness*(jump_u[1] - fDelta[1])/(*fArea);
				}
				else
				{
					/* new rates and error */
					Rates(fState, jump_u, fTraction, fdD, fdq, feq_active);

					/* check for constraints */
					t_free = feq_active[0] > 0.5; /* tangential opening */
					n_free = feq_active[1] > 0.5; /* normal opening */
					if (!t_free) fTraction[0] = traction_last[0] + fConstraintStiffness*(jump_u[0] - fDelta[0])/(*fArea);
					if (!n_free) fTraction[1] = traction_last[1] + fConstraintStiffness*(jump_u[1] - fDelta[1])/(*fArea);
					
					/* compute residual */
					fR[0] = -fdD[0]*feq_active[0];
					fR[1] = -fdD[1]*feq_active[1];
					fR[2] = -fdq[0]*feq_active[2];
					if (fabs((*fTimeStep)) > kSmall) {
						fR[0] += feq_active[0]*(jump_u[0] - fDelta[0])/(*fTimeStep);
						fR[1] += feq_active[1]*(jump_u[1] - fDelta[1])/(*fTimeStep);
						fR[2] += feq_active[2]*(fphi - phi_n)/(*fTimeStep);
					}

					error = fR.Magnitude();
				}
			}

			/* not converged */
			if (count >= max_iter)
				ExceptionT::BadJacobianDet(caller, "not converged after %d its: e_rel = %g > %g", max_iter, error/error0, frel_tol);

			/* incremental and integrated dissipation */
			double* T_n = state.Pointer(k_dex_T_n);
			fState[k_dex_incr_diss] = 0.5*(T_n[0] + fTraction[0])*(jump_u[0] - fDelta[0]) + 
			                          0.5*(T_n[1] + fTraction[1])*(jump_u[1] - fDelta[1]);
			fState[k_dex_dissipation] += fState[k_dex_incr_diss];

			/* update state variables */
			fDelta = jump_u;
			state = fState;

//TEMP
#ifdef DEBUG
fOut << "t_free = " << t_free << '\n';
fOut << "n_free = " << n_free << '\n';
fOut << "traction = " << fTraction.no_wrap() << '\n';
fOut << "state_out = " << state_temp.wrap(5) << '\n';
#endif
		}
	}	

	return fTraction;
}

/* potential stiffness */
const dMatrixT& InelasticDuctile_RP2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, 
	const dArrayT& sigma)
{
	const char caller[] = "InelasticDuctile_RP2DT::Stiffness";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
#endif

	/* not free */
	if (fabs(state[k_tied_flag] - kFreeNode) > kSmall || 
		state[k_dex_phi_s] > fphi_max || /* failed */
		fUpdateIterations > 1) /* subiterations */
	{
		/* no contribution to the stiffness */
		fStiffness = 0.0;
	}
	else /* surface is active */
	{
		/* copy state variables */
		fState = state;

		bool t_free = feq_active[0] > 0.5; /* tangential opening */
		bool n_free = feq_active[1] > 0.5; /* normal opening */

		/* something active */
		if (t_free || n_free)
		{
			/* compute Jacobian of local problem */
			Jacobian(state, jump_u, fTraction, fdq, fK);		
		
			/* all free */
			if (t_free && n_free)
			{
				/* pull values from Jacobian */
				fStiffness(0,0) = (*fTimeStep)*(fK(0,0) - (fK(0,2)*fK(2,0))/fK(2,2));
				fStiffness(0,1) = (*fTimeStep)*(fK(0,1) - (fK(0,2)*fK(2,1))/fK(2,2));
				fStiffness(1,0) = (*fTimeStep)*(fK(1,0) - (fK(1,2)*fK(2,0))/fK(2,2));
				fStiffness(1,1) = (*fTimeStep)*(fK(1,1) - (fK(1,2)*fK(2,1))/fK(2,2));
				fStiffness.Inverse();
			}
			else if (t_free) /* n constrained */
			{
				double A_tt = fK(0,0) - (fK(0,2)*fK(2,0))/fK(2,2);

				fStiffness(0,0) = 1.0/((*fTimeStep)*A_tt);
				fStiffness(0,1) = 0.0;
				fStiffness(1,0) = 0.0;		
				fStiffness(1,1) = 0.0;
			}
			else /* n_free and t constrained */
			{
				double A_nn = fK(1,1) - (fK(1,2)*fK(2,1))/fK(2,2);

				fStiffness(0,0) = 0.0;
				fStiffness(0,1) = 0.0;		
				fStiffness(1,0) = 0.0;
				fStiffness(1,1) = 1.0/((*fTimeStep)*A_nn);
			}
		}
		else /* fully constrained */
			fStiffness = 0.0;		
	}

	return fStiffness;
}

/* potential stiffness */
const dMatrixT& InelasticDuctile_RP2DT::Stiffness_penalty(const dArrayT& jump_u, const ArrayT<double>& state, 
	const dArrayT& sigma)
{
	const char caller[] = "InelasticDuctile_RP2DT::Stiffness_penalty";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
#endif

	/* not free */
	if (fabs(state[k_tied_flag] - kFreeNode) > kSmall || 
	    fabs(state[k_dex_phi_s] - fphi_max) < kSmall)
	{
		/* no contribution to the stiffness */
		fStiffness = 0.0;
	}
	else if (fUpdateIterations > 1) /* subiterations */
		fStiffness = 0.0;
	else /* surface is active */
	{
		/* copy state variables */
		fState = state;

		//TEMP
//		dArrayT state_temp;
//		state_temp.Alias(fState);
//		fOut << "state = " << state_temp.no_wrap() << '\n';

		bool t_free = feq_active[0] > 0.5; /* tangential opening */
		bool n_free = feq_active[1] > 0.5; /* normal opening */

		//TEMP
//		fOut << "t_free = " << t_free << '\n';
//		fOut << "n_free = " << n_free << '\n';

		/* something active */
		if (t_free || n_free)
		{
			/* compute Jacobian of local problem */
			Jacobian(state, jump_u, fTraction, fdq, fK);		
		
			/* all free */
			if (t_free && n_free)
			{
				/* pull values from Jacobian */
				fStiffness(0,0) = (*fTimeStep)*(fK(0,0) - (fK(0,2)*fK(2,0))/fK(2,2));
				fStiffness(0,1) = (*fTimeStep)*(fK(0,1) - (fK(0,2)*fK(2,1))/fK(2,2));
				fStiffness(1,0) = (*fTimeStep)*(fK(1,0) - (fK(1,2)*fK(2,0))/fK(2,2));
				fStiffness(1,1) = (*fTimeStep)*(fK(1,1) - (fK(1,2)*fK(2,1))/fK(2,2));
				fStiffness.Inverse();
			}
			else if (t_free) /* n constrained */
			{
				double A_tt = fK(0,0) - (fK(0,2)*fK(2,0))/fK(2,2);
				double A_tn = fK(0,1) - (fK(0,2)*fK(2,1))/fK(2,2);

				fStiffness(0,0) = 1.0/((*fTimeStep)*A_tt);
				fStiffness(0,1) =-fConstraintStiffness*A_tn/A_tt;
				fStiffness(1,0) = 0.0;		
				fStiffness(1,1) = fConstraintStiffness/(*fArea);
			}
			else /* n_free and t constrained */
			{
				double A_nt = fK(1,0) - (fK(1,2)*fK(2,0))/fK(2,2);
				double A_nn = fK(1,1) - (fK(1,2)*fK(2,1))/fK(2,2);

				fStiffness(0,0) = fConstraintStiffness/(*fArea);
				fStiffness(0,1) = 0.0;		
				fStiffness(1,0) =-fConstraintStiffness*A_nt/A_nn;		
				fStiffness(1,1) = 1.0/((*fTimeStep)*A_nn);
			}
		}
		else /* fully constrained */
		{
			fStiffness(0,0) = fConstraintStiffness/(*fArea);
			fStiffness(0,1) = 0.0;		
			fStiffness(1,0) = 0.0;		
			fStiffness(1,1) = fConstraintStiffness/(*fArea);
		}

		//TEMP
//		fOut << "K = \n" << fStiffness << endl;
	}

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT InelasticDuctile_RP2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch("InelasticDuctile_RP2DT::Status");
#endif

	if (state[k_dex_phi_s] < fphi_init)
		return Precritical;
	else if (state[k_dex_phi_s] < fphi_max)
		return Critical;
	else
		return Failed;
}

void InelasticDuctile_RP2DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    Inelastic ductile process zone 2D \n";
#endif
}

/* print parameters to the output stream */
void InelasticDuctile_RP2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Initial width of the process zone . . . . . . . = " << fw_0 << '\n';
	out << " Critical void volume fraction . . . . . . . . . = " << fphi_init << '\n';
	out << " Failure damage. . . . . . . . . . . . . . . . . = " << fphi_max << '\n';
	out << " Damage evolution exponent (n > 1) . . . . . . . = " << fdamage_exp << '\n';
	out << " Absolute convergence tolerance. . . . . . . . . = " << fabs_tol << '\n';
	out << " Relative convergence tolerance. . . . . . . . . = " << frel_tol << '\n';
	out << " Reversible damage . . . . . . . . . . . . . . . = " << ((fReversible) ? "true" : "false") << '\n';
	out << " Constraint stiffness. . . . . . . . . . . . . . = " << fConstraintStiffness << '\n';
	out << " BCJ model kinetic parameters:\n";
	out << "     temperature = " << fTemperature << '\n';
	out << "              C1 = " << fC1 << '\n';
	out << "              C2 = " << fC2 << '\n';
	out << "              C3 = " << fC3 << '\n';
	out << "              C4 = " << fC4 << '\n';
	out << "              C5 = " << fC5 << '\n';
	out << "              C6 = " << fC6 << '\n';
	out << "             C19 = " << fC19 << '\n';
	out << "             C20 = " << fC20 << '\n';
	out << "             C21 = " << fC21 << '\n';
	out << "           phi_0 = " << fphi_0 << '\n';
	out << " Rate-independent strain rate limit (f). . . . . = " << feps_0 << '\n';
	out << " Strain rate sensitivity (V) . . . . . . . . . . = " << fV << '\n';
	out << " Initial yield stress (Y). . . . . . . . . . . . = " << fY << '\n';
	out << " Number of element groups for bulk information . = " << iBulkGroups.Length() << '\n';
	iArrayT tmp;
	tmp.Alias(iBulkGroups); /* non-const alias */
	tmp++;
	out << " Group numbers : " << tmp.no_wrap() << '\n';
	tmp--;
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int InelasticDuctile_RP2DT::NumOutputVariables(void) const { return 3; }

void InelasticDuctile_RP2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(3);
	labels[0] = "kappa";
	labels[1] = "phi";
	labels[2] = "phi_max";
}

/* release condition depends on this bulk quantity */
int InelasticDuctile_RP2DT::NodalQuantityNeeded(void) const
{
	return kMaterialOutputCode; /* SolidElementT::iMaterialData to get the damage */
}

/*************************************************************************
 * Protected
 *************************************************************************/

void InelasticDuctile_RP2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) ExceptionT::GeneralFail("InelasticDuctile_RP2DT::ComputeOutput");
#endif	
	output[0] = state[k_dex_kappa];
	output[1] = state[k_dex_phi];
	output[2] = state[k_dex_phi_s];
}

bool InelasticDuctile_RP2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const InelasticDuctile_RP2DT* p_inelastic = dynamic_cast<const InelasticDuctile_RP2DT*>(&potential);
	return p_inelastic != NULL;
#endif
}

/* evaluate the rates */
void InelasticDuctile_RP2DT::Rates(const ArrayT<double>& q, const dArrayT& D, 
	const dArrayT& T, dArrayT& dD, dArrayT& dq, dArrayT& active)
{
#pragma unused(D)
#if __option(extended_errorcheck)
	const char caller[] = "InelasticDuctile_RP2DT::Rates";
	if (active.Length() != knumDOF + 1) ExceptionT::SizeMismatch(caller);
#endif

	/* state variables */
	double kappa = q[k_dex_kappa];
	double   phi = q[k_dex_phi];
	double phi_s = q[k_dex_phi_s];

	double k1 = 1.0 - phi;
	double k2 = 1.0 - phi_s;
	double k21 = pow(k2, fdamage_exp);

	/* tangential traction */
	double T_t = T[0];
	double sign_T_t = 0.0;
	if (T_t > 0.0) sign_T_t = 1.0;
	else if (T_t < 0.0) sign_T_t = -1.0;

	/* effective traction */
	double T_eff = T[1] + fabs(T_t);
	double sign_T_eff = 0.0;
	if (T_eff > 0.0) sign_T_eff = 1.0;
	else if (T_eff < 0.0) sign_T_eff = -1.0;

	/* sinh flow arguments */
	double   arg_dq = (fabs(T_eff)/k21 - (kappa + fY))/fV;
	double arg_dD_t = (fabs(T_t)/k21 - (kappa + fY))/fV;

	/* damage rates */
	dq[0] = phi*feps_0*sign_T_eff*sinh(arg_dq);
	dq[1] = (dq[0] > 0.0) ? dq[0] : 0.0;
	
	/* opening rates */
	dD[0] = phi*fw_0*feps_0*sign_T_t*sinh(arg_dD_t);
	dD[1] = fw_0*dq[0]/(k1*k1);
}

void InelasticDuctile_RP2DT::Jacobian(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
	const dArrayT& dq, dMatrixT& K)
{
#pragma unused(D)

	/* state variables */
	double kappa = q[k_dex_kappa];
	double   phi = q[k_dex_phi];
	double phi_s = q[k_dex_phi_s];

	double k1 = 1.0 - phi;
	double k2 = 1.0 - phi_s;
	double k21 = pow(k2, fdamage_exp);
	double k22 = k21*k2;

	/* tangential traction */
	double T_t = T[0];
	double sign_T_t = 0.0;
	if (T_t > 0.0) sign_T_t = 1.0;
	else if (T_t < 0.0) sign_T_t = -1.0;

	/* effective traction */
	double T_eff = T[1] + fabs(T_t);
	double sign_T_eff = 0.0;
	if (T_eff > 0.0) sign_T_eff = 1.0;
	else if (T_eff < 0.0) sign_T_eff = -1.0;

	/* flow arguments */
	double   arg_dq = (fabs(T_eff)/k21 - (kappa + fY))/fV;
	double arg_dD_t = (fabs(T_t)/k21 - (kappa + fY))/fV;

	double k3 = feps_0*cosh(arg_dq);
	double k4 = fw_0*feps_0*cosh(arg_dD_t);

	double k5 = feps_0*sinh(arg_dq);
	double k6 = fw_0*feps_0*sinh(arg_dD_t);

	/* irreversible damage - phi_s is changing */
	if (phi_s < fphi_max && (fReversible || dq[0] > 0.0))
	{
		K(0,0) = phi*sign_T_t*sign_T_t*k4/(fV*k21);
		K(0,1) = 0.0;
		K(0,2) = sign_T_t*k6 + fdamage_exp*phi*sign_T_t*fabs(T_t)*k4/(fV*k22);

		K(2,0) = phi*sign_T_t*sign_T_eff*sign_T_eff*k3/(fV*k21);
		K(2,1) = phi*sign_T_eff*sign_T_eff*k3/(fV*k21);
		K(2,2) = sign_T_eff*k5 + fdamage_exp*phi*sign_T_eff*fabs(T_eff)*k3/(fV*k22);
	} 
	else /* phi_s not changing */
	{
		K(0,0) = phi*sign_T_t*sign_T_t*k4/(fV*k21);
		K(0,1) = 0.0;
		K(0,2) = sign_T_t*k6;

		K(2,0) = phi*sign_T_t*sign_T_eff*sign_T_eff*k3/(fV*k21);
		K(2,1) = phi*sign_T_eff*sign_T_eff*k3/(fV*k21);
		K(2,2) = sign_T_eff*k5;
	}

	/* d_delta_dot_n */
	K(1,0) = fw_0*K(2,0)/(k1*k1);
	K(1,1) = fw_0*K(2,1)/(k1*k1);
	K(1,2) = fw_0*(K(2,2) + 2.0*dq[0]/k1)/(k1*k1);

	/* because: R_q = q_dot - 1/dt (q_n+1 - q_n) */ 
	if (fabs((*fTimeStep)) > kSmall) K(2,2) -= 1.0/(*fTimeStep);
}

/* reduce local the active equations */
void InelasticDuctile_RP2DT::Assemble(const dArrayT& R, const dMatrixT& K, const dArrayT& active)
{
	/* dimension check */
	if (active.Length() != 3) ExceptionT::GeneralFail("InelasticDuctile_RP2DT::Assemble");

	/* count active equations */
	int dim = 0;
	if (active[0] > 0.5) dim++;
	if (active[1] > 0.5) dim++;
	if (active[2] > 0.5) dim++;

	/* dimension working space */
	fR_man.SetLength(dim, false);
	fK_man.SetDimensions(dim);

	int row_dex = 0;
	for (int row = 0; row < 3; row++)
		if (active[row] > 0.5)
		{	
			fR_temp[row_dex] = R[row];
	
			int col_dex = 0;
			for (int col = 0; col < 3; col++)
				if (active[col] > 0.5)
				{
					fK_temp(row_dex, col_dex) = K(row, col);
				
					/* next col */
					col_dex++;				
				}
	
			/* next row */
			row_dex++;
		}
}
