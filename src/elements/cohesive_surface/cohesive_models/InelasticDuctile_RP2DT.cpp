/* $Id: InelasticDuctile_RP2DT.cpp,v 1.1 2003-08-08 16:03:18 paklein Exp $  */
#include "InelasticDuctile_RP2DT.h"
#include "ifstreamT.h"
#include "dArrayT.h"

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
		+ 1 /* tied node flag (last slot) */;
        
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

/* values in the BCJ output array - defined in BCJHypoIsoDamageYC3D.cpp */
const int kBCJ_kappa_dex = 0;
const int   kBCJ_phi_dex = 3;

/* code for material output */
const int kMaterialOutputCode = 6;

/* local iteration tolerances */
const double abs_tol = 1.0e-10;
const double rel_tol = 1.0e-12;
const int   max_iter = 25;
const double phi_max = 0.999;

/* constructor */
InelasticDuctile_RP2DT::InelasticDuctile_RP2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step),
	fw_0(-1.0),
	feps_0(-1.0),
	fphi_init(-1.0),
	fdD(2),
	
	/* work space */
	fR(knumDOF + 1), /* the traction and phi */
	fK(knumDOF + 1),

	fdD_active(knumDOF),
	fdq_active(knumInternal),

	/* state variable data */
	fState(kNumState),
	fDelta(2, fState.Pointer(k_dex_Delta_n)),
	fkappa(fState[k_dex_kappa]),
	fphi(fState[k_dex_phi]),
	fphi_s(fState[k_dex_phi_s]),
	fdissipation(fState[k_dex_dissipation]),
	fdq(2, fState.Pointer(k_dex_dq))
{
	const char caller[] = "InelasticDuctile_RP2DT::InelasticDuctile_RP2DT";

	/* reset memory in traction return value */
	fTraction.Set(2, fState.Pointer(k_dex_T_n));

	/* traction potential parameters */
	int reverse = -1;
	in >> fw_0; if (fw_0 < 0.0) ExceptionT::BadInputValue(caller);
	in >> fphi_init; if (fphi_init < 0.0) ExceptionT::BadInputValue(caller);
	in >> fkappa_scale; if (fkappa_scale < 0.0) ExceptionT::BadInputValue(caller);
	in >> reverse; if (reverse != 0 && reverse != 1) ExceptionT::BadInputValue(caller);
	fReversible = bool(reverse);
	in >> fConstraintStiffness; if (fConstraintStiffness < 0.0) ExceptionT::BadInputValue(caller);

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
	   >> fC21;

	/* BCJ model derived parameters */
	fV = fC1*exp(-fC2/fTemperature);
	fY = (fC3/(fC21 + exp(-fC4/fTemperature)))*0.5*(1.0 + tanh(fC19*(fC20 - fTemperature)));
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
}

/* return the number of state variables needed by the model */
int InelasticDuctile_RP2DT::NumStateVariables(void) const { return kNumState; }

/* location in state variable array of the state flag */
int InelasticDuctile_RP2DT::TiedStatusPosition(void) const { return k_tied_flag; }

/* initialize the state variable array. By default, initialization
 * involves only setting the array to zero. */
void InelasticDuctile_RP2DT::InitStateVariables(ArrayT<double>& state)
{
	state = 0.0;

#pragma message("InelasticDuctile_RP2DT::InitStateVariables: using temporary values")
	state[k_dex_phi]   = 0.05;
	state[k_dex_phi_s] = 0.05;
	state[k_dex_kappa] = 1.0;
	
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
	if (fTimeStep < 0.0) ExceptionT::BadInputValue(caller, "expecting non-zero time increment: %g", fTimeStep);
#endif

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
			state[k_dex_kappa] = fkappa_scale*sigma[kBCJ_kappa_dex];
		}
	}
	else /* calculate tractions */
	{
		/* copy state variables */
		fState = state;
	
		/* integrate traction and state variables */
		if (qIntegrate)
		{
			double kappa = state[k_dex_kappa];
			double phi_s = state[k_dex_phi_s];

			/* supply initial guess (for T_t) */
			if (fabs(fTraction[0]) < kSmall)
				fTraction[0] = 1.0001*(kappa + fY)*(1.0 - phi_s); // just a little more than the resistance

			/* supply initial guess (for T_n) */
			if (fabs(fTraction[1]) < kSmall)
				fTraction[1] = 1.01*(kappa + fY)*(1.0 - phi_s); // just a little more than the resistance
		
			/* compute rates */
			double& phi_n = state[k_dex_phi];
			Rates(fState, jump_u, fTraction, fdD, fdq, fdD_active, fdq_active);
			
			//TEMP?
			//LimitRates(jump_u, fdD, fState[k_dex_phi], fdq);
			
			fR[0] = -fdD[0];
			fR[1] = -fdD[1];
			fR[2] = -fdq[0];
			if (fabs(fTimeStep) > kSmall) {
				fR[0] += (jump_u[0] - fDelta[0])/fTimeStep;
				fR[1] += (jump_u[1] - fDelta[1])/fTimeStep;
				fR[2] += (fphi - phi_n)/fTimeStep;
			}

			double error, error0;
			error = error0 = fR.Magnitude();

			int count = 0;
			double phi_s_in = fphi_s;
			while (count++ < max_iter && error > abs_tol && error/error0 > rel_tol) 
			{
				/* compute Jacobian */
				Jacobian(fState, jump_u, fTraction, fdq, fK);		

				/* solve update */
				fK.LinearSolve(fR);

				//TEMP? - scale update vector to limit phi
				if (fR[2] > 0.05)
				{
					fR *= 0.05/fR[2];
				}
				else if (fR[2] < -0.05)
				{
					fR *= -0.05/fR[2];
				}

				/* updates */
				fTraction[0] += fR[0];
				fTraction[1] += fR[1];				
				fphi += fR[2];
				if (fReversible || fphi > phi_s_in) 
					fphi_s = fphi;
				else
					fphi_s = phi_s_in;
		
				/* max damage */
				if (fphi_s > phi_max)
				{
					fphi_s = phi_max;

					/* no traction */
					fTraction[0] = 0.0;
					fTraction[1] = 0.0;
		
					error = 0.0;
				}
				else
				{
					/* new rates and error */
					Rates(fState, jump_u, fTraction, fdD, fdq, fdD_active, fdq_active);
					
					//TEMP?
					//LimitRates(jump_u, fdD, fState[k_dex_phi], fdq);
					
					fR[0] = -fdD[0];
					fR[1] = -fdD[1];
					fR[2] = -fdq[0];
					if (fabs(fTimeStep) > kSmall) {
						fR[0] += (jump_u[0] - fDelta[0])/fTimeStep;
						fR[1] += (jump_u[1] - fDelta[1])/fTimeStep;
						fR[2] += (fphi - phi_n)/fTimeStep;
					}

					error = fR.Magnitude();
				}
			}

			/* not converged */
			if (count >= max_iter)
				ExceptionT::BadJacobianDet(caller, "not converged after %d its.", max_iter);

			/* incremental and integrated dissipation */
			double* T_n = state.Pointer(k_dex_T_n);
			fState[k_dex_incr_diss] = 0.5*(T_n[0] + fTraction[0])*(jump_u[0] - fDelta[0]) + 
			                          0.5*(T_n[1] + fTraction[1])*(jump_u[1] - fDelta[1]);
			fState[k_dex_dissipation] += fState[k_dex_incr_diss];

			/* update state variables */
			fDelta = jump_u;
			state = fState;
		}
	}	

#if 0
dArrayT tmp;
tmp.Alias(state);
cout << "F:   T = " << fTraction.no_wrap() << '\n';
cout << "F: fdq = " << fdq.no_wrap() << '\n';
cout << "F: jmp = " << jump_u.no_wrap() << '\n';
cout << "F:  st = \n" << tmp << '\n';
#endif

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
	if (fabs(state[k_tied_flag] - kFreeNode) > kSmall)
	{
		/* no contribution to the stiffness */
		fStiffness = 0.0;
	}
	else /* surface is active */
	{
		/* copy state variables */
		fState = state;

		/* compute Jacobian of local problem */
		Jacobian(state, jump_u, fTraction, fdq, fK);		

		/* pull values from Jacobian */
		fStiffness(0,0) = fTimeStep*(fK(0,0) - (fK(0,2)*fK(2,0))/fK(2,2));
		fStiffness(0,1) = fTimeStep*(fK(0,1) - (fK(0,2)*fK(2,1))/fK(2,2));
		fStiffness(1,0) = fTimeStep*(fK(1,0) - (fK(1,2)*fK(2,0))/fK(2,2));
		fStiffness(1,1) = fTimeStep*(fK(1,1) - (fK(1,2)*fK(2,1))/fK(2,2));
		fStiffness.Inverse();
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
	else if (state[k_dex_phi_s] < phi_max)
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
	out << " Reversible damage . . . . . . . . . . . . . . . = " << ((fReversible) ? "true" : "false") << '\n';
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
	const dArrayT& T, dArrayT& dD, dArrayT& dq, ArrayT<bool>& dD_active, ArrayT<bool>& dq_active)
{
#pragma unused(D)
#if __option(extended_errorcheck)
	const char caller[] = "InelasticDuctile_RP2DT::Rates";
	if (dD_active.Length() != 2) ExceptionT::SizeMismatch(caller);
	if (dq_active.Length() != 1) ExceptionT::SizeMismatch(caller);
#endif

	/* state variables */
	double kappa = q[k_dex_kappa];
	double   phi = q[k_dex_phi];
	double phi_s = q[k_dex_phi_s];

	double k1 = 1.0 - phi;
	double k2 = 1.0 - phi_s;

	double kY = (kappa + fY);
	double kV = fV;

	/* tangential traction */
	double T_t = T[0];
	double sign_T_t = 0.0;
	if (T_t > 0.0) sign_T_t = 1.0;
	else if (T_t < 0.0) sign_T_t = -1.0;

	/* effective traction */
	double T_eff = T[1] + fabs(T[0]);
	double sign_T_eff = 0.0;
	if (T_eff > 0.0) sign_T_eff = 1.0;
	else if (T_eff < 0.0) sign_T_eff = -1.0;

	/* sinh flow arguments */
	double   arg_dq = (fabs(T_eff)/k2 - (kappa + fY))/fV;
	double arg_dD_t = (fabs(T_t)/k2 - (kappa + fY))/fV;

	/* damage rates */
	dq[0] = feps_0*sign_T_eff*sinh(arg_dq);
	dq[1] = (dq[0] > 0.0) ? dq[0] : 0.0;
	
	/* opening rates */
	dD[0] = fw_0*feps_0*sign_T_t*sinh(arg_dD_t);
	dD[1] = fw_0*dq[0]/(k1*k1);
	
	/* what's active */
	dq_active[0] = (arg_dq > 0.0) ? true : false;
	dD_active[0] = (arg_dD_t > 0.0) ? true : false;
	dD_active[1] = dq_active[0];
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

	double kY = (kappa + fY)*(1.0 - phi_s);
	double kV = fV*(1.0 - phi_s);

	double sign_T_t = 0.0;
	if (T[0] > 0.0) sign_T_t = 1.0;
	else if (T[0] < 0.0) sign_T_t = -1.0;

	double k3 = feps_0*cosh((fabs(T[0]) + T[1] - kY)/kV);
	double k4 = fw_0*feps_0*cosh(sign_T_t*(fabs(T[0]) - kY)/kV);

	/* irreversible damage - phi_s is changing */
	if (phi_s < phi_max && (fReversible || dq[0] > 0.0))
	{
		K(0,0) = k4/kV;
		K(0,1) = 0.0;
		K(0,2) = sign_T_t*k4*((kappa + kY)/kV + (fabs(T[0]) - kY)/(kV*k2));

		K(2,0) = sign_T_t*k3/kV;
		K(2,1) = k3/kV;
		K(2,2) = k3*((kappa + kY)/kV + (fabs(T[0]) + T[1] - kY)/(kV*k2));
	} 
	else /* phi_s not changing */
	{
		K(0,0) = fw_0*feps_0*cosh(T[0]/(kappa*k2))/(kappa*k2);
		K(0,1) = 0.0;
		K(0,2) = 0.0;

		K(2,0) = sign_T_t*k3/kV;
		K(2,1) = k3/kV;
		K(2,2) = 0.0;
	}

	/* d_delta_dot_n */
	K(1,0) = fw_0*K(2,0)/(k1*k1);
	K(1,1) = fw_0*K(2,1)/(k1*k1);
	K(1,2) = fw_0*(K(2,2) + 2.0*dq[0]/k1)/(k1*k1);

	/* because: R_q = q_dot - 1/dt (q_n+1 - q_n) */ 
	if (fabs(fTimeStep) > kSmall) K(2,2) -= 1.0/fTimeStep;
}
