/* $Id: InelasticDuctile2DT.cpp,v 1.1 2003-01-24 18:46:24 paklein Exp $  */
//DEVELOPMENT
#include "InelasticDuctile2DT.h"
#include "ifstreamT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;

/* state variable information */
const int kNumState = 
    knumDOF /* \Delta_n */
  + knumDOF /* T_n */
        + 1 /* \kappa */
        + 1 /* \phi_n */
        + 1 /* \phi^*_n */
        + 1 /* \int \mathbf{T} \cdot d\boldsymbol{\Delta} */;
        
/* indicies in state variable array */
const int         k_dex_T_n = 0;
const int         k_dex_phi = 2;
const int       k_dex_phi_s = 3;
const int k_dex_dissipation = 4;
const int       k_dex_kappa = 5;
const int     k_dex_Delta_n = 6;

/* local iteration tolerances */
const double abs_tol = 1.0e-10;
const double rel_tol = 1.0e-12;
const int   max_iter = 15;

/* constructor */
InelasticDuctile2DT::InelasticDuctile2DT(ifstreamT& in, const double& time_step): 
	SurfacePotentialT(knumDOF),
	fTimeStep(time_step),
	fw_0(-1.0),
	feps_0(-1.0),
	fphi_init(-1.0),
	fdq(2), fdD(2),
	
	/* residual */
	fR(knumDOF + 1), /* the traction and phi */
	fK(knumDOF + 1),

	/* state variable data */
	fState(kNumState),
	fDelta(2, fState.Pointer(k_dex_Delta_n)),
	fTraction(2, fState.Pointer(k_dex_T_n)),
	fkappa(fState[k_dex_kappa]),
	fphi(fState[k_dex_phi]),
	fphi_s(fState[k_dex_phi_s]),
	fdissipation(fState[k_dex_dissipation])
{
	const char caller[] = "InelasticDuctile2DT::InelasticDuctile2DT";

	/* traction potential parameters */
	in >> fw_0; if (fw_0 < 0.0) ExceptionT::BadInputValue(caller);
	in >> feps_0; if (feps_0 < 0.0) ExceptionT::BadInputValue(caller);
	in >> fphi_init; if (fphi_init < 0.0) ExceptionT::BadInputValue(caller);
}

/* return the number of state variables needed by the model */
int InelasticDuctile2DT::NumStateVariables(void) const { return kNumState; }

/* surface potential */ 
double InelasticDuctile2DT::FractureEnergy(const ArrayT<double>& state) 
{
	return state[k_dex_dissipation]; 
}

double InelasticDuctile2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
#pragma unused(jump_u)
#pragma unused(state)
	return 0.0;
}
	
/* traction vector given displacement jump vector */	
const dArrayT& InelasticDuctile2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma)
{
	const char caller[] = "InelasticDuctile2DT::Traction";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
	if (fTimeStep <= 0.0) ExceptionT::BadInputValue(caller, "expecting positive time increment: %g", fTimeStep);
#endif

	/* values from t_n */
	dArrayT D_n(knumDOF, state.Pointer(k_dex_Delta_n));
	dArrayT T_n(knumDOF, state.Pointer(k_dex_T_n));
	double& phi_n = state[k_dex_phi];

	/* copy state variables */
	fState = state;

	/* compute rates */
	Rates(fState, jump_u, fTraction, fdD, fdq);
	fR[0] = (jump_u[0] - fDelta[0])/fTimeStep - fdD[0];
	fR[1] = (jump_u[1] - fDelta[1])/fTimeStep - fdD[1];
	fR[2] = (fphi - phi_n)/fTimeStep - fdq[0];

	double error, error0;
	error = error0 = fR.Magnitude();
	int count = 0;
	while (count++ < max_iter && error > abs_tol && error/error0 > rel_tol) 
	{
		/* compute Jacobian */
		Jacobian(fState, jump_u, fTraction, fdq, fK);		

//TEMP
#if 0
cout << "iteration, error = " << setw(kIntWidth) << count
                              << setw(kDoubleWidth) << error << '\n';
#endif
//TEMP

#if 0
//TEMP
cout << "iteration =  " << count << '\n';
cout << "dD = \n" << fdD << '\n';
cout << "dq = \n" << fdq << '\n';
cout << "D = \n" << jump_u << '\n';
cout << "T = \n" << fTraction << '\n';
cout << "R = \n" << fR << '\n';
cout << "K = \n" << fK << endl;
//TEMP
#endif

		/* solve update */
		fK.LinearSolve(fR);
		
		/* updates */
		fTraction[0] += fR[0];
		fTraction[1] += fR[1];
		fphi += fR[2];
//		if (fdq[0] > 0.0) fphi_s += fR[2];
		fphi_s += fR[2];
	
		/* new rates and error */
		Rates(fState, jump_u, fTraction, fdD, fdq);
		fR[0] = (jump_u[0] - fDelta[0])/fTimeStep - fdD[0];
		fR[1] = (jump_u[1] - fDelta[1])/fTimeStep - fdD[1];
		fR[2] = (fphi - phi_n)/fTimeStep - fdq[0];
		error = fR.Magnitude();
	}

	/* not converged */
	if (count == max_iter)
		ExceptionT::BadJacobianDet(caller, "not converged after %d its.", max_iter);

	/* update state variables */
	fDelta = jump_u;
	state = fState;

	return fTraction;
}

/* potential stiffness */
const dMatrixT& InelasticDuctile2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
	const char caller[] = "InelasticDuctile2DT::Stiffness";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::GeneralFail(caller);
#endif

//TEMP
ExceptionT::Stop(caller, "not implemented");

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT InelasticDuctile2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch("InelasticDuctile2DT::Status");
#endif

	if (state[k_dex_phi_s] < fphi_init)
		return Precritical;
	else if (state[k_dex_phi_s] < 0.99)
		return Critical;
	else
		return Failed;
}

void InelasticDuctile2DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    Inelastic ductile process zone 2D \n";
#endif
}

/* print parameters to the output stream */
void InelasticDuctile2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Initial width of the process zone . . . . . . . = " << fw_0 << '\n';
	out << " Rate-independent strain rate limit. . . . . . . = " << feps_0 << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int InelasticDuctile2DT::NumOutputVariables(void) const { return 3; }

void InelasticDuctile2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(3);
	labels[0] = "kappa";
	labels[1] = "phi";
	labels[2] = "phi_max";
}

/*************************************************************************
 * Protected
 *************************************************************************/

void InelasticDuctile2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) ExceptionT::GeneralFail("InelasticDuctile2DT::ComputeOutput");
#endif	
	output[0] = state[k_dex_kappa];
	output[1] = state[k_dex_phi];
	output[2] = state[k_dex_phi_s];
}

bool InelasticDuctile2DT::CompatibleOutput(const SurfacePotentialT& potential) const
{
#ifdef __NO_RTTI__
#pragma unused(potential)
	return false;
#else
	const InelasticDuctile2DT* p_inelastic = dynamic_cast<const InelasticDuctile2DT*>(&potential);
	return p_inelastic != NULL;
#endif
}

/* evaluate the rates */
void InelasticDuctile2DT::Rates(const ArrayT<double>& q, const dArrayT& D, 
	const dArrayT& T, dArrayT& dD, dArrayT& dq)
{
	/* state variables */
	double kappa = q[k_dex_kappa];
	double   phi = q[k_dex_phi];
	double phi_s = q[k_dex_phi_s];

	double T_mag = sqrt(T[0]*T[0] + T[1]*T[1]);
	double T_dot_D = D[0]*T[0] + D[1]*T[1];
	double sgn_T_dot_D = (T_dot_D < 0.0) ? -1.0 : 1.0;

	dq[0] = feps_0*sinh(sgn_T_dot_D*T_mag/(kappa*(1 - phi_s)));
	dq[1] = (dq[0] > 0.0) ? dq[0] : 0.0;
	
	dD[0] = fw_0*feps_0*sinh(T[0]/(kappa*(1 - phi_s)));
	dD[1] = fw_0*dq[0]/((1.0 - phi)*(1.0 - phi));
}

void InelasticDuctile2DT::Jacobian(const ArrayT<double>& q, const dArrayT& D, const dArrayT& T,
	const dArrayT& dq, dMatrixT& K)
{
	/* state variables */
	double kappa = q[k_dex_kappa];
	double   phi = q[k_dex_phi];
	double phi_s = q[k_dex_phi_s];

	double k1 = 1.0 - phi;
	double k2 = 1.0 - phi_s;

	double T_mag = sqrt(T[0]*T[0] + T[1]*T[1]);
	double T_dot_D = D[0]*T[0] + D[1]*T[1];
	double sgn_T_dot_D = (T_dot_D < 0.0) ? -1.0 : 1.0;

	double k3 = feps_0*cosh(sgn_T_dot_D*T_mag/(kappa*k2))/(kappa*k2);

	/* irreversible damage - phi_s is changing */
	if (true || dq[0] > 0.0) 
	{
		K(0,0) = fw_0*feps_0*cosh(T[0]/(kappa*k2))/(kappa*k2);
		K(0,1) = 0.0;
		K(0,2) = K(0,0)*T[0]/k2;

		K(2,0) = sgn_T_dot_D*T[0]*k3/T_mag; 
		K(2,1) = sgn_T_dot_D*T[1]*k3/T_mag;
		K(2,2) = k3*sgn_T_dot_D*T_mag/k2 - 1.0/fTimeStep;	
	} 
	else /* phi_s not changing */
	{
		K(0,0) = fw_0*feps_0*cosh(T[0]/(kappa*k2))/(kappa*k2);
		K(0,1) = 0.0;
		K(0,2) = 0.0;

		K(2,0) = T[0]*k3/T_mag; 
		K(2,1) = T[1]*k3/T_mag;
		K(2,2) = 0.0;
	}

	/* d_delta_dot_n */
	K(1,0) = fw_0*K(2,0)/(k1*k1);
	K(1,1) = fw_0*K(2,1)/(k1*k1);
	K(1,2) = fw_0*2.0*dq[0]/(k1*k1*k1) + fw_0*k3*sgn_T_dot_D*T_mag/(k2*k2);
}
