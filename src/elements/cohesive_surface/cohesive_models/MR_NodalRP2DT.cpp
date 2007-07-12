/* $Id: MR_NodalRP2DT.cpp,v 1.16 2007-07-12 21:26:27 skyu Exp $  */
#include "MR_NodalRP2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "dArrayT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* class parameters */
const int   knumDOF = 2;

/* state variable information */
const int kNumState =
	  knumDOF 	/* \bT */
	+ knumDOF 	/* \bDelta */
	+ 1 		/* \epsilon^p_n */
	+ 1 		/* \epsilon^p_s */
	+ 1 		/* \chi */
	+ 1 		/* c */
	+ 1 		/* \phi */
	+ 1 		/* \psi */
	+ 1 		/* yield */
	+ 1 		/* norm of residual */
	+ knumDOF 	/* \bu^p */
	+ 1 		/* plastic multiplier */
	+ 1 		/* plasticity flag */
	+ 1 		/* iteration number */
	+ 1 		/* quantity A? */
	+ 1 		/* quantity B? */
	+ 1 		/* status flag */
	+ 1 		/* tied node flag (last slot) */
	+ (knumDOF+2);	/* 0.0/1.0 for active local equations (\bT, up_t and up_n) */

/* indices in state variable array */
const int	         k_T_t = 0;
const int	         k_T_n = 1;
const int	     k_Delta_t = 2;
const int	     k_Delta_n = 3;
const int	         k_esp = 4;
const int	         k_enp = 5;
const int	        k_fchi = 6;
const int	          k_fc = 7;
const int	        k_fphi = 8;
const int	        k_fpsi = 9;
const int	       k_yield = 10;
const int	  k_norm_resid = 11;
const int	        k_up_t = 12;
const int	        k_up_n = 13;
const int	  k_plast_mult = 14;
const int	     k_plastic = 15;
const int	   k_iteration = 16;
const int	           k_A = 17;
const int	           k_B = 18;
const int	 k_status_flag = 19;
const int	   k_tied_flag = 20;
const int	k_active_flags = 21;

/* code for material output */
const int kMaterialOutputCode = 6;

/* local iteration tolerances */
const int max_iter = 25;

//TEMP - enable local debugging code
//#undef DEBUG

/* constructor */
MR_NodalRP2DT::MR_NodalRP2DT(void): 
	SurfacePotentialT(knumDOF),
	fTimeStep(NULL),
	fArea(NULL),
	fIteration(NULL),
	
	fFixedIterationCount(0),
	fUpdateIterations(1),

	/* state variable data */
	fState(kNumState),
	feq_active(knumDOF + 2, fState.Pointer(k_active_flags))
{
	const char caller[] = "MR_NodalRP2DT::MR_NodalRP2DT";

	/* reset memory in traction return value */
	fTraction.Alias(2, fState.Pointer(k_T_t));
	
	SetName("nodal-rigid-plastic_MR_RP2D");
}

/* create an alias to the traction vector within the state variable array */
void MR_NodalRP2DT::GetTraction(const ArrayT<double>& state, dArrayT& traction)
{
	traction.Alias(knumDOF, state.Pointer(k_T_t));
}

/* make alias to the active flags */
void MR_NodalRP2DT::GetActiveFlags(const ArrayT<double>& state, dArrayT& active)
{
	active.Alias(knumDOF, state.Pointer(k_active_flags));
}

void MR_NodalRP2DT::RigidQ(const ArrayT<double>& state, ArrayT<bool>& rigid) const
{
	// not sure this will work, RAR
	
	/* check for yielding */
	double yield = state[k_yield];

	/* evaluate constraints */
	rigid[0] = (yield < 0.0);
	rigid[1] = (yield < 0.0);
}

/* return the number of state variables needed by the model */
int MR_NodalRP2DT::NumStateVariables(void) const { return kNumState; }

/* location in state variable array of the state flag */
int MR_NodalRP2DT::TiedStatusPosition(void) const { return k_tied_flag; }

/* describe the parameters needed by the interface */
void MR_NodalRP2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	/* model parameters */
	ParameterT Gf_I(fGf_I, "Gf_I");
	Gf_I.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Gf_I);
	
	ParameterT Gf_II(fGf_II, "Gf_II");
	Gf_II.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Gf_II);
	
	ParameterT chi_p(fchi_p, "chi_p");
	chi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi_p);
	
	ParameterT chi_r(fchi_r, "chi_r");
	chi_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi_r);
	
	ParameterT c_p(fc_p, "c_p");
	c_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c_p);
	
	ParameterT c_r(fc_r, "c_r");
	c_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c_r);
	
	ParameterT phi_p(fphi_p, "phi_p");
	phi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_p);
	
	ParameterT phi_r(fphi_r, "phi_r");
	phi_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_r);
	
	ParameterT psi_p(fpsi_p, "psi_p");
	psi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(psi_p);

	ParameterT alpha_chi(falpha_chi, "alpha_chi");
	alpha_chi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_chi);
	
	ParameterT alpha_c(falpha_c, "alpha_c");
	alpha_c.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_c);
	
	ParameterT alpha_phi(falpha_phi, "alpha_phi");
	alpha_phi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_phi);
	
	ParameterT alpha_psi(falpha_psi, "alpha_psi");
	alpha_psi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_psi);
	
	ParameterT Tol_1(fTol_1, "Tol_1");
	Tol_1.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Tol_1);
	
	ParameterT Tol_2(fTol_2, "Tol_2");
	Tol_2.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Tol_2);
}

/* information about subordinate parameter lists */
void MR_NodalRP2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* MR_NodalRP2DT::NewSub(const StringT& name) const
{
	return SurfacePotentialT::NewSub(name);
}

/* accept parameter list */
void MR_NodalRP2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);
	
	/* extract parameters */
	fGf_I = list.GetParameter("Gf_I");
	fGf_II = list.GetParameter("Gf_II");
	fchi_p = list.GetParameter("chi_p");
	fchi_r = list.GetParameter("chi_r");
	fc_p = list.GetParameter("c_p");
	fc_r = list.GetParameter("c_r");
	fphi_p = list.GetParameter("phi_p");
	fphi_r = list.GetParameter("phi_r");
	fpsi_p = list.GetParameter("psi_p");
	falpha_chi = list.GetParameter("alpha_chi");
	falpha_c = list.GetParameter("alpha_c");
	falpha_phi = list.GetParameter("alpha_phi");
	falpha_psi = list.GetParameter("alpha_psi");
	fTol_1 = list.GetParameter("Tol_1");
	fTol_2 = list.GetParameter("Tol_2");
	
	/* calculate limits on fc_p and fc_r */
	double const1 = tan(fphi_p)*fchi_p;
	if (fc_p < const1) fc_p = const1;
	const1 = tan(fphi_r)*fchi_r;
	if (fc_r < const1) fc_r = const1;
	
	/* setup output file and format */
	outputPrecision = 10;
	outputFileWidth = outputPrecision + 8;
	mr_rp_2d_out.open("mr_rp_2d.info");
}

/* initialize the state variable array. By default, initialization
 * involves only setting the array to zero. */
void MR_NodalRP2DT::InitStateVariables(ArrayT<double>& state)
{
	int num_state = NumStateVariables();
	if (state.Length() != num_state) {
		cout << "\n MR_NodalRP2DT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
		throw ExceptionT::kSizeMismatch;	
	}

	/* clear */
	if (num_state > 0) state = 0.0;
	
	/* Initializing internal state variables */
	double enp  = state[k_enp];
	double esp  = state[k_esp];
	double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
	double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
	double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
	double ftan_psi = (tan(fpsi_p))*exp(-falpha_psi*esp);	
	state[k_fchi]  = fchi;
	state[k_fc]  = fc ;
	state[k_fphi]  = ftan_phi;
	state[k_fpsi]  = ftan_psi;
	state[k_norm_resid] = 0.;
	state[k_tied_flag] = kTiedNode;
	//state[k_tied_flag] = kFreeNode;
}

double MR_NodalRP2DT::signof(double r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

/* Value of the Yield Function */ 
double MR_NodalRP2DT::YFValue(const ArrayT<double>& state)
{
	return state[k_yield];
}

/** dissipated energy. Total amount of energy dissipated reaching
	 * the current state. */
double MR_NodalRP2DT::FractureEnergy(const ArrayT<double>& state)
{
	double FE = state[k_T_n]*state[k_up_n] + state[k_T_t]*state[k_up_t];
	return FE;
}

/** potential energy */
double MR_NodalRP2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
	double PT = (jump_u[0]*state[k_T_t] + jump_u[1]*state[k_T_n]);
	return PT;
}

	
/* traction vector given displacement jump vector */	
const dArrayT& MR_NodalRP2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
	const char caller[] = "MR_NodalRP2DT::Traction";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
	if ((*fTimeStep) < 0.0) ExceptionT::BadInputValue(caller, "expecting non-zero time increment: %g", (*fTimeStep));
#endif

	/* check pointers */
	if (!fArea) ExceptionT::GeneralFail(caller, "area pointer not set");
	if (!(*fTimeStep)) ExceptionT::GeneralFail(caller, "time step pointer not set");

	/*
	if (state[k_tied_flag] == kTiedNode)
	{
		fTraction = 0.;
		return fTraction;
	}
	else
	{
	//nodes are either untied or will be next step
	}
	*/

	/*
	//this condition may be obsolete, RAR
	if (state[k_tied_flag] == kReleaseNextStep)
	{	// nodes will be untied next step. Store stress and move along
		// qIntegrate is guaranteed to be true
		if (!qIntegrate)
			ExceptionT::GeneralFail("MR_NodalRP2DT::Traction","nodes freed and !qIntegrate");
		state[k_T_n] = sigma[2];
		state[k_T_t] = sigma[1];
		state[k_tied_flag] = kFirstFreeStep;
		fTraction = 0.;
			
		return fTraction;
	} 

	//this condition may be obsolete, RAR
	if (state[k_tied_flag] == kFirstFreeStep) // First timestep with free nodes
	{
		if (qIntegrate) state[k_tied_flag] = kFreeNode;
	}
	*/

	/*
	if (fabs(jump_u[0]) < kSmall && fabs(jump_u[1]) < kSmall)
	{
		fTraction = 0.;
		return fTraction;
	}
	*/
	
	bool update_traction = true;

	/* value of Parameter fUpdateIterations doesn't be obtained from the input
	if (fUpdateIterations > 1)
	{
		if (!fIteration) ExceptionT::GeneralFail(caller, "iteration pointer required for %d subiterations", fUpdateIterations);
		int mod = *fIteration % fUpdateIterations;
		update_traction = (mod == 0);
	}
	*/

	/* integrate traction and state variables */
	if (qIntegrate && update_traction)
	{
		int i; int j; int kk; int iplastic;

		dMatrixT AA(6,6); dMatrixT KE(2,2); dMatrixT KE_Inv(2,2);
		dMatrixT I_mat(4,4); dMatrixT CMAT(6,6); dMatrixT A_qq(4,4);
		dMatrixT A_uu(2,2); dMatrixT A_uq(2,4); dMatrixT A_qu(4,2);
		dMatrixT ZMAT(2,4); dMatrixT ZMATP(4,2); dMatrixT dQdSig2(2,2);
		dMatrixT dqbardq(4,4); dMatrixT dQdSigdq(2,4);
		dMatrixT dqbardSig(4,2); dMatrixT AA_inv(6,6);
		dMatrixT X(6,1); dMatrixT Y(6,1);

		dArrayT up(2); dArrayT dup(2); dArrayT dSig(2); dArrayT qn(4);
		dArrayT qo(4); dArrayT Rvec(6); dArrayT Cvec(6); dArrayT upo(2);
		dArrayT R(6); dArrayT Rmod(6); dArrayT Sig(2); dArrayT Sig_I(2);
		dArrayT dQdSig(2); dArrayT dfdq(4); dArrayT qbar(4);
		dArrayT R2(6); dArrayT V_sig(2); dArrayT V_q(4);
		dArrayT dfdSig(2); dArrayT dq(4);

		double ff; double bott; double topp; double dlam; double dlam2;
		double normr; double normflow; double normdup;

		/* initialize the necessary vectors */
		I_mat = 0.;
		ZMAT = 0.; ZMATP = 0.; dlam = 0.; dlam2 = 0.; normr = 0.;
		    
		up[0] = jump_u[0];
		up[1] = jump_u[1];
		dup[0] = up[0] - state[k_up_t];
		dup[1] = up[1] - state[k_up_n];
		upo[0] = state[k_up_t];
		upo[1] = state[k_up_n];
		Sig_I[0] = state[k_T_t];
		Sig_I[1] = state[k_T_n];

		qn[0] = state[k_fchi];
		qn[1] = state[k_fc];
		qn[2] = state[k_fphi];
		qn[3] = state[k_fpsi];

		for (i = 0; i<=3; ++i)
		{
			qo[i] = qn[i];
			I_mat(i,i) = 1.;
		}
	    
		Sig = Sig_I;
		dQdSig_f(Sig, qn, dQdSig);
		qbar_f(Sig, qn, qbar);
	    
		/* first estimate of plastic consistency parameter */
		normflow = dQdSig.Magnitude();
		normdup  = dup.Magnitude();
		dlam     = normdup;
		dlam    /= normflow;
	    
		if (normdup <=kSmall)
		{
			fTraction[0] = state[k_T_t];
			fTraction[1] = state[k_T_n];

			return fTraction;
		}
	    
		/* calculate residuals */    
		for (i = 0; i<=1; ++i) 
		{
			R[i]  = upo[i];
			R[i] -= up[i];
			R[i] += dlam*dQdSig[i];
		}

		for (i = 0; i<=3; ++i) 
		{
			R[i+2]  = qo[i];
			R[i+2] -= qn[i];
			R[i+2] += dlam*qbar[i];
		}
		normr = R.Magnitude();
	    
		/* Local Iteration */
		kk = 0;
		// iplastic = 1;
		Yield_f(Sig, qn, ff); //check for yield and correct the tied flag

		if (ff >= 0.0)
		{
			state[k_tied_flag] = kFreeNode;
			iplastic = 1;
		}
		else
		{
			state[k_tied_flag] = kTiedNode;
			iplastic = 0;
		}

		//begin iteration loop
		while (ff > fTol_1 || ff < 0.0 || normr > fTol_2)
		//while (ff > fTol_1 || normr > fTol_2)
		{
			//check for the local iteration
			if (kk <= 10)
			{
				mr_rp_2d_out << setw(outputFileWidth) << "yield_f = ----------" << "     "
					<< setw(outputFileWidth) << "norm_R = ----------"
					<< endl;

				mr_rp_2d_out << setw(outputFileWidth) << "global iteration # = " << (*fIteration) << "     "
					<< setw(outputFileWidth) << "local iteration # = " << kk
					<< endl;

				mr_rp_2d_out << setw(outputFileWidth) << "yield_f = " << ff
						<< setw(outputFileWidth) << "norm_R = " << normr
						<< endl;

				mr_rp_2d_out << setw(outputFileWidth) << "Tt = " << Sig[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "Tn = " << Sig[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "jump_u[0] = " << jump_u[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "jump_u[1] = " << jump_u[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "up_t = " << up[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "up_n = " << up[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "dlam = " << dlam << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qn[0] = " << qn[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qn[1] = " << qn[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qn[2] = " << qn[2] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qn[3] = " << qn[3] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qbar[0] = " << qbar[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qbar[1] = " << qbar[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qbar[2] = " << qbar[2] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "qbar[3] = " << qbar[3] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "R[0] = " << R[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "R[1] = " << R[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "R[2] = " << R[2] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "R[3] = " << R[3] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "R[4] = " << R[4] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "R[5] = " << R[5] << endl;

				// Check the stiffness at each local iteration
				Stiffness(jump_u, state, sigma);

				mr_rp_2d_out << setw(outputFileWidth) << "KEP(0,0) = " << fStiffness[0] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "KEP(0,1) = " << fStiffness[1] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "KEP(1,0) = " << fStiffness[2] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "KEP(1,1) = " << fStiffness[3] << endl;

				mr_rp_2d_out << setw(outputFileWidth) << "yield_f = ----------" << "     "
					<< setw(outputFileWidth) << "norm_R = ----------"
					<< endl;
			}
			else
			{
				mr_rp_2d_out << setw(outputFileWidth) << "yield_f = " << ff
					<< setw(outputFileWidth) << "norm_R = " << normr
					<< endl;
			}
			//End

			if (kk > 500) {
				ExceptionT::GeneralFail("MR_NodalRP2DT::Traction","Too Many Iterations");
			}

			dfdSig_f(Sig, qn, dfdSig);
			dQdSig_f(Sig, qn, dQdSig);
			qbar_f(Sig, qn, qbar);
			dfdq_f(Sig, qn, dfdq);
			dQdSig2_f(qn, dQdSig2);
			dQdSigdq_f(Sig, qn, A_uq);
			dqbardSig_f(Sig, qn, A_qu);
			dqbardq_f(Sig, qn, A_qq);

			for (i = 0; i<=5; ++i)
			{
				for (j = 0; j<=5; ++j)
				{
					if (i<=1 & j<=1)
					{
						AA_inv(i,j)  = dQdSig2(i,j);
						AA_inv(i,j) *= dlam;
					}
					if (i<=1 & j>1)
					{
						AA_inv(i,j)  = A_uq(i,j-2);
						AA_inv(i,j) *= dlam;
					}
					if (i>1 & j<=1)
					{
						AA_inv(i,j)  = A_qu(i-2,j);
						AA_inv(i,j) *= dlam;
					} 
					if (i>1 & j>1) 
					{
						AA_inv(i,j)  = I_mat(i-2,j-2);
						AA_inv(i,j) *= -1.;
						AA_inv(i,j) += dlam*A_qq(i-2,j-2);
					} 
				}
			}
	        
			AA.Inverse(AA_inv);
	        
			V_sig = dfdSig;
			V_q = dfdq;

			for (i = 0; i<=5; ++i) 
			{
				if (i<=1)
				{
					Rvec[i] = V_sig[i];
					Cvec[i] = dQdSig[i];
				}

				if (i > 1)
				{
					Rvec[i] = V_q[i-2];
					Cvec[i] = qbar[i-2];
				}
			}

			Yield_f(Sig, qn, ff);
			dArrayT tmpVec(6);
			AA.Multx(R,tmpVec);
			topp = ff;
			topp -= dArrayT::Dot(Rvec,tmpVec);        
			AA.Multx(Cvec,tmpVec);
			bott = dArrayT::Dot(Rvec,tmpVec); 		
			dlam2 = topp/bott;

			for (i = 0; i<=5; ++i)
			{
				if (i<=1) Rmod[i] = dQdSig[i];
				if (i >1) Rmod[i] = qbar[i-2];
			}

			Rmod *= dlam2;
			R2 = R;
			R2 += Rmod;
			AA.Multx(R2,X);
			Y = 0.;
			Y -= X;

			for (i = 0; i<=5; ++i)
			{
				if (i<=1) dSig[i] = Y[i];
				if (i > 1) dq[i-2] = Y[i];
			}
	        
			/*  Update stresses and internal variables */
			Sig += dSig;
			qn  += dq;
			dlam = dlam + dlam2;
			kk = kk + 1;
	        
			/*  Calculation of Yield Function and Residuals for next iteration check */
			Yield_f(Sig, qn, ff);
			dQdSig_f(Sig, qn, dQdSig);
			qbar_f(Sig, qn, qbar);

			if (ff >= 0.0)
			{
				state[k_tied_flag] = kFreeNode;
				iplastic = 1;
			}
			else
			{
				state[k_tied_flag] = kTiedNode;
				iplastic = 0;
			}

			for (i = 0; i<=1; ++i)
			{
				R[i]  = upo[i];
				R[i] -= up[i];
				R[i] += dlam*dQdSig[i];
			}

			for (i = 0; i<=3; ++i)
			{
				R[i+2]  = qo[i];
				R[i+2] -= qn[i];
				R[i+2] += dlam*qbar[i];
			}
			normr = R.Magnitude();
	        
		} //end iteration loop

		/* update the state variables after convergence is achieved */
		state[k_T_t] = Sig[0];
		state[k_T_n] = Sig[1];
		fTraction[0] = state[k_T_t];
		fTraction[1] = state[k_T_n];
		state[k_Delta_t] = jump_u[0];
		state[k_Delta_n] = jump_u[1];
		// state[k_Delta_t] = up[0];
		// state[k_Delta_n] = up[1];
		state[k_up_t] = up[0];
		state[k_up_n] = up[1];
		state[k_fchi] = qn[0];
		state[k_fc]   = qn[1];
		state[k_fphi] = qn[2];
		state[k_fpsi] = qn[3];
		state[k_yield] = ff;
		state[k_plast_mult] = dlam;
		state[k_plastic] = double(iplastic);
		state[k_norm_resid] = normr;
		dQdSig_f(Sig, qn, dQdSig);
		state[k_A] = Sig[0]*dQdSig[0];
		state[k_A] += (Sig[1] + fabs(Sig[1]))*dQdSig[1]/2.;
		state[k_A] /=fGf_I;
		state[k_A] *=dlam;
		state[k_B]  = signof(Sig[0]);
		state[k_B] -= signof(Sig[0])*fabs(Sig[1]*qn[2]);
		state[k_B] *= dQdSig[0];
		state[k_B] /= fGf_II;
		state[k_B] *=dlam;
		state[k_iteration] = double(kk);

		// upadte the state variable of tied flag
		Yield_f(Sig, qn, ff);
		
		if (ff >= 0.0)
		{
			state[k_tied_flag] = kFreeNode;
		}
		else
		{
			state[k_tied_flag] = kTiedNode;
		}

		//check for yield and norm_R after convergence is achieved
		mr_rp_2d_out << setw(outputFileWidth) << "     " << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "yield_f = **********" << "     "
				<< setw(outputFileWidth) << "norm_R = **********"
				<< endl;

		mr_rp_2d_out << setw(outputFileWidth) << "yield_f is converged " << "     "
				<< setw(outputFileWidth) << "norm_R is converged "
				<< endl;

		mr_rp_2d_out << setw(outputFileWidth) << "global iteration # = " << (*fIteration) << "     "
				<< setw(outputFileWidth) << "local iteration # = " << kk
				<< endl;

		mr_rp_2d_out << setw(outputFileWidth) << "yield_f = " << ff
				<< setw(outputFileWidth) << "norm_R = " << normr
				<< endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_T_t] = " << state[k_T_t] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_T_n] = " << state[k_T_n] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_Delta_t] = " << state[k_Delta_t] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_Delta_n] = " << state[k_Delta_n] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_up_t] = " << state[k_up_t] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_up_n] = " << state[k_up_n] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_fchi] = " << state[k_fchi] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_fc] = " << state[k_fc] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_fphi] = " << state[k_fphi] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "state[k_fpsi] = " << state[k_fpsi] << endl;

		// check for the stiffness after convergence is achieved
		Stiffness(jump_u, state, sigma);
		mr_rp_2d_out << setw(outputFileWidth) << "KEP(0,0) = " << fStiffness[0] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "KEP(0,1) = " << fStiffness[1] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "KEP(1,0) = " << fStiffness[2] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "KEP(1,1) = " << fStiffness[3] << endl;

		mr_rp_2d_out << setw(outputFileWidth) << "yield_f = **********" << "     "
				<< setw(outputFileWidth) << "norm_R = **********"
				<< endl;

		mr_rp_2d_out << setw(outputFileWidth) << "     " << endl;
		//End

		return fTraction;
	}

	else
	{
		dArrayT Sig(2); dArrayT qn(4);
		double ff;

		Sig[0] = state[k_T_t];
		Sig[1] = state[k_T_n];
		qn[0] = state[k_fchi];
		qn[1] = state[k_fc];
		qn[2] = state[k_fphi];
		qn[3] = state[k_fpsi];

		// upadte the state variable of tied flag
		Yield_f(Sig, qn, ff);

		if (ff >= 0.0)
		{
			state[k_tied_flag] = kFreeNode;
		}
		else
		{
			state[k_tied_flag] = kTiedNode;
		}

		// Nothing to do for tractions
		fTraction[0] = state[k_T_t];
		fTraction[1] = state[k_T_n];

		return fTraction;
	}
}


double& MR_NodalRP2DT::Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff)
{
	double tmp1, tmp11, tmp22, tmp3, tmp31, tmp32, tmp4, tmp5;
  
	tmp1   = qn[1];
	tmp11  = Sig[1];
	tmp11 *= qn[2];
	tmp1  -= tmp11;
  
	tmp22  = Sig[0];
	tmp22 *= Sig[0];
  
	tmp3   = qn[1];
	tmp31  = qn[0];
	tmp31 *= qn[2];
	tmp3  -= tmp31;
	tmp32  = tmp3;
	tmp32 *= tmp3;

	tmp4  = tmp22;
	tmp4 += tmp32;
  
	tmp5 = sqrt(tmp4);
  
	ff  = tmp5;
	ff -= tmp1;
  
	return ff;
}

/* calculation of qbar_f */
dArrayT& MR_NodalRP2DT::qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar)
{
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
	double DQDT = 2.*Sig[0];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
   
	qbar[0] = A1*B1*DQDN + A1*B2*DQDT;
	qbar[1] = A2*B3*DQDT;
	qbar[2] = A3*B3*DQDT;
	qbar[3] = A4*B3*DQDT;

	return qbar;
}


/* calculation of dQdSig2_f */
dMatrixT& MR_NodalRP2DT::dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2)
{
	dQdSig2(0,0) = 2.;
	dQdSig2(1,1) = -2.*qn[3]*qn[3];
	dQdSig2(0,1) = 0.;
	dQdSig2(1,0) = 0.;
  
	return dQdSig2;
}

/* calculation of dfdSig_f */
dArrayT& MR_NodalRP2DT::dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
	double Shear = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[2])*(qn[1] - qn[0]*qn[2]);
	dfdSig[0] = Sig[0]/sqrt(Shear);
	dfdSig[1] = qn[2];
  
	return dfdSig;
}

/* calculation of dQdSig_f */
dArrayT& MR_NodalRP2DT::dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
	dQdSig[0] = 2.*Sig[0];
	dQdSig[1] = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
  
	return dQdSig;
}


/* calculation of dfdq_f */
dArrayT& MR_NodalRP2DT::dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
	double Shear = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[2])*(qn[1] - qn[0]*qn[2]);
	double zeta = (qn[1] - qn[0]*qn[2])/sqrt(Shear);
	dfdq[0] = -qn[2]*zeta;
	dfdq[1] = zeta - 1.;
	dfdq[2] = -qn[0]*zeta + Sig[1];
	dfdq[3] = 0.;
  
	return dfdq;
}

/* calculation of dQdSigdq_f */
dMatrixT& MR_NodalRP2DT::dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
	dQdSigdq(0,0) = 0.;
	dQdSigdq(0,1) = 0.;
	dQdSigdq(0,2) = 0.;
	dQdSigdq(0,3) = 0.;
	dQdSigdq(1,0) = 0.;
	dQdSigdq(1,1) = 2.*qn[3];
	dQdSigdq(1,2) = 0.;
	dQdSigdq(1,3) = 2.*qn[1] - 4.*Sig[1]*qn[3];
  
	return dQdSigdq;
}

/* calculation of dqbardSig_f */
dMatrixT& MR_NodalRP2DT::dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
	double DQDT = 2.*Sig[0];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))/fGf_II/2.;
	double DB3_DTt = 1./fGf_II;
	double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])/fGf_II;
	double DQDN2 = -2.*qn[3]*qn[3];
	double DQDT2 = 2.;
	double DQDTN = 0.;
	double DQDNT = 0.;
	double SN = signof(Sig[1]);
	double DB1DN = (SN + fabs(SN))/2./fGf_I;
   
	dqbardSig(0,0) = A1*B2*DQDT2 + A1*DQDT/fGf_I;
	dqbardSig(0,1) = A1*B1*DQDN2 + A1*DQDN*DB1DN;
	dqbardSig(1,0) = A2*B3*DQDT2 + A2*DQDT*DB3_DTt;
	dqbardSig(1,1) = A2*DQDT*DB3_DTn;
	dqbardSig(2,0) = A3*B3*DQDT2 + A3*DQDT*DB3_DTt;
	dqbardSig(2,1) = A3*DQDT*DB3_DTn;
	dqbardSig(3,0) = A4*B3*DQDT2 + A4*DQDT*DB3_DTt;
	dqbardSig(3,1) = A4*DQDT*DB3_DTn;
 
	return dqbardSig;
}
  
/* calculation of dqbardq_f */
dMatrixT& MR_NodalRP2DT::dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
	double DQDT = 2.*Sig[0];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))/fGf_II/2.;
	double DB3_DTt = 1./fGf_II;
	double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])/fGf_II;
	double DQDN2 = -2.*qn[3]*qn[3];
	double DQDT2 = 2.;
	double DQDTN = 0.;
	double DQDNT = 0.;
	double SN = signof(Sig[1]);
	double DB1DN = (SN + fabs(SN))/2./fGf_I;
   
	dqbardq(0,0) = -falpha_chi*(B1*DQDN + B2*DQDT);
	dqbardq(0,1) =  A1*B1*(2.*qn[3]);
	dqbardq(0,2) = 0.;
	dqbardq(0,3) =  A1*B1*(2.*qn[1] - 4.*Sig[1]*qn[3]);
	dqbardq(1,0) = 0.;
	dqbardq(1,1) = -falpha_c*B3*DQDT;
	dqbardq(1,2) = A2*DQDT*DB3_DTanphi;
	dqbardq(1,3) = 0.;
	dqbardq(2,0) = 0.;
	dqbardq(2,1) = 0.;
	dqbardq(2,2) = -falpha_phi*B3*DQDT + A3*DQDT*DB3_DTanphi;
	dqbardq(2,3) = 0.;
	dqbardq(3,0) = 0.;
	dqbardq(3,1) = 0.;
	dqbardq(3,2) = A4*DQDT*DB3_DTanphi;
	dqbardq(3,3) = -falpha_psi*B3*DQDT;

	return dqbardq;
}

/* potential stiffness */
const dMatrixT& MR_NodalRP2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state,
	const dArrayT& sigma)
{
	const char caller[] = "MR_NodalRP2DT::Stiffness";

#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) ExceptionT::SizeMismatch(caller);
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch(caller);
#endif

	int i, j;

	dMatrixT AA(6,6), I_mat(4,4), CMAT(6,6),AA_inv(6,6),
	         A_qq(4,4), A_uu(2,2), A_uq(2,4), A_qu(4,2), ZMAT(2,4),
	         ZMATP(4,2), dQdSig2(2,2), dqdbar(4,4), dqbardSig(4,2),
	         dQdSigdq(2,4), KP(2,2), KP2(2,2), KEP(2,2);
	         
	dMatrixT I_m(2,2), Rmat(2,2), R_Inv(2,2), KE(2,2), KE_Inv(2,2),
	         Ch(4,4), Ch_Inv(4,4), KE1(4,2), KE2(2,2), KE3(2,4);

	dArrayT  u(2), up(2), du(2), dup(2), qn(4), qo(4), Rvec(6),Cvec(6),
	         R(6), Rmod(6), Sig(2), Sig_I(2), dQdSig(2), dfdq(4), qbar(4),
	         R2(6), X(6), V_sig(2), V_q(4), dfdSig(2), K1(2), K2(2);

	double bott, dlam, ff;

	fStiffness[1] = fStiffness[2] = 0.;
	I_m(0,0) = 1.; I_m(0,1) =0.; I_m(1,0) = 0.; I_m(1,1) = 1.;
	I_mat = 0.;
	qn[0] = state[k_fchi];
	qn[1] = state[k_fc];
	qn[2] = state[k_fphi];
	qn[3] = state[k_fpsi];

	for (i = 0; i<=3; ++i) 
	{
		I_mat(i,i) = 1.;
	}
    
	Sig[0] = state[k_T_t];
	Sig[1] = state[k_T_n];

	// not free
	
	//value of Parameter fUpdateIterations doesn't be obtained from the input
	// if (fabs(state[k_tied_flag] - kFreeNode) > kSmall ||
	//	jump_u[0] < kSmall && jump_u[1] < kSmall ||
	//	fUpdateIterations > 1) //subiterations

	if (fabs(state[k_tied_flag] - kFreeNode) > kSmall ||
		jump_u[0] < kSmall && jump_u[1] < kSmall)
	{
		//no contribution to the stiffness
		fStiffness = 0.;

		return fStiffness;
	}

	//check for the plastic flag
	
	if (state[k_plastic] == 0.)
	{
		//no contribution to the stiffness
		// fStiffness[0] = 1.e20;
		// fStiffness[3] = 1.e20;
		fStiffness[0] = 0.;
		fStiffness[3] = 0.;
	}
	else if (state[k_plastic] == 1.)
	{
		dlam = state[k_plast_mult];
		dQdSig2_f(qn, dQdSig2);
		dqbardSig_f(Sig, qn, A_qu);
		dqbardq_f(Sig, qn, A_qq);
		dQdSigdq_f(Sig, qn, A_uq);
		Ch  = A_qq;
		Ch *= -dlam;
		Ch += I_mat;
		Ch_Inv.Inverse(Ch);
		KE1.MultAB(Ch_Inv,A_qu);
		KE.MultAB(A_uq,KE1);
		KE *= state[k_plast_mult];
		KE *= state[k_plast_mult];
		// KE = 0.;
		KE2  = dQdSig2;
		KE2 *= state[k_plast_mult];
		KE  += KE2;
	        
		KE_Inv.Inverse(KE);
	     
		for (i = 0; i<=5; ++i) 
		{
			for (j = 0; j<=5; ++j)
			{
				if (i<=1 & j<=1)
				{
					AA_inv(i,j)  = 0.;
					AA_inv(i,j) += dlam*dQdSig2(i,j);
				}
				if (i<=1 & j>1)
				{
					AA_inv(i,j)  = A_uq(i,j-2);
					AA_inv(i,j) *= dlam;
				} 
				if (i>1 & j<=1)
				{
					AA_inv(i,j)  = A_qu(i-2,j);
					AA_inv(i,j) *= dlam;
				} 
				if (i>1 & j>1)
				{
					AA_inv(i,j)  = I_mat(i-2,j-2);
					AA_inv(i,j) *= -1.;
					AA_inv(i,j) += dlam*A_qq(i-2,j-2);
				} 
			}
		}
		AA.Inverse(AA_inv);
	
		dfdSig_f(Sig, qn, dfdSig);
		V_sig = dfdSig;
		dfdq_f(Sig,qn, dfdq);
		V_q = dfdq;
		dQdSig_f(Sig, qn, dQdSig);
		qbar_f(Sig, qn, qbar);

		for (i = 0; i<=5; ++i)
		{
			if (i<=1)
			{
				Rvec[i] = V_sig[i];
				Cvec[i] = dQdSig[i];
			}
			if (i>1)
			{
				Rvec[i] = V_q[i-2];
				Cvec[i] = qbar[i-2];
			}
		}

		dArrayT tmpVec(6), Vvec(2), dVec(2);
		AA.Multx(Cvec,tmpVec);
		bott = dArrayT::Dot(Rvec,tmpVec);
            
		for (i = 0; i<=1; ++i)
		{
			Vvec[i] = 0.;
			for (j = 0; j<=5; ++j) Vvec[i] += Rvec[j]*AA(j,i);
		}
            
		for (i = 0; i<=1; ++i)
		{
			for (j = 0; j<=1; ++j) KP(i,j) = dQdSig[i]*Vvec[j];
		}
            
		KE3.MultAB(A_uq, Ch_Inv);
		KE3.Multx(qbar,dVec);

		for (i = 0; i<=1; ++i)
		{
			for (j = 0; j<=1; ++j) KP2(i,j) = dVec[i]*Vvec[j];
		}
	        
		KP2 *= state[k_plast_mult];
		KP  += KP2;
		KP  /= -bott;
		KP  += I_m;
		KEP.MultAB(KE_Inv, KP);
	/*
 	int i, j;

	dMatrixT KEP(2,2), KEP_Inv(2,2);

	dArrayT qn(4), Sig(2), dQdSig(2), dfdq(4), qbar(4), dfdSig(2);

	double ff, HP;

	fStiffness[1] = fStiffness[2] = 0.;
	HP = 0.;
	Sig[0] = state[k_T_t];
	Sig[1] = state[k_T_n];
	qn[0] = state[k_fchi];
	qn[1] = state[k_fc];
	qn[2] = state[k_fphi];
	qn[3] = state[k_fpsi];

	//check for the yield and tied flag
	Yield_f(Sig, qn, ff);

	mr_rp_2d_out << setw(outputFileWidth) << "yield_f = " << ff << "     "
				<< setw(outputFileWidth) << "state[k_tied_flag] = " << state[k_tied_flag]
				<< endl;

	// not free
	if (fabs(state[k_tied_flag] - kFreeNode) > kSmall ||
		jump_u[0] < kSmall && jump_u[1] < kSmall)
	{
		//no contribution to the stiffness
		fStiffness = 0.;

		return fStiffness;
	}

	if (state[k_plastic] == 0.)
	{
		//no contribution to the stiffness
		fStiffness[0] = 1.e20;
		fStiffness[3] = 1.e20;
	}
	else if (state[k_plastic] == 1.)
	{
		dfdSig_f(Sig, qn, dfdSig);
		dQdSig_f(Sig, qn, dQdSig);
		qbar_f(Sig, qn, qbar);
		dfdq_f(Sig, qn, dfdq);

		HP -= dArrayT::Dot(dfdq,qbar);

		KEP.Outer(dQdSig,dfdSig);
		KEP_Inv.Inverse(KEP);

		for (i = 0; i <= 1; ++i){
			for (j = 0; j<=1; ++j){
				KEP(i,j) = KEP_Inv(i,j)*HP;
			}
		}
	*/
		fStiffness[0] = KEP(0,0);
		fStiffness[1] = KEP(0,1);
		fStiffness[2] = KEP(1,0);
		fStiffness[3] = KEP(1,1);
	}

	return fStiffness;
}

/* surface status */
SurfacePotentialT::StatusT MR_NodalRP2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) ExceptionT::SizeMismatch("MR_NodalRP2DT::Status");
#endif
	
	if (state[k_yield] < 0.0)
		return Precritical;
	else if (state[k_yield] > 0.0)
		return Critical;
	else
		return Critical;

	/* else if ((jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])>100000.)
		return Failed; */
}

void MR_NodalRP2DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    Rigid plastic MR 2D \n";
#endif
}

/* print parameters to the output stream */
void MR_NodalRP2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Mode_I Fracture Energy            . . . . . . . = " << fGf_I     << '\n';
	out << " Mode_II Fracture Energy            . . .  . . . = " << fGf_II << '\n';
	out << " Peak Cohesion                 . . . . . . . . . = " << fc_p    << '\n';
	out << " Residual Cohesion             . . . . . . . . . = " << fc_r    << '\n';
	out << " Peak Tensile Strength   . . . . . . . . . . . . = " << fchi_p << '\n';
	out << " Residual Tensile Strength . . . . . . . . . . . = " << fchi_r << '\n';
	out << " Peak Friction Angle         . . . . . . . . . . = " << fphi_p   << '\n';
	out << " Critical State Friction Angle       . . . . . . = " << fphi_r  << '\n';
	out << " Peak Dilation Angle.. . . . . . . . . . . . . . = " << fpsi_p   << '\n';
	out << " Coefficient of Tensile Strength Degradation ..  = " << falpha_chi   << '\n';
	out << " Coefficient of Cohesion Degradation. .. . . . . = " << falpha_c   << '\n';
	out << " Coefficient for Frictional Angle Degradation .  = " << falpha_phi   << '\n';
	out << " Coefficient for Dilation Angle Degradation  . . = " << falpha_psi   << '\n';
	out << " Error Tolerance for Yield Function. . . . . . . = " << fTol_1 << '\n';
	out << " Error Tolerance for Residual    . . . . . . . . = " << fTol_2 << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int MR_NodalRP2DT::NumOutputVariables(void) const { return 8; }

void MR_NodalRP2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(8);
	labels[0] = "up_t";
	labels[1] = "up_n";
	labels[2] = "Chi";
	labels[3] = "Cohesion";
	labels[4] = "Friction Angle";
	labels[5] = "Yield Function Value";
	labels[6] = "Norm of residuals";
	labels[7] = "No. of Iterations";
}

/* release condition depends on this quantity */
int MR_NodalRP2DT::NodalQuantityNeeded(void) const
{
	//what is this?, RAR
	return 2;
}


/*************************************************************************
 * Protected
 *************************************************************************/

void MR_NodalRP2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) ExceptionT::GeneralFail("MR_NodalRP2DT::ComputeOutput");
#endif	
	output[0] = state[k_up_t];
	output[1] = state[k_up_n];
	output[2] = state[k_fchi];
	output[3] = state[k_fc];
	//output[4] = state[k_fphi];
	output[4] = atan(state[k_fphi]);  // state[k_fphi] = ftan_phi
	output[5] = state[k_yield];
	output[6] = state[k_norm_resid];
	output[7] = state[k_iteration];
}

