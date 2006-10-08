/*$Id: MR_RP2DT.cpp,v 1.29 2006-10-08 19:22:19 regueiro Exp $*/
/* created by manzari*/
/* Rigid Plastic Cohesive Model for Geomaterials*/
#include "MR_RP2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nMatrixT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 2;
const int nTiedFlag = 8*knumDOF+1;

/* constructor */
MR_RP2DT::MR_RP2DT(void): SurfacePotentialT(knumDOF),
	TiedPotentialBaseT()
{
	const char caller[] = "MR_RP2DT::MR_RP2DT";
	
	SetName("rigid-plastic_MR_RP2D");

#if 0	
	/* Elastic and Fracrure Energy parameters */
	in >> fGf_I; if (fGf_I < 0) throw ExceptionT::kBadInputValue;
	in >> fGf_II; if (fGf_II < 0) throw ExceptionT::kBadInputValue;

	/* Inelastic Response initiation parameters */
	in >> fchi_p; if (fchi_p < 0) throw ExceptionT::kBadInputValue;
	in >> fchi_r; if (fchi_r < 0) throw ExceptionT::kBadInputValue;
	in >> fc_p; if (fc_p < 0) throw ExceptionT::kBadInputValue;
	in >> fc_r; if (fc_r < 0) throw ExceptionT::kBadInputValue;
    in >> fphi_p; if (fphi_p < 0) throw ExceptionT::kBadInputValue;
	in >> fphi_r; if (fphi_r < 0) throw ExceptionT::kBadInputValue;
	in >> fpsi_p; if (fpsi_p < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_chi; if (falpha_chi < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_c; if (falpha_c < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_phi; if (falpha_phi < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_psi; if (falpha_psi < 0) throw ExceptionT::kBadInputValue;
	in >> fTol_1; if (fTol_1 < 0) throw ExceptionT::kBadInputValue;
	in >> fTol_2; if (fTol_2 < 0) throw ExceptionT::kBadInputValue;
#endif	
}

/* return the number of state variables needed by the model */
int MR_RP2DT::NumStateVariables(void) const { return 8*knumDOF +2; }

/* describe the parameters needed by the interface */
void MR_RP2DT::DefineParameters(ParameterListT& list) const
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
void MR_RP2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* MR_RP2DT::NewSub(const StringT& name) const
{
	return SurfacePotentialT::NewSub(name);
}

/* accept parameter list */
void MR_RP2DT::TakeParameterList(const ParameterListT& list)
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
}


/* location in state variable array of the state flag */
int MR_RP2DT::TiedStatusPosition(void) const { return nTiedFlag; }

/* initialize the state variable array */
void MR_RP2DT::InitStateVariables(ArrayT<double>& state)
{
	int num_state = NumStateVariables();
	if (state.Length() != num_state) {
#ifndef _SIERRA_TEST_	
		cout << "\n SurfacePotentialT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
#endif
		throw ExceptionT::kSizeMismatch;	
	}

	/* clear */
	if (num_state > 0) state = 0.0;
	
	/* Initializing internal state variables */
	double enp  = state[14];
	double esp  = state[15];
	double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
	double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
	double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
	double ftan_psi = (tan(fpsi_p))*exp(-falpha_psi*esp);	
	state[6]  = fchi;
	state[7]  = fc ;
    state[8]  = ftan_phi;
    state[9]  = ftan_psi;
    state[13] = 0.;
    state[nTiedFlag] = kTiedNode;
}

/* Value of the Yield Function */ 
double MR_RP2DT::YFValue(const ArrayT<double>& state) 
{
	return state[10]; 
}	

/** dissipated energy. Total amount of energy dissipated reaching
	 * the current state. */
double MR_RP2DT::FractureEnergy(const ArrayT<double>& state)
{
    double FE = state[0]*state[4] + state[1]*state[5];
    return FE;
}

/** potential energy */
double MR_RP2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
    double PT = (jump_u[0]*state[0] + jump_u[1]*state[1]);
    return PT;
}

/* surface status */
SurfacePotentialT::StatusT MR_RP2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
    if (state[10]<0.){
       int StatusT = 0;
    }
    if (state[10]>0) {
       int StatusT = 1;
    }
    if ((jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])>100000.) {
       int StatusT = 2;
    }
      
    /*return SurfacePotentialT::Precritical;*/
    return StatusT();
}

void MR_RP2DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    MR process zone 2D \n";
#endif
}

		
/* traction vector given displacement jump vector */	

const dArrayT& MR_RP2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

if (state[nTiedFlag] == kTiedNode)
{  
	fTraction = 0.;		
	return fTraction;
}
else 
{  	/* nodes are either untied or will be next step */

	if (state[nTiedFlag] == kReleaseNextStep) 
	{	/* nodes will be untied next step. Store stress and move along */
		/* qIntegrate is guaranteed to be true */
		if (!qIntegrate)
			ExceptionT::GeneralFail("MR_RP2DT::Traction","nodes freed and !qIntegrate");
		state[0] = sigma[2];
		state[1] = sigma[1];
		state[nTiedFlag] = kFirstFreeStep;
		fTraction = 0.;
			
		return fTraction;
	} 

	if (state[nTiedFlag] == kFirstFreeStep) /* First timestep with free nodes */
	{  
		if (qIntegrate)
			state[nTiedFlag] = kFreeNode;
	}
	
	if (fabs(jump_u[0]) < kSmall && fabs(jump_u[1]) < kSmall)
	{
		fTraction = 0.;
		return fTraction;
	}

	if (!qIntegrate) 
	{   // Nothing to do here
	    fTraction[0] = state[0];
	    fTraction[1] = state[1];
	    return fTraction;
	}
	else
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
		    
		for (i = 0; i<=1; ++i) {
	       up[i] = jump_u[i];
	       dup[i] = up[i] - state[i+2];
	       upo[i] = state[i+4];
	       Sig_I[i] = state[i];
	    }
		    
	    for (i = 0; i<=3; ++i) {
	        qn[i] = state[i+6];
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
	    
	    if (normdup <=kSmall) {
	    fTraction[0] = state[0];
	    fTraction[1] = state[1];
	    return fTraction;
	    }
	    
		/* calculate residuals */    
	    for (i = 0; i<=1; ++i) {
	          R[i]  = upo[i];
	          R[i] -= up[i];
	          R[i] += dlam*dQdSig[i];
	    }
	        for (i = 0; i<=3; ++i) {
	          R[i+2]  = qo[i];
	          R[i+2] -= qn[i];
	          R[i+2] += dlam*qbar[i];
	    }
	    normr = R.Magnitude();
	    
	    
		/* Local Iteration */
	    kk = 0;
	    iplastic = 1;
	    Yield_f(Sig, qn, ff);
	    while (ff > fTol_1 | normr > fTol_2) {
	        if (kk > 500) {
	        	ExceptionT::GeneralFail("MR_RP2DT::Traction","Too Many Iterations");
	        }
	        
	        dfdSig_f(Sig, qn, dfdSig);
	        dQdSig_f(Sig, qn, dQdSig);
	        qbar_f(Sig, qn, qbar);
	        dfdq_f(Sig, qn, dfdq);
	        dQdSig2_f(qn, dQdSig2);
	        dQdSigdq_f(Sig, qn, A_uq);
	        dqbardSig_f(Sig, qn, A_qu);
	        dqbardq_f(Sig, qn, A_qq);
	        
	        for (i = 0; i<=5; ++i) {
	          for (j = 0; j<=5; ++j) {
	            if (i<=1 & j<=1){
	             AA_inv(i,j)  = dQdSig2(i,j);
	             AA_inv(i,j) *= dlam;
	            }
	            if (i<=1 & j>1){
	             AA_inv(i,j)  = A_uq(i,j-2);
	             AA_inv(i,j) *= dlam;
	            } 
	            if(i>1 & j<=1){
	             AA_inv(i,j)  = A_qu(i-2,j);
	             AA_inv(i,j) *= dlam;
	            } 
	            if(i>1 & j >1) {
	             AA_inv(i,j)  = I_mat(i-2,j-2);
	             AA_inv(i,j) *= -1.;
	             AA_inv(i,j) += dlam*A_qq(i-2,j-2);
	            } 
	          }
	        }
	        
	        AA.Inverse(AA_inv);
	        
	        V_sig = dfdSig;
	        V_q = dfdq;
	        for (i = 0; i<=5; ++i) {
	            if (i<=1){
	             Rvec[i] = V_sig[i];
	             Cvec[i] = dQdSig[i];
	            }
	            if (i > 1){
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
	         for (i = 0; i<=5; ++i) {
	            if (i<=1){
	             Rmod[i] = dQdSig[i];
	            }
	            if (i >1){
	             Rmod[i] = qbar[i-2];
	            }
	        }
	        Rmod *= dlam2;
	        R2 = R;
	        R2 += Rmod;
	        AA.Multx(R2,X);
	        Y = 0.;
	        Y -= X;
	        for (i = 0; i<=5; ++i) {
	            if (i<=1) {
	             dSig[i] = Y[i];
	            }
	            if (i > 1) {
	             dq[i-2] = Y[i];
	            }
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
	        
	        for (i = 0; i<=1; ++i) {
	          R[i]  = upo[i];
	          R[i] -= up[i];
	          R[i] += dlam*dQdSig[i];
	        }
	        for (i = 0; i<=3; ++i) {
	          R[i+2]  = qo[i];
	          R[i+2] -= qn[i];
	          R[i+2] += dlam*qbar[i];
	        }
	        normr = R.Magnitude();
	        
	      }
	      
		/* update the state variables after convergence is achieved */
	    state[0] = Sig[0];
	    state[1] = Sig[1];     
		fTraction[0] = state[0];
		fTraction[1] = state[1];
		state[2] = jump_u[0];
		state[3] = jump_u[1];
		state[4] = up[0];
		state[5] = up[1];
		state[6] = qn[0];
		state[7] = qn[1];
		state[8] = qn[2];
		state[9] = qn[3];
		state[10] = ff;
		state[11] = dlam;
		state[12] = double(iplastic);
		state[13] = normr;
		dQdSig_f(Sig, qn, dQdSig);
		state[14] = Sig[0]*dQdSig[0];
		state[14] += (Sig[1] + fabs(Sig[1]))*dQdSig[1]/2.;
		state[14] /=fGf_I;
		state[14] *=dlam;
		state[15]  = signof(Sig[0]);
		state[15] -= signof(Sig[0])*fabs(Sig[1]*qn[2]);
		state[15] *= dQdSig[0];
		state[15] /= fGf_II;
		state[15] *=dlam;
		state[16] = double(kk);

		return fTraction;
		}
	}
}

/* calculation of Yield_f */

double& MR_RP2DT::Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff)
{
  double tmp1, tmp11, tmp22, tmp3, tmp31, tmp32, tmp4, tmp5;
  
  tmp1   = qn[1];
  tmp11  = Sig[1];
  tmp11 *= qn[2];
  tmp1  -= tmp11;
  
  tmp22  = Sig[0];
  tmp22 *= Sig[0];
  
  tmp3  = qn[1];
  tmp31 = qn[0];
  tmp31 *= qn[2];
  tmp3 -= tmp31;
  tmp32 = tmp3;
  tmp32 *=tmp3;
  
  tmp4  = tmp22;
  tmp4 += tmp32;
  
  tmp5 = sqrt(tmp4);
  
  ff = tmp5;
  ff -=tmp1;
  
  return ff;
}


/* calculation of qbar_f */

dArrayT& MR_RP2DT::qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar)
{
 
   double A1 = -falpha_chi*(qn[0] - fchi_r);
   double B1 = (Sig[1]+fabs(Sig[1]))/2./fGf_I;
   double B2 = Sig[0]/fGf_I;
   double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
   double DQDT = 2.*Sig[0];
   double A2 = -falpha_c*(qn[1] - fc_r);
   double TNA = (Sig[1]-fabs(Sig[1]))/2.;
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

dMatrixT& MR_RP2DT::dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2)
{
  dQdSig2(0,0) = 2.;
  dQdSig2(1,1) = -2.*qn[3]*qn[3];
  dQdSig2(0,1) = 0.;
  dQdSig2(1,0) = 0.;
  
  return dQdSig2;
}

/* calculation of dfdSig_f */

dArrayT& MR_RP2DT::dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
  double Shear = Sig[0]*Sig[0]+(qn[1]-qn[0]*qn[2])*(qn[1]-qn[0]*qn[2]);
  dfdSig[0] = Sig[0]/sqrt(Shear);
  dfdSig[1] = qn[2];
  
  return dfdSig;
}

/* calculation of dQdSig_f */

dArrayT& MR_RP2DT::dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
  dQdSig[0] = 2.*Sig[0];
  dQdSig[1] = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
  
  return dQdSig;
}


/* calculation of dfdq_f */

dArrayT& MR_RP2DT::dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
  double Shear = Sig[0]*Sig[0]+(qn[1]-qn[0]*qn[2])*(qn[1]-qn[0]*qn[2]);
  double zeta = (qn[1] - qn[0]*qn[2])/sqrt(Shear);
  dfdq[0] = -qn[2]*zeta;
  dfdq[1] = zeta - 1.;
  dfdq[2] = -qn[0]*zeta + Sig[1];
  dfdq[3] = 0.;
  
  return dfdq;
}

/* calculation of dQdSigdq_f */

dMatrixT& MR_RP2DT::dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
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

dMatrixT& MR_RP2DT::dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{

   double A1 = -falpha_chi*(qn[0] - fchi_r);
   double B1 = (Sig[1]+fabs(Sig[1]))/2./fGf_I;
   double B2 = Sig[0]/fGf_I;
   double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
   double DQDT = 2.*Sig[0];
   double A2 = -falpha_c*(qn[1] - fc_r);
   double TNA = (Sig[1]-fabs(Sig[1]))/2.;
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
   double DB1DN = (SN +fabs(SN))/2./fGf_I;
   
   dqbardSig(0,0) = A1*B2*DQDT2 + A1*DQDT/fGf_I;
   dqbardSig(0,1) = A1*B1*DQDN2 + A1*DQDN*DB1DN;
   dqbardSig(1,0) = A2*B3*DQDT2 + A2*DQDT*DB3_DTt;
   dqbardSig(1,1) = A2*DQDT*DB3_DTn;
   dqbardSig(2,0) = A3*B3*DQDT2 + A3*DQDT*DB3_DTt;
   dqbardSig(2,1) = A3*DQDT*DB3_DTn;
   dqbardSig(3,0) = A3*B3*DQDT2 + A4*DQDT*DB3_DTt;
   dqbardSig(3,1) = A4*DQDT*DB3_DTn;
 
    return dqbardSig;
}
  
/* calculation of dqbardq_f */

dMatrixT& MR_RP2DT::dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{

   double A1 = -falpha_chi*(qn[0] - fchi_r);
   double B1 = (Sig[1]+fabs(Sig[1]))/2./fGf_I;
   double B2 = Sig[0]/fGf_I;
   double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
   double DQDT = 2.*Sig[0];
   double A2 = -falpha_c*(qn[1] - fc_r);
   double TNA = (Sig[1]-fabs(Sig[1]))/2.;
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
   double DB1DN = (SN +fabs(SN))/2./fGf_I;
   
   dqbardq(0,0) = -falpha_chi*(B1*DQDN + B2*DQDT);
   dqbardq(0,1) =  A1*B1*(2.*qn[3]);
   dqbardq(0,2) = 0.;
   dqbardq(0,3) =  A1*B1*(2.*qn[1]-4.*Sig[1]*qn[3]);   
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


/* rigidplastic consistent tangent operator*/
const dMatrixT& MR_RP2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	if (state[nTiedFlag] != kFreeNode || jump_u[0] < kSmall && jump_u[1] < kSmall)
	{
		fStiffness = 0.;
		return fStiffness;
	}

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
         
double bott, dlam;
	
	fStiffness[1] = fStiffness[2] = 0.;
	I_m(0,0) = 1.; I_m(0,1) =0.; I_m(1,0) = 0.; I_m(1,1) = 1.;
	I_mat = 0.;
	for (i = 0; i<=3; ++i) {
      qn[i] = state[i+6];
      I_mat(i,i) = 1.;
    }
    
	Sig[0] = state[0];
	Sig[1] = state[1];
    

	if (state[12] == 0.) 
	{
	    fStiffness[3] = 1.e20;
	    fStiffness[0] = 1.e20;
	}
	else 
	  	if (state[12] == 1.) 
	  	{
	  	    dlam = state[11];
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
	        KE *= state[11];
	        KE *= state[11];
	        /*KE = 0.;*/
	        KE2 = dQdSig2;
	        KE2 *=state[11];
	        KE += KE2;
	        
	        KE_Inv.Inverse(KE);
	     
            for (i = 0; i<=5; ++i) {
              for (j = 0; j<=5; ++j) {
                if (i<=1 & j<=1){
                 AA_inv(i,j)  = 0.;
                 AA_inv(i,j) += dlam*dQdSig2(i,j);
                }
                if (i<=1 & j>1){
                  AA_inv(i,j) = A_uq(i,j-2);
                  AA_inv(i,j) *= dlam;
                } 
                if(i>1 & j<=1){
                  AA_inv(i,j) = A_qu(i-2,j);
                  AA_inv(i,j) *= dlam;
                } 
                if(i>1 & j >1) {
                  AA_inv(i,j)  = I_mat(i-2,j-2);
                  AA_inv(i,j)  *= -1.; 
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
            for (i = 0; i<=5; ++i) {
              if (i<=1) {
                Rvec[i] = V_sig[i];
                Cvec[i] = dQdSig[i];
              }
              if (i>1) {
                Rvec[i] = V_q[i-2];
                Cvec[i] = qbar[i-2];
              }
            }
            dArrayT tmpVec(6), Vvec(2), dVec(2);
            AA.Multx(Cvec,tmpVec);
            bott = dArrayT::Dot(Rvec,tmpVec);
            
            for (i = 0; i<=1; ++i) {
                  Vvec[i] = 0.;
	   		    for (j = 0; j<=5; ++j) {
	   		      Vvec[i] += Rvec[j]*AA(j,i);
                }
	        }
            
            for (i = 0; i<=1; ++i) {
	   		    for (j = 0; j<=1; ++j) {
	   		      KP(i,j) = dQdSig[i]*Vvec[j];
                }
	        }
            
            KE3.MultAB(A_uq, Ch_Inv);
            KE3.Multx(qbar,dVec);
            for (i = 0; i<=1; ++i) {
	   		    for (j = 0; j<=1; ++j) {
	   		      KP2(i,j) = dVec[i]*Vvec[j];
                }
	        }
	        
	        KP2 *= state[11];
	        KP += KP2;
	        KP /= -bott;
            KP += I_m;
            KEP.MultAB(KE_Inv, KP);
 
	   		fStiffness[0] = KEP(0,0);
	   		fStiffness[1] = KEP(0,1);
	   		fStiffness[2] = KEP(1,0);
	   		fStiffness[3] = KEP(1,1);
	       }

	return fStiffness;

}


/* print parameters to the output stream */
void MR_RP2DT::Print(ostream& out) const
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
int MR_RP2DT::NumOutputVariables(void) const { return 8; }

void MR_RP2DT::OutputLabels(ArrayT<StringT>& labels) const
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

void MR_RP2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	
	output[0] = state[4];
	output[1] = state[5];;
	output[2] = state[6];
	output[3] = state[7];
	output[4] = state[8];
	output[5] = state[10];
	output[6] = state[13];
	output[7] = state[16];
}

bool MR_RP2DT::NeedsNodalInfo(void) const
{ 
	return true; 
}

int MR_RP2DT::NodalQuantityNeeded(void) const
{ 
	return 2; 
}

bool MR_RP2DT::NodesMayRetie(void) const
{	
	return false;
}

double MR_RP2DT::signof(double r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

bool MR_RP2DT::InitiationQ(const nArrayT<double>& sigma) const
{

  double tmp1, tmp11, tmp2, tmp3, tmp31, tmp32, tmp4, tmp5, ff;
  dArrayT Sig(2);
  
  double enp = 0.;
  double esp = 0.;
  double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
  double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
  double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
  
  Sig[0] = sigma[2];
  Sig[1] = sigma[1];

  tmp11  = Sig[1];
  tmp11 *= ftan_phi;
  tmp1   = fc;
  tmp1  -= tmp11;
  
  tmp2  = Sig[0];
  tmp2 *= Sig[0];
  
  tmp31  = fchi;
  tmp31  *= ftan_phi;
  tmp3 = fc; 
  tmp3 -= tmp31;
  tmp32 = tmp3;
  tmp32 *=tmp3;
  
  tmp4 = tmp2;
  tmp4 += tmp32;
  
  tmp5 = sqrt(tmp4);
  
  ff = tmp5;
  ff -=tmp1;
  
  return ff >= 0.;
}