#include "MR2DT.h"

#include <iostream.h>
#include <math.h>

#include "ExceptionT.h"
#include "fstreamT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nMatrixT.h"
#include "utils.h"

/* class parameters */

using namespace Tahoe;

const int knumDOF = 2;

/* constructor */
MR2DT::MR2DT(ifstreamT& in): SurfacePotentialT(knumDOF)
{
	/* Elastic and Fracrure Energy parameters */
	in >> fE_n; if (fE_n < 0) throw ExceptionT::kBadInputValue;
	in >> fE_t; if (fE_t < 0) throw ExceptionT::kBadInputValue;
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
	double esp = 0.;
	double enp = 0.;
	fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
	fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
	fphi = fphi_r + (fphi_p - fphi_r)*exp(-falpha_phi*esp);
	fpsi = (fpsi_p)*exp(-falpha_psi*esp);
}

/* return the number of state variables needed by the model */
int MR2DT::NumStateVariables(void) const { return 6*knumDOF; }

/* initialize the state variable array */
void MR2DT::InitStateVariables(ArrayT<double>& state)
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
	
	state[1] = fchi;
	state[2] = fc ;
    state[3] = tan(fphi);
    state[4] = tan(fpsi);	
}

/* Value of the Yield Function */ 
double MR2DT::YFValue(const ArrayT<double>& state) 
{
	return state[12]; 
}	

	
/* traction vector given displacement jump vector */	
const dArrayT& MR2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

int i; int j; int kk;

dMatrixT AA(6,6); dMatrixT KE(2,2); dMatrixT KE_Inv(2,2); dMatrixT I_mat(4,4); 
dMatrixT CMAT(4,4); dMatrixT A_qq(4,4); dMatrixT A_uu(2,2); dMatrixT A_uq(2,4);
dMatrixT A_qu(4,2); dMatrixT ZMAT(2,4); dMatrixT ZMATP(4,2);
dMatrixT dQdSig2(2,2); dMatrixT dqbardq(4,4); dMatrixT dQdSigdq(2,4);
dMatrixT dqbardSig(4,2); dMatrixT AA_inv(6,6);

/*dMatrixT AA[6][6]; dMatrixT KE[2][2]; dMatrixT KE_Inv[2][2]; dMatrixT I_mat[4][4]; 
dMatrixT CMAT[4][4]; dMatrixT A_qq[4][4]; dMatrixT A_uu[2][2]; dMatrixT A_uq[2][4];
dMatrixT A_qu[4][2]; dMatrixT ZMAT[4][4]; dMatrixT ZMATP[4][4];
dMatrixT dqdSig2[2][2]; dMatrixT dqdbar[4][4];*/

dArrayT u(2); dArrayT up(2); dArrayT du(2); dArrayT dup(2); dArrayT qn(4);
dArrayT qo(4); dArrayT Rvec(6); dArrayT Cvec(6); dArrayT upo(2);
dArrayT R(6); dArrayT Rmod(6); dArrayT Sig(2); dArrayT Sig_I(2);
dArrayT dQdSig(2); dArrayT dfdq(4); dArrayT qbar(4);
dArrayT R2(6); dMatrixT X(6,1); dArrayT V_sig(2); dArrayT V_q(4); 
dArrayT dfdSig(2); dArrayT dq(4);

/*dArrayT u[2]; dArrayT up[2];dArrayT du[2]; dArrayT dup[2];
dArrayT qn[4];dArrayT qo[4]; dArrayT Rvec[6]; dArrayT Cvec[6]; 
dArrayT upo[2]; dArrayT R[6]; dArrayT Rmod[6]; dArrayT Sig[2];
dArrayT Sig_I[2]; dArrayT dqdSig[2]; dArrayT dfdq[4]; dArrayT qbar[4];
dArrayT R2[6]; dArrayT X[6]; dArrayT V_sig[2]; dArrayT V_q[4]; 
dArrayT dfdSig[2];*/

double ff; double dlam; double dlam2; double normr;

/* Calculate incremental jumps and initialize the neecessary vectors */
    for (i = 0; i<=1; ++i) {
       u[i+1] = jump_u[i];
       du[i+1] = u[i+1] - state[i+2];
       up[i+1] = state[i+4];
       upo[i+1] = up[i+1];
       Sig_I[i+1] = state[i+10];
       KE(i+1,i+1) = fE_t;
       KE(i+2,i+2) = fE_n;
       KE(i+1, i+2) = 0.;
       KE(i+2,i+1) = 0.;
       
    }

      KE.Inverse(KE_Inv);
    for (i = 1; i<=4; ++i) {
        qn[i] = state[i+5];
        qo[i] = qn[i];
        I_mat[i] = 1.;
    }
     
    for (i = 1; i<=2; ++i) {
       for (j = 1; j<=2; ++i) {   
         Sig[i] = (KE(i,j))*(u[j] - up[j]) + Sig_I[i];
       } 
    }
    
    int iplastic;
    dlam = 0.; dlam2 = 0.; normr = 0.;
/* Check the yield function */
    ff = Sig[1]*Sig[1] - (qn[2] - Sig[2]*qn[3])*(qn[2] - Sig[2]*qn[3]) + (qn[2] - qn[1]*qn[3])*(qn[2] - qn[1]*qn[3]);
    if (ff <0.) {
      int iplastic = 0;
      state[1] = ff;
      normr = 0.;
      state[2] = normr;
    }  
    else {
      kk = 0;
      iplastic = 1;
      while (ff > fTol_1 | normr > fTol_2)
        if (kk > 500) {
        	ExceptionT::GeneralFail("MR2DT::Traction","Too Many Iterations");
        }
        for (i = 0; i<=1; ++i) {
          for (j = 1; j<=2; ++i) { 
           Sig[i] = KE(i,j)*(u[j] - up[j]) + Sig_I[i];
          }
        }
        ff = Sig[1]*Sig[1] - (qn[2] - Sig[2]*qn[3])*(qn[2] - Sig[2]*qn[3]) + (qn[2] - qn[1]*qn[3])*(qn[2] - qn[1]*qn[3]);
        dQdSig[1] = 2.*Sig[1];
        dQdSig[2] = 2.*qn[4]*qn[2] - Sig[2]*qn[4];
        qbar_f(Sig, qn, qbar);
        for (i = 1; i<=2; ++i) {
          R[i] = -up[1] + upo[i] + dlam*dQdSig[i];
        }
        for (i = 1; i<=4; ++i) {
          R[i+2] = -qn[i] + qo[i] + dlam*qbar[i];
        }
        normr = R.Magnitude();
        dQdSig2_f(qn,dQdSig2);
        dQdSigdq_f(Sig, qn, A_uq);
        dqbardSig_f(Sig, qn, A_qu);
        dqbardq_f(Sig, qn, A_qq);
        for (i = 1; i<=6; ++i) {
          for (j = 0; j<=6; ++i) {
            if (i<=2 & j<=2)
             AA_inv(i,j) = KE_Inv(i,j) + dlam*dQdSig2(i,j);
            if (i<=2 & j>2)
             AA_inv(i,j) = A_qu(i,j-2);
            if(i>2 & j<=2)
             AA_inv(i,j) = A_uq(i-2,j);
            if(i>2 & j >2)
             AA_inv(i,j) = -I_mat(i-2,j-2) + dlam*A_qq(i-2,j-2);
           }
        }
        AA_inv.Inverse(AA);
        dfdSig_f(Sig, qn, dfdSig);
        V_sig = dfdSig;
        dfdq_f(Sig, qn, dfdq);
        V_q = dfdq;
        for (i = 1; i<=6; ++i) {
            if (i<=2)
             Rvec[i] = V_sig[i];
             Cvec[i] = dQdSig[i];
            if (i > 2)
             Rvec[i] = V_q[i-2];
             Cvec[i] = qbar[i-2];
        }
        dArrayT tmpVec(6);
        AA.MultTx(R,tmpVec);
        double topp = ff - Rvec.Dot(Rvec,tmpVec);        
        AA.MultTx(Cvec,tmpVec);
        double bott = Rvec.Dot(Rvec,tmpVec); 		
        double dlam2 = topp/bott;
        for (i = 1; i<=6; ++i) {
          for (j = 0; j<=6; ++i) {
            if (i<=2 & j<=2)
             CMAT(i,j) = KE_Inv(i,j);
            if (i<=2 & j>2)
             CMAT(i,j) = ZMAT(i,j-2);
            if(i>2 & j<=2)
             CMAT(i,j) = ZMATP(i-2,j);
            if(i>2 & j >2)
             CMAT(i,j) = -I_mat(i-2,j-2);
           }
        }
         for (i = 1; i<=6; ++i) {
            if (i<=2)
             Rmod[i] = dQdSig[i];
            if (i > 2)
             Rmod[i] = qbar[i-2];
        }
        Rmod *= dlam2;
        R2 = R;
        R2 += Rmod;
        AA.MultTx(R2,X);
        CMAT.MultTx(X,X);
        for (i = 1; i<=6; ++i) {
            if (i<=2)
             dup[i] = X[i];
            if (i > 2)
             dq[i] = X[i];
        }
        up += dup;
        qn += dq;
        dlam = dlam + dlam2;
        kk = kk + 1;
    }
    state[0] = Sig[1];
    state[1] = Sig[2];     
	fTraction[0] = state[0];
	fTraction[1] = state[1];
	state[2] = jump_u[0];
	state[3] = jump_u[1];
	state[4] = up[1];
	state[5] = up[2];
	state[6] = qn[1];
	state[7] = qn[2];
	state[8] = qn[3];
	state[9] = qn[4];
	state[10] = ff;
	state[11] = dlam;
	state[12] = double(iplastic);

	return fTraction;
}

/* calculation of qbar_f */

dArrayT& MR2DT::qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar)
{
 
   dMatrixT AMAT(4,2); dMatrixT BMAT(2,2); dArrayT DQDSIG(2);
   double A1 = -falpha_chi*(qn[1] - fchi_r);
   double B1 = (Sig[2]+abs(Sig[2]))/2./fGf_I;
   double B2 = Sig[1]/fGf_I;
   double DQDN = 2.*qn[4]*(qn[2] - Sig[2]*qn[4]);
   double DQDT = 2.*Sig[2];
   double A2 = -falpha_c*(qn[2] - fc_r);
   double TNA = (Sig[2]-abs(Sig[2]))/2.;
   double B3 = (Sig[2] - abs(TNA*qn[3])*sign(Sig[1]))/fGf_II;
   double A3 = -falpha_phi*(qn[3] - tan(fphi_r));
   double A4 = -falpha_psi*qn[4];
   double DB3_DTn = -qn[3]*sign(Sig[1])*sign(TNA)*(1. - sign(Sig[2]))/fGf_II/2.;
   double DB3_DTt = 1./fGf_II;
   double DB3_DTanphi = -abs(TNA)*sign(Sig[1])/fGf_II;
   double DQDN2 = -2.*qn[4]*qn[4];
   double DQDT2 = 2.;
   double DQDTN = 0.;
   double DQDNT = 0.;
   double DB1DN = B1;
   DQDSIG[1] = DQDT;
   DQDSIG[2] = DQDN;
   AMAT(1,1) = A1;
   AMAT(1,2) = 0.;
   AMAT(2,1) = 0.;
   AMAT(2,2) = A2;
   AMAT(3,1) = 0.;
   AMAT(3,2) = A3;
   AMAT(4,1) = 0.;
   AMAT(4,2) = A4;
   BMAT(1,1) = B1;
   BMAT(1,2) = B2;
   BMAT(2,1) = 0.;
   BMAT(2,2) = B3;
       
   AMAT.MultAB(AMAT, BMAT);
   AMAT.MultTx(DQDSIG,qbar);
   return qbar;
 }


/* calculation of dQdSig2_f */

dMatrixT& MR2DT::dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2)
{
  dQdSig2(1,1) = 2.;
  dQdSig2(2,2) = -2.*qn[4]*qn[4];
  dQdSig2(1,2) = 0.;
  dQdSig2(2,1) = 0.;
  
  return dQdSig2;
}

/* calculation of dfdSig_f */

dArrayT& MR2DT::dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
  dfdSig[1] = 2.*Sig[2];
  dfdSig[2] = 2.*qn[3]*qn[2] - Sig[2]*qn[3];
  
  return dfdSig;
}


/* calculation of dfdq_f */

dArrayT& MR2DT::dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
  dfdq[1] = -2.*qn[3]*(qn[2]-qn[1]*qn[3]);
  dfdq[2] = 2.*(Sig[2] - qn[1]*qn[3]);
  dfdq[3] = 2.*Sig[2]*(qn[2] - Sig[2]*qn[3]) - 2.*qn[1]*(qn[2]-qn[1]*qn[3]);
  dfdq[4] = 0.;
  
  return dfdq;
}

/* calculation of dqdSigdq_f */

dMatrixT& MR2DT::dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
  dQdSigdq(1,1) = 0.;
  dQdSigdq(1,2) = 0.;
  dQdSigdq(1,3) = 0.;
  dQdSigdq(1,4) = 0.;
  dQdSigdq(2,1) = 0.;
  dQdSigdq(2,2) = 2.*qn[4];
  dQdSigdq(2,3) = 0.;
  dQdSigdq(2,4) = 2.*qn[2] - 4.*Sig[2]*qn[4];
  
  
  return dQdSigdq;
}

/* calculation of dqbardSig_f */

dMatrixT& MR2DT::dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{

   double A1 = -falpha_chi*(qn[1] - fchi_r);
   double B1 = (Sig[2]+abs(Sig[2]))/2./fGf_I;
   double B2 = Sig[1]/fGf_I;
   double DQDN = 2.*qn[4]*(qn[2] - Sig[2]*qn[4]);
   double DQDT = 2.*Sig[2];
   double A2 = -falpha_c*(qn[2] - fc_r);
   double TNA = (Sig[2]-abs(Sig[2]))/2.;
   double B3 = (Sig[2] - abs(TNA*qn[3])*sign(Sig[1]))/fGf_II;
   double A3 = -falpha_phi*(qn[3] - tan(fphi_r));
   double A4 = -falpha_psi*qn[4];
   double DB3_DTn = -qn[3]*sign(Sig[1])*sign(TNA)*(1. - sign(Sig[2]))/fGf_II/2.;
   double DB3_DTt = 1./fGf_II;
   double DB3_DTanphi = -abs(TNA)*sign(Sig[1])/fGf_II;
   double DQDN2 = -2.*qn[4]*qn[4];
   double DQDT2 = 2.;
   double DQDTN = 0.;
   double DQDNT = 0.;
   double DB1DN = B1;
   
   dqbardSig(1,1) = A1*B2*DQDT2 + A1*DQDT/fGf_I;
   dqbardSig(1,2) = A1*B1*DQDN2 + A1*DQDN*DB1DN;
   dqbardSig(2,1) = A2*B3*DQDT2 + A2*DQDT*DB3_DTt;
   dqbardSig(2,2) = A2*DQDT*DB3_DTn;
   dqbardSig(3,1) = A3*B3*DQDT2 + A3*DQDT*DB3_DTt;
   dqbardSig(3,2) = A3*DQDT*DB3_DTt;
   dqbardSig(4,1) = A3*B3*DQDT2 + A4*DQDT*DB3_DTt;
   dqbardSig(4,2) = A4*DQDT*DB3_DTn;
 
    return dqbardSig;
}
  
/* calculation of dqbardq_f */

dMatrixT& MR2DT::dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{

   double A1 = -falpha_chi*(qn[1] - fchi_r);
   double B1 = (Sig[2]+abs(Sig[2]))/2./fGf_I;
   double B2 = Sig[1]/fGf_I;
   double DQDN = 2.*qn[4]*(qn[2] - Sig[2]*qn[4]);
   double DQDT = 2.*Sig[2];
   double A2 = -falpha_c*(qn[2] - fc_r);
   double TNA = (Sig[2]-abs(Sig[2]))/2.;
   double B3 = (Sig[2] - abs(TNA*qn[3])*sign(Sig[1]))/fGf_II;
   double A3 = -falpha_phi*(qn[3] - tan(fphi_r));
   double A4 = -falpha_psi*qn[4];
   double DB3_DTn = -qn[3]*sign(Sig[1])*sign(TNA)*(1. - sign(Sig[2]))/fGf_II/2.;
   double DB3_DTt = 1./fGf_II;
   double DB3_DTanphi = -abs(TNA)*sign(Sig[1])/fGf_II;
   double DQDN2 = -2.*qn[4]*qn[4];
   double DQDT2 = 2.;
   double DQDTN = 0.;
   double DQDNT = 0.;
   double DB1DN = B1;
   
   dqbardq(1,1) = -falpha_chi*(B1*DQDN + B2*DQDT);
   dqbardq(1,2) =  A1*B1*(2.*qn[4]);
   dqbardq(1,3) = 0.;
   dqbardq(1,4) =  A1*B1*(2.*qn[2]-4.*Sig[2]*qn[4]);
   
   dqbardq(2,1) = 0.;
   dqbardq(2,2) = -falpha_c*B3*DQDT;
   dqbardq(2,3) = A2*DQDT*DB3_DTanphi;
   dqbardq(2,4) = 0.;
   dqbardq(3,1) = 0.;
   dqbardq(3,2) = 0.;
   dqbardq(3,3) = -falpha_phi*B3*DQDT + A3*DQDT*DB3_DTanphi;
   dqbardq(3,4) = 0.;
   dqbardq(4,1) = 0.;
   dqbardq(4,2) = 0.;
   dqbardq(4,3) = A4*DQDT*DB3_DTanphi;
   dqbardq(4,4) = -falpha_psi*B3*DQDT;
   
    return dqbardq;
}


/* elastoplastic consistent tangent operator*/
const dMatrixT& MR2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

int i, j;

dMatrixT AA(6,6), I_mat(4,4), CMAT(4,4),AA_inv(6,6), 
         A_qq(4,4), A_uu(2,2), A_uq(2,4), A_qu(4,2), ZMAT(2,4),
         ZMATP(4,2), dQdSig2(2,2), dqdbar(4,4), KP(2,2), KEP(2,2);
dMatrixT I_m(2,2), Rmat(2,2), R_Inv(2,2), KE(2,2), KE_Inv(2,2), 
         KEA(2,2), KEA_Inv(2,2);
dArrayT  u(2), up(2), du(2), dup(2), qn(4), qo(4), Rvec(6),Cvec(6),
         R(6), Rmod(6), Sig(2), Sig_I(2), dQdSig(2), dfdq(4), qbar(4),
         R2(6), X(6), V_sig(2), V_q(4), dfdSig(2), K1(2), K2(2);
	
	fStiffness[1] = fStiffness[2] = 0.;
	I_m(1,1) = 1.;
	I_m(2,2) = 2.;
	Sig[1] = sigma[0];
	Sig[2] = sigma[1];
	KE(1,1) = fE_t;
    KE(2,2) = fE_n;
    KE_Inv.Inverse(KE);
    for (i = 1; i<=4; ++i) {
      qn[i] = state[i+5];
      I_mat(i,i) = 1.;
    }

	if (state[12] < 1) 
	{
	    fStiffness[3] = fE_n;
	    fStiffness[0] = fE_t;
	}
	else 
	  	if (state[12] > 0) 
	  	{
	  	    dQdSig2_f(qn,dQdSig2);
	  	    Rmat = dQdSig2;
	  	    Rmat *= state[11];
	   		Rmat += I_mat;
	   		R_Inv.Inverse(Rmat);
	   		KEA.MultAB(R_Inv, KE);
	   		KEA_Inv.Inverse(KEA);
            for (i = 1; i<=6; ++i) {
             for (j = 0; j<=6; ++i) {
               if (i<=2 & j<=2)
                 AA_inv(i,j) = KE_Inv(i,j);
               if (i<=2 & j>2)
                 AA_inv(i,j) = 0.;
               if(i>2 & j<=2)
                 AA_inv(i,j) = 0.;
               if(i>2 & j >2)
                 AA_inv(i,j) = -I_mat(i-2,j-2);
             }
            }
            AA_inv.Inverse(AA);
            dfdSig_f(Sig, qn, dfdSig);
            V_sig = dfdSig;
            dfdq_f(Sig,qn, dfdq);
            V_q = dfdq;
            dQdSig[1] = 2.*Sig[1];
            dQdSig[2] = 2.*qn[4]*qn[2] - Sig[2]*qn[4];
            for (i = 1; i<=6; ++i) {
              if (i<=2)
                Rvec[i] = V_sig[i];
                Cvec[i] = dQdSig[i];
              if (i > 2)
                Rvec[i] = V_q[i-2];
                Cvec[i] = qbar[i-2];
            }
            dArrayT tmpVec(6);
            AA.MultTx(Cvec,tmpVec);
            double bott = Rvec.Dot(Rvec,tmpVec);
            KEA.MultTx(dQdSig, K1);
            KEA.MultTx(dfdSig, K2);
            for (i = 1; i<=2; ++i) {
	   		    for (j = 1; j<=2; ++i) {
	   		      KP(i,j) = K1[i]*K2[j];
                }
	        } 
            KP /=bott;
            KEP = KEA;
            KEP -= KP;
	   		for (i = 1; i<=2; ++i) {
	   		    for (j = 1; j<=2; ++i) {
	   		      fStiffness[0] = KEP(1,1);
	   		      fStiffness[1] = KEP(1,2);
	   		      fStiffness[2] = KEP(2,1);
	   		      fStiffness[3] = KEP(2,2);
                }
	        }
	       }
	      

	return fStiffness;

}


/* print parameters to the output stream */
void MR2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Elastic tangential stiffness. . . . . . . . . . = " << fE_t << '\n';
	out << " Elastic Normal stiffness . . .. . . . . . . . . = " << fE_n     << '\n';
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
int MR2DT::NumOutputVariables(void) const { return 6; }

void MR2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(6);
	labels[0] = "up_t";
	labels[1] = "up_n";
	labels[2] = "Chi";
	labels[3] = "Cohesion";
	labels[4] = "Friction Angle";
	labels[5] = "Dilation Angle";
}

void MR2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	
	output[0] = (jump_u[0]-state[2]);
	output[1] = (jump_u[1]-state[3]);;
	output[2] = state[6];
	output[3] = state[7];
	output[4] = state[8];
	output[5] = state[9];
	
}


bool MR2DT::NeedsNodalInfo(void) { return false; }

int MR2DT::NodalQuantityNeeded(void) 
{ 
        return 2; 
}



void MR2DT::SetElementGroupsNeeded(iArrayT& iGroups) 
{	
	iGroups[0] = 1;
}

