/* $Id: MRSSNLHardT.cpp,v 1.8 2005-10-31 18:03:12 kyonten Exp $ */
/* created: Majid T. Manari (04/16/2003)              */

/* Interface for a nonassociative, small strain,      */
/* pressure dependent plasticity model with nonlinear */ 
/* isotropic hardening/softening.                     */

#include "MRSSNLHardT.h"
#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

/* class constants */

using namespace Tahoe;

const int    kNumInternal = 28; // number of internal state variables
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;

/* constructor */
MRSSNLHardT::MRSSNLHardT(int num_ip, double mu, double lambda):
	fNumIP(num_ip),
	fmu(mu),
	flambda(lambda),
	fkappa(flambda + (2.0/3.0*fmu)),
	fMeanStress(0.0)
{
	SetName("MR_SS_nonlinear_hardening");
}

const dSymMatrixT& MRSSNLHardT::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* remove plastic strain */
	if (element.IsAllocated()) 
	{
		/* load internal variables */
		LoadData(element, ip);

		/* compute elastic strain */
		/*fElasticStrain.DiffOf(totalstrain, fPlasticStrain);*/
		fElasticStrain = totalstrain;
	
		return fElasticStrain;
	}	
	/* no plastic strain */
	else	
		return totalstrain;
}

/* return correction to stress vector computed by mapping the
 * stress back to the yield surface, if needed */
const dSymMatrixT& MRSSNLHardT::StressCorrection(
      const dSymMatrixT& trialstrain, ElementCardT& element, int ip)
{

  	int i; int j; int kk, PLI;

    dMatrixT AA(10,10); dMatrixT KE(6,6); dMatrixT KE_Inv(6,6); dMatrixT I_mat(4,4); 
    dMatrixT CMAT(10,10); dMatrixT A_qq(4,4); dMatrixT A_uu(6,6); dMatrixT A_uq(6,4);
    dMatrixT A_qu(4,6); dMatrixT ZMAT(6,4); dMatrixT ZMATP(4,6);
    dMatrixT dQdSig2(6,6); dMatrixT dqbardq(4,4); dMatrixT dQdSigdq(6,4);
    dMatrixT dqbardSig(4,6); dMatrixT AA_inv(10,10);

    dArrayT u(6); dArrayT up(6); dArrayT du(6); dArrayT dup(6); dArrayT qn(4);
    dArrayT qo(4); dArrayT Rvec(10); dArrayT Cvec(10); dArrayT upo(6);
    dArrayT R(10); dArrayT Rmod(10); dArrayT Sig(6); dArrayT Sig_I(6);
    dArrayT dQdSig(6); dArrayT dfdq(4); dArrayT qbar(4);
    dArrayT R2(10); dMatrixT X(10,1); dArrayT V_sig(6); dArrayT V_q(4); 
    dArrayT dfdSig(6); dArrayT dq(4); dArrayT Y(10); dArrayT state(28);
    dArrayT Sig_trial(6);

    double ff; double bott; double topp; double dlam; double dlam2; double normr;
    
    KE = 0.;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
    I_mat = 0.;
    ZMAT = 0.; ZMATP = 0.;
      
    u[0] = trialstrain(0,0);
	u[1] = trialstrain(1,1);
	u[2] = trialstrain(2,2);
	u[3] = trialstrain(1,2);
	u[4] = trialstrain(0,2);
	u[5] = trialstrain(0,1);
	
	PLI = PlasticLoading(trialstrain, element, ip);
	
	if (PlasticLoading(trialstrain, element, ip) && 
	    element.IsAllocated()) {
	  LoadData(element, ip);
	  for (i = 0; i<=27; ++i) {
       state[i] = fInternal[i];   
      }
	}
    
    if (!PlasticLoading(trialstrain, element, ip) && 
	    element.IsAllocated())
	{
		/* initialize element data */
		double enp  = 0.;
        double esp  = 0.;
        double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
        double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
        double ftan_psi = (tan(fphi_p))*exp(-falpha_psi*esp);
        state = 0.;
        state[18] = fchi;
        state[19] = fc;
        state[20] = ftan_phi;
        state[21] = ftan_psi;
	}
	
	if (!PlasticLoading(trialstrain, element, ip) && 
	    !element.IsAllocated())
	{
		/* initialize element data */
		double enp  = 0.;
        double esp  = 0.;
        fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
        double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
        double ftan_psi = (tan(fphi_p))*exp(-falpha_psi*esp);
        state = 0.;
        state[18] = fchi;
        state[19] = fc;
        state[20] = ftan_phi;
        state[21] = ftan_psi;
	}
	
	/* initialize in the case of first plastic loading*/
	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, element, ip) && 
	    !element.IsAllocated())
	{
		/* new plastic element */
		AllocateElement(element); 
		PlasticLoading(trialstrain, element, ip); 
		/* initialize element data */
		double enp  = 0.;
        double esp  = 0.;
        fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
        double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
        double ftan_psi = (tan(fphi_p))*exp(-falpha_psi*esp);
        state = 0.;
        state[18] = fchi;
        state[19] = fc;
        state[20] = ftan_phi;
        state[21] = ftan_psi;
	}
	
	/* Calculate incremental strains and initialize the neecessary vectors */
    for (i = 0; i<=5; ++i) {
       du[i] = u[i] - state[i+6];
       up[i] = state[i+12];
       upo[i] = up[i];
       Sig_I[i] = 0.;
    }
    
    KE_Inv.Inverse(KE);
    
    for (i = 0; i<=3; ++i) {
        qn[i] = state[i+18];
        qo[i] = qn[i];
        I_mat(i,i) = 1.;
    }
     
    Sig = Sig_I;
    dArrayT ue(6), Sig_e(6);
    ue = u;
    ue -= up;
    KE.MultTx(ue,Sig_e);
    Sig +=Sig_e;
    Sig_trial = Sig;
    
    int iplastic;
    dlam = 0.; dlam2 = 0.; normr = 0.;
    
/* Check the yield function */
     
    Yield_f(Sig, qn, ff);
    if (ff <kYieldTol) {
      iplastic = 0;
      state[22] = ff;
      state[27] = ff;
      normr = 0.;
      state[25] = normr;
      kk = 0.;
    }
      
    else {
      state[27] = ff;
      kk = 0;
      iplastic = 1;
    
      while (ff > fTol_1 | normr > fTol_2) {
        if (kk > 500) {
        	ExceptionT::GeneralFail("MRSSNLHardT::StressCorrection","Too Many Iterations");
        }
        
        Sig = Sig_I;
        ue = u;
        ue -= up;
        KE.Multx(ue,Sig_e);
        Sig +=Sig_e;
        
        Yield_f(Sig, qn, ff);
        dQdSig_f(Sig, qn, dQdSig);
        
        qbar_f(Sig, qn, qbar);
        for (i = 0; i<=5; ++i) {
          R[i] = upo[i];
          R[i] -=up[i];
          R[i] +=dlam*dQdSig[i];
        }
        for (i = 0; i<=3; ++i) {
          R[i+6] = qo[i];
          R[i+6] -=qn[i];
          R[i+6] +=dlam*qbar[i];
        }
                
        normr = R.Magnitude();
        dQdSig2_f(qn,dQdSig2);
        dQdSigdq_f(Sig, qn, A_uq);
        dqbardSig_f(Sig, qn, A_qu);
        dqbardq_f(Sig, qn, A_qq);
       
        for (i = 0; i<=9; ++i) {
          for (j = 0; j<=9; ++j) {
            if (i<=5 & j<=5){
             AA_inv(i,j)  = KE_Inv(i,j);
             AA_inv(i,j) += dlam*dQdSig2(i,j);
            }
            if (i<=5 & j>5){
             AA_inv(i,j) = A_uq(i,j-6);
             AA_inv(i,j) *= dlam;
            } 
            if(i>5 & j<=5){
             AA_inv(i,j) = A_qu(i-6,j);
             AA_inv(i,j) *= dlam;
            } 
            if(i>5 & j >5) {
             AA_inv(i,j)  = I_mat(i-6,j-6);
             AA_inv(i,j)  *= -1.; 
             AA_inv(i,j) += dlam*A_qq(i-6,j-6);
            } 
          }
        }
        AA.Inverse(AA_inv);
        dfdSig_f(Sig, qn, dfdSig);
        V_sig = dfdSig;
        dfdq_f(Sig, qn, dfdq);
        V_q = dfdq;
        for (i = 0; i<=9; ++i) {
            if (i<=5){
             Rvec[i] = V_sig[i];
             Cvec[i] = dQdSig[i];
            }
            if (i > 5){
             Rvec[i] = V_q[i-6];
             Cvec[i] = qbar[i-6];
            }
        }
        dArrayT tmpVec(10);
        AA.Multx(R,tmpVec);
        topp = ff;
        topp -= dArrayT::Dot(Rvec,tmpVec);
        AA.Multx(Cvec,tmpVec);
        bott = dArrayT::Dot(Rvec,tmpVec); 		
        dlam2 = topp/bott;
        for (i = 0; i<=9; ++i) {
          for (j = 0; j<=9; ++j) {
            if (i<=5 & j<=5){
             CMAT(i,j) = KE_Inv(i,j);
            }
            if (i<=5 & j>5) {
             CMAT(i,j) = ZMAT(i,j-6);
            }
            if(i>5 & j<=5) {
             CMAT(i,j) = ZMATP(i-6,j);
            }
            if(i>5 & j >5) {
             CMAT(i,j) = -I_mat(i-6,j-6);
            }
          }
        }
        for (i = 0; i<=9; ++i) {
            if (i<=5){
             Rmod[i] = dQdSig[i];
            }
            if (i >5){
             Rmod[i] = qbar[i-6];
            }
        }
        Rmod *= dlam2;
        R2 = R;
        R2 += Rmod;
        AA.Multx(R2,X);
        CMAT.Multx(X,Y);
        for (i = 0; i<=9; ++i) {
            if (i<=5) {
             dup[i] = Y[i];
            }
            if (i > 5) {
             dq[i-6] = Y[i];
            }
        }
        up += dup;
        qn += dq;
        dlam = dlam + dlam2;
        kk = kk + 1;
      }
    }
    state[0] = Sig[0];
    state[1] = Sig[1];
    state[2] = Sig[2];
    state[3] = Sig[3];
    state[4] = Sig[4];
    state[5] = Sig[5];     
	state[6] = trialstrain(0,0);
	state[7] = trialstrain(1,1);
	state[8] = trialstrain(2,2);
	state[9] = trialstrain(1,2);
	state[10] = trialstrain(0,2);
	state[11] = trialstrain(0,1);
	state[12] = up[0];
	state[13] = up[1];
	state[14] = up[2];
	state[15] = up[3];
	state[16] = up[4];
	state[17] = up[5];
	state[18] = qn[0];
	state[19] = qn[1];
	state[20] = qn[2];
	state[21] = qn[3];
	state[22] = ff;
	state[23] = dlam;
	state[24] = double(iplastic);
	state[25] = normr;
	state[26] = double(kk);
	for (i = 0; i<=5; ++i) {
	      fStressCorr[i] = state[i];
    }
	if (iplastic>0) {
	   for (i =0; i<=27; ++i) {
		  fInternal[i] = state[i];
	   }
	   for (i = 0; i<=5; ++i) {
	      fPlasticStrain[i] = state[i+12];
       }
	}		
 return fStressCorr;
}

/*
 * Returns the value of the yield function given the
 * stress vector and state variables, where alpha
 * represents isotropic hardening.
 */
double& MRSSNLHardT::Yield_f(const dArrayT& Sig, 
			const dArrayT& qn, double& ff)
{
  double kTemp1, kTemp2, kTemp3, kTemp4;
  double fc, fchi, ffriction, fpress;
  dMatrixT devstress(3,3);
  
  fpress  = Sig[0]+Sig[1]+Sig[2];
  fpress /=3.;
  devstress(0,0) = Sig[0] - fpress;
  devstress(1,1) = Sig[1] - fpress;
  devstress(2,2) = Sig[2] - fpress;
  devstress(1,2) = Sig[3];
  devstress(0,2) = Sig[4];
  devstress(0,1) = Sig[5];
  devstress(2,1) = Sig[3];
  devstress(2,0) = Sig[4];
  devstress(1,0) = Sig[5];

  fc = qn[1];
  ffriction = qn[2];
  fchi = qn[0];
  ff = dMatrixT::Dot(devstress,devstress);
  ff /= 2.;
  kTemp2  = (fc - ffriction*fpress);
  kTemp1  = kTemp2;
  kTemp1 *= kTemp2;
  ff  -= kTemp1;
  kTemp3  = (fc - ffriction*fchi);
  kTemp4  = kTemp3;
  kTemp4 *= kTemp3;
  ff  += kTemp4;
  return  ff;
}


/* calculation of qbar_f */

dArrayT& MRSSNLHardT::qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, B2dQdS, B3dQdS;
   dMatrixT Sig_Dev(3,3), B2, B3, dQdS;
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   Sig_Dev(0,0) = Sig[0] - Sig_p;
   Sig_Dev(1,1) = Sig[1] - Sig_p;
   Sig_Dev(2,2) = Sig[2] - Sig_p;
   Sig_Dev(1,2) = Sig[3];
   Sig_Dev(0,2) = Sig[4];
   Sig_Dev(0,1) = Sig[5];
   Sig_Dev(2,1) = Sig[3];
   Sig_Dev(2,0) = Sig[4];
   Sig_Dev(1,0) = Sig[5];
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   dQdS = Sig_Dev;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   B2dQdS = dMatrixT::Dot(B2,dQdS);
   B3dQdS = dMatrixT::Dot(B3,dQdS);
      
   qbar[0]  = A1*B1*dQdP; 
   qbar[0] += A1*B2dQdS;
   qbar[1]  = B3dQdS;
   qbar[1]  *=A2;
   qbar[2]  = B3dQdS;
   qbar[2]  *=A3;
   qbar[3]  = B3dQdS;
   qbar[3]  *=A4;
   return qbar;
 }


/* calculation of dQdSig2_f */

dMatrixT& MRSSNLHardT::dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2)
{
  int i, j;
  double Fac;
  dMatrixT II_mat(6,6), I_mat(6,6);
  
  II_mat = 0.;
  II_mat(0,0) = II_mat(1,1) = II_mat(2,2) = 1.; 
  II_mat(3,3) = II_mat(4,4) = II_mat(5,5) = 1.;
  
  I_mat = 0.;
  for (i = 0; i<=2; ++i) {
     for (j = 0; j<=2; ++j) {
      I_mat(i,j) = 1.;
     }
  }
  
  Fac = 2./3.;
  Fac *= qn[3]*qn[3];
  Fac += 1.;
  Fac /= 3.;
  I_mat *= Fac;
  
  
  dQdSig2  = II_mat;
  dQdSig2 -= I_mat;
  
  return dQdSig2;
}

/* calculation of dfdSig_f */

dArrayT& MRSSNLHardT::dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
  double Sig_p, temp;
  int i; 
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   dfdSig[0] = Sig[0] - Sig_p;
   dfdSig[1] = Sig[1] - Sig_p;
   dfdSig[2] = Sig[2] - Sig_p;
   dfdSig[3] = Sig[3];
   dfdSig[4] = Sig[4];
   dfdSig[5] = Sig[5];

   temp = 2./3.;
   temp  *= qn[2];
   temp *= (qn[1] - Sig_p*qn[2]);
   for (i = 0; i<=2; ++i) {
      dfdSig[i] += temp;
   }
  
  return dfdSig;
}

/* calculation of dQdSig_f */

dArrayT& MRSSNLHardT::dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
  double Sig_p, temp;
  int i;
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   dQdSig[0] = Sig[0] - Sig_p;
   dQdSig[1] = Sig[1] - Sig_p;
   dQdSig[2] = Sig[2] - Sig_p;
   dQdSig[3] = Sig[3];
   dQdSig[4] = Sig[4];
   dQdSig[5] = Sig[5];

   temp = 2./3.;
   temp  *= qn[3];
   temp *= (qn[1] - Sig_p*qn[3]);
   for (i = 0; i<=2; ++i) {
      dQdSig[i] += temp;
   }
  
  return dQdSig;
}


/* calculation of dfdq_f */

dArrayT& MRSSNLHardT::dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
  double Sig_p;
  Sig_p = (Sig[0] + Sig[1] + Sig[2])/3.;
  dfdq[0] = -2.*qn[2]*(qn[1]-qn[0]*qn[2]);
  dfdq[1] = 2.*(Sig_p - qn[0])*qn[2];
  dfdq[2] = 2.*Sig_p*(qn[1] - Sig_p*qn[2]) - 2.*qn[0]*(qn[1]-qn[0]*qn[2]);
  dfdq[3] = 0.;
  
  return dfdq;
}

/* calculation of dQdSigdq_f */

dMatrixT& MRSSNLHardT::dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
  double Sig_p;
  Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
  dQdSigdq = 0.;
  dQdSigdq(0,1) = 2.*qn[3]/3.;
  dQdSigdq(1,1) = 2.*qn[3]/3.;
  dQdSigdq(2,1) = 2.*qn[3]/3.; 
  dQdSigdq(0,3) = (2.*qn[1] - 4.*Sig_p*qn[3])/3.;
  dQdSigdq(1,3) = (2.*qn[1] - 4.*Sig_p*qn[3])/3.;
  dQdSigdq(2,3) = (2.*qn[1] - 4.*Sig_p*qn[3])/3.;
    
  return dQdSigdq;
}

/* calculation of dqbardSig_f */

dMatrixT& MRSSNLHardT::dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, d2QdP2, dB1dP, SN;
   dMatrixT dhchi_dSig(3,3), dhc_dSig(3,3), dhtanphi_dSig(3,3), dhtanpsi_dSig(3,3);
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3), dB2dS_dQdS(3,3), dB3dS_dQdS(3,3);
   dMatrixT I_mat(3,3); dMatrixT tempmat(3,3);
   
   I_mat = 0.;
   I_mat(0,0) = I_mat(1,1) = I_mat(2,2) = 1.;
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   Sig_Dev(0,0) = Sig[0] - Sig_p;
   Sig_Dev(1,1) = Sig[1] - Sig_p;
   Sig_Dev(2,2) = Sig[2] - Sig_p;
   Sig_Dev(1,2) = Sig[3];
   Sig_Dev(0,2) = Sig[4];
   Sig_Dev(0,1) = Sig[5];
   Sig_Dev(2,1) = Sig[3];
   Sig_Dev(2,0) = Sig[4];
   Sig_Dev(1,0) = Sig[5];
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p + fabs(Sig_p))/2./fGf_I;
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   dQdS = Sig_Dev;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   
   d2QdP2      =  -2.*qn[3]*qn[3];
   dB2dS_dQdS  = Sig_Dev;
   dB2dS_dQdS /= fGf_I;
   dB3dS_dQdS  = Sig_Dev;
   dB3dS_dQdS /= fGf_II;
   SN = signof(Sig_p);
   dB1dP = (SN +fabs(SN))/2./fGf_I;
   
   dhchi_dSig  = I_mat;
   dhchi_dSig *= (A1*B1*d2QdP2+A1*dQdP*dB1dP)/3.;
   tempmat =  dB2dS_dQdS; 
   tempmat += B2; 
   tempmat *= A1;
   dhchi_dSig += tempmat;
   dhc_dSig   = B3;
   dhc_dSig += dB3dS_dQdS;
   dhc_dSig  *= A2;
   dhtanphi_dSig  = B3;
   dhtanphi_dSig += dB3dS_dQdS;
   dhtanpsi_dSig *= A3;
   dhtanpsi_dSig  = B3;
   dhtanpsi_dSig += dB3dS_dQdS;
   dhtanpsi_dSig *= A4;
   
   dqbardSig(0,0) = dhchi_dSig(0,0);
   dqbardSig(0,1) = dhchi_dSig(1,1);
   dqbardSig(0,2) = dhchi_dSig(2,2);
   dqbardSig(0,3) = dhchi_dSig(1,2);
   dqbardSig(0,4) = dhchi_dSig(0,2);
   dqbardSig(0,5) = dhchi_dSig(0,1);
   dqbardSig(1,0) = dhc_dSig(0,0);
   dqbardSig(1,1) = dhc_dSig(1,1);
   dqbardSig(1,2) = dhc_dSig(2,2);
   dqbardSig(1,3) = dhc_dSig(1,2);
   dqbardSig(1,4) = dhc_dSig(0,2);
   dqbardSig(1,5) = dhc_dSig(0,1);
   dqbardSig(2,0) = dhtanphi_dSig(0,0);
   dqbardSig(2,1) = dhtanphi_dSig(1,1);
   dqbardSig(2,2) = dhtanphi_dSig(2,2);
   dqbardSig(2,3) = dhtanphi_dSig(1,2);
   dqbardSig(2,4) = dhtanphi_dSig(0,2);
   dqbardSig(2,5) = dhtanphi_dSig(0,1);
   dqbardSig(3,0) = dhtanpsi_dSig(0,0);
   dqbardSig(3,1) = dhtanpsi_dSig(1,1);
   dqbardSig(3,2) = dhtanpsi_dSig(2,2);
   dqbardSig(3,3) = dhtanpsi_dSig(1,2);
   dqbardSig(3,4) = dhtanpsi_dSig(0,2);
   dqbardSig(3,5) = dhtanpsi_dSig(0,1);
   
    return dqbardSig;
}
  
/* calculation of dqbardq_f */

dMatrixT& MRSSNLHardT::dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, B2dQdS, B3dQdS;
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3);
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   Sig_Dev(0,0) = Sig[0] - Sig_p;
   Sig_Dev(1,1) = Sig[1] - Sig_p;
   Sig_Dev(2,2) = Sig[2] - Sig_p;
   Sig_Dev(1,2) = Sig[3];
   Sig_Dev(0,2) = Sig[4];
   Sig_Dev(0,1) = Sig[5];
   Sig_Dev(2,1) = Sig[3];
   Sig_Dev(2,0) = Sig[4];
   Sig_Dev(1,0) = Sig[5];
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   dQdS = Sig_Dev;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   B2dQdS = dMatrixT::Dot(B2,dQdS);
   B3dQdS = dMatrixT::Dot(B3,dQdS);
   
   dqbardq(0,0) = -falpha_chi*(B1*dQdP + B2dQdS);
   dqbardq(0,1) =  A1*B1*(2.*qn[3]);
   dqbardq(0,2) = 0.;
   dqbardq(0,3) =  A1*B1*(2.*qn[1]-4.*Sig_p*qn[3]);   
   dqbardq(1,0) = 0.;
   dqbardq(1,1) = -falpha_c*B3dQdS;
   dqbardq(1,2) = 0.;
   dqbardq(1,3) = 0.;
   dqbardq(2,0) = 0.;
   dqbardq(2,1) = 0.;
   dqbardq(2,2) = -falpha_phi*B3dQdS;
   dqbardq(2,3) = 0.;
   dqbardq(3,0) = 0.;
   dqbardq(3,1) = 0.;
   dqbardq(3,2) = 0.;
   dqbardq(3,3) = -falpha_psi*B3dQdS;
   
    return dqbardq;
}

/* return the consistent elstoplastic moduli 
 *
 * Note: Return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& MRSSNLHardT::Moduli(const ElementCardT& element, 
	int ip)
{
	    
	 int i; int j;
	 double bott, dlam;
     dMatrixT AA(10,10); dMatrixT KE(6,6); dMatrixT KE_Inv(6,6); dMatrixT I_mat(4,4); 
     dMatrixT CMAT(10,10); dMatrixT A_qq(4,4); dMatrixT A_uu(6,6); dMatrixT A_uq(6,4);
     dMatrixT A_qu(4,6); dMatrixT ZMAT(6,4); dMatrixT ZMATP(4,6), I_m(6,6);
     dMatrixT Rmat(6,6), dQdSig2(6,6); dMatrixT dqbardq(4,4); dMatrixT dQdSigdq(6,4);
     dMatrixT dqbardSig(4,6); dMatrixT AA_inv(10,10), R_Inv(6,6), KEA(6,6), KEA_Inv(6,6);
     dMatrixT KP(6,6); dMatrixT KP2(6,6); dMatrixT KEP(6,6); dMatrixT KES(6,6);
     dMatrixT KES_Inv(6,6);
    
     dArrayT state(28);
        
     dMatrixT Ch(4,4), Ch_Inv(4,4), KE1(4,6), KE2(6,6), KE3(6,4);
         
     dArrayT  u(6), up(6), du(6), dup(6), qn(4), qo(4), Rvec(10), Cvec(10),
              R(10), Rmod(10), Sig(6), Sig_I(6), dQdSig(6), dfdq(4), qbar(4),
              R2(10), X(10), V_sig(6), V_q(4), dfdSig(6), K1(6), K2(6);

	KE = 0.;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(2,1) = KE(1,0) = KE(0,1) = KE(2,0) = KE(0,2) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
	if(!element.IsAllocated()) {
	  	fModuli = KE;
	  	return fModuli;
	}
	
    I_mat = 0.; I_m = 0.;
    ZMAT = 0.; ZMATP = 0.;
    
    /* load internal state variables */
    if(element.IsAllocated()) {
	  	LoadData(element,ip);
	  	for (i =0; i<=27; ++i) {
		  state[i] = fInternal[i];
		}
	}
	  	
    for (i = 0; i<=5; ++i) {
       Sig[i] = state[i];
    }
    
    KE_Inv.Inverse(KE);
    
    for (i = 0; i<=3; ++i) {
      qn[i] = state[i+18];
      I_mat(i,i) = 1.;
    }
    for (i = 0; i<=5; ++i) {
      I_m(i,i) = 1.;
    }

	if (state[24] == 0.) 
	{
	    fModuli = KE;
	    fModuli.CopySymmetric();
	}
	else 
	  	if (state[24] == 1.) 
	  	{
	  	    dlam = state[23];
	  	    dQdSig2_f(qn, dQdSig2);
	        dqbardSig_f(Sig, qn, A_qu);
	        dqbardq_f(Sig, qn, A_qq);
	        dQdSigdq_f(Sig, qn, A_uq);
	        Ch  = A_qq;
	        Ch *= -dlam;
	        Ch += I_mat;
	        Ch_Inv.Inverse(Ch);
	        KE1.MultAB(Ch_Inv,A_qu);
	        KES.MultAB(A_uq,KE1);
	        KES *= state[23];
	        KES *= state[23];
	        /*KES = 0.;*/
	        KE2 = dQdSig2;
	        KE2 *=state[23];
	        KES += KE2;
	        KES += KE;
	        
	        KES_Inv.Inverse(KES);
	     	
            for (i = 0; i<=9; ++i) {
              for (j = 0; j<=9; ++j) {
                if (i<=5 & j<=5){
                 AA_inv(i,j)  = KE_Inv(i,j);
                 AA_inv(i,j) += dlam*dQdSig2(i,j);
                }
                if (i<=5 & j>5){
                  AA_inv(i,j) = A_uq(i,j-6);
                  AA_inv(i,j) *= dlam;
                } 
                if(i>5 & j<=5){
                  AA_inv(i,j) = A_qu(i-6,j);
                  AA_inv(i,j) *= dlam;
                } 
                if(i>5 & j >5) {
                  AA_inv(i,j)  = I_mat(i-6,j-6);
                  AA_inv(i,j)  *= -1.; 
                  AA_inv(i,j) += dlam*A_qq(i-6,j-6);
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
            for (i = 0; i<=9; ++i) {
              if (i<=5) {
                Rvec[i] = V_sig[i];
                Cvec[i] = dQdSig[i];
              }
              if (i>5) {
                Rvec[i] = V_q[i-6];
                Cvec[i] = qbar[i-6];
              }
            }
            dArrayT tmpVec(10), Vvec(6), dVec(6);
            AA.Multx(Cvec,tmpVec);
            bott = dArrayT::Dot(Rvec,tmpVec);
            
            for (i = 0; i<=5; ++i) {
                  Vvec[i] = 0.;
	   		    for (j = 0; j<=9; ++j) {
	   		      Vvec[i] += Rvec[j]*AA(j,i);
                }
	        }
            
            for (i = 0; i<=5; ++i) {
	   		    for (j = 0; j<=5; ++j) {
	   		      KP(i,j) = dQdSig[i]*Vvec[j];
                }
	        }
            
            KE3.MultAB(A_uq, Ch_Inv);
            KE3.Multx(qbar,dVec);
            for (i = 0; i<=5; ++i) {
	   		    for (j = 0; j<=5; ++j) {
	   		      KP2(i,j) = dVec[i]*Vvec[j];
                }
	        }
	        
	        KP2 *= state[11];
	        KP += KP2;
	        KP /= -bott;
            KP += I_m;
            KEP.MultAB(KES_Inv, KP);
	   		fModuli = KEP;
	   		
	       }
	return fModuli;
}


/* return the correction to modulus Cep~, checking for discontinuous
 *   bifurcation */


const dMatrixT& MRSSNLHardT::ModuliPerfPlas(const ElementCardT& element, 
	int ip)
{
	/* initialize */
	fModuliPerfPlas = 0.0;

	if (element.IsAllocated() && 
	   (element.IntegerData())[ip] == kIsPlastic)
	{

	}

	return fModuliPerfPlas;
}	

 	 	
/* return a pointer to a new plastic element object constructed with
 * the data from element */
void MRSSNLHardT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += kNumInternal*fNumIP;        //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0
}


/* accept parameter list */
void MRSSNLHardT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MRPrimitiveT::TakeParameterList(list);

	/* dimension work space */
	fElasticStrain.Dimension(kNSD);
	fStressCorr.Dimension(kNSD);
	fModuli.Dimension(dSymMatrixT::NumValues(kNSD));
	fModuliPerfPlas.Dimension(dSymMatrixT::NumValues(kNSD));
	fDevStress.Dimension(kNSD);
	fDevStrain.Dimension(kNSD); 
	fTensorTemp.Dimension(dSymMatrixT::NumValues(kNSD));
	IdentityTensor2.Dimension(kNSD);
	One.Dimension(kNSD);
    
	/* initialize constant tensor */
	One.Identity();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* element level data */
void MRSSNLHardT::Update(ElementCardT& element)
{
	/* get flags */
	iArrayT& Flags = element.IntegerData();

	/* check if reset state */
	if (Flags[0] == kReset)
	{
		Flags = kIsElastic; //don't update again
		return; 
	}

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
		if (Flags[ip] == kIsPlastic) /* plastic update */
		{
			/* do not repeat if called again. */
			Flags[ip] = kIsElastic;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* fetch element data */
			LoadData(element, ip);
		}
}

/* resets to the last converged solution */
void MRSSNLHardT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* load element data for the specified integration point */
void MRSSNLHardT::LoadData(const ElementCardT& element, int ip)
{
	/* check */
	if (!element.IsAllocated()) 
	    ExceptionT::GeneralFail("MRSSNLHardT::LoadData","The element should have been allocated");
	/* fetch arrays */
	 const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Alias(        dim, &d_array[           dex]);
	/*fUnitNorm.Set(        kNSD, &d_array[  offset + dex]); */    
	fInternal.Alias(kNumInternal, &d_array[2*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int MRSSNLHardT::PlasticLoading(const dSymMatrixT& trialstrain, 
	  ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated()) 
		return( YieldCondition(DeviatoricStress(trialstrain,element),
			       MeanStress(trialstrain,element)) > kYieldTol );
        /* already plastic */
	else 
	{
	/* get flags */
	 iArrayT& Flags = element.IntegerData();
		
	/* load internal variables */
	LoadData(element, ip);

		/* plastic */
		if (fInternal[kplastic] > 0.5)
		{		
			/* set flag */
			Flags[ip] = kIsPlastic;
	
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
		    Flags[ip] = kIsElastic; //removed to avoid restting 7/01
			
			return 0;
		}
	}
}	

/* Computes the stress corresponding to the given element
 * and elastic strain.  The function returns a reference to the
 * stress in fDevStress */
dSymMatrixT& MRSSNLHardT::DeviatoricStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain);

	/* compute deviatoric elastic stress */
	fDevStress.SetToScaled(2.0*fmu,fDevStrain);

	return fDevStress;
}

/* computes the hydrostatic (mean) stress */
double MRSSNLHardT::MeanStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

  fMeanStress = fkappa*trialstrain.Trace();
  return fMeanStress;
}

double MRSSNLHardT::signof(double& r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}