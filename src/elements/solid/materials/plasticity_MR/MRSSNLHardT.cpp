/* $Id: MRSSNLHardT.cpp,v 1.15 2006-02-09 15:12:19 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003)              */

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

const int    kNumInternal = 30; //28; // number of internal state variables
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
        //fElasticStrain.DiffOf(totalstrain, fPlasticStrain);
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

    int kk, iplastic;
    double ff, bott, topp, dlam, dlam2, normr;

    /* allocate matrices */
    dMatrixT AA(10), AA_inv(10), CMAT(10), KE(6), KE_Inv(6); 
    dMatrixT A_qq(4), A_uu(6), A_uq(6,4), A_qu(4,6);
    dMatrixT dQdSig2(6), dQdSigdq(6,4);
    dMatrixT dqbardq(4), dqbardSig(4,6);
    
    /* allocate reduced index vector of symmetric matrices */
    dSymMatrixT u(3), up(3), du(3), dup(3), upo(3), ue(3);
    dSymMatrixT Sig(3), Sig_I(3), Sig_trial(3), Sig_e(3);
    
    /* allocate vectors */ 
    dArrayT qn(4), qo(4), dq(4), qbar(4);
    dArrayT Rvec(10), Cvec(10), R(10), Rmod(10);
    dArrayT dQdSig(6), dfdSig(6), dfdq(4);
    dArrayT X(10), Y(10), state(30);
    
    /* elastic moduli tensor */
    KE = 0.;
    KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
    KE(1,2) = KE(0,1) = KE(0,2) = flambda;
    KE(2,1) = KE(1,0) = KE(2,0) = flambda;
    KE(5,5) = KE(4,4) = KE(3,3) = fmu;
    
    /* initialize element data */
    state = 0.;
    //double enp  = 0.;
    //double esp  = 0.;
    double enp  = state[28];
    double esp  = state[29];
    double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
    double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
    double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
    double ftan_psi = (tan(fphi_p))*exp(-falpha_psi*esp);
    state[18] = fchi;
    state[19] = fc;
    state[20] = ftan_phi;
    state[21] = ftan_psi;
    
    /* initialize in the case of first plastic loading */
    /* check consistency and initialize plastic element */
    if (PlasticLoading(trialstrain, element, ip) && 
        element.IsAllocated()) {
      LoadData(element, ip);
      
       /* fetch internal variables */
      state.CopyIn(0, fInternal);
    }
    
    if (PlasticLoading(trialstrain, element, ip) && 
        !element.IsAllocated())
    {
        /* new plastic element */
        AllocateElement(element);
        
        /* initialize element data */ 
        PlasticLoading(trialstrain, element, ip); 
    }
    
    /* calculate incremental strains and initialize the neecessary vectors */
    u = trialstrain;	
    for (int i = 0; i < 6; i++) 
       du[i] = u[i] - state[i+6];
       
    up.CopyPart(0, state, 12, up.Length());
    upo = up;
    qn.CopyPart(0, state, 18, qn.Length());
    qo = qn;
    
    KE_Inv.Inverse(KE);
  
    Sig_I = 0.;
    Sig = Sig_I;
    ue.DiffOf(u, up);
    KE.MultTx(ue, Sig_e);
    Sig += Sig_e;
    Sig_trial = Sig;
    
    dlam = 0.; dlam2 = 0.; normr = 0.;
    
/* check the yield function */
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
    
      while (ff > fTol_1 || normr > fTol_2) {
        if (kk > 500)
            ExceptionT::GeneralFail("MRSSNLHardT::StressCorrection","Too Many Iterations");
        
        Sig = Sig_I;
        ue.DiffOf(u, up);
        KE.Multx(ue, Sig_e);
        Sig += Sig_e;
        
        /* check yield condition */
        Yield_f(Sig, qn, ff);
        
        /* residuals for plastic strain and internal variables */
        dQdSig_f(Sig, qn, dQdSig);
        qbar_f(Sig, qn, qbar);
        for (int i = 0; i < 6; i++) {
          R[i]  = upo[i];
          R[i] -= up[i];
          R[i] += dlam*dQdSig[i];
        }
        for (int i = 0; i < 4; i++) {
          R[i+6]  = qo[i];
          R[i+6] -= qn[i];
          R[i+6] += dlam*qbar[i];
        }
        
        /* L2 norm of the residual vector */
        normr = R.Magnitude();
        
        /* form AA_inv matrix */
        dQdSig2_f(qn, dQdSig2);
        dQdSigdq_f(Sig, qn, dQdSigdq);
        dqbardSig_f(Sig, qn, dqbardSig);
        dqbardq_f(Sig, qn, dqbardq);
        A_uu.SetToScaled(dlam, dQdSig2);
        A_uu += KE_Inv;
        A_uq.SetToScaled(dlam, dQdSigdq);
        A_qu.SetToScaled(dlam, dqbardSig);
        A_qq.SetToScaled(dlam, dqbardq);
        A_qq.PlusIdentity(-1.0);
        AA_inv = 0.0;
        AA_inv.AddBlock(0,           0,           A_uu);
        AA_inv.AddBlock(0,           A_uu.Cols(), A_uq);
        AA_inv.AddBlock(A_uu.Rows(), 0,           A_qu);
        AA_inv.AddBlock(A_uu.Rows(), A_uu.Cols(), A_qq);
        
        /* calculate dlam2 */
        AA.Inverse(AA_inv);
        dfdSig_f(Sig, qn, dfdSig);
        dfdq_f(Sig, qn, dfdq);
        Rvec.CopyIn(0, dfdSig);
        Rvec.CopyIn(dfdSig.Length(), dfdq);
        Cvec.CopyIn(0, dQdSig);
        Cvec.CopyIn(dQdSig.Length(), qbar);

        dArrayT tmpVec(10);
        AA.Multx(R,tmpVec);
        topp = ff - dArrayT::Dot(Rvec,tmpVec);
        AA.Multx(Cvec,tmpVec);
        bott = dArrayT::Dot(Rvec,tmpVec);       
        dlam2 = topp/bott;
        if (kk == 0 && dlam2 < 0.0 ) {
        	//cout << "dlam2 =" << dlam2 << endl << endl;
        	//cout << "topp = " << topp << endl << endl;
        	cout << "bott =" << bott << endl << endl;
        	//cout << "tmpVec =" << tmpVec << endl << endl;
        	//cout << "Rvec =" << Rvec << endl << endl;
        	//cout <<"dQdSig =" << dQdSig << endl << endl;
        	//cout <<"qbar =" << qbar << endl << endl;
        }
        
        /* calculate dup and dq */
        dMatrixT I_mat(4,4);
        I_mat = 0.0;
     	CMAT = 0.0; 
        I_mat.PlusIdentity(-1.0);
        CMAT.AddBlock(0, 0, KE_Inv);
        CMAT.AddBlock(KE_Inv.Rows(), KE_Inv.Cols(), I_mat);
        Rmod.CopyIn(0, dQdSig);
        Rmod.CopyIn(dQdSig.Length(), qbar);
        Rmod *= dlam2;
        Rmod += R;
        AA.Multx(Rmod, X);
        CMAT.Multx(X,Y);
        
        /* separate dup and dq */
        dup.CopyPart(0, Y, 0, dup.Length());
        dq.CopyPart(0, Y, dup.Length(), dq.Length());
        
        /* update plastic strains and internal variables */
        up += dup;
        qn += dq;
        dlam += dlam2;
        cout << "k = " << kk << " dlam2 = " << dlam2
             << " dlam = " << dlam << endl << endl;
        kk++;
      } // while (ff > fTol_1 || normr > fTol_2)
    } //  if (ff <kYieldTol)
    
     /* update state variables */
    state.CopyIn(0, Sig);
    state.CopyIn(Sig.Length(), trialstrain);       
    state.CopyIn(12, up);
    state.CopyIn(18, qn);
    state[22] = ff;
    state[23] = dlam;
    state[24] = double(iplastic);
    state[25] = normr;
    state[26] = double(kk);
    
    /* calculate and store enp and esp */
    dQdSig_f(Sig, qn, dQdSig);
    dSymMatrixT Sig_Dev(3), dQdS(3);
    double Sig_p = Sig.Trace()/3.0;
   	double B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   	double dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   	Sig_Dev.Deviatoric(Sig);
   	dQdS = Sig_Dev;
   	state[28] =  dMatrixT::Dot(Sig_Dev,dQdS)/fGf_I;
   	state[28] += B1*dQdP;
   	state[28] *= dlam;
   	state[29] = dMatrixT::Dot(Sig_Dev,dQdS)/fGf_II;
   	state[29] *= dlam;
    
    fStressCorr = Sig;
          
    if (iplastic > 0) {
       fInternal.CopyIn(0, state);
       fPlasticStrain = up;
    }                
 return fStressCorr;
}

/*
 * returns the value of the yield function given the
 * stress vector and state variables, where alpha
 * represents isotropic hardening.
 */
void MRSSNLHardT::Yield_f(const dSymMatrixT& Sig, 
            const dArrayT& qn, double& ff)
{
  dSymMatrixT devstress(3);
  double kTemp1, kTemp2, kTemp3, kTemp4;
  double fchi = qn[0];
  double fc = qn[1];
  double ffriction = qn[2]; 
  double fpress = Sig.Trace()/3.0;
  
  devstress.Deviatoric(Sig);
  ff = (devstress.ScalarProduct())/2.0;
  kTemp2  = (fc - ffriction*fpress);
  kTemp1  = kTemp2;
  kTemp1 *= kTemp2;
  ff  -= kTemp1;
  kTemp3  = (fc - ffriction*fchi);
  kTemp4  = kTemp3;
  kTemp4 *= kTemp3;
  ff  += kTemp4;
}

/* calculation of dfdSig_f */
void MRSSNLHardT::dfdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
   double fc = qn[1];
   double ftan_phi = qn[2]; 
   double Sig_p = Sig.Trace()/3.0;
   dfdSig[0] = Sig[0] - Sig_p;
   dfdSig[1] = Sig[1] - Sig_p;
   dfdSig[2] = Sig[2] - Sig_p;
   dfdSig[3] = Sig[3];
   dfdSig[4] = Sig[4];
   dfdSig[5] = Sig[5];

   double temp = 2./3.;
   temp  *= ftan_phi;
   temp *= (fc - Sig_p*ftan_phi);
   
   for (int i = 0; i < 3; i++) 
      dfdSig[i] += temp;
}

/* calculation of dfdq_f */
void MRSSNLHardT::dfdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
  double fchi = qn[0];
  double fc = qn[1];
  double ftan_phi = qn[2];
  double Sig_p = Sig.Trace()/3.0;
   
  dfdq[0] = -2.*ftan_phi*(fc-fchi*ftan_phi);
  dfdq[1] = 2.*(Sig_p - fchi)*ftan_phi;
  dfdq[2] = 2.*Sig_p*(fc - Sig_p*ftan_phi) - 2.*fchi*(fc-fchi*ftan_phi);
  dfdq[3] = 0.;
}

/* calculation of dQdSig_f */
void MRSSNLHardT::dQdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
   double fc = qn[1];
   double ftan_psi = qn[3];
   double Sig_p = Sig.Trace()/3.0;
   dQdSig[0] = Sig[0] - Sig_p;
   dQdSig[1] = Sig[1] - Sig_p;
   dQdSig[2] = Sig[2] - Sig_p;
   dQdSig[3] = Sig[3];
   dQdSig[4] = Sig[4];
   dQdSig[5] = Sig[5];

   double temp = 2./3.;
   temp  *= ftan_psi;
   temp *= (fc - Sig_p*ftan_psi);
   
   for (int i = 0; i < 3; i++)
      dQdSig[i] += temp;
}

/* calculation of dQdSig2_f */
void MRSSNLHardT::dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2)
{
  double ftan_psi = qn[3];
  dMatrixT I_mat(6);
  
  I_mat = 0.; 
  dQdSig2 = 0.;
  for (int i = 0; i < 3; i++)
  	for (int j = 0; j < 3; j++)
    	I_mat(i,j) = 1.;
  
  double Fac = 2./3.;
  Fac *= ftan_psi*ftan_psi;
  Fac += 1.;
  Fac /= 3.;
  I_mat *= Fac;
  
  dQdSig2.PlusIdentity();
  dQdSig2 -= I_mat;
}

/* calculation of dQdSigdq_f */
void MRSSNLHardT::dQdSigdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
  double fc = qn[1];
  double ftan_psi = qn[3];
  double Sig_p = Sig.Trace()/3.0;
  
  dQdSigdq = 0.;
  dQdSigdq(0,1) = 2.*ftan_psi/3.;
  dQdSigdq(1,1) = dQdSigdq(2,1) = dQdSigdq(0,1);
  dQdSigdq(0,3) = (2.*fc - 4.*Sig_p*ftan_psi)/3.;
  dQdSigdq(1,3) = dQdSigdq(2,3) = dQdSigdq(0,3); 
}

/* calculation of qbar_f */
void MRSSNLHardT::qbar_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& qbar)
{
   dSymMatrixT Sig_Dev(3), B2(3), B3(3), dQdS(3);
   
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A1 = -falpha_chi*(fchi - fchi_r);
   double A2 = -falpha_c*(fc - fc_r);
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   double Sig_p = Sig.Trace()/3.0;
   double B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   double dQdP = 2.*ftan_psi*(fc - Sig_p*ftan_psi);
   
   Sig_Dev.Deviatoric(Sig);
   dQdS = Sig_Dev;
   B2.SetToScaled(1.0/fGf_I, Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
      
   qbar[0]  = A1*B1*dQdP; 
   qbar[0] += A1*dMatrixT::Dot(B2,dQdS);
   qbar[1]  = A2*dMatrixT::Dot(B3,dQdS);
   qbar[2]  = A3*dMatrixT::Dot(B3,dQdS);
   qbar[3]  = A4*dMatrixT::Dot(B3,dQdS);
 }
 
/* calculation of dqbardSig_f */
void MRSSNLHardT::dqbardSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
   dSymMatrixT dhchi_dSig(3), dhc_dSig(3), dhtanphi_dSig(3), dhtanpsi_dSig(3);
   dSymMatrixT Sig_Dev(3), B2(3), B3(3), dQdS(3), dB2dS_dQdS(3), dB3dS_dQdS(3);
   dSymMatrixT tempmat(3);
   
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A1 = -falpha_chi*(fchi - fchi_r);
   double A2 = -falpha_c*(fc - fc_r);
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   double Sig_p = Sig.Trace()/3.0;
   double SN = signof(Sig_p);
   double B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   double dB1dP = (SN +fabs(SN))/2./fGf_I;
   double dQdP = 2.*ftan_psi*(fc - Sig_p*ftan_psi);
   double d2QdP2 =  -2.*ftan_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   dQdS = Sig_Dev;
   B2.SetToScaled(1.0/fGf_I, Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   dB2dS_dQdS.SetToScaled(1.0/fGf_I, Sig_Dev);
   dB3dS_dQdS.SetToScaled(1.0/fGf_II, Sig_Dev);
   
   dhchi_dSig = 0.0;
   dhchi_dSig.PlusIdentity();
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
   
   dqbardSig(0,0) = dhchi_dSig[0];
   dqbardSig(0,1) = dhchi_dSig[1];
   dqbardSig(0,2) = dhchi_dSig[2];
   dqbardSig(0,3) = dhchi_dSig[3];
   dqbardSig(0,4) = dhchi_dSig[4];
   dqbardSig(0,5) = dhchi_dSig[5];
   dqbardSig(1,0) = dhc_dSig[0];
   dqbardSig(1,1) = dhc_dSig[1];
   dqbardSig(1,2) = dhc_dSig[2];
   dqbardSig(1,3) = dhc_dSig[3];
   dqbardSig(1,4) = dhc_dSig[4];
   dqbardSig(1,5) = dhc_dSig[5];
   dqbardSig(2,0) = dhtanphi_dSig[0];
   dqbardSig(2,1) = dhtanphi_dSig[1];
   dqbardSig(2,2) = dhtanphi_dSig[2];
   dqbardSig(2,3) = dhtanphi_dSig[3];
   dqbardSig(2,4) = dhtanphi_dSig[4];
   dqbardSig(2,5) = dhtanphi_dSig[5];
   dqbardSig(3,0) = dhtanpsi_dSig[0];
   dqbardSig(3,1) = dhtanpsi_dSig[1];
   dqbardSig(3,2) = dhtanpsi_dSig[2];
   dqbardSig(3,3) = dhtanpsi_dSig[3];
   dqbardSig(3,4) = dhtanpsi_dSig[4];
   dqbardSig(3,5) = dhtanpsi_dSig[5];
}
  
/* calculation of dqbardq_f */
void MRSSNLHardT::dqbardq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{
   dSymMatrixT Sig_Dev(3), B2(3), B3(3), dQdS(3);
   
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A1 = -falpha_chi*(fchi - fchi_r);
   double A2 = -falpha_c*(fc - fc_r);
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   double Sig_p = Sig.Trace()/3.0;
   double B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   double dQdP = 2.*ftan_psi*(fc - Sig_p*ftan_psi);
   
   Sig_Dev.Deviatoric(Sig);
   dQdS = Sig_Dev;
   B2.SetToScaled(1.0/fGf_I, Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   double B2dQdS = dMatrixT::Dot(B2,dQdS);
   double B3dQdS = dMatrixT::Dot(B3,dQdS);
   
   dqbardq = 0.0;
   dqbardq(0,0) = -falpha_chi*(B1*dQdP + B2dQdS);
   dqbardq(0,1) =  A1*B1*(2.*ftan_psi);
   dqbardq(0,3) =  A1*B1*(2.*fc-4.*Sig_p*ftan_psi);   
   dqbardq(1,1) = -falpha_c*B3dQdS;
   dqbardq(2,2) = -falpha_phi*B3dQdS;
   dqbardq(3,3) = -falpha_psi*B3dQdS;
}

/* return the consistent elstoplastic moduli 
 *
 * Note: Return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& MRSSNLHardT::Moduli(const ElementCardT& element, 
    int ip)
{
    double bott, dlam;
    
    /* allocate matrices */
    dMatrixT AA(10), AA_inv(10), KE(6), KE_Inv(6), CMAT(10); 
    dMatrixT A_qq(4), A_uu(6), A_uq(6,4), A_qu(4,6), A_qu_trans(6,4);
    dMatrixT Rmat(6), dQdSig2(6), dqbardq(4), dQdSigdq(6,4);
    dMatrixT dqbardSig(4,6), R_Inv(6), KEA(6), KEA_Inv(6);
    dMatrixT KP(6), KP2(6), KEP(6), KES(6), KES_Inv(6);   
    dMatrixT Ch(4), Ch_Inv(4), KE1(4,6), KE2(6), KE3(6,4);
         
    /* allocate reduced index vector of symmetric matrices */
    dSymMatrixT Sig(3), Sig_I(3);
    
    /* allocate vectors */
    dArrayT  Rvec(10), Cvec(10), R(10), Rmod(10); 
    dArrayT  qn(4), qo(4), dQdSig(6), dfdq(4), qbar(4);
    dArrayT  R2(10), X(10), dfdSig(6), K1(6), K2(6);
    dArrayT  state(30);

     /* elastic moduli tensor */
    KE = 0.;
    KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
    KE(1,2) = KE(2,1) = KE(1,0) = KE(0,1) = KE(2,0) = KE(0,2) = flambda;
    KE(5,5) = KE(4,4) = KE(3,3) = fmu;
    
    if(!element.IsAllocated()) {
    	fModuli = KE;
        return fModuli;
    }
    
    /* load internal state variables */
    if(!element.IsAllocated()) {
    	LoadData(element,ip);
        state.CopyIn(0, fInternal);
    }
              
    Sig.CopyPart(0, state, 0, Sig.Length());          
    qn.CopyPart(0, state, 18, qn.Length());
    KE_Inv.Inverse(KE);

    if (state[24] == 0.) 
    {
        fModuli = KE;
        fModuli.CopySymmetric();
    }
    else 
        if (state[24] == 1.) 
        {
        	dlam = state[23];
        	
        	/* calculate the first part of Cep */
            dQdSig2_f(qn, dQdSig2);
            dqbardSig_f(Sig, qn, A_qu);
            dqbardq_f(Sig, qn, A_qq);
            dQdSigdq_f(Sig, qn, A_uq);
            Ch.SetToScaled(-dlam, A_qq);
            Ch.PlusIdentity();
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
                 
            /* form AA_inv matrix */
        	A_uu.SetToScaled(dlam, dQdSig2);
        	A_uu += KE_Inv;
        	A_uq.SetToScaled(dlam, dQdSigdq);
        	A_qu.SetToScaled(dlam, dqbardSig);
        	A_qq.SetToScaled(dlam, dqbardq);
        	A_qq.PlusIdentity(-1.0);
        	AA_inv = 0.0;
        	AA_inv.AddBlock(0,           0,           A_uu);
        	AA_inv.AddBlock(0,           A_uu.Cols(), A_uq);
        	AA_inv.AddBlock(A_uu.Rows(), 0,           A_qu);
        	AA_inv.AddBlock(A_uu.Rows(), A_uu.Cols(), A_qq);
            AA.Inverse(AA_inv);
    
            /* calculate second part of Cep */
            dArrayT tmpVec(10), Vvec(6), Vvec2(6), dVec(6);
        	dfdSig_f(Sig, qn, dfdSig);
        	dfdq_f(Sig,qn, dfdq);
        	dQdSig_f(Sig, qn, dQdSig);
        	qbar_f(Sig, qn, qbar);
       	 	Rvec.CopyIn(0, dfdSig);
        	Rvec.CopyIn(dfdSig.Length(), dfdq);
        	Cvec.CopyIn(0, dQdSig);
        	Cvec.CopyIn(dQdSig.Length(), qbar);
        	AA.Multx(Cvec, tmpVec);
        	bott = dArrayT::Dot(Rvec, tmpVec); /* H (scalar) */
        	A_uu.Multx(dfdSig, Vvec);
        	A_qu_trans.Transpose(A_qu);
        	A_qu_trans.Multx(dfdq, Vvec2);
        	Vvec += Vvec2;      /* H (scalar) */
        	KP.Outer(dQdSig, Vvec);     
        	KE3.MultAB(dQdSigdq, Ch_Inv);
        	KE3.Multx(qbar, dVec);
        	KP2.Outer(dVec, Vvec);
        	KP2 *= dlam;
        	KP += KP2;
        	KP /= -bott;
        	KP.PlusIdentity();
        
        	/* calculate Cep */
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

/* computes the stress corresponding to the given element
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
