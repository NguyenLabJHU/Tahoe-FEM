/* $Id: MRSSNLHardT.cpp,v 1.11 2005-12-24 16:19:17 kyonten Exp $ */
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

const int    kNumInternal = 28; // number of internal state variables
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;
const int    kNSTR        = dSymMatrixT::NumValues(kNSD);
const double ratio23      = 2.0/3.0;

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

  	int kk, iplastic;
  	double dlam = 0.0;
  	double dlam2 = 0.0;
  	double ff, bott, topp, normr;
  	
	/* allocate matrices */
    dMatrixT KE(6), AA(10), AA_inv(10), KE_Inv(6), CMAT(10);  
    dMatrixT A_uu(6), A_uq(6,4), A_qu(4,6), A_qq(4);
    dMatrixT dQdSig2(6), dQdSigdq(6,4), dqbardq(4), dqbardSig(4,6);
    
	/* allocate reduced index vector of symmetric matrices */
    dSymMatrixT u(3), up(3), upo(3), du(3), dup(3), ue(3);
    dSymMatrixT Sig(3), Sig_I(3), Sig_trial(3), Sig_e(3);;
    
    /* allocate vectors */
    dArrayT Rvec(10), Cvec(10), R(10), Rmod(10), X(10), Y(10);
    dArrayT qo(4), qn(4), dq(4), R_up(6), R_q(4); 
    dArrayT dfdSig(6), dfdq(4), dQdSig(6), qbar(4);  
    dArrayT state(28);
	
	/* elastic moduli tensor */
	KE = 0.0;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
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
    
	/* initialize in the case of first plastic loading*/
	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, element, ip) && 
	    !element.IsAllocated())
	{
		/* new plastic element */
		AllocateElement(element);
		
		/* initialize element data */ 
		PlasticLoading(trialstrain, element, ip); 
		
		/* fetch internal variables */
		state.CopyIn(0, fInternal);
		state[18] = fchi;
    	state[19] = fc;
    	state[20] = ftan_phi;
    	state[21] = ftan_psi;
	}
    
	/* calculate incremental strains and initialize the necessary vectors */
	u = trialstrain;
    for (int i = 0; i < 6; i++)
       du[i] = u[i] - state[i+6];
    up.CopyPart(0, state, 12, up.Length());
    upo = up;
    qn.CopyPart(0, state, 18, qn.Length());
    qo = qn;
    
    /* calculate stress */
    Sig_I = 0.; 
    Sig = Sig_I;
    ue.DiffOf(u, up);
    KE.MultTx(ue, Sig_e);
    Sig += Sig_e;
    Sig_trial = Sig;
    KE_Inv.Inverse(KE);
    
/* check the yield function */
    Yield_f(Sig, qn, ff);
    if (ff < kYieldTol) {
      iplastic = kIsElastic;
      state[22] = ff;
      state[27] = ff;
      normr = 0.;
      state[25] = normr;
      kk = 0.;
    }
      
    else {
      state[27] = ff;
      kk = 0;
      iplastic = kIsPlastic;
      bool TolExceeded = true;
       
      while (TolExceeded) {
        if (kk > 500)
        	ExceptionT::GeneralFail("MRSSNLHardT::StressCorrection","Too Many Iterations");
        
        /* calculate stress */
        Sig = Sig_I;
        ue.DiffOf(u, up);
        KE.Multx(ue, Sig_e);
        Sig += Sig_e;
        
        /* check yield condition */
        Yield_f(Sig, qn, ff);
        
        /* residuals for plastic strain and internal variables */
        dQdSig_f(Sig, qn, dQdSig);
        qbar_f(Sig, qn, qbar);
        R_up.SetToScaled(dlam, dQdSig); 
        R_up += upo;
        R_up -= up;
        R_q.SetToScaled(dlam, qbar);
        R_q += qo;
        R_q -= qn;
        R.CopyIn(0, R_up);
        R.CopyIn(R_up.Length(), R_q);
        /*
        cout << "kk = " << kk << endl;
        cout << "up " << endl;
        cout << up << endl << endl;
        cout << "upo " << endl;
        cout << upo << endl << endl;
        */       
        /* L2 norms of the residual vectors */
        normr = R.Magnitude();
        double norm_up = R_up.Magnitude();
        double norm_q = R_q.Magnitude();
        
        /* exit the loop if ff < fTol_1 && normr < fTol_2 */
        //cout << "k=" << kk << "   ff=" << ff << "      norm=" << normr << endl;
        if (ff < fTol_1 && normr < fTol_2) TolExceeded = false;
        
        /* print on screen to check which conditions are met */
        //if (ff < fTol_1) cout << "ff < fTol_1 is satisfied" << endl;
        //if (normr < fTol_2) cout << "normr < fTol_2 is satisfied" << endl;
        
        /* check residuals of plastic strain and internal variables separately */
		/*
        cout << "k=" << kk << "   ff=" << ff << "      norm_up=" << norm_up
             << "      norm_q=" << norm_q << endl;
        if(ff < fTol_1 && norm_up < fTol_2 && norm_q < fTol_2)
        	TolExceeded = false;
        if (ff < fTol_1) cout << "ff < fTol_1 is satisfied" << endl;
        if (norm_up < fTol_2) cout << "norm_up < fTol_2 is satisfied" << endl;
        if (norm_q < fTol_2) cout << "norm_q < fTol_2 is satisfied" << endl;
        */
        
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
        A_qq -= Identity4x4;
        AA_inv = 0.0;
        AA_inv.AddBlock(0,           0,           A_uu);
        AA_inv.AddBlock(0,           A_uu.Cols(), A_uq);
        AA_inv.AddBlock(A_uu.Rows(), 0,           A_qu);
        AA_inv.AddBlock(A_uu.Rows(), A_uu.Cols(), A_qq);
        
        /* calculate dlam2 */
        dArrayT tmpVec(10); /* work space */
        dfdSig_f(Sig, qn, dfdSig);
        dfdq_f(Sig, qn, dfdq);
        Rvec.CopyIn(0, dfdSig);
        Rvec.CopyIn(dfdSig.Length(), dfdq);
        Cvec.CopyIn(0, dQdSig);
        Cvec.CopyIn(dQdSig.Length(), qbar);
        AA.Inverse(AA_inv);
        AA.Multx(R, tmpVec);
        topp = ff;
        topp -= dArrayT::Dot(Rvec, tmpVec);
        AA.Multx(Cvec, tmpVec);
        bott = dArrayT::Dot(Cvec, tmpVec);
        dlam2 = topp/bott;
        //cout << "dlam2 = " << dlam2 << endl;
        //cout << "sign of topp " << signof(topp) << endl;
        //cout << "sign of bott " << signof(bott) << endl;		
        //cout << "k = " << kk << " dlam = " << dlam << endl;
        
        /* calculate dup and dq */
        dMatrixT I_mat(4,4);
        CMAT = 0.0; 
        I_mat.SetToScaled(-1.0, Identity4x4);
        CMAT.AddBlock(0, 0, KE_Inv);
        CMAT.AddBlock(KE_Inv.Rows(), KE_Inv.Cols(), I_mat);
        Rmod.CopyIn(0, dQdSig);
        Rmod.CopyIn(dQdSig.Length(), qbar);
        Rmod *= dlam2;
        Rmod += R;
        AA.Multx(Rmod, X);
        CMAT.Multx(X, Y);
        dup.CopyPart(0, Y, 0, dup.Length());
        dq.CopyPart(0, Y, dup.Length(), dq.Length());
        
        /* update state variables and plastic multiplier */
        //upo = up; // previous up
        //qo = qn; // previous qn
        up += dup;
        qn += dq;
        dlam += dlam2;
        
        /* print dlam to screen */
        cout << "dlam = " << dlam << endl;
        kk++;
      }
    }
    
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
	
	fStressCorr = Sig;
	      
	if (iplastic == kIsPlastic) {
	   fInternal.CopyIn(0, state);
	   fPlasticStrain = up;
	}		
 return fStressCorr;
}

/*
 * Returns the value of the yield function given the
 * stress vector and state variables
 */
void MRSSNLHardT::Yield_f(const dSymMatrixT& Sig, 
			const dArrayT& qn, double& ff)
{
  dSymMatrixT Sig_Dev(3);
  double fchi = qn[0];
  double fc = qn[1];
  double ffriction = qn[2]; 
  double fpress = Sig.Trace()/3.0;
  
  Sig_Dev.Deviatoric(Sig);
  ff = Sig_Dev.Invariant2();
  ff -= pow((fc - ffriction*fpress), 2);
  ff += pow((fc - ffriction*fchi), 2);
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

   double temp  = ratio23*ftan_phi;
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

   double temp  = ratio23*ftan_psi;
   temp *= (fc - Sig_p*ftan_psi);
   for (int i = 0; i < 3; i++) 
      dQdSig[i] += temp;
}

/* calculation of dQdSig2_f */
void MRSSNLHardT::dQdSig2_f(const dArrayT& qn, dMatrixT& dQdSig2)
{
  double ftan_psi = qn[3];
  dMatrixT I_mat(6,6);
  
  I_mat = 0.;
  for (int i = 0; i < 3; i++) 
  {
     for (int j = 0; j < 3; j++) 
     	I_mat(i,j) = 1.;
  }
  
  double Fac = ratio23*ftan_psi*ftan_psi;
  Fac += 1.;
  Fac /= 3.;
  I_mat *= Fac;
  
  dQdSig2 = Identity6x6;
  dQdSig2 -= I_mat;
}

/* calculation of dQdSigdq_f */
void MRSSNLHardT::dQdSigdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
  double fc = qn[1];
  double ftan_psi = qn[3];
  double Sig_p = Sig.Trace()/3.0;
  
  dQdSigdq = 0.;
  dQdSigdq(1,0) = ratio23*ftan_psi;
  dQdSigdq(1,1) = dQdSigdq(1,2) = dQdSigdq(1,0);
  dQdSigdq(3,0) = ratio23*(fc - 2.*Sig_p*ftan_psi);
  dQdSigdq(3,1) = dQdSigdq(3,2) = dQdSigdq(3,0); 
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
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdS = Sig_Dev;
   B3 = Sig_Dev;
   B3 /= fGf_II; 
      
   qbar[0]  = A1*B1*dQdP; 
   qbar[0] += A1*dMatrixT::Dot(B2,dQdS);
   qbar[1]  = dMatrixT::Dot(B3,dQdS);
   qbar[1]  *=A2;
   qbar[2]  = dMatrixT::Dot(B3,dQdS);
   qbar[2]  *=A3;
   qbar[3]  = dMatrixT::Dot(B3,dQdS);
   qbar[3]  *=A4;
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
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdS = Sig_Dev;
   B3 = Sig_Dev;
   B3 /= fGf_II;
   dB2dS_dQdS  = Sig_Dev;
   dB2dS_dQdS /= fGf_I;
   dB3dS_dQdS  = Sig_Dev;
   dB3dS_dQdS /= fGf_II;
   dhchi_dSig  = Identity3x3;
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
   
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdS = Sig_Dev;
   B3 = Sig_Dev;
   B3 /= fGf_II;
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
     dMatrixT KE(6), AA(10), AA_inv(10), CMAT(10), KE_Inv(6); 
     dMatrixT A_uu(6), A_uq(6,4), A_qu(4,6), A_qu_trans(6,4), A_qq(4);
     dMatrixT dQdSig2(6), dqbardq(4), dQdSigdq(6,4), dqbardSig(4,6);
     dMatrixT KP(6), KP2(6), KEP(6), KES(6), KES_Inv(6);
     dMatrixT Ch(4), Ch_Inv(4), KE1(4,6), KE2(6), KE3(6,4);
     
     dSymMatrixT Sig(3), Sig_I(3);     
     dArrayT dfdSig(6), dfdq(4), dQdSig(6), qbar(4);
     dArrayT qn(4), Rvec(10), Cvec(10);
	
    /* elastic moduli tensor */
	KE = 0.0;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
    if(element.IsAllocated() && (element.IntegerData())[ip] == kIsPlastic) {
	  	/* load internal state variables */
	  	LoadData(element, ip);
	  	Sig.CopyPart(0, fInternal, 0, Sig.Length());  	
    	qn.CopyPart(0, fInternal, 18, qn.Length());
		dlam = fInternal[23];
		KE_Inv.Inverse(KE);
		
		/* calculate the first part of Cep */
	  	dQdSig2_f(qn, dQdSig2);
	    dqbardSig_f(Sig, qn, dqbardSig);
	    dqbardq_f(Sig, qn, dqbardq);
	    dQdSigdq_f(Sig, qn, dQdSigdq);
	    Ch.SetToScaled(-dlam, dqbardq);
	    Ch += Identity4x4;
	    Ch_Inv.Inverse(Ch);
	    KE1.MultAB(Ch_Inv, dqbardSig);
	    KES.MultAB(dQdSigdq, KE1);
	    KES *= pow(dlam, 2);
	    /*KES = 0.;*/
	    KE2.SetToScaled(dlam, dQdSig2);
	    KES += KE2;
	    KES += KE_Inv;
	    KES_Inv.Inverse(KES);
	    
        /* form AA_inv matrix */
        A_uu.SetToScaled(dlam, dQdSig2);
        A_uu += KE_Inv;
        A_uq.SetToScaled(dlam, dQdSigdq);
        A_qu.SetToScaled(dlam, dqbardSig);
        A_qq.SetToScaled(dlam, dqbardq);
        A_qq -= Identity4x4;
        AA_inv = 0.0;
        AA_inv.AddBlock(0,           0,           A_uu);
        AA_inv.AddBlock(0,           A_uu.Cols(), A_uq);
        AA_inv.AddBlock(A_uu.Rows(), 0,           A_qu);
        AA_inv.AddBlock(A_uu.Rows(), A_uu.Cols(), A_qq);
	
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
        Vvec += Vvec2;       /* V (vector) */
        KP.Outer(dQdSig, Vvec);     
        KE3.MultAB(dQdSigdq, Ch_Inv);
        KE3.Multx(qbar, dVec);
        KP2.Outer(dVec, Vvec);
	    KP2 *= dlam;
	    KP += KP2;
	    KP /= -bott;
        KP += Identity6x6;
        
        /* calculate Cep */
        KEP.MultAB(KES_Inv, KP);
	    fModuli = KEP;
	    return fModuli;
	}
	else {
		fModuli = KE;
	  	return fModuli;
	}
	
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
	fModuli.Dimension(kNSTR);
	fModuliPerfPlas.Dimension(kNSTR);
	fDevStress.Dimension(kNSD);
	fDevStrain.Dimension(kNSD); 
	fTensorTemp.Dimension(kNSTR);
	IdentityTensor2.Dimension(kNSD);
	Identity3x3.Dimension(kNSD); 
	Identity4x4.Dimension(kNSD+1);
	Identity6x6.Dimension(kNSTR);
    
	/* initialize constant tensor */
	Identity3x3.Identity();
	Identity4x4.Identity(); 
	Identity6x6.Identity();
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
		if (fInternal[kplastic] == kIsPlastic)
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