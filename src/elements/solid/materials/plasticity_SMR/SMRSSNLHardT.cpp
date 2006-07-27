/* $Id: SMRSSNLHardT.cpp,v 1.2 2006-07-27 21:56:03 kyonten Exp $ */
/* created: Karma Yonten */

/* Interface for a nonassociative, small strain,      */
/* pressure dependent plasticity model with nonlinear */ 
/* isotropic hardening/softening.                     */

#include "SMRSSNLHardT.h"
#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

/* class constants */

using namespace Tahoe;

const int    kNumInternal = 31; // number of internal state variables
const double kYieldTol    = 1.0e-8;
const int    kNSD         = 3;
const int    kNSTR        = dSymMatrixT::NumValues(kNSD);
const double ratio32      = 3.0/2.0;

/* constructor */
SMRSSNLHardT::SMRSSNLHardT(int num_ip, double mu, double lambda):
	fNumIP(num_ip),
	fmu(mu),
	flambda(lambda),
	fkappa(flambda + (2.0/3.0*fmu)),
	fMeanStress(0.0)
{
	SetName("SMR_SS_nonlinear_hardening");
}

const dSymMatrixT& SMRSSNLHardT::ElasticStrain(const dSymMatrixT& totalstrain, 
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
const dSymMatrixT& SMRSSNLHardT::StressCorrection(
      const dSymMatrixT& totalstrain_curr,
      ElementCardT& element, int ip)
{
  	int kk, iplastic;
    double ff, dlam, dlam2, normr;
  	
	/* allocate matrices */
    dMatrixT KE(6), AA(10), AA_inv(10), KE_Inv(6), CMAT(10);  
    dMatrixT Auu_inv(6), Auq_inv(6,4), Aqu_inv(4,6), Aqq_inv(4);
    dMatrixT dQdSig2(6), dQdSigdq(6,4), dqbardq(4), dqbardSig(4,6);
    
	/* allocate reduced index vector of symmetric matrices */
    dSymMatrixT u(3), up(3), upo(3), du(3), dup(3), ue(3);
    dSymMatrixT Sig(3), Sig_I(3), Sig_trial(3), Sig_e(3);
    
    /* allocate vectors */
    dArrayT Rvec(10), Cvec(10), R(10), Rmod(10), X(10), Y(10);
    dArrayT qo(4), qn(4), dq(4), R_up(6), R_q(4), R2(10); 
    dArrayT dfdSig(6), dfdq(4), dQdSig(6), qbar(4);  
    dArrayT state(31);
	
	/* elastic moduli tensor */
	KE = 0.0;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
        
	//if(ip == 0)
		//cout << "ip "<< ip << endl; // ip #
		
	/* case 1: element/ip elastic */
	if (!PlasticLoading(totalstrain_curr, element, ip) && 
	    !element.IsAllocated())
	{
		/* initialize element data */
        double esp  = 0.;
        state = 0.;
        state[ktanphi] = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
        state[ktanpsi] = (tan(fphi_p))*exp(-falpha_psi*esp);
        //if(ip==0) cout << "elastic step " << endl;
	}
	
	/* case 2: element/ip plastic for the first time */
	if (PlasticLoading(totalstrain_curr, element, ip) && 
	    !element.IsAllocated())
	{
		/* new plastic element */
		AllocateElement(element);
		
		/* initialize element data */ 
		PlasticLoading(totalstrain_curr, element, ip); 
		//if(ip==0) cout << "first plastic step " << endl;
	}
	
	/* case 3: element/ip plastic currently and/or previously */
	// for other ips in the plastic element initialize the internal variables
	// if same ip of the element is encountered the second time, use the 
	// previous values of the internal variables/state variables
	if (PlasticLoading(totalstrain_curr, element, ip) && 
        element.IsAllocated()) {
     	
     	LoadData(element, ip);
     	
     	/* fetch internal variables */
     	state.CopyIn(0, fInternal);
     	
        //if(ip==0) cout << "plastic step " << endl;
    }
    
    /* case 4: element/ip previously plastic is unloaded */
	// set internal variables to previous values
	if (!PlasticLoading(totalstrain_curr, element, ip) && 
	    element.IsAllocated())
	{
        LoadData(element, ip);
        
       	/* fetch internal variables */
     	state.CopyIn(0, fInternal);
     	//if(ip==0) cout << "unloading step" << endl;
	}
	
	/* initialize and copy the necessary vectors */
    up.CopyPart(0, state, 12, up.Length());
    upo = up;
    qn.CopyPart(0, state, 18, qn.Length());
    qo = qn;
    
    dup = 0.; dq = 0.;
    dlam = 0.; dlam2 = 0.; normr = 0.;
    
    /* calculate stress */
    ue.DiffOf(totalstrain_curr, up);
    //KE.Multx(ue, Sig_e);
    Sig_e.A_ijkl_B_kl(KE, ue); 
    Sig = Sig_e; 
    KE_Inv.Inverse(KE);
  
/* check the yield function */
    ff = Yield_f(Sig, qn);
    if (ff < kYieldTol) {
      iplastic = 0;
      state[27] = ff;
      normr = 0.;
      state[25] = normr;
      kk = 0;
    } 
    else {
      state[27] = ff;
      kk = 0;
      iplastic = 1;
      bool NotConverged = true;
      while (NotConverged) 
      {
        if (kk > 100)
        	ExceptionT::GeneralFail("SMRSSNLHardT::StressCorrection","Too Many Iterations");
        
        /* calculate stress */
        ue.DiffOf(totalstrain_curr, up);
    	//KE.Multx(ue, Sig_e);
    	Sig_e.A_ijkl_B_kl(KE, ue); 
    	Sig = Sig_e;
        
        /* check yield condition */
        ff = Yield_f(Sig, qn);
 
        /* residuals for plastic strain and internal variables */
        dQdSig_f(Sig, qn, dQdSig);
        qbar_f(Sig, qn, qbar);

        for (int i = 0; i < 6; i++) {
          R[i]  = upo[i]-up[i];
          R[i] += dlam*dQdSig[i];
        }
        for (int i = 0; i < 4; i++) {
          R[i+6]  = qo[i]-qn[i]; 
          R[i+6] += dlam*qbar[i]; 
        }
            
        /* L2 norms of the residual vectors */
        normr = R.Magnitude();
        
        /* form AA_inv matrix  and calculate AA */
        dQdSig2_f(Sig, qn, dQdSig2);
        dQdSigdq_f(dQdSigdq);
        dqbardSig_f(Sig, qn, dqbardSig);
        dqbardq_f(Sig, qn, dqbardq);
        Auu_inv.SetToScaled(dlam, dQdSig2);
        Auu_inv += KE_Inv;
        Auq_inv.SetToScaled(dlam, dQdSigdq);
        Aqu_inv.SetToScaled(dlam, dqbardSig);
        Aqq_inv.SetToScaled(dlam, dqbardq);
        Aqq_inv.PlusIdentity(-1.0);
        AA_inv = 0.0;
        AA_inv.AddBlock(0,           0,           Auu_inv);
        AA_inv.AddBlock(0,           Auu_inv.Cols(), Auq_inv);
        AA_inv.AddBlock(Auu_inv.Rows(), 0,           Aqu_inv);
        AA_inv.AddBlock(Auu_inv.Rows(), Auu_inv.Cols(), Aqq_inv);
       	AA.Inverse(AA_inv);
       	
        /* calculate dlam2 */
        dArrayT tmpVec(10); 
        dfdSig_f(Sig, qn, dfdSig);
        dfdq_f(Sig, dfdq);
        Rvec.CopyIn(0, dfdSig);
        Rvec.CopyIn(dfdSig.Length(), dfdq);
        Cvec.CopyIn(0, dQdSig);
        Cvec.CopyIn(dQdSig.Length(), qbar);
        
        /* include contribution of all off diagonal terms 
         * in reduced vector of  symmetric matrices 
         * dfdSig, dQdSig, up */
        R2=R;
        for(int i = 0; i < 3; i++)
        {
        	Rvec[i+3] = Rvec[i+3]*2.0;
        	Cvec[i+3] = Cvec[i+3]*2.0; 
        	R2[i+3] = R2[i+3]*2.0;
        }
        
        AA.Multx(R2, tmpVec);
        double topp = (ff - dArrayT::Dot(Rvec, tmpVec));
        AA.Multx(Cvec, tmpVec);
        double bott = dArrayT::Dot(Rvec, tmpVec);
        dlam2 = topp/bott;
        
        /* calculate dup and dq */
        dMatrixT I_mat(4);
        CMAT = 0.0;  
        I_mat.Identity(-1.0);
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
        
        /* update internal variables and plastic multiplier */
        up += dup;
        qn += dq;
        dlam += dlam2;                                                                                                                                                                                                                                                                                                    
        /*
        if(ip == 0 ) { 
        	cout << "ip="<< ip <<" kk=" << kk << " ff=" << ff << " normr=" << normr
             	 << " dlam2=" << dlam2 << " dlam=" << dlam << endl << endl;
        }
        */
        kk++;
        
        /* exit while loop if solution has converged */
        /* convergence criteria: stresses brought back to the yield 
         * surface (i.e. ff ~= 0) and residuals of plastic strains 
         * and internal variables are very small */
        if (abs(ff) < fTol_1 && normr < fTol_2) {
        	NotConverged = false;
        }
        
      } // while (NotConverged)
    } // if (ff < kYieldTol)
    
    /* update state variables */
    state.CopyIn(0, Sig);
    state.CopyIn(Sig.Length(), totalstrain_curr); 	   
	state.CopyIn(12, up);
	state.CopyIn(18, qn);
	state[22] = ff;
	state[23] = dlam;
	state[24] = double(iplastic);
	state[25] = normr;
	state[26] = double(kk);
	
	fStressCorr = Sig;
	if (iplastic == 1) {
   		/* indicator during the first elastic to plastic transition
   		 * use qn from the fInternal(i.e. plastic) and not elastic qn in second plastic step */ 
   		state[30] = 1.; 
   		
	   	fInternal.CopyIn(0, state);
	   	fPlasticStrain = up;
	}	
	
 return fStressCorr;
}

/*
 * returns the value of the yield function given the
 * stress vector and state variables
 */
double SMRSSNLHardT::Yield_f(const dSymMatrixT& Sig, 
			const dArrayT& qn)
{
  dSymMatrixT Sig_Dev(3);
  double ffriction = qn[2]; 
  double fpress = Sig.Trace()/3.0;
  
  Sig_Dev.Deviatoric(Sig);
  double temp  = (Sig_Dev.ScalarProduct())/2.0;
  double ff = sqrt(3.*temp) + (ffriction*fpress) - fc;
  return ff;
}

/* calculation of dfdSig_f */
void SMRSSNLHardT::dfdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
   dSymMatrixT Sig_Dev(3);
   double ftan_phi = qn[2]; 
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double deno = 2.*sqrt(3.*temp)/3.;
   Sig_Dev /= deno;
   dfdSig=0.0;
   for (int i = 0; i < 6; i++)
   		if (i < 3)  
   			dfdSig[i] = Sig_Dev[i] + (ftan_phi/3.0);
   		else 
   			dfdSig[i] = Sig_Dev[i];
}

/* calculation of dfdq_f */
void SMRSSNLHardT::dfdq_f(const dSymMatrixT& Sig, dArrayT& dfdq)
{
   double Sig_p = Sig.Trace()/3.0;
   dfdq = 0.0;
   dfdq[2] = Sig_p;
}

/* calculation of dQdSig_f */
void SMRSSNLHardT::dQdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
   dSymMatrixT Sig_Dev(3);
   double ftan_psi = qn[3];
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double deno = 2.*sqrt(3.*temp)/3.;
   Sig_Dev /= deno;
   dQdSig=0.0;
   for (int i = 0; i < 6; i++)
   		if (i < 3)  
   			dQdSig[i] = Sig_Dev[i] + (ftan_psi/3.0);
   		else 
   			dQdSig[i] = Sig_Dev[i];
}

/* calculation of dQdSig2_f */
void SMRSSNLHardT::dQdSig2_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSig2)
{
  dSymMatrixT Sig_Dev(3); 
  dMatrixT I_mat(6), I(3), matrix2(6);
  double ftan_psi = qn[3];
  
  Sig_Dev.Deviatoric(Sig);
  double temp  = (Sig_Dev.ScalarProduct())/2.0;
  double q = sqrt(3.*temp);
  I_mat = 0.;
  for (int i = 0; i < 3; i++) 
  {
     for (int j = 0; j < 3; j++) 
     	I_mat(i,j) = 1.;
  }
  dQdSig2=0.0;
  I_mat /= 3.0;
  dQdSig2.Identity(); 
  dQdSig2 -= I_mat;
  dQdSig2 /= q; 
  matrix2.Outer(Sig_Dev,Sig_Dev);
  matrix2 /= (2.*q*q*q/3.);
  dQdSig2 -= matrix2;
  dQdSig2 /= 2./3.;
}

/* calculation of dQdSigdq_f */
void SMRSSNLHardT::dQdSigdq_f(dMatrixT& dQdSigdq)
{
  dQdSigdq = 0.;
  for (int i = 0; i < 3; i++)
  	dQdSigdq(i,3) = 1.0/3.0;
}

/* calculation of qbar_f */
void SMRSSNLHardT::qbar_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& qbar)
{
   dSymMatrixT Sig_Dev(3), B2(3), B3(3), dQdS(3); 
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double deno = 2.*sqrt(3.*temp)/3.;
   dQdS.SetToScaled(1.0/deno,Sig_Dev);
   B3.SetToScaled(1.0/fGf_II,Sig_Dev);
   qbar = 0.0;   
   qbar[2]  = A3*B3.ScalarProduct(dQdS); 
   qbar[3]  = A4*B3.ScalarProduct(dQdS); 
 }
 
/* calculation of dqbardSig_f */
void SMRSSNLHardT::dqbardSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
   dSymMatrixT dhtanphi_dSig(3), dhtanpsi_dSig(3);
   dSymMatrixT Sig_Dev(3), B2(3), B3(3), dQdS(3), dB3dS_dQdS(3), B3_dQdQ_dSdS(3);
   dMatrixT tempmat(6);
   
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double q = sqrt(3.*temp);
   double deno = 2.*q/3.;
   dQdS.SetToScaled(1.0/deno,Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   tempmat.Outer(Sig_Dev,Sig_Dev);
   //tempmat.Multx(Sig_Dev, B3_dQdQ_dSdS);
   Contract4To2(tempmat,Sig_Dev,B3_dQdQ_dSdS);
   B3_dQdQ_dSdS /= -2.*q*q*q/3.;
   B3_dQdQ_dSdS.AddScaled(1./q, Sig_Dev);
   B3_dQdQ_dSdS /= (2.*fGf_II)/3.;
   dB3dS_dQdS.SetToScaled(1.0/fGf_II, dQdS);
   dhtanphi_dSig  = B3_dQdQ_dSdS;
   dhtanphi_dSig += dB3dS_dQdS;
   dhtanphi_dSig *= A3;
   dhtanpsi_dSig  = B3_dQdQ_dSdS;
   dhtanpsi_dSig += dB3dS_dQdS;
   dhtanpsi_dSig *= A4;
   
   dqbardSig = 0.0;
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
void SMRSSNLHardT::dqbardq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{
   dSymMatrixT Sig_Dev(3), B2(3), B3(3), dQdS(3);
   
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double deno = 2.*sqrt(3.*temp)/3.;
   dQdS.SetToScaled(1.0/deno,Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   double B3dQdS = B3.ScalarProduct(dQdS); 
   dqbardq = 0.0;
   dqbardq(2,2)  = -falpha_phi*B3dQdS;
   dqbardq(3,3)  = -falpha_psi*B3dQdS;
}

/* return the consistent elastoplastic moduli 
 *
 * Note: Return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& SMRSSNLHardT::Moduli(const ElementCardT& element, 
	int ip)
{
	 double dlam;
	 
	  /* allocate matrices */
     dMatrixT KE(6), AA(10), AA_inv(10), CMAT(10), KE_Inv(6); 
     dMatrixT Auu_inv(6), Auq_inv(6,4), Aqu_inv(4,6), Aqq_inv(4);
     dMatrixT Auu_Aqu(10,6);
     dMatrixT dQdSig2(6), dqbardq(4), dQdSigdq(6,4), dqbardSig(4,6);
     dMatrixT KP(6), KP2(6), KEP(6), KES(6), KES_Inv(6);
     dMatrixT Ch(4), Ch_Inv(4), KE2(6), KE4(4,6);
     
     /* allocate reduced index vector of symmetric matrices */
     dSymMatrixT Sig(3);  
     
     /* allocate vectors */   
     dArrayT dfdSig(6), dfdq(4), dQdSig(6), qbar(4);
     dArrayT qn(4), Rvec(10), Rvec2(10), Cvec(10);
     
    /* elastic moduli tensor */
	KE = 0.0;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
    if(element.IsAllocated() && (element.IntegerData())[ip] == kIsPlastic) 
    {
	  	/* load internal state variables */
	  	LoadData(element, ip);
	  	Sig.CopyPart(0, fInternal, 0, Sig.Length());
    	qn.CopyPart(0, fInternal, 18, qn.Length());
		KE_Inv.Inverse(KE);
		
		if(fInternal[kplastic] == 0.) 
			fModuli = KE;
		else 
		{
			/* calculate the first part of Cep */
			dlam = fInternal[kdlambda];
	  		dQdSig2_f(Sig, qn, dQdSig2);
	    	dqbardSig_f(Sig, qn, dqbardSig);
	    	dqbardq_f(Sig, qn, dqbardq);
	    	dQdSigdq_f(dQdSigdq);
	    	Ch.SetToScaled(-dlam, dqbardq);
	    	Ch.PlusIdentity();
	    	Ch_Inv.Inverse(Ch);
	    	KES.MultABC(dQdSigdq,Ch_Inv,dqbardSig);
	    	KES *= (dlam*dlam);
	    	KE2.SetToScaled(dlam, dQdSig2);
	    	KES += KE2;
	    	KES += KE_Inv;
	    	KES_Inv.Inverse(KES);
	    
        	/* form AA_inv matrix and calculate AA */
        	Auu_inv.SetToScaled(dlam, dQdSig2);
        	Auu_inv += KE_Inv;
        	Auq_inv.SetToScaled(dlam, dQdSigdq);
        	Aqu_inv.SetToScaled(dlam, dqbardSig);
        	Aqq_inv.SetToScaled(dlam, dqbardq);
        	Aqq_inv.PlusIdentity(-1.0);
        	AA_inv = 0.0;
        	AA_inv.AddBlock(0,           0,           Auu_inv);
        	AA_inv.AddBlock(0,           Auu_inv.Cols(), Auq_inv);
        	AA_inv.AddBlock(Auu_inv.Rows(), 0,           Aqu_inv);
        	AA_inv.AddBlock(Auu_inv.Rows(), Auu_inv.Cols(), Aqq_inv);
        	AA.Inverse(AA_inv);
	
        	/* calculate second part of Cep */
        	dArrayT tmpVec(10), Vvec(6);
        	dMatrixT Rvec_t(1,10),Vvec_t(1,6);
        	dfdSig_f(Sig, qn, dfdSig);
        	dfdq_f(Sig,dfdq);
        	dQdSig_f(Sig, qn, dQdSig);
        	qbar_f(Sig, qn, qbar);
        	Rvec.CopyIn(0, dfdSig);
        	Rvec.CopyIn(dfdSig.Length(), dfdq);
        	Cvec.CopyIn(0, dQdSig);
        	Cvec.CopyIn(dQdSig.Length(), qbar);
        	
        	/* include contribution of all off diagonal terms 
             * in reduced vector of  symmetric matrices 
             * dfdSig and dQdSig  */
            Rvec2 = Rvec;
            for(int i = 0; i < 3; i++)
            {
        	    Rvec2[i+3] = Rvec2[i+3]*2.0;
        	    Cvec[i+3] = Cvec[i+3]*2.0; 
            }
        	AA.Multx(Cvec, tmpVec);
        	double H = dArrayT::Dot(Rvec2, tmpVec); /* H (scalar) */
        	
        	/* collect Auu and Aqu into Auu_Aqu */
        	for (int i = 0; i < 10; i++)
        		for (int j = 0; j < 6; j++) 
        			Auu_Aqu(i,j) = AA(i,j);
            for (int i=0; i<10; i++) //transpose Rvec
        	   Rvec_t(0,i)=Rvec[i];
            Vvec_t.MultAB(Rvec_t,Auu_Aqu); /* V (vector) */ 
            for (int i=0; i<6; i++) 
        	   Vvec[i]=Vvec_t(0,i);     
        	KP.Outer(dQdSig,Vvec);  
        	KE4.Outer(qbar,Vvec);
        	KP2.MultABC(dQdSigdq,Ch_Inv,KE4);
	    	KP2 *= dlam;
	    	KP += KP2;
	    	KP /= -H;
        	KP.PlusIdentity();
        
        	/* calculate Cep */
        	KEP.MultAB(KES_Inv, KP);
        	
        	/* continuum jacobian, i.e., "inconsistent" tangent operator */
   			dMatrixT dfmat(6,1),dQmat(6,1);
   			for (int i=0; i<6; i++) {
   			   if (i>2) {
   			      dfdSig[i] *= 2.;
   			      dfdSig[i] *= 2.;
   			   }
   			   dfmat(i,0) = dfdSig[i];
   			   dQmat(i,0) = dQdSig[i];
   			}
   			
   			KP2.MultABCT(KE,dQmat,dfmat);
   			KP.MultAB(KP2,KE);
   			double bott=KE.MultmBn(dQdSig,dfdSig);
   			bott-=dArrayT::Dot(dfdq,qbar);
   			KP/=bott;
   			//KEP.DiffOf(KE,KP); //uncomment to activate continuum jacobian
   			
	    	fModuli = KEP;
	    }
	}
	else
		fModuli = KE;
		
	//fModuli = KE;  // uncomment to use constant stiffness
	return fModuli;
}

/* return the correction to modulus Cep~, checking for discontinuous
 *   bifurcation */
const dMatrixT& SMRSSNLHardT::ModuliPerfPlas(const ElementCardT& element, 
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
void SMRSSNLHardT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	//d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += kNumInternal*fNumIP;        //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0
}

/* accept parameter list */
void SMRSSNLHardT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SMRPrimitiveT::TakeParameterList(list);

	/* dimension work space */
	fElasticStrain.Dimension(kNSD);
	fStressCorr.Dimension(kNSD);
	fModuli.Dimension(kNSTR);
	fModuliPerfPlas.Dimension(kNSTR);
	fDevStress.Dimension(kNSD);
	fDevStrain.Dimension(kNSD); 
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* element level data */
void SMRSSNLHardT::Update(ElementCardT& element)
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
	// updating done in SMRSSNLHardT::StressCorrection
	//for (int ip = 0; ip < fNumIP; ip++)
		//if (Flags[ip] == kIsPlastic) /* plastic update */
		//{
			/* do not repeat if called again. */
			//Flags[ip] = kIsElastic;
			/* NOTE: ComputeOutput writes the updated internal variables
			 *       for output even during iteration output, which is
			 *       called before UpdateHistory */

			/* fetch element data */
			//LoadData(element, ip);
		//}
}

/* resets to the last converged solution */
void SMRSSNLHardT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* load element data for the specified integration point */
void SMRSSNLHardT::LoadData(const ElementCardT& element, int ip)
{
	/* check */
	if (!element.IsAllocated()) 
	    ExceptionT::GeneralFail("SMRSSNLHardT::LoadData","The element should have been allocated");
	/* fetch arrays */
	 const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Alias(        dim, &d_array[           dex]);
	/*fUnitNorm.Set(        kNSD, &d_array[  offset + dex]); */    
	fInternal.Alias(kNumInternal, &d_array[offset + ip*kNumInternal]); //2*offset if fUnitNorm
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int SMRSSNLHardT::PlasticLoading(const dSymMatrixT& trialstrain, 
	  ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated()) {
	
		double esp  = 0.;
        double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
		return(YieldCondition(DeviatoricStress(trialstrain,element),
			       MeanStress(trialstrain,element),ftan_phi,fc) > kYieldTol );
	}
        /* already plastic */
	else 
	{
	/* get flags */
	 iArrayT& Flags = element.IntegerData();
		
	/* load internal variables */
	LoadData(element, ip);
	
	if(fInternal[30] == 0.) {
		double esp = 0.;
    	fInternal[ktanphi] = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
		fInternal[ktanpsi] = (tan(fphi_p))*exp(-falpha_psi*esp);
	}
	
	dSymMatrixT elasticstrain(3);
	elasticstrain.DiffOf(trialstrain, fPlasticStrain);
	fInternal[kftrial] = YieldCondition(DeviatoricStress(elasticstrain,element),
			             MeanStress(elasticstrain,element),fInternal[ktanphi],fc);

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{		
			/* set flag */
			Flags[ip] = kIsPlastic;
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
		    Flags[ip] = kIsElastic; //removed to avoid resetting 7/01
			return 0;
		}
	}
}	

/* Computes the stress corresponding to the given element
 * and elastic strain.  The function returns a reference to the
 * stress in fDevStress */
dSymMatrixT& SMRSSNLHardT::DeviatoricStress(const dSymMatrixT& trialstrain,
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
double SMRSSNLHardT::MeanStress(const dSymMatrixT& trialstrain,
	const ElementCardT& element)
{
#pragma unused(element)

  fMeanStress = fkappa*trialstrain.Trace();
  return fMeanStress;
}

double SMRSSNLHardT::signof(double& r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

/* off-diagonal terms in reduced symmetric matrix multiplied by 2
   before dot-product operation
*/
double SMRSSNLHardT::DotProduct2(const dArrayT& vec1,const dArrayT& vec2)
{
	dArrayT A(vec1.Length()), B(vec2.Length());
	A = vec1;
	B = vec2;
	for (int i=0; i<3; i++) {
	   A[i+3]=vec1[i+3] *2.;
	   B[i+3]=vec2[i+3] *2.;
	}
	return dArrayT::Dot(A,B);
}

void SMRSSNLHardT::Contract4To2(const dMatrixT& mat,const dSymMatrixT& vec,
                                      dSymMatrixT& res)
{
	dSymMatrixT tmp(3);
	tmp = vec;
	for (int i=0; i<3; i++) {
	   tmp[i+3]=vec[i+3] *2.;
	}
	mat.Multx(tmp, res);
}