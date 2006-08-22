/* $Id: SMRSSNLHardT.cpp,v 1.4 2006-08-22 14:37:11 kyonten Exp $ */
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

using namespace Tahoe;

/* class constants */
const int    kNumInternal = 26; // number of internal state variables
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
  	int kk,iplastic;
	double ff,dlam,dlam2,normr;
  	
	/* allocate matrices */
    dMatrixT KE(6),KE_Inv(6),AA(8),AA_inv(8),Auu_inv(6),Auq_inv(6,2),
             Aqu_inv(2,6),Aqq_inv(2),dQdSig2(6),dQdSigdq(6,2),
             dqbardq(2),dqbardSig(2,6);
    
	/* allocate reduced index vector of symmetric matrices */
    dSymMatrixT u(3),up(3),upo(3),du(3),dup(3),ue(3),Sig(3),dSig(3),
                Sig_e(3),dfdSig(3),dQdSig(3); 
    
    /* allocate vectors */
    dArrayT Rvec(8),Cvec(8),R(8),Rmod(8),X(8),qo(2),qn(2),dq(2),R2(8), 
            dfdq(2),qbar(2),state(26);
	
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
    qn.CopyPart(0, state, 18, qn.Length());
    upo = up; qo = qn;
    dup = 0.; dq = 0.;
    dlam = 0.; dlam2 = 0.; normr = 0.;
    
    /* calculate stress */
    ue.DiffOf(totalstrain_curr, up);
    Sig_e.A_ijkl_B_kl(KE, ue);
    //KE.Multx(ue, Sig_e); 
    Sig = Sig_e; 
    KE_Inv.Inverse(KE);
   
/* check the yield function */
    ff = Yield_f(Sig, qn);
    if (ff < fTol_1) {
      iplastic = 0;
      kk = 0;
    } 
    else {
      state[20] = ff;
      kk = 0;
      iplastic = 1;
      bool NotConverged = true;
      while (NotConverged) 
      {
        if (kk > 20)
        	ExceptionT::GeneralFail("SMRSSNLHardT::StressCorrection","Too Many Iterations");
        
        /* calculate stress */
        ue.DiffOf(totalstrain_curr, up);
    	Sig_e.A_ijkl_B_kl(KE, ue);
    	//KE.Multx(ue, Sig_e);  
    	Sig = Sig_e;
        
        /* check yield condition */
        ff = Yield_f(Sig, qn);
 
        /* residuals for plastic strain and internal variables */
        dQdSig_f(Sig, qn, dQdSig);
        qbar_f(Sig, qn, qbar);

        for (int i=0; i<6; i++) {
          R[i]  = upo[i]-up[i];
          R[i] += dlam*dQdSig[i];
        }
        for (int i=0; i<2; i++) {
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
        dfdSig_f(Sig, qn, dfdSig);
        dfdq_f(Sig, dfdq);
        Rvec.CopyIn(0, dfdSig);
        Rvec.CopyIn(dfdSig.Length(), dfdq);
        Cvec.CopyIn(0, dQdSig);
        Cvec.CopyIn(dQdSig.Length(), qbar);
        
        /* include contribution of all off diagonal terms 
         * in reduced vectors of symmetric matrices for 
         * double contraction 
         * dfdSig, dQdSig, up */
        R2=R;
        for(int i = 0; i < 3; i++)
        {
        	Rvec[i+3] = Rvec[i+3]*2.;
        	Cvec[i+3] = Cvec[i+3]*2.; 
        	R2[i+3] = R[i+3]*2.;
        }
        
        double topp = ff - AA.MultmBn(Rvec, R2);
        double bott = AA.MultmBn(Rvec, Cvec);
        dlam2 = topp/bott;
        
        /* calculate dup and dq */
        Rmod.CopyIn(0, dQdSig);
        Rmod.CopyIn(dQdSig.Length(), qbar);
        Rmod *= dlam2;
        Rmod += R;
        Rmod *= -1.; 
        AA.Multx(Rmod, X);
        dSig.CopyPart(0, X, 0, dSig.Length());
        dq.CopyPart(0, X, dSig.Length(), dq.Length());
        
        /* update internal variables and plastic multiplier */
        KE_Inv.Multx(dSig,dup);
        dup *= -1.;
        up += dup;
        qn += dq;
        dlam += dlam2;
                                                                                                                                                                                                                                                                                                            
        /*
        // check the local convergence
        if(ip == 0 ) { 
        	cout << "kk=" << kk << " ff=" << ff << " normr=" << normr
             	 << " dlam2=" << dlam2 << " dlam=" << dlam << endl << endl;
        }
        */
        kk++;
        
        /* exit while loop if solution has converged */
        /* convergence criteria: stresses brought back to the yield 
         * surface (i.e. ff ~= 0) and residuals of plastic strains 
         * and internal variables are very small */
        if (fabs(ff) < fTol_1 && normr < fTol_2) {
        	NotConverged = false;
        }
        
      } // while (NotConverged)
    } // if (ff < fTol_1)
    
    /* update state variables */
    state.CopyIn(0, Sig);
    state.CopyIn(Sig.Length(), totalstrain_curr); 	   
	state.CopyIn(12, up);
	state.CopyIn(18, qn);
	state[20] = ff;
	state[21] = dlam;
	state[22] = double(iplastic);
	state[23] = normr;
	state[24] = double(kk);
	
	fStressCorr = Sig;
	if (iplastic == 1) {
   		/* indicator during the first elastic to plastic transition
   		 * use qn from the fInternal(i.e. plastic) and not elastic qn in second plastic step */ 
   		state[25] = 1.; 
   		
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
  double ftan_phi = qn[0]; 
  double p = Sig.Trace()/3.0;
  
  Sig_Dev.Deviatoric(Sig);
  double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
  double ff = q + (ftan_phi*p) - fc;
  return ff;
}

/* calculation of dfdSig_f */
void SMRSSNLHardT::dfdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
   dSymMatrixT Sig_Dev(3);
   double ftan_phi = qn[0]; 
   
   Sig_Dev.Deviatoric(Sig);
   double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
   Sig_Dev *= ratio32/q;
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
   double p = Sig.Trace()/3.0;
   dfdq = 0.0;
   dfdq[0] = p;
}

/* calculation of dQdSig_f */
void SMRSSNLHardT::dQdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
   dSymMatrixT Sig_Dev(3);
   double ftan_psi = qn[1];
   
   Sig_Dev.Deviatoric(Sig);
   double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
   Sig_Dev *= ratio32/q;
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
  dSymMatrixT Sig_Dev(3), Ones(3); 
  dMatrixT I_mat(6),I(3),matrix2(6);
  double ftan_psi = qn[1];
  
  Sig_Dev.Deviatoric(Sig);
  double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
  Ones.Identity();
  I_mat.Outer(Ones,Ones);
  I_mat /= 3.0;
  dQdSig2=0.0;
  dQdSig2.Identity(); 
  dQdSig2 -= I_mat;
  dQdSig2 /= q; 
  matrix2.Outer(Sig_Dev,Sig_Dev);
  matrix2 *= ratio32/(q*q*q);
  dQdSig2 -= matrix2;
  dQdSig2 *= ratio32;
}

/* calculation of dQdSigdq_f */
void SMRSSNLHardT::dQdSigdq_f(dMatrixT& dQdSigdq)
{
  dQdSigdq = 0.;
  for (int i = 0; i < 3; i++)
  	dQdSigdq(i,1) = 1.0/3.0;
}

/* calculation of qbar_f */
void SMRSSNLHardT::qbar_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& qbar)
{
   dSymMatrixT Sig_Dev(3),B2(3),B3(3),dQdS(3); 
   double ftan_phi = qn[0];
   double ftan_psi = qn[1]; 
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
   dQdS.SetToScaled(ratio32/q,Sig_Dev);
   B3.SetToScaled(1.0/fGf_II,Sig_Dev);
   qbar = 0.0;   
   qbar[0]  = A3*B3.ScalarProduct(dQdS); 
   qbar[1]  = A4*B3.ScalarProduct(dQdS); 
 }
 
/* calculation of dqbardSig_f */
void SMRSSNLHardT::dqbardSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
   dSymMatrixT dhtanphi_dSig(3),dhtanpsi_dSig(3),Sig_Dev(3), 
               B2(3),B3(3),dQdS(3),dB3dS_dQdS(3),B3_dQdQ_dSdS(3);
   dMatrixT tempmat(6);
   
   double ftan_phi = qn[0];
   double ftan_psi = qn[1]; 
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
   dQdS.SetToScaled(ratio32/q,Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   tempmat.Outer(Sig_Dev,Sig_Dev);
   //tempmat.Multx(Sig_Dev, B3_dQdQ_dSdS);
   Contract4To2(tempmat,Sig_Dev,B3_dQdQ_dSdS);
   B3_dQdQ_dSdS *= -ratio32/(q*q*q);
   B3_dQdQ_dSdS.AddScaled(1./q, Sig_Dev);
   B3_dQdQ_dSdS *= ratio32/fGf_II;
   dB3dS_dQdS.SetToScaled(1.0/fGf_II, dQdS);
   dhtanphi_dSig  = B3_dQdQ_dSdS;
   dhtanphi_dSig += dB3dS_dQdS;
   dhtanphi_dSig *= A3;
   dhtanpsi_dSig  = B3_dQdQ_dSdS;
   dhtanpsi_dSig += dB3dS_dQdS;
   dhtanpsi_dSig *= A4;
   
   dqbardSig = 0.0;
   dqbardSig(0,0) = dhtanphi_dSig[0];
   dqbardSig(0,1) = dhtanphi_dSig[1];
   dqbardSig(0,2) = dhtanphi_dSig[2];
   dqbardSig(0,3) = dhtanphi_dSig[3];
   dqbardSig(0,4) = dhtanphi_dSig[4];
   dqbardSig(0,5) = dhtanphi_dSig[5];
   dqbardSig(1,0) = dhtanpsi_dSig[0];
   dqbardSig(1,1) = dhtanpsi_dSig[1];
   dqbardSig(1,2) = dhtanpsi_dSig[2];
   dqbardSig(1,3) = dhtanpsi_dSig[3];
   dqbardSig(1,4) = dhtanpsi_dSig[4];
   dqbardSig(1,5) = dhtanpsi_dSig[5];
}
  
/* calculation of dqbardq_f */
void SMRSSNLHardT::dqbardq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{
   dSymMatrixT Sig_Dev(3),B2(3),B3(3),dQdS(3);
   
   double ftan_phi = qn[0];
   double ftan_psi = qn[1]; 
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double q  = sqrt((Sig_Dev.ScalarProduct())*ratio32);
   dQdS.SetToScaled(ratio32/q,Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   double B3dQdS = B3.ScalarProduct(dQdS); 
   dqbardq = 0.0;
   dqbardq(0,0)  = -falpha_phi*B3dQdS;
   dqbardq(1,1)  = -falpha_psi*B3dQdS;
}

/* return the consistent elastoplastic moduli 
 *
 * Note: Return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& SMRSSNLHardT::Moduli(const ElementCardT& element, 
	int ip)
{
	  /* allocate matrices */
     dMatrixT KE(6),AA(8),AA_inv(8),TT(8),KE_Inv(6), 
              Auu_inv(6),Auq_inv(6,2),Aqu_inv(2,6),Aqq_inv(2),
              dQdSig2(6),dqbardq(2),dQdSigdq(6,2),dqbardSig(2,6),
              KP(6),KEP(6);
     
     /* allocate reduced index vector of symmetric matrices */
     dSymMatrixT Sig(3),dfdSig(3),dQdSig(3);  
     
     /* allocate vectors */   
     dArrayT dfdq(2),qbar(2),qn(2),Rvec(8),Cvec(8),
             vec1(8),vec2(8);
     
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
    	double dlam = fInternal[kdlambda];
		KE_Inv.Inverse(KE);
		
		
	    /* calculate C_EPC */
        /* form AA_inv matrix and calculate AA */
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
	
        dfdSig_f(Sig, qn, dfdSig);
        dfdq_f(Sig, dfdq);
        dQdSig_f(Sig, qn, dQdSig);
        qbar_f(Sig, qn, qbar);
        	
        /* include contribution of all off diagonal terms 
         * in reduced vectors of symmetric matrices 
         * dfdSig and dQdSig  */
        for (int i=0; i<3; i++) {
   			 dfdSig[i+3] *= 2.;
   			 dQdSig[i+3] *= 2.;
   		}
   			
        Rvec.CopyIn(0, dfdSig);
        Rvec.CopyIn(dfdSig.Length(), dfdq);
        Cvec.CopyIn(0, dQdSig);
        Cvec.CopyIn(dQdSig.Length(), qbar);
        double H = AA.MultmBn(Rvec, Cvec); /* H (scalar) */
        AA.Multx(Cvec, vec1);
        AA.MultTx(Rvec, vec2);
        TT.Outer(vec1, vec2);
        TT /= H;
        AA -= TT;
        for (int i=0; i<6; i++)
        	for (int j=0; j<6; j++)
        	    KEP(i,j) = AA(i,j);
        	
        /*************************************************************/
        /* continuum jacobian, i.e., "inconsistent" tangent operator */
   		dArrayT KE_dfdSig(6),KE_dQdSig(6);
   		KE.Multx(dfdSig, KE_dfdSig);
   		KE.Multx(dQdSig, KE_dQdSig);
   		KP.Outer(KE_dfdSig, KE_dQdSig);
   		double bott = KE.MultmBn(dfdSig,dQdSig);
   		bott -= dArrayT::Dot(dfdq,qbar);
   		KP /= bott;
   		//KEP.DiffOf(KE,KP); //uncomment to activate continuum jacobian
   		/*************************************************************/
   			
	    fModuli = KEP;
	    
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
			       MeanStress(trialstrain,element),ftan_phi,fc) > fTol_1 );
	}
        /* already plastic */
	else 
	{
	/* get flags */
	 iArrayT& Flags = element.IntegerData();
		
	/* load internal variables */
	LoadData(element, ip);
	
	if(fInternal[25] == 0.) {
		double esp = 0.;
    	fInternal[ktanphi] = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
		fInternal[ktanpsi] = (tan(fphi_p))*exp(-falpha_psi*esp);
	}
	
	dSymMatrixT elasticstrain(3);
	elasticstrain.DiffOf(trialstrain, fPlasticStrain);
	fInternal[kftrial] = YieldCondition(DeviatoricStress(elasticstrain,element),
			             MeanStress(elasticstrain,element),fInternal[ktanphi],fc);

		/* plastic */
		if (fInternal[kftrial] > fTol_1)
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

/* computes the stress corresponding to the given element
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

double SMRSSNLHardT::signof(double r)
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
