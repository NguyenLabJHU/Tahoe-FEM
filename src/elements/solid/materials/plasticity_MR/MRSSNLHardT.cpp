/* $Id: MRSSNLHardT.cpp,v 1.16 2006-08-22 14:39:17 kyonten Exp $ */
/* created: Karma Yonten */

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

const int    kNumInternal = 31; // number of internal state variables
const int    kNSD         = 3;
const int    kNSTR        = dSymMatrixT::NumValues(kNSD);
const double ratio32      = 3.0/2.0;

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
      const dSymMatrixT& totalstrain_curr, ElementCardT& element, int ip)
{
  	int kk,iplastic;
    double ff,dlam,dlam2,normr;
  	
	/* allocate matrices */
    dMatrixT KE(6),AA(10),AA_inv(10),KE_Inv(6),CMAT(10),  
             Auu_inv(6),Auq_inv(6,4),Aqu_inv(4,6),Aqq_inv(4),
             dQdSig2(6),dQdSigdq(6,4),dqbardq(4),dqbardSig(4,6);
    
	/* allocate reduced index vector of symmetric matrices */
    dSymMatrixT u(3),up(3),upo(3),du(3),dup(3),ue(3),
                Sig(3),dSig(3),Sig_e(3);
    
    /* allocate vectors */
    dArrayT Rvec(10),Cvec(10),R(10),Rmod(10),X(10),Y(10),
            qo(4),qn(4),dq(4),R_up(6),R_q(4),R2(10), 
            dfdSig(6),dfdq(4),dQdSig(6),qbar(4),state(31);
	
	/* elastic moduli tensor */
	KE = 0.0;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
        
	if(ip == 0)
		//cout << "ip "<< ip << endl; // ip #
		
	/* case 1: element/ip elastic */
	if (!PlasticLoading(totalstrain_curr, element, ip) && 
	    !element.IsAllocated())
	{
		/* initialize element data */
        double enp  = 0.;
        double esp  = 0.;
        state = 0.;
        state[kchi] = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        state[kc] = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
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
    //KE.MultTx(ue, Sig_e);
    Sig_e.A_ijkl_B_kl(KE, ue); 
    Sig = Sig_e; 
    KE_Inv.Inverse(KE); 
  
/* check the yield function */
    ff = Yield_f(Sig, qn);
    if (ff < fTol_1) {
      iplastic = 0;
      kk = 0;
    } 
    else {
      state[27] = ff;
      kk = 0;
      iplastic = 1;
      bool NotConverged = true;
      while (NotConverged) 
      {
        if (kk > 20)
        	ExceptionT::GeneralFail("MRSSNLHardT::StressCorrection","Too Many Iterations");
        
        /* calculate stress */
        ue.DiffOf(totalstrain_curr, up);
    	//KE.MultTx(ue, Sig_e);
    	Sig_e.A_ijkl_B_kl(KE,ue);
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
        dQdSigdq_f(Sig, qn, dQdSigdq);
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
        dfdq_f(Sig, qn, dfdq);
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
        if (abs(ff) < fTol_1 && normr < fTol_2) {
        	NotConverged = false;
        }
      } // while (NotConverged)
    } // if (ff < fTol_1)
    
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
double MRSSNLHardT::Yield_f(const dSymMatrixT& Sig, 
			const dArrayT& qn)
{
  dSymMatrixT Sig_Dev(3);
  double fchi = qn[0];
  double fc = qn[1];
  double ffriction = qn[2]; 
  double fpress = Sig.Trace()/3.0;
  
  Sig_Dev.Deviatoric(Sig);
  double temp  = (Sig_Dev.ScalarProduct())/2.0;
  double temp2 = fc - ffriction*fchi;
  double temp3 = temp2 * temp2;
  temp += temp3;
  double ff = sqrt(temp)-(fc - ffriction*fpress); 
  return ff;
}

/* calculation of dfdSig_f */
void MRSSNLHardT::dfdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
   dSymMatrixT Sig_Dev(3);
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2]; 
   double Sig_p = Sig.Trace()/3.0;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double temp2 = fc - ftan_phi*fchi;
   double temp3 = temp2 * temp2;
   temp += temp3;
   double deno = 2.0*sqrt(temp);
   Sig_Dev /= deno;
   dfdSig=0.0;
   for (int i = 0; i < 6; i++) 
   		if (i < 3)  
   			dfdSig[i] = Sig_Dev[i] + (ftan_phi/3.0);
   		else 
   			dfdSig[i] = Sig_Dev[i];
}

/* calculation of dfdq_f */
void MRSSNLHardT::dfdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
   dSymMatrixT Sig_Dev(3);
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double Sig_p = Sig.Trace()/3.0;
   
   double const1 = fc-fchi*ftan_phi; 
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   temp += (const1*const1);
   double deno = sqrt(temp);
   dfdq=0.0;
   dfdq[0] = -ftan_phi*const1/deno;
   dfdq[1] = const1/deno - 1.0;
   dfdq[2] = -fchi*const1/deno + Sig_p;
}

/* calculation of dQdSig_f */
void MRSSNLHardT::dQdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
   dSymMatrixT Sig_Dev(3);
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_psi = qn[3];
   double Sig_p = Sig.Trace()/3.0;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double temp2 = fc - ftan_psi*fchi;
   double temp3 = temp2 * temp2;
   temp += temp3;
   double deno = 2.0*sqrt(temp);
   Sig_Dev /= deno;
   dQdSig=0.0;
   for (int i = 0; i < 6; i++)
   		if (i < 3)  
   			dQdSig[i] = Sig_Dev[i] + (ftan_psi/3.0);
   		else 
   			dQdSig[i] = Sig_Dev[i];
}

/* calculation of dQdSig2_f */
void MRSSNLHardT::dQdSig2_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSig2)
{
  dSymMatrixT Sig_Dev(3),Ones(3); 
  dMatrixT I_mat(6),matrix2(6);
  double fchi = qn[0];
  double fc = qn[1];
  double ftan_psi = qn[3];
  double Sig_p = Sig.Trace()/3.0;
  
  Sig_Dev.Deviatoric(Sig);
  double temp  = (Sig_Dev.ScalarProduct())/2.0;
  double temp2 = fc - ftan_psi*fchi;
  double temp3 = temp2 * temp2;
  temp += temp3;
  double deno = 2.0*sqrt(temp);
  Ones.Identity();
  I_mat.Outer(Ones,Ones);
  I_mat /= 3.0;
  dQdSig2 = 0.0;
  dQdSig2.PlusIdentity();
  dQdSig2 -= I_mat;
  dQdSig2 /= deno;   
  matrix2.Outer(Sig_Dev,Sig_Dev);
  matrix2 /= 4.0*pow(temp, ratio32);
  
  dQdSig2 -= matrix2;
}

/* calculation of dQdSigdq_f */
void MRSSNLHardT::dQdSigdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
  dSymMatrixT Sig_Dev(3);
  dMatrixT col_1(6,1),col_2(6,1),col_4(6,1);
  double fchi = qn[0];
  double fc = qn[1];
  double ftan_psi = qn[3];
  double Sig_p = Sig.Trace()/3.0;
  double const1 = fc - ftan_psi*fchi;
  
  Sig_Dev.Deviatoric(Sig);
  double temp  = (Sig_Dev.ScalarProduct())/2.0;
  temp += (const1*const1);
  const1 /= 2.0*pow(temp, ratio32); 
       
  col_1.SetToScaled(ftan_psi*const1, Sig_Dev);
  col_2.SetToScaled(-1.0*const1, Sig_Dev);
  col_4.SetToScaled(fchi*const1, Sig_Dev);
  for (int i=0; i<3; i++)
     col_4(i,0)/=3.;     
  dQdSigdq = 0.;
  dQdSigdq.AddBlock(0, 0, col_1);
  dQdSigdq.AddBlock(0, 1, col_2); 
  dQdSigdq.AddBlock(0, 3, col_4);
}

/* calculation of qbar_f */
void MRSSNLHardT::qbar_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& qbar)
{
   dSymMatrixT Sig_Dev(3),B2(3),B3(3),dQdS(3); 

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
   double dQdP = ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double temp2 = fc - ftan_psi*fchi;
   double temp3 = temp2 * temp2;
   temp += temp3;
   double deno = 2.0*sqrt(temp);
   B2 = Sig_Dev;
   B2 /= fGf_I;
   dQdS = Sig_Dev;
   dQdS /= deno;
   B3 = Sig_Dev;
   B3 /= fGf_II; 
   qbar=0.0;   
   qbar[0]  = A1*B1*dQdP; 
   qbar[0] += A1*B2.ScalarProduct(dQdS); 
   qbar[1]  = A2*B3.ScalarProduct(dQdS); 
   qbar[2]  = A3*B3.ScalarProduct(dQdS); 
   qbar[3]  = A4*B3.ScalarProduct(dQdS); 
 }
 
/* calculation of dqbardSig_f */
void MRSSNLHardT::dqbardSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
   dSymMatrixT dhchi_dSig(3),dhc_dSig(3),dhtanphi_dSig(3),dhtanpsi_dSig(3),
               Sig_Dev(3),dQdS(3),II_S(3),SxS_S(3),tmp1(3),tmp2(3);
   dMatrixT SxS(6,6),II(6,6);
   
   II=0.;
   for (int i=0; i<6; i++)
       II(i,i)=1.;
               
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
   double dQdP = ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double temp2 = fc - ftan_psi*fchi;
   double temp3 = temp2 * temp2;
   temp += temp3;
   double deno1 = 2.0*sqrt(temp);
   double deno2 = 4.0*pow(temp, ratio32);
   dQdS = Sig_Dev;
   dQdS /= deno1;
   
   II.Multx(Sig_Dev,II_S);
   II_S /= deno1/2.;
   SxS.Outer(Sig_Dev,Sig_Dev);
   //SxS.Multx(Sig_Dev,SxS_S);
   Contract4To2(SxS,Sig_Dev,SxS_S);
   SxS_S /= deno2;
   
   tmp1 = II_S;
   tmp1 -= SxS_S;
   tmp2 = tmp1;
   tmp1 /= fGf_I;
   tmp2 /= fGf_II;
   
   dhchi_dSig = 0.0;
   dhchi_dSig = tmp1;
   for (int i=0; i<3; i++)
       dhchi_dSig[i] += dQdP*dB1dP/3.;
   dhchi_dSig *= A1; 
   dhc_dSig = tmp2;
   dhc_dSig *= A2;
   dhtanphi_dSig = tmp2;
   dhtanphi_dSig *= A3;
   dhtanpsi_dSig = tmp2;
   dhtanpsi_dSig *= A4;
  
   dqbardSig=0.0;
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
   dSymMatrixT Sig_Dev(3),B2(3),B3(3),dQdS(3);
   
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
   double dQdP = ftan_psi;
   double const1 = fc - ftan_psi*fchi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   temp += (const1*const1);
   double deno = 2.0*sqrt(temp);
   const1 /= 2.0*pow(temp, ratio32);
   dQdS = Sig_Dev;
   dQdS /= deno;
   B2.SetToScaled(1.0/fGf_I, Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   double B2S = B2.ScalarProduct(Sig_Dev); 
   double B3S = B3.ScalarProduct(Sig_Dev); 
   double B2dQdS = B2.ScalarProduct(dQdS); 
   double B3dQdS = B3.ScalarProduct(dQdS); 
   dqbardq = 0.0;
   dqbardq(0,0)  = -falpha_chi*(B1*dQdP + B2dQdS);
   dqbardq(0,0) += A1*ftan_psi*const1*B2S;
   dqbardq(0,1)  =  -A1*const1*B2S;
   dqbardq(0,3)  =  A1*B1+A1*fchi*const1*B2S;
   dqbardq(1,0)  = A3*ftan_psi*const1*B3S;   
   dqbardq(1,1)  = -falpha_c*B3dQdS;
   dqbardq(1,1) -= A2*const1*B3S;
   dqbardq(1,3)  = A2*fchi*const1*B3S;
   dqbardq(2,0)  = A3*ftan_psi*const1*B3S;
   dqbardq(2,1)  = -A3*const1*B3S;
   dqbardq(2,2)  = -falpha_phi*B3dQdS;
   dqbardq(2,3)  = A3*fchi*const1*B3S;
   dqbardq(3,0)  = A4*ftan_psi*const1*B3S;
   dqbardq(3,1)  = -A4*const1*B3S;
   dqbardq(3,3)  = -falpha_psi*B3dQdS;
   dqbardq(3,3) += A4*fchi*const1*B3S;
}

/* return the consistent elastoplastic moduli 
 *
 * Note: Return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& MRSSNLHardT::Moduli(const ElementCardT& element, 
	int ip)
{
	  /* allocate matrices */
     dMatrixT KE(6),AA(10),AA_inv(10),TT(10),KE_Inv(6), 
              Auu_inv(6),Auq_inv(6,4),Aqu_inv(4,6),Aqq_inv(4),
              dQdSig2(6),dqbardq(4),dQdSigdq(6,4),dqbardSig(4,6),
              KP(6),KEP(6);
     
     /* allocate reduced index vector of symmetric matrices */
     dSymMatrixT Sig(3),dfdSig(3),dQdSig(3);  
     
     /* allocate vectors */   
     dArrayT dfdq(4),qbar(4),qn(4),Rvec(10),Cvec(10),
             vec1(10),vec2(10);
     
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
        dQdSigdq_f(Sig, qn, dQdSigdq);
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
        dfdq_f(Sig, qn, dfdq);
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
	//d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
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
	// updating done in MRSSNLHardT::StressCorrection
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
	fInternal.Alias(kNumInternal, &d_array[offset + ip*kNumInternal]); //2*offset if fUnitNorm
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int MRSSNLHardT::PlasticLoading(const dSymMatrixT& trialstrain, 
	  ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated()) {
	    double enp  = 0.;
        double esp  = 0.;
        double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
        double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
        double ftan_psi = (tan(fphi_p))*exp(-falpha_psi*esp);
		return( YieldCondition(DeviatoricStress(trialstrain,element),
			       MeanStress(trialstrain,element),fchi,fc,ftan_phi) > fTol_1 );
	}
        /* already plastic */
	else 
	{
	/* get flags */
	 iArrayT& Flags = element.IntegerData();
		
	/* load internal variables */
	LoadData(element, ip);
	
    if(fInternal[30] == 0.) {
		double enp  = 0.;
		double esp = 0.;
    	fInternal[kchi] = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        fInternal[kc]   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
    	fInternal[ktanphi] = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
		fInternal[ktanpsi] = (tan(fphi_p))*exp(-falpha_psi*esp);
	}
	
	dSymMatrixT elasticstrain(3);
	elasticstrain.DiffOf(trialstrain, fPlasticStrain);
	fInternal[kftrial] = YieldCondition(DeviatoricStress(elasticstrain,element),
			             MeanStress(elasticstrain,element),fInternal[kchi],
			             fInternal[kc],fInternal[ktanphi]);
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

double MRSSNLHardT::signof(double r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

void MRSSNLHardT::Contract4To2(const dMatrixT& mat,const dSymMatrixT& vec,
                                      dSymMatrixT& res)
{
	dSymMatrixT tmp(3);
	tmp = vec;
	for (int i=0; i<3; i++) {
	   tmp[i+3]=vec[i+3] *2.;
	}
	mat.Multx(tmp, res);
}