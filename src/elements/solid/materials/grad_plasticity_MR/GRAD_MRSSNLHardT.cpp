/* $Id: GRAD_MRSSNLHardT.cpp,v 1.17 2005-09-06 15:02:51 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   Gradient Enhanced MR Model
*/

/* interface for a nonassociative, small strain,      */
/* pressure dependent gradient plasticity model       */
/* with nonlinear isotropic hardening/softening.      */

#include "GRAD_MRSSNLHardT.h"
#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "LAdMatrixT.h" //for solving Ax=b linear system of equations

/* class constants */

using namespace Tahoe;

const int    kNumInternal = 40; // number of internal state variables
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;
//const int    kNSTR        = dSymMatrixT::NumValues(kNSD);

/* constructor */
GRAD_MRSSNLHardT::GRAD_MRSSNLHardT(int num_ip, double mu, double lambda):
    fNumIP(num_ip),
	fmu(mu),
	flambda(lambda),
	fkappa(flambda + (2.0/3.0*fmu)),
	fmu_ast(flse_s*flse_s*mu),
	flambda_ast(flse_v*flse_v*fkappa-2.0/3.0*fmu_ast),
	fkappa_ast(flambda_ast + (2.0/3.0*fmu_ast))
{
	SetName("GRAD_MR_SS_nonlinear_hardening");
}
const dSymMatrixT& GRAD_MRSSNLHardT::ElasticStrain(const dSymMatrixT& totalstrain, 
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

const dSymMatrixT& GRAD_MRSSNLHardT::LapElasticStrain(const dSymMatrixT& lap_totalstrain, 
	const ElementCardT& element, int ip) 
	
{
	/* remove plastic strain */
	if (element.IsAllocated()) 
	{
		/* load internal variables */
		LoadData(element, ip);

	    /* compute gradient elastic strain */
		/*fLapElasticStrain.DiffOf(lap_totalstrain, fLapPlasticStrain);*/
		fLapElasticStrain = lap_totalstrain;
		return fLapElasticStrain; 
	}	
	/* no plastic strain */
	else	
		return lap_totalstrain; 
}

/* return correction to stress vector computed by mapping the
 * stress back to the yield surface, if needed */
const dSymMatrixT& GRAD_MRSSNLHardT::StressCorrection(const dSymMatrixT& trialstrain, 
                  const dSymMatrixT& lap_trialstrain, const dArrayT& triallambda, const dArrayT& lap_triallambda,
                  ElementCardT& element, int ip)
{	
  	int kk, PLI;
  	double ff;

    /* define and allocate matrices */
    dMatrixT KE(6,6); dMatrixT KE_AST(6,6);  
    dMatrixT dhdSig(4,6); dMatrixT dhdq(4,4); dMatrixT dhdm(4,6);
    dMatrixT dgdSig(4,6); dMatrixT dgdq(4,4);
    dMatrixT dmdSig(6,6); dMatrixT dmdq(6,4);
    dMatrixT dRSig_dSig(6,6); dMatrixT dRSig_dq(6,4); 
    dMatrixT dRq_dSig(4,6); dMatrixT dRq_dq(4,4);  dMatrixT Y(6,6); 
    LAdMatrixT dLHS(10); // 10x10: set up matrix A for solving Ax=b system
     
    /* define and allocate vectors */
    dArrayT u(6); dArrayT up(6); dArrayT du(6); dArrayT dup(6); dArrayT upo(6); 
    dArrayT lap_u(6); dArrayT lap_up(6); dArrayT lap_du(6); dArrayT lap_dup(6); 
    dArrayT lap_upo(6);  
    dArrayT qn(4); dArrayT qo(4);
    dArrayT Sig(6); dArrayT dSig(6); dArrayT Sig_I(6);
    dArrayT mm(6); dArrayT rr(4); dArrayT nn(6); 
    dArrayT dq(4); dArrayT hh(4); dArrayT gg(4); 
    dArrayT state(40); dArrayT Sig_trial(6);
    dArrayT dRHS(10); dArrayT ls(4);
    dArrayT RSig(6); dArrayT Rq(4);
    
    /* initialize */
    dRHS = 0.0; dLHS = 0.0;
    KE = 0.; KE_AST = 0.;
    fIVFlag = 0; // elastic or plastic but zero lambda value
    fIniInternal = 0.; 
    
	/* length scales parameters */
	ls[0] = flse_v; // lse_v: pore space length scale (elastic)  
	ls[1] = flse_s; // lse_s: grain size length scale (elastic)
	ls[2] = flsp_v; // lsp_v: pore space length scale (plastic)   
	ls[3] = flsp_s; // lsp_s: grain size length scale (plastic) 
	
	/* C and C_AST matrices */
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
	KE_AST(2,2) = KE_AST(1,1) = KE_AST(0,0) = flambda_ast + 2.0*fmu_ast;
	KE_AST(1,2) = KE_AST(0,1) = KE_AST(0,2) = flambda_ast;
	KE_AST(2,1) = KE_AST(1,0) = KE_AST(2,0) = flambda_ast;
	KE_AST(5,5) = KE_AST(4,4) = KE_AST(3,3) = fmu_ast;
	 
    /* take out the symmetric part of the strain tensor
       and laplacian of strain tensor   */
	for (int ij = 0; ij < dSymMatrixT::NumValues(kNSD); ij++)
	{
		int i, j;
		dSymMatrixT::ExpandIndex(kNSD, ij, i, j);
		u[ij] = trialstrain(i,j);
		lap_u[ij] = lap_trialstrain(i,j);
	} 
        
	PLI = PlasticLoading(trialstrain, lap_trialstrain, element, ip);
	
	if (PlasticLoading(trialstrain, lap_trialstrain, element, ip) && 
	    element.IsAllocated()) 
    {
	  LoadData(element, ip);
	  state.CopyIn(0, fInternal);
	}
    
	if (!PlasticLoading(trialstrain, lap_trialstrain, element, ip) && 
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
        state[30] = fchi;
        state[31] = fc;
        state[32] = ftan_phi;
        state[33] = ftan_psi;
	}
	
	if (!PlasticLoading(trialstrain, lap_trialstrain, element, ip) && 
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
        state[30] = fchi;
        state[31] = fc;
        state[32] = ftan_phi;
        state[33] = ftan_psi;
	}
	
	/* initialize in the case of first plastic loading  */
	/* check consistency and initialize plastic element */
	if (PlasticLoading(trialstrain, lap_trialstrain, element, ip) && 
	    !element.IsAllocated())
	{
		/* new plastic element */
		AllocateElement(element); 
		PlasticLoading(trialstrain, lap_trialstrain, element, ip); 
		/* initialize element data */
		double enp  = 0.;
        double esp  = 0.;
        fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
        double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
        double ftan_psi = (tan(fphi_p))*exp(-falpha_psi*esp);
        state = 0.;
        state[30] = fchi;
        state[31] = fc;
        state[32] = ftan_phi;
        state[33] = ftan_psi;
	}
	
	/* calculate incremental strains and initialize the necessary vectors */
    for (int i = 0; i < 6; i++) 
    {
       du[i] = u[i] - state[i+6];
       lap_du[i] = lap_u[i] - state[i+12]; //laplacian of du
       up[i] = state[i+18];	
       lap_up[i] = state[i+24];	//laplacian of up
       upo[i] = up[i];
       lap_upo[i] = lap_up[i];
       Sig_I[i] = 0.;
       if (i < 4) {
       		qn[i] = state[i+30];
        	qo[i] = qn[i];
       }
    }
     
    Sig = Sig_I;
    dArrayT ue(6), lap_ue(6), Sig_e(6), lap_Sig_e(6); 
    ue = u;
    ue -= up;
    lap_ue = lap_u;
    lap_ue -= lap_up;
    KE.Multx(ue, Sig_e);
    KE_AST.Multx(lap_ue, lap_Sig_e);
    Sig += Sig_e; 
    Sig -= lap_Sig_e;
    //cout << endl << "sig_e =" << endl << Sig_e << endl;
    //cout << endl << "lap_sig_e =" << endl << lap_Sig_e << endl;
    Sig_trial = Sig;
    
    int iplastic; 
    double dlam = triallambda[0]; double lap_dlam = lap_triallambda[0];
 
/* check the yield function */
	//cout << endl << "ip = " << ip << endl;
    //cout << endl << "sig = " << endl << Sig << endl;
    //cout << endl << "qn = " << endl << qn << endl;
    
    Yield_f(Sig, qn, ff);
    
    /* if yield condition is satisfied, but lambda is zero
       skip plastic part, and send back current stresses and 
       internal variables to the global level
    */ 
    if (ff < kYieldTol) 
    {
      iplastic = 0;
      state[34] = ff;
      state[39] = ff;
      kk = 0;
    }
      
    else  
    {
      cout << endl << " ip yield condition satisfied!!! " << ff << endl;
      state[39] = ff;
      iplastic = 1;
      if (dlam > kYieldTol) fIVFlag = 1;
      kk = 0;
      while (ff > fTol_1) 
      {
        if (kk > 500) 
        	ExceptionT::GeneralFail("GRAD_MRSSNLHardT::StressCorrection","Too Many Iterations");
        
        /* terminate iteration since lambda = 0 will run into convergence problem */ 
        if (fIVFlag == 0) break;
         
        Sig = Sig_I;
        ue = u;
        ue -= up;
        lap_ue = lap_u;
        lap_ue -= lap_up;
        KE.Multx(ue, Sig_e);
        KE_AST.Multx(lap_ue, lap_Sig_e);
        Sig += Sig_e; 
        Sig -= lap_Sig_e;
        
        Yield_f(Sig, qn, ff);  //yield function
        
        /* calculate all the necessary derivatives */
        m_f(Sig, qn, mm); dmdSig_f(qn, dmdSig); dmdq_f(Sig, qn, dmdq);
        h_f(Sig, qn, hh); dhdSig_f(Sig, qn, dhdSig); dhdq_f(Sig, qn, dhdq); dhdm_f(Sig, qn, dhdm);  
        g_f(Sig, qn, ls, gg); dgdSig_f(Sig, qn, ls, dgdSig); dgdq_f(Sig, qn, ls, dgdq);   
        
        /* calculate R_Sig_Sig and R_Sig_q matrices */
        dMatrixT RSigMat1(6,6), RSigMat2(6,6); /* work space */ 
        RSigMat1.SetToScaled(dlam, KE);       
        RSigMat2.SetToScaled(lap_dlam, KE_AST);
        RSigMat1 -= RSigMat2;  
        dRSig_dSig.MultAB(RSigMat1, dmdSig);
        dRSig_dSig += IdentityMatrix6;
         
        dRSig_dq.MultAB(RSigMat1, dmdq);
        
        /* calculate R_q_Sig and R_q_q matrices */
        dMatrixT RqMat1(4,6), RqMat2(4,4); /* work space */
        dRq_dSig.MultAB(dhdm, dmdSig);
        dRq_dSig += dhdSig;
        dRq_dSig *= -dlam;
        RqMat1.SetToScaled(lap_dlam, dgdSig);
        dRq_dSig += RqMat1; 
        
        dRq_dq.MultAB(dhdm, dmdq);
        dRq_dq += dhdq;
        dRq_dq *= -dlam;
        RqMat2.SetToScaled(lap_dlam, dgdq);
        dRq_dq += RqMat2;
        dRq_dq += IdentityMatrix4;
         
        /* calculate R_Sig vector */
         /* work space */
        dArrayT RSigTemp1(6), RSigTemp2(6), RSigTemp3(6), RSigTemp4(6);
        dArrayT RSigTemp5(6), RSigTemp6(6);
        KE.Multx(du, RSigTemp1);
        KE_AST.Multx(lap_du, RSigTemp2);
        RSig.DiffOf(RSigTemp2, RSigTemp1); // RSigTemp2 - RSigTemp1
        KE.Multx(mm, RSigTemp3);  
        KE_AST.Multx(mm, RSigTemp4);
        RSigTemp5.SetToScaled(dlam, RSigTemp3); // dlam*RSigTemp3;
        RSig += RSigTemp5;				   // dlam2 == dlam??
        RSigTemp6.SetToScaled(lap_dlam, RSigTemp4); // lap_dlam*RSigTemp4;
        RSig -= RSigTemp6;
        
        /* calculate R_q vector */						//lap_dlam2 == lap_dlam?? 
        dArrayT Rq_temp(4); /* work space */
        Rq.SetToScaled(lap_dlam, gg);  //lap_dlam*gg;
        Rq_temp.SetToScaled(dlam, hh);  //dlam2 == dlam??
        Rq -= Rq_temp; 
        
        /* check the matrices are correct before assemblying */
        /*if (kk == 500) {
        	cout << "dRSig_dSig =" << endl << dRSig_dSig << endl;
        	cout << "dRSig_dq =" << endl << dRSig_dq << endl;
        	cout << "dRq_dSig =" << endl << dRq_dSig << endl;
        	cout << "dRq_dq =" << endl << dRq_dq << endl;
        	cout << "RSig =" << endl << RSig << endl;
        	cout << "Rq =" << endl << Rq << endl;
        }*/
        
        /* form dLHS (10,10) matrix */
        dLHS.AddBlock(0,                 0,                 dRSig_dSig);
		dLHS.AddBlock(0,                 dRSig_dSig.Cols(), dRSig_dq);
		dLHS.AddBlock(dRSig_dSig.Rows(), 0,                 dRq_dSig);
		dLHS.AddBlock(dRSig_dSig.Rows(), dRSig_dSig.Cols(), dRq_dq);
		
        /* form dRHS (10) vector */
		dRHS.CopyIn(0, RSig);
		dRHS.CopyIn(RSig.Length(), Rq);
		
		/* check to see dLHS and dRHS are assembled correctly */
		/*if (kk == 500) {
			cout << "dLHS = " << endl << dLHS << endl;
			cout << "dRHS = " << endl << dRHS << endl;
		}*/
        
        /* solve the linear sys. of eqs. */
        dLHS.LinearSolve(dRHS); 
        
        /* extract stress vector and internal variables from dRHS */
        for (int i = 0; i < dRHS.Length(); i++) 
        {
            if (i < 6) 
            	dSig[i] = dRHS[i];
            else 
            	dq[i-6] = dRHS[i];
        }
        
        /* check the solutions */
        /*if (kk == 500) {
        	cout << "dSig = " << endl << dSig << endl;
        	cout << "dq = " << endl << dq << endl;
        }*/
        
       /* update plastic strain, gradient plastic strain 
       	  and internal variables 						*/
        m_f(dSig, qn, mm);
        dup.SetToScaled(dlam, mm);
        up += dup;
        lap_up.SetToScaled(lap_dlam, mm);  //lap_dlam2 == lap_dlam??
        lap_up += lap_dup;
        Sig += dSig;   // stress automatically updated when up & lap_up are updated??
        qn += dq;
        kk = kk + 1;
      }	// while (ff > fTol_1) 
    }	// if (ff <kYieldTol)
    
    /* update state variables */
    state.CopyIn(0, Sig); 
	state[6] = trialstrain(0,0);
	state[7] = trialstrain(1,1);
	state[8] = trialstrain(2,2);
	state[9] = trialstrain(1,2);
	state[10] = trialstrain(0,2);
	state[11] = trialstrain(0,1);
	state[12] = lap_trialstrain(0,0);
	state[13] = lap_trialstrain(1,1);
	state[14] = lap_trialstrain(2,2);
	state[15] = lap_trialstrain(1,2);
	state[16] = lap_trialstrain(0,2);
	state[17] = lap_trialstrain(0,1);
	state.CopyIn(18, up);
	state.CopyIn(24, lap_up);
	state.CopyIn(30, qn);
	state[34] = ff; 
	state[35] = dlam;
	state[36] = lap_dlam;
	state[37] = double(iplastic);
	state[38] = double(kk);
 	
	fYield = state[34];
	for (int i = 0; i < 6; i++) 
		fStressCorr[i] = state[i];
	
	for (int i = 0; i < 4; i++) 
		fIniInternal[i] = state[i+30];
    
	if (iplastic > 0) 
	{
	   for (int i = 0; i < 39; i++) { 
	   		fInternal[i] = state[i];
	   		
	   		/* collect plastic strain and its laplacian */
	   		if (i > 17 && i < 24)
	      		fPlasticStrain[i-18] = state[i];
	      	else if (i > 23 && i < 30)
	        	fLapPlasticStrain[i-24] = state[i];
	   }
	}
				
 return fStressCorr;
}	//end StressCorrection

/*
 * returns the value of the yield function given the
 * stress vector and state variables, where alpha
 * represents isotropic hardening.
 */
double& GRAD_MRSSNLHardT::Yield_f(const dArrayT& Sig, 
			const dArrayT& qn, double& ff)
{
  double kTemp1, kTemp2, kTemp3, kTemp4;
  double fc, fchi, ffriction, fpress;
  dMatrixT devstress(3,3);
  devstress_f(Sig, devstress);
  fpress  = Sig[0]+Sig[1]+Sig[2];
  fpress /=3.;
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

/* returns the deviatoric stress tensor given the stress vector */
dMatrixT& GRAD_MRSSNLHardT::devstress_f(const dArrayT& Sig, dMatrixT& devstress)
{
	double fpress;
  	fpress  = Sig[0]+Sig[1]+Sig[2];
  	fpress /=3.;
  	devstress(0,0) = Sig[0] - fpress;
  	devstress(1,1) = Sig[1] - fpress;
  	devstress(2,2) = Sig[2] - fpress;
  	devstress(1,2) = devstress(2,1) = Sig[3];
  	devstress(0,2) = devstress(2,0) = Sig[4];
  	devstress(0,1) = devstress(1,0) = Sig[5];
  	return devstress;
}


/* calculation of h_f vector */
dArrayT& GRAD_MRSSNLHardT::h_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& hh)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, B2dQdS, B3dQdS;
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3); //need allocation of matrix dimension
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
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
      
   hh[0]  = A1*B1*dQdP; 
   hh[0] += A1*B2dQdS;
   hh[1]  = B3dQdS;
   hh[1] *= A2;
   hh[2]  = B3dQdS;
   hh[2] *= A3;
   hh[3]  = B3dQdS;
   hh[3] *= A4;
   return hh;
 }
 
/* calculation of g_f vector */
dArrayT& GRAD_MRSSNLHardT::g_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& ls, dArrayT& gg)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, B2dQdS, B3dQdS;
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3); // need allocation of matrix dimension
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   B2 = Sig_Dev;
   B2 /= fGf_I;
   B1 *= ls[2]*ls[2];
   B2 *= ls[3]*ls[3];
   dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   dQdS = Sig_Dev;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   B3 *= ls[3]*ls[3]; 
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   B2dQdS = dMatrixT::Dot(B2,dQdS);
   B3dQdS = dMatrixT::Dot(B3,dQdS);
      
   gg[0]  = A1*B1*dQdP; 
   gg[0] += A1*B2dQdS;
   gg[1]  = B3dQdS;
   gg[1] *= A2;
   gg[2]  = B3dQdS;
   gg[2] *= A3;
   gg[3]  = B3dQdS;
   gg[3] *= A4;
   return gg;
 }

/* calculation of dQdSig2_f or dmdSig_f */
dMatrixT& GRAD_MRSSNLHardT::dmdSig_f(const dArrayT& qn, dMatrixT& dmdSig)
{
  double Fac;
  dMatrixT I_mat(6,6);
  
  I_mat = 0.;
  for (int i = 0; i < 3; i++) 
  {
     for (int j = 0; j < 3; j++) 
     	I_mat(i,j) = 1.;
  }
  
  Fac = 2./3.;
  Fac *= qn[3]*qn[3];
  Fac += 1.;
  Fac /= 3.;
  I_mat *= Fac;
  
  dmdSig  = IdentityMatrix6;
  dmdSig -= I_mat;
  return dmdSig;
}

/* calculation of dfdSig_f or n_f*/
dArrayT& GRAD_MRSSNLHardT::n_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& nn)
{
  double Sig_p, temp; 
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   nn[0] = Sig[0] - Sig_p;
   nn[1] = Sig[1] - Sig_p;
   nn[2] = Sig[2] - Sig_p;
   nn[3] = Sig[3];
   nn[4] = Sig[4];
   nn[5] = Sig[5];

   temp = 2./3.;
   temp  *= qn[2];
   temp *= (qn[1] - Sig_p*qn[2]);
   for (int i = 0; i < 3; i++) 
   		nn[i] += temp;
  return nn;
}

/* calculation of dQdSig_f or m_f*/
dArrayT& GRAD_MRSSNLHardT::m_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& mm)
{
  double Sig_p, temp;
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   mm[0] = Sig[0] - Sig_p;
   mm[1] = Sig[1] - Sig_p;
   mm[2] = Sig[2] - Sig_p;
   mm[3] = Sig[3];
   mm[4] = Sig[4];
   mm[5] = Sig[5];

   temp = 2./3.;
   temp  *= qn[3];
   temp *= (qn[1] - Sig_p*qn[3]);
   for (int i = 0; i < 3; i++) 
      mm[i] += temp;
  return mm;
}


/* calculation of dfdq_f or r_f*/
dArrayT& GRAD_MRSSNLHardT::r_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& rr)
{
  double Sig_p;
  Sig_p = (Sig[0] + Sig[1] + Sig[2])/3.;
  rr[0] = -2.*qn[2]*(qn[1]-qn[0]*qn[2]);
  rr[1] = 2.*(Sig_p - qn[0])*qn[2];
  rr[2] = 2.*Sig_p*(qn[1] - Sig_p*qn[2]) - 2.*qn[0]*(qn[1]-qn[0]*qn[2]);
  rr[3] = 0.;
  return rr;
}

/* calculation of dQdSigdq_f or dmdq_f */
dMatrixT& GRAD_MRSSNLHardT::dmdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dmdq)
{
  double Sig_p;
  Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
  dmdq = 0.;
  dmdq(1,0) = 2.*qn[3]/3.;
  dmdq(1,1) = 2.*qn[3]/3.;
  dmdq(1,2) = 2.*qn[3]/3.; 
  dmdq(3,0) = (2.*qn[1] - 4.*Sig_p*qn[3])/3.;
  dmdq(3,1) = (2.*qn[1] - 4.*Sig_p*qn[3])/3.;
  dmdq(3,2) = (2.*qn[1] - 4.*Sig_p*qn[3])/3.;  
  return dmdq;
}

/* calculation of dhdSig_f */
dMatrixT& GRAD_MRSSNLHardT::dhdSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdSig)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, d2QdP2, dB1dP, SN;
   dMatrixT dhchi_dSig(3,3), dhc_dSig(3,3), dhtanphi_dSig(3,3), dhtanpsi_dSig(3,3);
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3), dB2dS_dQdS(3,3), dB3dS_dQdS(3,3);
   dMatrixT tempmat(3,3);
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
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
   
   dhchi_dSig  = IdentityMatrix3;
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
   
   dhdSig(0,0) = dhchi_dSig(0,0);
   dhdSig(0,1) = dhchi_dSig(1,1);
   dhdSig(0,2) = dhchi_dSig(2,2);
   dhdSig(0,3) = dhchi_dSig(1,2);
   dhdSig(0,4) = dhchi_dSig(0,2);
   dhdSig(0,5) = dhchi_dSig(0,1);
   dhdSig(1,0) = dhc_dSig(0,0);
   dhdSig(1,1) = dhc_dSig(1,1);
   dhdSig(1,2) = dhc_dSig(2,2);
   dhdSig(1,3) = dhc_dSig(1,2);
   dhdSig(1,4) = dhc_dSig(0,2);
   dhdSig(1,5) = dhc_dSig(0,1);
   dhdSig(2,0) = dhtanphi_dSig(0,0);
   dhdSig(2,1) = dhtanphi_dSig(1,1);
   dhdSig(2,2) = dhtanphi_dSig(2,2);
   dhdSig(2,3) = dhtanphi_dSig(1,2);
   dhdSig(2,4) = dhtanphi_dSig(0,2);
   dhdSig(2,5) = dhtanphi_dSig(0,1);
   dhdSig(3,0) = dhtanpsi_dSig(0,0);
   dhdSig(3,1) = dhtanpsi_dSig(1,1);
   dhdSig(3,2) = dhtanpsi_dSig(2,2);
   dhdSig(3,3) = dhtanpsi_dSig(1,2);
   dhdSig(3,4) = dhtanpsi_dSig(0,2);
   dhdSig(3,5) = dhtanpsi_dSig(0,1);
   return dhdSig;
}

/* calculation of dgdSig_f */
dMatrixT& GRAD_MRSSNLHardT::dgdSig_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& ls, dMatrixT& dgdSig)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, d2QdP2, dB1dP, SN;
   dMatrixT dgchi_dSig(3,3), dgc_dSig(3,3), dgtanphi_dSig(3,3), dgtanpsi_dSig(3,3);
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3), dB2dS_dQdS(3,3), dB3dS_dQdS(3,3);
   dMatrixT tempmat(3,3);
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p + fabs(Sig_p))/2./fGf_I;
   B1 *= ls[2]*ls[2];
   B2 = Sig_Dev;
   B2 /= fGf_I;
   B2 *= ls[3]*ls[3];
   dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   dQdS = Sig_Dev;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   B3 *= ls[3]*ls[3];
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   
   d2QdP2      =  -2.*qn[3]*qn[3];
   dB2dS_dQdS  = Sig_Dev;
   dB2dS_dQdS /= fGf_I;
   dB2dS_dQdS *= ls[2]*ls[2]; 
   dB3dS_dQdS  = Sig_Dev;
   dB3dS_dQdS /= fGf_II;
   dB3dS_dQdS *= ls[3]*ls[3];
   SN = signof(Sig_p);
   dB1dP = (SN +fabs(SN))/2./fGf_I;
   dB1dP *= ls[2]*ls[2];
   
   dgchi_dSig  = IdentityMatrix3;
   dgchi_dSig *= (A1*B1*d2QdP2+A1*dQdP*dB1dP)/3.;
   tempmat =  dB2dS_dQdS; 
   tempmat += B2;  
   tempmat *= A1;
   dgchi_dSig += tempmat;
   dgc_dSig   = B3;
   dgc_dSig += dB3dS_dQdS;
   dgc_dSig  *= A2;
   dgtanphi_dSig  = B3;
   dgtanphi_dSig += dB3dS_dQdS;
   dgtanpsi_dSig *= A3;
   dgtanpsi_dSig  = B3;
   dgtanpsi_dSig += dB3dS_dQdS;
   dgtanpsi_dSig *= A4;
   
   dgdSig(0,0) = dgchi_dSig(0,0);
   dgdSig(0,1) = dgchi_dSig(1,1);
   dgdSig(0,2) = dgchi_dSig(2,2);
   dgdSig(0,3) = dgchi_dSig(1,2);
   dgdSig(0,4) = dgchi_dSig(0,2);
   dgdSig(0,5) = dgchi_dSig(0,1);
   dgdSig(1,0) = dgc_dSig(0,0);
   dgdSig(1,1) = dgc_dSig(1,1);
   dgdSig(1,2) = dgc_dSig(2,2);
   dgdSig(1,3) = dgc_dSig(1,2);
   dgdSig(1,4) = dgc_dSig(0,2);
   dgdSig(1,5) = dgc_dSig(0,1);
   dgdSig(2,0) = dgtanphi_dSig(0,0);
   dgdSig(2,1) = dgtanphi_dSig(1,1);
   dgdSig(2,2) = dgtanphi_dSig(2,2);
   dgdSig(2,3) = dgtanphi_dSig(1,2);
   dgdSig(2,4) = dgtanphi_dSig(0,2);
   dgdSig(2,5) = dgtanphi_dSig(0,1);
   dgdSig(3,0) = dgtanpsi_dSig(0,0);
   dgdSig(3,1) = dgtanpsi_dSig(1,1);
   dgdSig(3,2) = dgtanpsi_dSig(2,2);
   dgdSig(3,3) = dgtanpsi_dSig(1,2);
   dgdSig(3,4) = dgtanpsi_dSig(0,2);
   dgdSig(3,5) = dgtanpsi_dSig(0,1);
   return dgdSig;
}
  
/* calculation of dhdq_f or dqbardq_f*/
dMatrixT& GRAD_MRSSNLHardT::dhdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdq)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, B2dQdS, B3dQdS;
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3);
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
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
   
   dhdq(0,0) = -falpha_chi*(B1*dQdP + B2dQdS);
   dhdq(0,1) =  A1*B1*(2.*qn[3]);
   dhdq(0,2) = 0.;
   dhdq(0,3) =  A1*B1*(2.*qn[1]-4.*Sig_p*qn[3]);   
   dhdq(1,0) = 0.;
   dhdq(1,1) = -falpha_c*B3dQdS;
   dhdq(1,2) = 0.;
   dhdq(1,3) = 0.;
   dhdq(2,0) = 0.;
   dhdq(2,1) = 0.;
   dhdq(2,2) = -falpha_phi*B3dQdS;
   dhdq(2,3) = 0.;
   dhdq(3,0) = 0.;
   dhdq(3,1) = 0.;
   dhdq(3,2) = 0.;
   dhdq(3,3) = -falpha_psi*B3dQdS;
   return dhdq;
}

/* calculation of dgdq_f */
dMatrixT& GRAD_MRSSNLHardT::dgdq_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& ls, dMatrixT& dgdq)
{
   double Sig_p, A1, B1, A2, A3, A4, dQdP, B2dQdS, B3dQdS;
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3);
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   B2 = Sig_Dev;
   B2 /= fGf_I;
   B1 *= ls[2]*ls[2];
   B2 *= ls[3]*ls[3];
   dQdP = 2.*qn[3]*(qn[1] - Sig_p*qn[3]);
   dQdS = Sig_Dev;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   B3 *= ls[3]*ls[3];
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   B2dQdS = dMatrixT::Dot(B2,dQdS);
   B3dQdS = dMatrixT::Dot(B3,dQdS);
   
   dgdq(0,0) = -falpha_chi*(B1*dQdP + B2dQdS);
   dgdq(0,1) =  A1*B1*(2.*qn[3]);
   dgdq(0,2) = 0.;
   dgdq(0,3) =  A1*B1*(2.*qn[1]-4.*Sig_p*qn[3]);   
   dgdq(1,0) = 0.;
   dgdq(1,1) = -falpha_c*B3dQdS;
   dgdq(1,2) = 0.;
   dgdq(1,3) = 0.;
   dgdq(2,0) = 0.;
   dgdq(2,1) = 0.;
   dgdq(2,2) = -falpha_phi*B3dQdS;
   dgdq(2,3) = 0.;
   dgdq(3,0) = 0.;
   dgdq(3,1) = 0.;
   dgdq(3,2) = 0.;
   dgdq(3,3) = -falpha_psi*B3dQdS;
   return dgdq;
}

/* calculation of dhdm_f */
dMatrixT& GRAD_MRSSNLHardT::dhdm_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dhdm)
{
   double Sig_p, A1, B1, A2, A3, A4;
   dMatrixT Sig_Dev(3,3), B2(3,3), B3(3,3), dQdS(3,3);
   dMatrixT dhchi_dm(3,3), dhc_dm(3,3); 
   dMatrixT dhtanphi_dm(3,3), dhtanpsi_dm(3,3);
   dMatrixT tempMat(3,3);
   
   Sig_p = (Sig[0]+Sig[1]+Sig[2])/3.0;
   devstress_f(Sig, Sig_Dev);
   
   A1 = -falpha_chi*(qn[0] - fchi_r);
   B1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   B2 = Sig_Dev;
   B2 /= fGf_I;
   A2 = -falpha_c*(qn[1] - fc_r);
   B3 = Sig_Dev;
   B3 /= fGf_II;
   A3 = -falpha_phi*(qn[2] - tan(fphi_r));
   A4 = -falpha_psi*qn[3];
   dhchi_dm  = IdentityMatrix3;
   dhchi_dm *= (A1*B1)/3.;
   tempMat = B2;
   tempMat *= A1;
   dhchi_dm += tempMat;
   dhc_dm = B3;
   dhc_dm *= A2;
   dhtanphi_dm = B3;
   dhtanphi_dm *= A3;
   dhtanpsi_dm = B3;
   dhtanpsi_dm *= A4;
   
   dhdm(0,0) = dhchi_dm(0,0);
   dhdm(0,1) = dhchi_dm(1,1);
   dhdm(0,2) = dhchi_dm(2,2);
   dhdm(0,3) = dhchi_dm(1,2);
   dhdm(0,4) = dhchi_dm(0,2);
   dhdm(0,5) = dhchi_dm(0,1);
   dhdm(1,0) = dhc_dm(0,0);
   dhdm(1,1) = dhc_dm(1,1);
   dhdm(1,2) = dhc_dm(2,2);
   dhdm(1,3) = dhc_dm(1,2);
   dhdm(1,4) = dhc_dm(0,2);
   dhdm(1,5) = dhc_dm(0,1);
   dhdm(2,0) = dhtanphi_dm(0,0);
   dhdm(2,1) = dhtanphi_dm(1,1);
   dhdm(2,2) = dhtanphi_dm(2,2);
   dhdm(2,3) = dhtanphi_dm(1,2);
   dhdm(2,4) = dhtanphi_dm(0,2);
   dhdm(2,5) = dhtanphi_dm(0,1);
   dhdm(3,0) = dhtanpsi_dm(0,0);
   dhdm(3,1) = dhtanpsi_dm(1,1);
   dhdm(3,2) = dhtanpsi_dm(2,2);
   dhdm(3,3) = dhtanpsi_dm(1,2);
   dhdm(3,4) = dhtanpsi_dm(0,2);
   dhdm(3,5) = dhtanpsi_dm(0,1);
   return dhdm;
}

/* collect the C matrices (C_UU1, C_UU2, C_ULambda1, C_ULambda2, C_LambdaU1, 
 * C_LamdaU2; C_LambdaLambda1,and C_LambdaLambda2) 
 *
 * Note: return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& GRAD_MRSSNLHardT::Moduli(const ElementCardT& element, 
	int ip)
{

    /* define and allocate matrices */
    dMatrixT KE(6,6); dMatrixT KE_AST(6,6);
    dMatrixT KE_UU1(6,6); dMatrixT KE_UU2(6,6);
    dMatrixT KE_ULambda1(6,1); dMatrixT KE_ULambda2(6,1); //6x1
    dMatrixT KE_LambdaU1(1,6); dMatrixT KE_LambdaU2(1,6); //1x6
    dMatrixT KE_LambdaLambda1(1,1); dMatrixT KE_LambdaLambda2(1,1);  
    dMatrixT dhdSig(4,6); dMatrixT dhdq(4,4); dMatrixT dhdm(4,6);
    dMatrixT dgdSig(4,6); dMatrixT dgdq(4,4);
    dMatrixT dmdSig(6,6); dMatrixT dmdq(6,4);
    dMatrixT dRSig_dSig(6,6); dMatrixT dRSig_dq(6,4); 
    dMatrixT dRq_dSig(4,6); dMatrixT dRq_dq(4,4); dMatrixT dRq_dq_Inv(4,4); 
    dMatrixT dRR(10,10); dMatrixT Y(6,6); dMatrixT Y_Inv(6,6); 
     
    /* define and allocate vectors */
    dArrayT u(6); dArrayT up(6); dArrayT du(6); dArrayT dup(6); dArrayT upo(6); 
    dArrayT lap_u(6); dArrayT lap_up(6); dArrayT lap_du(6); dArrayT lap_dup(6); 
    dArrayT lap_upo(6);  
    dArrayT qn(4); dArrayT qo(4);
    dArrayT Sig(6); dArrayT Sig_I(6);
    dArrayT RSig(6), Rq(4);
    dArrayT mm(6); dArrayT rr(4); dArrayT nn(6); 
    dArrayT dq(4); dArrayT hh(4); dArrayT gg(4); 
    dArrayT state(40); dArrayT Sig_trial(6);
    dArrayT R(10); dArrayT ls(4);
    
    double dlam, lap_dlam;
    
	/* initialize */
	KE = 0.; KE_AST=0.;
	KE_UU1 =0.; KE_UU2 = 0.;
	KE_ULambda1 = 0.; KE_ULambda2 = 0.;
	KE_LambdaU1 = 0.; KE_LambdaU2 = 0.;
	KE_LambdaLambda1 = 0.; KE_LambdaLambda2 = 0.; 
	
	ls[0] = flse_v; // lse_v: pore space length scale (elastic)  
	ls[1] = flse_s; // lse_s: grain size length scale (elastic)
	ls[2] = flsp_v; // lsp_v: pore space length scale (plastic)   
	ls[3] = flsp_s; // lsp_s: grain size length scale (plastic)
	
	/* C and C_AST matrices */
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
	
	KE_AST(2,2) = KE_AST(1,1) = KE_AST(0,0) = flambda_ast + 2.0*fmu_ast;
	KE_AST(1,2) = KE_AST(0,1) = KE_AST(0,2) = flambda_ast;
	KE_AST(2,1) = KE_AST(1,0) = KE_AST(2,0) = flambda_ast;
	KE_AST(5,5) = KE_AST(4,4) = KE_AST(3,3) = fmu_ast; 
		
	if(!element.IsAllocated()) 
	{
	  	KE_UU1 = KE;
	  	KE_UU2.SetToScaled(-1., KE_AST);
	  	
	  	fModuli_UU1 = KE_UU1;
		fModuli_UU2 = KE_UU2;
		fModuli_ULam1 = KE_ULambda1;
		fModuli_ULam2 = KE_ULambda2;
		fModuli_LamU1 = KE_LambdaU1;
		fModuli_LamU2 = KE_LambdaU2;
		fModuli_LamLam1 = KE_LambdaLambda1;
		fModuli_LamLam2 = KE_LambdaLambda2;
	  	fModuli = 0.0;     
	  	return fModuli;
	}
    
    /* load internal state variables */
    if(!element.IsAllocated()) 
    {
	  	LoadData(element,ip);
	  	for (int i = 0; i < 39; i++) 
			state[i] = fInternal[i];
	}
	  	
    for (int i = 0; i < 6; i++) 
    {
       Sig[i] = state[i];
       if (i < 4)
       {
       		qn[i] = state[i+30];
       }
    }
	
	if (state[37] == 0) 
	{
	    KE_UU1 = KE;
	  	KE_UU2.SetToScaled(-1., KE_AST);
	  	
	  	fModuli_UU1 = KE_UU1;
		fModuli_UU2 = KE_UU2;
		fModuli_ULam1 = KE_ULambda1;
		fModuli_ULam2 = KE_ULambda2;
		fModuli_LamU1 = KE_LambdaU1;
		fModuli_LamU2 = KE_LambdaU2;
		fModuli_LamLam1 = KE_LambdaLambda1;
		fModuli_LamLam2 = KE_LambdaLambda2; 
	  	fModuli = 0.0;     
	}
	else  
	{
		dlam = state[35];
	  	lap_dlam = state[36];
	 
	  	/* compute the derivatives */
	  	m_f(Sig, qn, mm); dmdSig_f(qn, dmdSig); dmdq_f(Sig, qn, dmdq); 
       	h_f(Sig, qn, hh); dhdSig_f(Sig, qn, dhdSig); dhdq_f(Sig, qn, dhdq); dhdm_f(Sig, qn, dhdm); 
        g_f(Sig, qn, ls, gg); dgdSig_f(Sig, qn, ls, dgdSig); dgdq_f(Sig, qn, ls, dgdq);   
        	
        /* work space */
        dMatrixT tempMat1(6,6), tempMat2(6,6);  
        dMatrixT tempMat3(4,6), tempMat4(4,4);
        tempMat1.SetToScaled(dlam, KE); 
        tempMat2.SetToScaled(lap_dlam, KE_AST);
        tempMat1 -= tempMat2;
        dRSig_dSig.MultAB(tempMat1, dmdSig);
        dRSig_dSig += IdentityMatrix6; 
        dRSig_dq.MultAB(tempMat1, dmdq);
        
        dRq_dSig.MultAB(dhdm, dmdSig);
        dRq_dSig += dhdSig;
        dRq_dSig *= -dlam;
        tempMat3.SetToScaled(lap_dlam, dgdSig);
        dRq_dSig += tempMat3; 
        
        dRq_dq.MultAB(dhdm, dmdq);
        dRq_dq += dhdq;
        dRq_dq *= -dlam;
        tempMat4.SetToScaled(lap_dlam, dgdq);
        dRq_dq += tempMat4;
        dRq_dq += IdentityMatrix4;
         
        /* work space */
        dArrayT RSig1(6), RSig2(6), RSig3(6), RSig4(6);
        dArrayT temp1(6), temp2(6), temp3(4);
        	
        KE.Multx(du, RSig1);
       	KE_AST.Multx(lap_du, RSig2);
       	KE.Multx(mm, RSig3);  
        KE_AST.Multx(mm, RSig4);
        RSig.DiffOf(RSig2, RSig1); //RSig2-RSig1
        temp1.SetToScaled(dlam, RSig3);  //dlam2 == dlam??
        RSig += temp1;
        temp2.SetToScaled(lap_dlam, RSig4); // lap_dlam2 == dlam??	
        RSig -= temp2;
        
        Rq.SetToScaled(lap_dlam, gg);   //lap_dlam2/lap_dlam??
        temp3.SetToScaled(dlam, hh);     //dlam2/dlam??
        Rq -= temp3; 
        
        /* calculate KE_UU1 and KE_UU2 matrices */
        dMatrixT RRq_dqdSig(4,6);
        dRq_dq_Inv.Inverse(dRq_dq);
        RRq_dqdSig.MultAB(dRq_dq_Inv, dRq_dSig);
        Y.MultAB(dRSig_dq, RRq_dqdSig);
        Y -= dRSig_dSig;
        Y_Inv.Inverse(Y);
        KE_UU1.MultAB(Y_Inv, KE);
	  	KE_UU2.MultAB(Y_Inv, KE_AST);
        	
        /* calculate KE_ULambda1 and KE_ULambda2 matrices */
        /* work space */
        dArrayT termA1(6), termA2i(4), termA2(6);
        dArrayT termB1(6), termB2i(4), termB2(6);
        	
        KE.Multx(mm, termA1);
        dRq_dq_Inv.Multx(hh, termA2i);
        dRSig_dq.Multx(termA2i, termA2);
        termA2 += termA1;
        KE_AST.Multx(mm, termB1);
        dRq_dq_Inv.Multx(gg, termB2i);
        dRSig_dq.Multx(termB2i, termB2);
        termB2 += termB1;
	  	Y_Inv.Multx(termA2, KE_ULambda1);
	  	KE_ULambda1 *= -1.;
	  	Y_Inv.Multx(termB2, KE_ULambda2);
	  	    
	  	//cout << "Y inverse =" << endl << Y_Inv << endl;
	  	    
	  	/* calculate KE_LambdaU1 and KE_LambdaU2 matrices */
	  	/* work space */
	  	dArrayT mul1(4),term2C(6),TT(6);
	  	dMatrixT dRq_dSigT(6,4);
	  	    
	  	dRq_dq_Inv.Multx(rr, mul1);
	  	dRq_dSigT.Transpose(dRq_dSig);
	  	dRq_dSigT.Multx(mul1, term2C);
	  	term2C *= -1.;
	  	term2C += nn;
	  	Y_Inv.Multx(term2C, TT);
	  	KE.Multx(TT, KE_LambdaU1);
	  	KE_AST.Multx(TT, KE_LambdaU2);
	  	KE_LambdaU2 *= -1.;
	  	    
	  	/* calculate KE_LambdaLambda1 and KE_LambdaLambda2 matrices */
	  	dArrayT termD2i(4), termE2i(4);
	  	double termD1, termD2, termE1, termE2;
	  	termD1 = termD2 = termE1 = termE2 = 0.;
	  	dRq_dq_Inv.Multx(hh, termD2i);
	  	dRq_dq_Inv.Multx(gg, termE2i);
	  	termD1 = dMatrixT::Dot(TT, termA2); //vector dot product 
	  	termE1 = dMatrixT::Dot(TT, termB2);
	  	termD2 = dMatrixT::Dot(rr, termD2i); //vector dot product 
	  	termE2 = dMatrixT::Dot(rr, termE2i);
	  	termD1 += termD2;
	  	termE1 += termE2;
	  	KE_LambdaLambda1 = -1.* termD1;
	  	KE_LambdaLambda2 = termE1; 
	
        fModuli_UU1 = KE_UU1;
		fModuli_UU2 = KE_UU2;
		fModuli_ULam1 = KE_ULambda1;
		fModuli_ULam2 = KE_ULambda2;
		fModuli_LamU1 = KE_LambdaU1;
		fModuli_LamU2 = KE_LambdaU2;
		fModuli_LamLam1 = KE_LambdaLambda1;
		fModuli_LamLam2 = KE_LambdaLambda2; 
	  	fModuli = 0.0;
	} // if (state[37] == 0.)
	return fModuli; 
}


/* return the correction to modulus Cep~, checking for discontinuous
 *   bifurcation */
const dMatrixT& GRAD_MRSSNLHardT::ModuliPerfPlas(const ElementCardT& element, 
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
void GRAD_MRSSNLHardT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fLapPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += kNumInternal*fNumIP;        //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
	
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0
}

/* accept parameter list */
void GRAD_MRSSNLHardT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	GRAD_MRPrimitiveT::TakeParameterList(list);

	/* dimension work space */
	fElasticStrain.Dimension(kNSD);
	fLapElasticStrain.Dimension(kNSD);
	fStressCorr.Dimension(kNSD);
	fModuli.Dimension(dSymMatrixT::NumValues(kNSD));
	fModuliPerfPlas.Dimension(dSymMatrixT::NumValues(kNSD));
	fModuli_UU1.Dimension(dSymMatrixT::NumValues(kNSD));
	fModuli_UU2.Dimension(dSymMatrixT::NumValues(kNSD));
	fModuli_ULam1.Dimension(dSymMatrixT::NumValues(kNSD),1);
	fModuli_ULam2.Dimension(dSymMatrixT::NumValues(kNSD),1);
	fModuli_LamU1.Dimension(1,dSymMatrixT::NumValues(kNSD));
	fModuli_LamU2.Dimension(1,dSymMatrixT::NumValues(kNSD));
	fModuli_LamLam1.Dimension(1,1);
	fModuli_LamLam2.Dimension(1,1);
	fDevStress.Dimension(kNSD);
	fLapDevStress.Dimension(kNSD);
	fDevStrain.Dimension(kNSD);
	fLapDevStrain.Dimension(kNSD);
	IdentityMatrix3.Dimension(3); 
	IdentityMatrix4.Dimension(4);
	IdentityMatrix6.Dimension(6);
	fIniInternal.Dimension(4);
    
	/* initialize constant matrices */
	IdentityMatrix3.Identity();
	IdentityMatrix4.Identity();
	IdentityMatrix6.Identity();
}
/***********************************************************************
 * Protected
 ***********************************************************************/

/* element level data */
void GRAD_MRSSNLHardT::Update(ElementCardT& element)
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
void GRAD_MRSSNLHardT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* load element data for the specified integration point */
void GRAD_MRSSNLHardT::LoadData(const ElementCardT& element, int ip)
{
	/* check */
	if (!element.IsAllocated()) 
	    ExceptionT::GeneralFail("GRAD_MRSSNLHardT::LoadData","The element should have been allocated");
	/* fetch arrays */
	 const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
	
	fPlasticStrain.Alias(dim, &d_array[dex]);
	fLapPlasticStrain.Alias(dim, &d_array[dex]);
	/*fUnitNorm.Set(kNSD, &d_array[  offset + dex]); */    
	fInternal.Alias(kNumInternal, &d_array[2*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int GRAD_MRSSNLHardT::PlasticLoading(const dSymMatrixT& trialstrain, 
	  const dSymMatrixT& lap_trialstrain, ElementCardT& element, int ip) 
{
	/* not yet plastic */
	if (!element.IsAllocated()) { 
		//cout << "!element.IsAllocated()" << endl;
		//double yield = YieldCondition(DeviatoricStress(trialstrain,lap_trialstrain,element),
			   //MeanStress(trialstrain,lap_trialstrain,element));
		//cout << "failure fun = " << yield;
		
		return( YieldCondition(DeviatoricStress(trialstrain,lap_trialstrain,element),
			   MeanStress(trialstrain,lap_trialstrain,element)) > kYieldTol );
        /* already plastic */
    }
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
dSymMatrixT& GRAD_MRSSNLHardT::DeviatoricStress(const dSymMatrixT& trialstrain,
    const dSymMatrixT& lap_trialstrain, const ElementCardT& element)
{
#pragma unused(element)

	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain); 
	fLapDevStrain.Deviatoric(lap_trialstrain);

	/* compute deviatoric elastic stress */
	fDevStress.SetToScaled(2.0*fmu,fDevStrain);
	fLapDevStress.SetToScaled(2.0*fmu_ast,fLapDevStrain); 
	fDevStress -= fLapDevStress;
	return fDevStress;
}

/* computes the hydrostatic (mean) stress */
double GRAD_MRSSNLHardT::MeanStress(const dSymMatrixT& trialstrain,
    const dSymMatrixT& lap_trialstrain, const ElementCardT& element)
{
#pragma unused(element)

  fMeanStress = fkappa*trialstrain.Trace(); 
  fLapMeanStress = fkappa_ast*lap_trialstrain.Trace(); 
  fMeanStress -= fLapMeanStress;  
  return fMeanStress;
}

double GRAD_MRSSNLHardT::signof(double& r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}