/* $Id: GRAD_MRSSNLHardT.cpp,v 1.32 2006-08-29 21:18:41 kyonten Exp $ */
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

const int    kNumInternal = 8; // number of internal state variables
const int    kNSD         = 3;
const int    kNSTR        = dSymMatrixT::NumValues(kNSD);
const double ratio23      = 2.0/3.0;
const double ratio32      = 3.0/2.0;

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
const dSymMatrixT& GRAD_MRSSNLHardT::StressCorrection(const dSymMatrixT& totalstrain, 
                  const dSymMatrixT& lap_totalstrain, const dArrayT& lambda, 
                  const dArrayT& lap_lambda, ElementCardT& element, int ip)
{	
  	/* elastic step */
  	fStressCorr = DeviatoricStress(totalstrain, lap_totalstrain, element);
    fStressCorr.PlusIdentity(MeanStress(totalstrain, lap_totalstrain, element));
    
    /* plastic step */
    if (PlasticLoading(totalstrain, lap_totalstrain, lambda, element, ip) && 
	    !element.IsAllocated())
	{
		/* new plastic element */
		AllocateElement(element);
		
		/* initialize element data */ 
		PlasticLoading(totalstrain, lap_totalstrain, lambda, element, ip); 
	}
	
	if (element.IsAllocated()) 
	{		
    	/* define and allocate matrices */ 
    	dMatrixT KE(6),KE_AST(6),dhdSig(4,6),dhdq(4),dhdm(4,6),
             	 dgdSig(4,6),dgdq(4),dmdSig(6),dmdq(6,4),
             	 dRSig_dSig(6),dRSig_dq(6,4),RSigq_qq(6,4), 
             	 dRq_dSig(4,6),dRq_dq(4),dRq_dq_Inv(4), 
             	 RRq_dqdSig(4,6),Y(6),Y_Inv(6);
    
    	/* reduced index vectors of symmetric matrices */
    	dSymMatrixT up(3),du(3),dup(3),upo(3),lap_up(3),lap_du(3), 
                	lap_dup(3),lap_upo(3),Sig(3),dSig(3),ue(3), 
                	lap_ue(3),Sig_e(3),lap_Sig_e(3); 
     
    	/* define and allocate vectors */ 
    	dArrayT qn(4),dq(4),qo(4),mm(6),rr(4),nn(6),hh(4),gg(4), 
            	ls(2),RSig(6),Rq(4);
    
    	/* elastic moduli tensor */
		KE = 0.0;
		KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
		KE(1,2) = KE(0,1) = KE(0,2) = flambda;
		KE(2,1) = KE(1,0) = KE(2,0) = flambda;
		KE(5,5) = KE(4,4) = KE(3,3) = fmu;
    
    	/* elastic moduli tensor with length scale effect */
    	KE_AST = 0.0;
		KE_AST(2,2) = KE_AST(1,1) = KE_AST(0,0) = flambda_ast + 2.0*fmu_ast;
		KE_AST(1,2) = KE_AST(0,1) = KE_AST(0,2) = flambda_ast;
		KE_AST(2,1) = KE_AST(1,0) = KE_AST(2,0) = flambda_ast;
		KE_AST(5,5) = KE_AST(4,4) = KE_AST(3,3) = fmu_ast;
	
		/* get plastic multiplier and its laplacian */ 
		double dlam = lambda[0]; 
    	double lap_dlam = lap_lambda[0];
    	
    	bool print = false;
    	if(print) {
    		cout << "strain" << endl;
    		cout << totalstrain << endl << endl;
    		cout << "lap strain" << endl;
    		cout << lap_totalstrain << endl << endl;
    		cout << "lambda = " << dlam << endl << endl;
    		cout << "lap lambda = " << lap_dlam << endl << endl;
    	}
	
		/* initialize the necessary vectors */
    	up = fPlasticStrain;	
    	lap_up = fLapPlasticStrain;	//laplacian of up
    	qn.CopyPart(0, fInternal, 0, qn.Length());
    	upo = up; lap_upo = lap_up;
    	qo = qn;
    
    	/* check the yield function */
    	ue.DiffOf(totalstrain, up);
    	lap_ue.DiffOf(lap_totalstrain, lap_up);
    	Sig = DeviatoricStress(ue, lap_ue, element);
    	Sig.PlusIdentity(MeanStress(ue, lap_ue, element));
    	double ff = fInternal[kftrial];  
    	
    	if (ff > fTol_1 && dlam > 0.)
    	{ /* local Newton iteration */
      		int kk = 0;
      		int max_iteration = 15; 
      		while (fabs(ff) > fTol_1) 
      		{
        		/* calculate stress */
        		ue.DiffOf(totalstrain, up);
    			lap_ue.DiffOf(lap_totalstrain, lap_up);
    			Sig = DeviatoricStress(ue, lap_ue, element);
    			Sig.PlusIdentity(MeanStress(ue, lap_ue, element));
        		
        		/* calculate yield condition */
        		ff = YieldCondition(DeviatoricStress(ue, lap_ue, element),
			         MeanStress(ue, lap_ue, element),qn[0],qn[1],qn[2]);
 
    			/* calculate all the necessary derivatives */
   				m_f(Sig, qn, mm); 
   				dmdSig_f(Sig, qn, dmdSig); 
   				dmdq_f(Sig, qn, dmdq);
    			h_f(Sig, qn, hh); 
    			dhdSig_f(Sig, qn, dhdSig); 
    			dhdq_f(Sig, qn, dhdq); 
    			dhdm_f(Sig, qn, dhdm);  
    			g_f(Sig, qn, gg); 
    			dgdSig_f(Sig, qn, dgdSig); 
    			dgdq_f(Sig, qn, dgdq);   
        
    			/* calculate R_Sig_Sig and R_Sig_q matrices */
    			dMatrixT RSigMat1(6), RSigMat2(6); /* work space */ 
    			RSigMat1.SetToScaled(dlam, KE);       
    			RSigMat2.SetToScaled(lap_dlam, KE_AST);
    			RSigMat1 -= RSigMat2;  
    			dRSig_dSig.MultAB(RSigMat1, dmdSig);
    			dRSig_dSig.PlusIdentity();
         
    			dRSig_dq.MultAB(RSigMat1, dmdq);
        
    			/* calculate R_q_Sig and R_q_q matrices */
    			dMatrixT RqMat1(4,6), RqMat2(4); /* work space */
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
    			dRq_dq.PlusIdentity();
         
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
        
   				/* solving for dSig and dq */
   				/* work space */
   				LAdMatrixT RR(10);
   				dArrayT RR_vec(10);
   				RR = 0.0;
   				RR.AddBlock(0,                 0,                 dRSig_dSig);
   				RR.AddBlock(0,                 dRSig_dSig.Cols(), dRSig_dq);
   				RR.AddBlock(dRSig_dSig.Rows(), 0,                 dRq_dSig);
   				RR.AddBlock(dRSig_dSig.Rows(), dRSig_dSig.Cols(), dRq_dq);
   				RR_vec.CopyIn(0, RSig);
   				RR_vec.CopyIn(RSig.Length(), Rq);
   				RR.LinearSolve(RR_vec);
   				dSig.CopyPart(0, RR_vec, 0, dSig.Length());
   				dq.CopyPart(0, RR_vec, dSig.Length(), dq.Length()); 
   		
    			/* update plastic strain, gradient plastic strain 
    			and internal variables 						*/
   				m_f(dSig, qn, mm);
    			dup.SetToScaled(dlam, mm); //KE_Inv.Multx(dSig,dup);
    			up += dup;
    			lap_up.SetToScaled(lap_dlam, mm);  //lap_dlam2 == lap_dlam??
    			lap_up += lap_dup;
    			Sig += dSig;   // stress automatically updated when up & lap_up are updated??
    			qn += dq;
    			
        		/* terminate iteration */
        		if (++kk == max_iteration)
        			ExceptionT::GeneralFail("MRSSNLHardT::StressCorrection","local iteration failed after %d iterations", max_iteration);
      		} // while (fabs(ff) > fTol_1)
    	} // if (ff > fTol_1)
    
    	/* update state variables */
		fInternal.CopyIn(0, qn);
		fInternal[4] = ff; 
		fInternal[5] = dlam;
		fInternal[6] = lap_dlam;
		fInternal[7] = 1.;
 	
		/* update plastic strain */
   		fPlasticStrain = up;
   		fLapPlasticStrain = lap_up;
   		
   		/* update stress */	
   		fStress = Sig; 
		fStressCorr = Sig;
		
	} // if (element.IsAllocated())		
 return fStressCorr;
}	//end StressCorrection

/* calculation of dfdSig_f or n_f*/
void GRAD_MRSSNLHardT::n_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdSig)
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

/* calculation of dfdq_f or r_f*/
void GRAD_MRSSNLHardT::r_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dfdq)
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

/* calculation of dQdSig_f or m_f*/
void GRAD_MRSSNLHardT::m_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& dQdSig)
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

/* calculation of dQdSig2_f or dmdSig_f */
void GRAD_MRSSNLHardT::dmdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSig2)
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

/* calculation of dQdSigdq_f or dmdq_f */
void GRAD_MRSSNLHardT::dmdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
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

/* calculation of qbar_f or h_f vector */
void GRAD_MRSSNLHardT::h_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& qbar)
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
 
/* calculation of dqbardSig_f or dhdSig_f */
void GRAD_MRSSNLHardT::dhdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
   dSymMatrixT dhchi_dSig(3),dhc_dSig(3),dhtanphi_dSig(3),dhtanpsi_dSig(3),
               Sig_Dev(3),dQdS(3),II_S(3),SxS_S(3),tmp1(3),tmp2(3);
   dMatrixT SxS(6,6),II(6,6);
   
   II.Identity();
               
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

/* calculation of dqbardq_f or dhdq_f*/
void GRAD_MRSSNLHardT::dhdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
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

/* calculation of dhdm_f */
void GRAD_MRSSNLHardT::dhdm_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dhdm)
{
   dSymMatrixT dhchi_dm(3),dhc_dm(3),dhtanphi_dm(3),dhtanpsi_dm(3),
               Sig_Dev(3),B2(3),B3(3);
               
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
   
   Sig_Dev.Deviatoric(Sig);
   B2.SetToScaled(1.0/fGf_I, Sig_Dev);
   B3.SetToScaled(1.0/fGf_II, Sig_Dev);
   dhchi_dm.SetToScaled(A1, B2);
   for (int i=0; i<3; i++)
       dhchi_dm[i] += A1*B1/3.; 
   dhc_dm.SetToScaled(A2, B3);
   dhtanphi_dm.SetToScaled(A3, B3);
   dhtanpsi_dm.SetToScaled(A4, B3);
  
   dhdm = 0.0;
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
}

/* calculation of g_f vector */
void GRAD_MRSSNLHardT::g_f(const dSymMatrixT& Sig, const dArrayT& qn, dArrayT& gg)
{
   dSymMatrixT Sig_Dev(3),BB2(3),BB3(3),dQdS(3); 

   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A1 = -falpha_chi*(fchi - fchi_r);
   double A2 = -falpha_c*(fc - fc_r);
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   double lspv2 = flsp_v*flsp_v;
   double lsps2 = flsp_s*flsp_s;
   double Sig_p = Sig.Trace()/3.0;
   double BB1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   BB1 *= lspv2; 
   double dQdP = ftan_psi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   double temp2 = fc - ftan_psi*fchi;
   double temp3 = temp2 * temp2;
   temp += temp3;
   double deno = 2.0*sqrt(temp);
   BB2.SetToScaled(lsps2/fGf_I, Sig_Dev);
   BB3.SetToScaled(lsps2/fGf_II, Sig_Dev);
   dQdS.SetToScaled(1./deno, Sig_Dev);
 
   gg=0.0;   
   gg[0]  = A1*BB1*dQdP; 
   gg[0] += A1*BB2.ScalarProduct(dQdS); 
   gg[1]  = A2*BB3.ScalarProduct(dQdS); 
   gg[2]  = A3*BB3.ScalarProduct(dQdS); 
   gg[3]  = A4*BB3.ScalarProduct(dQdS); 
 }

/* calculation of dgdSig_f */
void GRAD_MRSSNLHardT::dgdSig_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dgdSig)
{
   dSymMatrixT dgchi_dSig(3),dgc_dSig(3),dgtanphi_dSig(3),dgtanpsi_dSig(3),
               Sig_Dev(3),dQdS(3),II_S(3),SxS_S(3),tmp1(3),tmp2(3);
   dMatrixT SxS(6,6),II(6,6);
   
   II.Identity();
               
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A1 = -falpha_chi*(fchi - fchi_r);
   double A2 = -falpha_c*(fc - fc_r);
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   double lspv2 = flsp_v*flsp_v;
   double lsps2 = flsp_s*flsp_s;
   double Sig_p = Sig.Trace()/3.0;
   double SN = signof(Sig_p);
   double BB1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   BB1 *= lspv2;
   double dBB1dP = (SN +fabs(SN))/2./fGf_I;
   dBB1dP *= lspv2; 
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
   tmp1 *= (lsps2/fGf_I);
   tmp2 *= (lsps2/fGf_II);
   
   dgchi_dSig = 0.0;
   dgchi_dSig = tmp1;
   for (int i=0; i<3; i++)
       dgchi_dSig[i] += dQdP*dBB1dP/3.;
   dgchi_dSig *= A1; 
   dgc_dSig = tmp2;
   dgc_dSig *= A2;
   dgtanphi_dSig = tmp2;
   dgtanphi_dSig *= A3;
   dgtanpsi_dSig = tmp2;
   dgtanpsi_dSig *= A4;
  
   dgdSig = 0.;
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
}

/* calculation of dgdq_f */
void GRAD_MRSSNLHardT::dgdq_f(const dSymMatrixT& Sig, const dArrayT& qn, dMatrixT& dgdq)
{
   dSymMatrixT Sig_Dev(3),BB2(3),BB3(3),dQdS(3);
   
   double fchi = qn[0];
   double fc = qn[1];
   double ftan_phi = qn[2];
   double ftan_psi = qn[3]; 
   double A1 = -falpha_chi*(fchi - fchi_r);
   double A2 = -falpha_c*(fc - fc_r);
   double A3 = -falpha_phi*(ftan_phi - tan(fphi_r));
   double A4 = -falpha_psi*ftan_psi;
   double lspv2 = flsp_v*flsp_v;
   double lsps2 = flsp_s*flsp_s;
   double Sig_p = Sig.Trace()/3.0;
   double BB1 = (Sig_p+fabs(Sig_p))/2./fGf_I;
   BB1 *= lspv2;
   double dQdP = ftan_psi;
   double const1 = fc - ftan_psi*fchi;
   
   Sig_Dev.Deviatoric(Sig);
   double temp  = (Sig_Dev.ScalarProduct())/2.0;
   temp += (const1*const1);
   double deno = 2.0*sqrt(temp);
   const1 /= 2.0*pow(temp, ratio32);
   dQdS = Sig_Dev;
   dQdS /= deno;
   BB2.SetToScaled(lsps2/fGf_I, Sig_Dev);
   BB3.SetToScaled(lsps2/fGf_II, Sig_Dev);
   double B2S = BB2.ScalarProduct(Sig_Dev); 
   double B3S = BB3.ScalarProduct(Sig_Dev); 
   double B2dQdS = BB2.ScalarProduct(dQdS); 
   double B3dQdS = BB3.ScalarProduct(dQdS); 
   dgdq = 0.0;
   dgdq(0,0)  = -falpha_chi*(BB1*dQdP + B2dQdS);
   dgdq(0,0) += A1*ftan_psi*const1*B2S;
   dgdq(0,1)  =  -A1*const1*B2S;
   dgdq(0,3)  =  A1*BB1+A1*fchi*const1*B2S;
   dgdq(1,0)  = A3*ftan_psi*const1*B3S;   
   dgdq(1,1)  = -falpha_c*B3dQdS;
   dgdq(1,1) -= A2*const1*B3S;
   dgdq(1,3)  = A2*fchi*const1*B3S;
   dgdq(2,0)  = A3*ftan_psi*const1*B3S;
   dgdq(2,1)  = -A3*const1*B3S;
   dgdq(2,2)  = -falpha_phi*B3dQdS;
   dgdq(2,3)  = A3*fchi*const1*B3S;
   dgdq(3,0)  = A4*ftan_psi*const1*B3S;
   dgdq(3,1)  = -A4*const1*B3S;
   dgdq(3,3)  = -falpha_psi*B3dQdS;
   dgdq(3,3) += A4*fchi*const1*B3S;
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
    dMatrixT KE(6),KE_AST(6),KE_UU1(6),KE_UU2(6),
             KE_ULambda1(6,1),KE_ULambda2(6,1), //6x1
             KE_LambdaU1(1,6),KE_LambdaU2(1,6), //1x6
             KE_LambdaLambda1(1),KE_LambdaLambda2(1),  
             dhdSig(4,6),dhdq(4),dhdm(4,6),
             dgdSig(4,6),dgdq(4),dmdSig(6),dmdq(6,4),
             dRSig_dSig(6),dRSig_dq(6,4), 
             dRq_dSig(4,6),dRq_dq(4),dRq_dq_Inv(4),
             RRq_dqdSig(4,6),RSigq_qq(6,4), 
             dRR(10),Y(6),Y_Inv(6);
    dSymMatrixT Sig(3); 
     
    /* define and allocate vectors */   
    dArrayT qn(4),RSig(6),Rq(4),T(6),mm(6),rr(4),nn(6), 
            hh(4),gg(4),R(10),ls(2);
    
	/* elastic moduli tensor */
	KE = 0.0;
	KE(2,2) = KE(1,1) = KE(0,0) = flambda + 2.0*fmu;
	KE(1,2) = KE(0,1) = KE(0,2) = flambda;
	KE(2,1) = KE(1,0) = KE(2,0) = flambda;
	KE(5,5) = KE(4,4) = KE(3,3) = fmu;
    
    /* elastic moduli tensor with length scale effect */
    KE_AST = 0.0;
	KE_AST(2,2) = KE_AST(1,1) = KE_AST(0,0) = flambda_ast + 2.0*fmu_ast;
	KE_AST(1,2) = KE_AST(0,1) = KE_AST(0,2) = flambda_ast;
	KE_AST(2,1) = KE_AST(1,0) = KE_AST(2,0) = flambda_ast;
	KE_AST(5,5) = KE_AST(4,4) = KE_AST(3,3) = fmu_ast;
	
	/* initialize */
	KE_UU1 = 0.; 
	KE_UU2 = 0.;
	KE_ULambda1 = 0.; 
	KE_ULambda2 = 0.;
	KE_LambdaU1 = 0.; 
	KE_LambdaU2 = 0.;
	KE_LambdaLambda1 = 0.; 
	KE_LambdaLambda2 = 0.; 
		
    if(element.IsAllocated() && (element.IntegerData())[ip] == kIsPlastic) 
    {
	  	/* load internal state variables */
	  	LoadData(element, ip);
	  	Sig = fStress;
	  	qn.CopyPart(0, fInternal, 0, qn.Length());
	 
		double dlam = fInternal[klambda];
		double lap_dlam = fInternal[klaplambda];
	  	
		/* compute all the required derivatives */
		m_f(Sig, qn, mm);
		n_f(Sig, qn, nn); 
   		r_f(Sig, qn, rr);  
		dmdSig_f(Sig, qn, dmdSig); 
		dmdq_f(Sig, qn, dmdq); 
    	h_f(Sig, qn, hh); 
    	dhdSig_f(Sig, qn, dhdSig); 
    	dhdq_f(Sig, qn, dhdq); 
    	dhdm_f(Sig, qn, dhdm); 
    	g_f(Sig, qn, gg); 
    	dgdSig_f(Sig, qn, dgdSig); 
    	dgdq_f(Sig, qn, dgdq);  
        	
    	/* work space */
    	dMatrixT tempMat1(6,6), tempMat2(6,6);  
    	dMatrixT tempMat3(4,6), tempMat4(4,4);
        
    	/* dRSig_dSig matrix */
    	tempMat1.SetToScaled(dlam, KE); 
    	tempMat2.SetToScaled(lap_dlam, KE_AST);
    	tempMat1 -= tempMat2;
    	dRSig_dSig.MultAB(tempMat1, dmdSig);
    	dRSig_dSig.PlusIdentity();
         
    	/* dRSig_dq matrix */
    	dRSig_dq.MultAB(tempMat1, dmdq);
        
    	/* dRq_dSig matrix */
    	dRq_dSig.MultAB(dhdm, dmdSig);
    	dRq_dSig += dhdSig;
    	dRq_dSig *= -dlam;
    	tempMat3.SetToScaled(lap_dlam, dgdSig);
    	dRq_dSig += tempMat3; 
        
    	/* dRq_dq matrix */
    	dRq_dq.MultAB(dhdm, dmdq);
    	dRq_dq += dhdq;
    	dRq_dq *= -dlam;
    	tempMat4.SetToScaled(lap_dlam, dgdq);
    	dRq_dq += tempMat4;
    	dRq_dq.PlusIdentity();
        
    	/* Y and Y_Inv matrix */
    	dRq_dq_Inv.Inverse(dRq_dq);
    	RRq_dqdSig.MultAB(dRq_dq_Inv, dRq_dSig);
    	Y.MultAB(dRSig_dq, RRq_dqdSig);
    	Y *= -1.0;
    	Y += dRSig_dSig;
    	Y_Inv.Inverse(Y);
        
    	/* T vector */
    	dArrayT TT(6); /* work space */
    	RRq_dqdSig.MultTx(rr, TT);
    	TT *= -1.0;
    	TT += nn;
    	Y_Inv.MultTx(TT, T);
        
    	/* RSigq_qq vector */
    	RSigq_qq.MultAB(dRSig_dq, dRq_dq_Inv);
        
    	/* calculate KE_UU1 and KE_UU2 matrices */
    	KE_UU1.MultAB(Y_Inv, KE);
		KE_UU2.MultAB(Y_Inv, KE_AST);
		KE_UU2 *= -1.0;
        	
    	/* calculate KE_ULambda1 and KE_ULambda2 matrices */
    	/* work space */
    	dArrayT termA1(6), termA2(6);
    	dArrayT termB1(6), termB2(6);
        	
    	KE.Multx(mm, termA1);
    	RSigq_qq.Multx(hh, termA2);
    	termA1 += termA2;
    	Y_Inv.Multx(termA1, KE_ULambda1);
    	KE_ULambda1 *= -1.0;
        
    	KE_AST.Multx(mm, termB1);
    	RSigq_qq.Multx(gg, termB2);
    	termB1 += termB2;
		Y_Inv.Multx(termB1, KE_ULambda2);
	  	    
		/* calculate KE_LambdaU1 and KE_LambdaU2 matrices */
		KE_UU1.Multx(T, KE_LambdaU1);
		KE_UU2.Multx(T, KE_LambdaU2);
		KE_LambdaU2 *= -1.0;
	  	    
		/* calculate KE_LambdaLambda1 and KE_LambdaLambda2 matrices */
	 	/* work space */
		dArrayT termD1a(6), termD2a(4);
		dArrayT termE1a(6), termE2a(4); 
		double termD2, termE2;
		Y.Multx(KE_ULambda1, termD1a);
		KE_LambdaLambda1 = dMatrixT::Dot(T, termD1a); //vector dot product
		dRq_dq_Inv.Multx(hh, termD2a); 
		termD2 = dMatrixT::Dot(rr, termD2a);  
		KE_LambdaLambda1 += termD2;
	  	
		Y.Multx(KE_ULambda2, termE1a);
		KE_LambdaLambda2 = dMatrixT::Dot(T, termE1a);
		dRq_dq_Inv.Multx(gg, termE2a);
		termE2 = dMatrixT::Dot(rr, termE2a);
		KE_LambdaLambda2 -= termE2; 
	
    	fModuli_UU1 = KE_UU1;
		fModuli_UU2 = KE_UU2;
		fModuli_ULam1 = KE_ULambda1;
		fModuli_ULam2 = KE_ULambda2;
		fModuli_LamU1 = KE_LambdaU1;
		fModuli_LamU2 = KE_LambdaU2;
		fModuli_LamLam1 = KE_LambdaLambda1;
		fModuli_LamLam2 = KE_LambdaLambda2;
		fModuli = 0.0;
		
		//*************debug*********************//
	StringT file_name;  
	file_name = "C:/Documents and Settings/kyonten/My Documents/tahoe_xml/C_Matrices";
	//file_name = "C:/Documents and Settings/Administrator/My Documents/tahoe/C_Matrices";
	file_name.Append(".txt");
	ofstream output(file_name);
	if (!output) {
		cout << "Error opening output file" << endl;
	}
			
		/* print C matrices */
		for (int i = 0; i < mm.Length(); i++)
				if (mm[i] != 0.0) 
					output << "m ("<<i<<"): "<< mm[i] << endl;
		output << endl;
			
		for (int i = 0; i < nn.Length(); i++)
			if (nn[i] != 0.0)
				output << "n ("<<i<<"): "<< nn[i] << endl;
		output << endl;
		
		for (int i = 0; i < rr.Length(); i++)
			if (rr[i] != 0.0)
				output << "r ("<<i<<"): "<< rr[i] << endl;
		output << endl;
		
		for (int i = 0; i < hh.Length(); i++)
			if (hh[i] != 0.0)
				output << "h ("<<i<<"): "<< hh[i] << endl;
		output << endl;
		
		for (int i = 0; i < gg.Length(); i++)
			if (gg[i] != 0.0)
				output << "g ("<<i<<"): "<< gg[i] << endl;
		output << endl;
		
		for (int i = 0; i < T.Length(); i++)
			if (T[i] != 0.0)
				output << "T ("<<i<<"): "<< T[i] << endl;
		output << endl;
		
		for (int i = 0; i < Y.Rows(); i++)
		{
			for (int j = 0; j < Y.Cols(); j++) 
				if (Y(i,j) != 0.0) 
					output << "Y("<< i << ","<< j <<"): " << Y(i,j) << endl;
		}
		output << endl;	
		
		output << "*******CUU1******* " << endl;
		for (int i = 0; i < fModuli_UU1.Rows(); i++)
		{
			for (int j = 0; j < fModuli_UU1.Cols(); j++) 
				if (fModuli_UU1(i,j) != 0.0) 
					output << "CUU1("<< i << ","<< j <<"): " << fModuli_UU1(i,j) << endl;
		}
		output << endl;	
		
		output << "*******CUU2******* " << endl;
		for (int i = 0; i < fModuli_UU2.Rows(); i++)
		{
			for (int j = 0; j < fModuli_UU2.Cols(); j++) 
				if (fModuli_UU2(i,j) != 0.0) 
					output << "CUU2("<< i << ","<< j <<"): " << fModuli_UU2(i,j) << endl;
		}
		output << endl;	
			
		output << "*******CULambda1******* " << endl;
		for (int i = 0; i < fModuli_ULam1.Rows(); i++)
		{
			for (int j = 0; j < fModuli_ULam1.Cols(); j++) 
				if (fModuli_ULam1(i,j) != 0.0)  
					output << "CULambda1("<< i << ","<< j <<"): " << fModuli_ULam1(i,j) << endl;
		}
		output << endl;	
		
		output << "*******CULambda2******* " << endl;
		for (int i = 0; i < fModuli_ULam2.Rows(); i++)
		{
			for (int j = 0; j < fModuli_ULam2.Cols(); j++) 
				if (fModuli_ULam2(i,j) != 0.0)    
					output << "CULambda2("<< i << ","<< j <<"): " << fModuli_ULam2(i,j) << endl;	
		}
		output << endl;	
				
		output << "*******CLambdaLambda1******* " << endl;
		for (int i = 0; i < fModuli_LamLam1.Rows(); i++)
		{
			for (int j = 0; j < fModuli_LamLam1.Cols(); j++) 
				if (fModuli_LamLam1(i,j) != 0.0)  
					output << "CLambdaLambda1("<< i << ","<< j <<"): " << fModuli_LamLam1(i,j) << endl;
		}
		output << endl;	
		
		output << "*******CLambdaLambda2******* " << endl;
		for (int i = 0; i < fModuli_LamLam2.Rows(); i++)
		{
			for (int j = 0; j < fModuli_LamLam2.Cols(); j++) 
				if (fModuli_LamLam2(i,j) != 0.0)  
					output << "CLambdaLambda2("<< i << ","<< j <<"): " << fModuli_LamLam2(i,j) << endl;
		}
		output << endl;	
		
		output << "*******CLambdaU1******* " << endl;
		for (int i = 0; i < fModuli_LamU1.Rows(); i++)
		{
			for (int j = 0; j < fModuli_LamU1.Cols(); j++) 
				if (fModuli_LamU1(i,j) != 0.0)  
					output << "CLambdaU1("<< i << ","<< j <<"): " << fModuli_LamU1(i,j) << endl;
		}
		output << endl;	
		
		output << "*******CLambdaU2******* " << endl;
		for (int i = 0; i < fModuli_LamU2.Rows(); i++)
		{
			for (int j = 0; j < fModuli_LamU2.Cols(); j++) 
				if (fModuli_LamU2(i,j) != 0.0)  
					output << "CLambdaU2("<< i << ","<< j <<"): " << fModuli_LamU2(i,j) << endl;
		}
		output << endl;	
		output.close();
		//*************debug*********************//	
		return fModuli;
	}
	else 
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
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fStress
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
	fModuli.Dimension(kNSTR);
	fModuliPerfPlas.Dimension(kNSTR);
	fModuli_UU1.Dimension(kNSTR);
	fModuli_UU2.Dimension(kNSTR);
	fModuli_ULam1.Dimension(kNSTR, 1);
	fModuli_ULam2.Dimension(kNSTR, 1);
	fModuli_LamU1.Dimension(1, kNSTR);
	fModuli_LamU2.Dimension(1, kNSTR);
	fModuli_LamLam1.Dimension(1);
	fModuli_LamLam2.Dimension(1);
	fDevStress.Dimension(kNSD);
	fLapDevStress.Dimension(kNSD);
	fDevStrain.Dimension(kNSD);
	fLapDevStrain.Dimension(kNSD);
	fIniInternal.Dimension(4);
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
	fStress.Alias(dim, &d_array[  offset + dex]);     
	fInternal.Alias(kNumInternal, &d_array[3*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int GRAD_MRSSNLHardT::PlasticLoading(const dSymMatrixT& totalstrain, 
	  const dSymMatrixT& lap_totalstrain, const dArrayT& lambda, 
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
		
		/* collect initial internal variables */
		fIniInternal[0] = fchi;
		fIniInternal[1] = fc;
		fIniInternal[2] = ftan_phi;
		fIniInternal[3] = ftan_psi;
		
		/* check if yield condition is satisfied */
		double ff = YieldCondition(DeviatoricStress(totalstrain,lap_totalstrain,element),
			        MeanStress(totalstrain,lap_totalstrain,element),fchi,fc,ftan_phi);
		int iplastic;
		if (ff < fTol_1 && lambda[0] == 0.) /* elastic */
			iplastic = 0;
		else if (ff > fTol_1 && lambda[0] > 0.) /* plastic */
		    iplastic = 1;
		else
			ExceptionT::GeneralFail("GRAD_MRSSNLHardT::PlasticLoading","conditions for elastic/plastic state check not met");
		return iplastic;
		
        /* already plastic */
    }
	else 
	{
		/* get flags */
	 	iArrayT& Flags = element.IntegerData();
		
		/* load internal variables */
		LoadData(element, ip);
		
		/* first time plasticity is reached */
		if(fInternal[7] == 0.) {
			double enp  = 0.;
			double esp = 0.;
    		fInternal[kchi] = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
        	fInternal[kc]   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
    		fInternal[ktanphi] = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
			fInternal[ktanpsi] = (tan(fphi_p))*exp(-falpha_psi*esp);
			
		}
	
		/* calculate trial elastic strain */
		dSymMatrixT elasticstrain(3),lap_elasticstrain(3);
		elasticstrain.DiffOf(totalstrain, fPlasticStrain);
		lap_elasticstrain.DiffOf(lap_totalstrain, fLapPlasticStrain);
		fInternal[kftrial] = YieldCondition(DeviatoricStress(elasticstrain,lap_elasticstrain,element),
			   				MeanStress(elasticstrain,lap_elasticstrain,element),fInternal[kchi],
			             	fInternal[kc],fInternal[ktanphi]); 
		int iplastic; /* plasticity index */
		
		/* plastic */
		if (fInternal[kftrial] > fTol_1 && lambda[0] > 0.)
		{		
			/* set flag */
			Flags[ip] = kIsPlastic;
			iplastic = 1;
		}
		/* elastic */
		else if (fInternal[kftrial] < fTol_1 && lambda[0] == 0.)
		{
			/* set flag */
		    Flags[ip] = kIsElastic; //removed to avoid restting 7/01
			iplastic = 0;
		}
		else
			ExceptionT::GeneralFail("GRAD_MRSSNLHardT::PlasticLoading","conditions for elastic/plastic state check not met");
		return iplastic;
	}
}	

/* computes the stress corresponding to the given element
 * and elastic strain. The function returns a reference to the
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

double GRAD_MRSSNLHardT::signof(double r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

void GRAD_MRSSNLHardT::Contract4To2(const dMatrixT& mat,const dSymMatrixT& vec,
                                      dSymMatrixT& res)
{
	dSymMatrixT tmp(3);
	tmp = vec;
	for (int i=0; i<3; i++) {
	   tmp[i+3]=vec[i+3] *2.;
	}
	mat.Multx(tmp, res);
}