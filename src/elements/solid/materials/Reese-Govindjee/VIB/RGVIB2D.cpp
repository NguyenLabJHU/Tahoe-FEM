/* $Id: RGVIB2D.cpp,v 1.3 2003-03-25 06:30:24 thao Exp $ */
/* created: TDN (01/22/2001) */

#include <math.h>
#include <iostream.h>
#include <stdlib.h>

#include "RGVIB2D.h"

#include "fstreamT.h"
#include "toolboxConstants.h"
#include "VariViscT.h"
#include "C1FunctionT.h"

/* point generator */
#include "EvenSpacePtsT.h"

#include <iostream.h>
#include <math.h>

using namespace Tahoe;

const int kNumOutputVar = 3;
static const char* Labels[kNumOutputVar] = {"Jv","Je","dW_disp"};

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
RGVIB2D::RGVIB2D(ifstreamT& in, const FSMatSupportT& support): 
	RGBaseT(in, support),
	ViscVIB(in, 2, 2, 3),
	fCircle(NULL),
        fb(2), 
        fEigs(2), 
        fEigs_e(2), 
        ftau_E(2), 
        fDtauDep_E(2), 
        ftau_I(2), 
        fDtauDep_I(2), 
        fModMat(3), 
        fModulus(3), 
        fStress(2), 
	fiKAB(2,2),
	fconst(0.5)
{
  	/* point generator */
	fCircle = new EvenSpacePtsT(in);

	/* set tables */
	Construct();
		
}

/* destructor */
RGVIB2D::~RGVIB2D(void)
{
	delete fCircle;
}

/* print parameters */
void RGVIB2D::Print(ostream& out) const
{
	/* inherited */
	RGBaseT::Print(out);
	ViscVIB::Print(out);

	fCircle->Print(out);
}

/* print name */
void RGVIB2D::PrintName(ostream& out) const
{
	/* inherited */
	RGBaseT::PrintName(out);
	ViscVIB::PrintName(out);
	out << "    2D\n";

	/* integration rule */
	fCircle->PrintName(out);
}

int RGVIB2D::NumOutputVariables() const {return kNumOutputVar;}

void RGVIB2D::OutputLabels(ArrayT<StringT>& labels) const
{
  //allocates space for labels
        labels.Dimension(kNumOutputVar);

        //copy labels
        for (int i = 0; i< kNumOutputVar; i++)
                labels[i] = Labels[i];
}
/* class specific initializations */
void RGVIB2D::Initialize(void)
{
  /* initial modulus */

        fEigs = 1.0;
        ddWddE(fEigs, ftau_E, fDtauDep_E, Elastic);
        ddWddE(fEigs, ftau_I, fDtauDep_I, Inelastic);

        double lambda = fDtauDep_E(0,1)+ fDtauDep_I(0,1);
        double mu = 0.5*((fDtauDep_E(0,0) - fDtauDep_E(0,1))+
                        (fDtauDep_I(0,0) - fDtauDep_I(0,1)));

	IsotropicT::Set_PurePlaneStress_mu_lambda(mu, lambda);
}

double RGVIB2D::StrainEnergyDensity(void)
{
	/* Calculates eigenvalues of total stretch tensor */
	Compute_b(fb);
	fb.PrincipalValues(fEigs);

	/*eigenvalues of elastic stretch tensor*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	ComputeEigs_e(fEigs, fEigs_e, ftau_I, fDtauDep_I);
	Store(element, CurrIP());
	
	/*get bond lengths*/
	ComputeLengths(fEigs, Elastic);
	ComputeLengths(fEigs_e, Inelastic);
	
	/*Update Potential Table*/
	fPotential_E->MapFunction(fLengths_E,fU_E);
	fPotential_I->MapFunction(fLengths_I,fU_I);
	
	double energy = 0.0;
	double* pU_E = fU_E.Pointer();
	double* pU_I = fU_I.Pointer();
	double* pj = fjacobian.Pointer();
	
	for (int i=0; i<fLengths_E.Length(); i++)
		energy += (*pj++) * ((*pU_E++)+(*pU_I++));
	
	return(energy);
}
/* modulus */ 
const dMatrixT& RGVIB2D::c_ijkl(void) 
{ 
  
        /* stretch tensor */ 
        Compute_b(fb); 
 
        /* spectral decomposition */ 
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);    

        /*calculates principal stretch and Jacobian determinant*/          
        fEigs = fSpectralDecompSpat.Eigenvalues(); 
	double J = sqrt(fEigs.Product());

        /*load the viscoelastic principal stretches from history variables*/ 
        ElementCardT& element = CurrentElement(); 
        Load(element, CurrIP()); 
	
	/*calc elastic principal stretch*/
        fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);   
        dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues(); 
	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;

        ddWddE(fEigs, ftau_E, fDtauDep_E, Elastic); 
        dSymMatrixT gamAB_E = fDtauDep_E; 
	gamAB_E(0,0) -= 2.0*ftau_E[0];
	gamAB_E(1,1) -= 2.0*ftau_E[1];

	dMatrixT Calg(2);
        Calgorithm(fEigs, fEigs_e, ftau_I, fDtauDep_I, Calg);
	dMatrixT gamAB_I = Calg;
	gamAB_I(0,0) -= 2.0*ftau_I[0];
	gamAB_I(1,1) -= 2.0*ftau_I[1];

        /*axial*/ 
        fModulus = fSpectralDecompSpat.EigsToRank4(gamAB_E);       
        fModulus += fSpectralDecompSpat.NonSymEigsToRank4(gamAB_I); 
 
        /* shear terms */ 
        const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors(); 
        double dl, coeff; 
	double tau0 = ftau_E[0]+ftau_I[0]; 
	double tau1 = ftau_E[1]+ftau_I[1]; 
	
	double& l0 = fEigs[0]; 
	double& l1 = fEigs[1]; 
	
	dl = l0 - l1; 
	/* modulus coefficient */ 
	if (fabs(dl) > kSmall) 
	  coeff = (tau0*l1 - tau1*l0)/dl; 
	else 
	{
	  //	  cout << "\nHere";
	  coeff = 0.5*(gamAB_E(0,0) - gamAB_E(0,1) +  
		       gamAB_I(0,0) - gamAB_I(0,1))-tau0;           
	}
	MixedRank4_2D(eigenvectors[0], eigenvectors[1], fModMat); 
	fModulus.AddScaled(2.0*coeff, fModMat); 

	fModulus /= J;	
	return fModulus; 
} 
const dSymMatrixT& RGVIB2D::s_ij(void) 
{ 
        /* stretch tensor */ 
        Compute_b(fb); 
 
        /* spectral decomposition */ 
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);    

        /*calculates principal stretch*/          
        fEigs = fSpectralDecompSpat.Eigenvalues(); 

        /*calc jacobian*/ 
        double J = sqrt(fEigs.Product()) ; 
         
        /*load viscoelastic principal stretches from history variables*/ 
        ElementCardT& element = CurrentElement(); 
        Load(element, CurrIP()); 

	if (fFSMatSupport.RunState() == GlobalT::kFormRHS) 
	{ 
	  double Jvn = sqrt(fC_vn.Det()); 
      	  if (Jvn < 1.2) 
       	  { 
	    SpectralDecompT SpectralDecompTrial(NumSD());
	    dSymMatrixT b_tr(NumSD()); 

	    /*calculate trial state;*/ 
	    dSymMatrixT iCvn = fC_vn; 
	    iCvn.Inverse(); 
	    const dMatrixT& F = F_mechanical(); 
	    b_tr.MultQBQT(F,iCvn);
	    SpectralDecompTrial.SpectralDecomp_Jacobi(b_tr, false);      

	    /*set initial value of elastic principal stretches to trial values*/ 
	    fEigs_e = SpectralDecompTrial.Eigenvalues(); 
	    ComputeEigs_e(fEigs, fEigs_e, ftau_I, fDtauDep_I); 

	    /*update viscuous stretch tensor*/ 
	    Compute_C(fC_v); 
	    fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v,false); 
	    dArrayT Eigs_v = fEigs; 
	    Eigs_v /= fEigs_e; 
	    fC_v = fSpectralDecompRef.EigsToRank2(Eigs_v); 
       	  } 
       	  else 
      	  { 
       	    fEigs_e = 1.0; 
       	    Compute_C(fC_v); 
       	  } 
	  Store(element, CurrIP()); 
	}        
        else  
	{ 
	  fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);   
	  dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues(); 
	  
	  fEigs_e = fEigs; 
	  fEigs_e /= Eigs_v; 
	} 
        dWdE(fEigs_e, ftau_I, Inelastic); 
        dWdE(fEigs, ftau_E, Elastic);
        fStress = fSpectralDecompSpat.EigsToRank2(ftau_E); 
        fStress += fSpectralDecompSpat.EigsToRank2(ftau_I); 
	fStress /= J;
	//	cout << "\nfEigs_e: "<<fEigs_e;
	//	cout << "\nfEigs: "<<fEigs;
 

        return fStress; 
} 
/* material description */ 
const dMatrixT& RGVIB2D::C_IJKL(void) 
{ 
  /* deformation gradient */ 
        const dMatrixT& Fmat = F(); 
         
        /* transform */ 
        fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl())); 
        return fModulus;         
} 
 
const dSymMatrixT& RGVIB2D::S_IJ(void) 
{ 
  /* deformation gradient */ 
        const dMatrixT& Fmat = F(); 
         
        /* transform */ 
        fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij())); 
        return fStress; 
} 
 
void RGVIB2D::ComputeOutput(dArrayT& output) 
{ 
        /* stretch tensor */ 
        Compute_b(fb); 
 
        /* spectral decomposition */ 
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);    
        /*calculates principal stretch*/ 
         
        fEigs = fSpectralDecompSpat.Eigenvalues(); 
	double J = sqrt(fEigs.Product());

        /*load the viscoelastic principal stretches from history variables*/ 
        ElementCardT& element = CurrentElement(); 
        Load(element, CurrIP()); 
 
        fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);   
        dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues(); 
	fEigs_e = fEigs;
	fEigs_e /= Eigs_v;

        double Je = sqrt(fEigs_e.Product()); 
	double Jv = J/Je;

        output[0] = Jv; 
        output[1] = Je; 
         
        dWdE(fEigs_e, ftau_I, Inelastic); 
        dWdE(fEigs, ftau_E, Elastic); 

        dSymMatrixT dev_stress(NumSD()); 
        dev_stress = fSpectralDecompSpat.EigsToRank2(ftau_I); 
	dev_stress /= J;

        double mean_stress = fconst*dev_stress.Sum();  
        dev_stress.PlusIdentity(-mean_stress); 
 
	fietaS = fShearVisc->Function(Jv,Je);
	fietaB = fBulkVisc->Function(Jv,Je);
	fietaS=1.0/fietaS;
	fietaB=1.0/fietaB;
        double rate_visc_disp = 0.5*(dev_stress.ScalarProduct()*0.5*fietaS+ 
          mean_stress*mean_stress*fietaB); 
         
        output[2] = rate_visc_disp*fFSMatSupport.TimeStep(); 
}  
/***********************************************************************
 * Protected
 ***********************************************************************/
void RGVIB2D::dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress, 
	int etype)
{
	/*Compute bond lengths*/
	ComputeLengths(eigenstretch, etype);
	
	double* pdU;
	double* pl;	
	/*Assign pointers to appropriate data storage*/
	if (etype == Elastic)
	{
		fPotential_E->MapDFunction(fLengths_E, fdU_E);
		pdU = fdU_E.Pointer();
		pl = fLengths_E.Pointer();
	}
	else 
	{ 	
		fPotential_I->MapDFunction(fLengths_I, fdU_I);
		pdU = fdU_I.Pointer();
		pl = fLengths_I.Pointer();
	}
	
	int length = fLengths_E.Length();
	double* pj = fjacobian.Pointer();
	
	/*Initialize kernel pointers*/
	double* pz0 = fStressTable(0);
	double* pz1 = fStressTable(1);
	
	/*Initialize stress pointer*/
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	
	/*stretch*/
	double& l0 = eigenstretch[0];
	double& l1 = eigenstretch[1];
	for (int i=0; i<length; i++)
	{
		double sfactor = (*pj++)*(*pdU++)/(*pl++);
	
		s0 += sfactor* (*pz0++);
		s1 += sfactor* (*pz1++);
	}	
	s0 *= l0;
	s1 *= l1;
}

void RGVIB2D::ddWddE(const dArrayT& eigenstretch, dArrayT& eigenstress, 
	dSymMatrixT& eigenmodulus, int etype)
{
	/*Compute bond lengths*/
	ComputeLengths(eigenstretch, etype);
	
	double* pdU;
	double* pddU;
	double* pl;		
	/*Assign pointers to appropriate data storage*/
	if (etype == Elastic)
	{
		fPotential_E->MapDFunction(fLengths_E, fdU_E);
		fPotential_E->MapDDFunction(fLengths_E, fddU_E);
		pdU = fdU_E.Pointer();
		pddU = fddU_E.Pointer();
		pl = fLengths_E.Pointer();
	}
	else
	{ 	
		fPotential_I->MapDFunction(fLengths_I, fdU_I);
		fPotential_I->MapDDFunction(fLengths_I, fddU_I);
		pdU = fdU_I.Pointer();
		pddU = fddU_I.Pointer();
		pl = fLengths_I.Pointer();
	}

	int length = fLengths_E.Length();	
	double* pj = fjacobian.Pointer();
	
	/*Initialize kernel pointers*/
	double* pz0 = fStressTable(0);
	double* pz1 = fStressTable(1);
	double* pz0z0 = fModuliTable(0);
	double* pz1z1 = fModuliTable(1);
	double* pz0z1 = fModuliTable(2);
	
	/*Initialize stress pointer*/
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	
	/*modulus*/
	double& c00 = eigenmodulus(0,0) = 0.0;
	double& c11 = eigenmodulus(1,1) = 0.0;
	double& c01 = eigenmodulus(0,1) = 0.0;
	
	/*stretch*/
	double& l0 = eigenstretch[0];
	double& l1 = eigenstretch[1];
	for (int i=0; i<length; i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj)*((*pddU)/((*pl)*(*pl))-
			(*pdU)/((*pl)*(*pl)*(*pl)));
		
		pj++;
		pl++;
		pdU++;
		pddU++;
		
		s0 += sfactor* (*pz0++);
		s1 += sfactor* (*pz1++);
		
		c00 += cfactor* (*pz0z0++);
		c11 += cfactor* (*pz1z1++);
		c01 += cfactor* (*pz0z1++);
	}
	s0 *= l0;
	s1 *= l1;
	
	c00 *= l0*l0;
	c11 *= l1*l1;
	c01 *= l0*l1;
	
	c00 += 2.0*s0;
	c11 += 2.0*s1;
}

void RGVIB2D::Calgorithm(const dArrayT& eigenstretch, const dArrayT& eigenstretch_e, 
			 dArrayT& eigenstress,dSymMatrixT& eigenmodulus, dMatrixT& Calg)
{	
	/*jacobian determinants*/
	double Je = sqrt(eigenstretch_e[0]*eigenstretch_e[1]);
	double J = sqrt(eigenstretch[0]*eigenstretch[1]);
	double Jv = J/Je;

        /*get sigma and dsigma/dep_e*/
	ddWddE(eigenstretch_e, eigenstress, eigenmodulus, Inelastic);
	
	double s0 = eigenstress[0]/J;
	double s1 = eigenstress[1]/J;

	/*caculate means*/
	double sm = fconst*(s0+s1);

	fietaS = 1.0/fShearVisc->Function(Jv,Je);
	fietaB = 1.0/fBulkVisc->Function(Jv,Je);

        /*Calculate derivatives of viscosities*/
        double DietaSDep = fietaS*fShearVisc->DFuncDJv(Jv, Je)*Jv;
        double DietaBDep = fietaB*fBulkVisc->DFuncDJv(Jv, Je)*Jv;

	/*evaluate KAB^-1 where KAB = 1+dt D/Dep_e(sigA_Idev/nD+isostress/nV)*/
	ComputeiKAB(J, Jv, eigenstress, eigenmodulus);

	/*GAB= 1 + dt D/Dep(sigA_Idev/nD+isostress/nV+deta_de*stress)*/
	double dt = fFSMatSupport.TimeStep();

	dMatrixT GAB(2,2);
	GAB(0,0) = 1+0.5*dt*fietaS*(s0-sm)*(1+DietaSDep)+
	               fconst*dt*fietaB*sm*(1+DietaBDep);
	GAB(1,1) = 1+0.5*dt*fietaS*(s1-sm)*(1+DietaSDep)+
	               fconst*dt*fietaB*sm*(1+DietaBDep);

	GAB(0,1) = 0.5*dt*fietaS*(s0-sm)*(1+DietaSDep)+
	             fconst*dt*fietaB*sm*(1+DietaBDep);
	GAB(1,0) = 0.5*dt*fietaS*(s1-sm)*(1+DietaSDep)+
	             fconst*dt*fietaB*sm*(1+DietaBDep);

	dMatrixT DAB(2);
	DAB.MultAB(fiKAB,GAB);

	Calg(0,0) = eigenmodulus(0,0)*DAB(0,0)+eigenmodulus(0,1)*DAB(1,0);
	Calg(1,0) = eigenmodulus(1,0)*DAB(0,0)+eigenmodulus(1,1)*DAB(1,0);
	Calg(0,1) = eigenmodulus(0,0)*DAB(0,1)+eigenmodulus(0,1)*DAB(1,1);
	Calg(1,1) = eigenmodulus(1,0)*DAB(0,1)+eigenmodulus(1,1)*DAB(1,1);
}

void RGVIB2D::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			    dArrayT& eigenstress, dSymMatrixT& eigenmodulus) 
{		
	const double ctol = 1.00e-14;
		
	/*set references to principle stretches*/
	double& l0 = eigenstretch[0];
	double& l1 = eigenstretch[1];
	double J = sqrt(l0*l1);
	double iJ = 1.0/J;

	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double Je = sqrt(le0*le1);
	
	double ltr0 = le0;
	double ltr1 = le1;

	double lvn0 = l0/le0;
	double lvn1 = l1/le1;
	double Jv=J/Je;
	
	double tol;
	
	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(ltr0);
	double ep_tr1 = 0.5*log(ltr1);
	
	/*initial guess for ep_e*/
	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;
	
	int counter = 0;
	
	/*initializes principle viscous stretch*/
	do 
	{
	        counter ++;
		/*calculate stresses and moduli*/
		ddWddE(eigenstretch_e, eigenstress, eigenmodulus, Inelastic);
	    
		/*initialize references*/
		double s0 = eigenstress[0]*iJ;
		double s1 = eigenstress[1]*iJ;
	    
		/*caculate means*/
		double sm = fconst*(s0+s1);
	    
		/*caluclate inverse viscosities*/
		Je = sqrt(le0*le1);
		Jv = J/Je;
	    
		fietaS = 1.0/fShearVisc->Function(Jv, Je);
		fietaB = 1.0/fBulkVisc->Function(Jv, Je);
	    
		ComputeiKAB(J,Je,eigenstress,eigenmodulus);
	    
		/*calculate the residual*/
		double dt = fFSMatSupport.TimeStep();
		double res0 = ep_e0 + dt*(0.5*fietaS*(s0-sm)+fconst*fietaB*sm) - ep_tr0;
		double res1 = ep_e1 + dt*(0.5*fietaS*(s1-sm)+fconst*fietaB*sm) - ep_tr1;
		
		/*solve for the principal strain increments*/
		double dep_e0 = -fiKAB(0,0)*res0 - fiKAB(0,1)*res1;
		double dep_e1 = -fiKAB(1,0)*res0 - fiKAB(1,1)*res1;
	    
		/*updates principal elastic stretches*/ 
		ep_e0 += dep_e0;
		ep_e1 += dep_e1;

		le0 = exp(2*ep_e0);
		le1 = exp(2*ep_e1);
		//	cout << "\n depsilon1 "<< dep_e0;
		
		if (counter > 300)
		{
	               ep_e0 = 0;
		       ep_e1 = 0;
		       le0 = 1.0;
		       le1 = 1.0;
		       counter = 0;
		       cout << "\nReset";
		}
	    
		/*Check that the L2 norm of the residual is less 
		 *than tolerance*/
		tol = sqrt(res0*res0 + res1*res1);
	}while (tol>ctol); 
}

/***********************************************************************
 * Private
 ***********************************************************************/
/* Initialize angle tables */

void RGVIB2D::ComputeiKAB(const double& J, const double& Je, 
			  const dArrayT& eigenstress, const dSymMatrixT& eigenmodulus)
{	
        /*Calculate derivatives of viscosities*/
        double Jv = J/Je;
        double DietaSDep_e = fietaS*(fShearVisc->DFuncDJe(Jv, Je)*Je - 
			     fShearVisc->DFuncDJv(Jv, Je)*Jv);
        double DietaBDep_e = fietaB*(fBulkVisc->DFuncDJe(Jv, Je)*Je - 
			     fBulkVisc->DFuncDJv(Jv, Je)*Jv);
	
	/*principal components of stress and modulus*/
	double iJ = 1.0/J;	
	double s0 = eigenstress[0]*iJ;
	double s1 = eigenstress[1]*iJ;
	
	double c0 = eigenmodulus[0]*iJ;
	double c1 = eigenmodulus[1]*iJ;
	double c01 = eigenmodulus[2]*iJ;

	/*mean values*/
	double sm = fconst*(s0+s1);
	double cm0 = fconst*(c0 + c01);
	double cm1 = fconst*(c01 + c1);

	dMatrixT& KAB = fiKAB;
		
	/*calculates  KAB = 1+dt*D(sig_Idev/nD+isostress/nV)/Dep_e*/

	double dt = fFSMatSupport.TimeStep();
	KAB(0,0) = 1.0+0.5*fietaS*dt*(c0-cm0-DietaSDep_e*(s0-sm))+
	             fconst*fietaB*dt*(cm0 - DietaBDep_e*sm);

	KAB(1,1) = 1.0+0.5*fietaS*dt*(c1-cm1-DietaSDep_e*(s1-sm))+
	             fconst*fietaB*dt*(cm1 - DietaBDep_e*sm);

	KAB(0,1) = 0.5*fietaS*dt*(c01-cm1-DietaSDep_e*(s0-sm))+
	             fconst*fietaB*dt*(cm1 - DietaBDep_e*sm);

	KAB(1,0) = 0.5*fietaS*dt*(c01-cm0-DietaSDep_e*(s1-sm))+
	             fconst*fietaB*dt*(cm0 - DietaBDep_e*sm);

	/*	cout <<"\n k1: "<< 0.5*fietaS*fdt*(c0-cm0)+
		fconst*fietaB*fdt*cm0+
		0.5*fietaS*fdt*(c01-cm1)+
		fconst*fietaB*fdt*cm1;

		cout <<"\n :c0 "<< c0;
		cout <<"\n :s0 "<< s0;
		
		cout <<"\n k2: "<< 0.5*fietaS*fdt*DietaSDep_e*(s0-sm)+
		fconst*fietaB*fdt*DietaBDep_e*sm+
		0.5*fietaS*fdt*DietaSDep_e*(s0-sm)+
		fconst*fietaB*fdt*DietaBDep_e*sm; */

	/*inverts KAB*/
	fiKAB.Inverse(KAB);
}

void RGVIB2D::ComputeLengths(const dArrayT& eigenstretch, int etype)
{
	double& e0 = eigenstretch[0];
	double& e1 = eigenstretch[1];
	
	double* pl;
	if(etype == Elastic)
		pl=fLengths_E.Pointer();
	else
		pl=fLengths_I.Pointer();

	/*sets pointers to bond directional vectors*/
	double* zo = fStressTable(0);
	double* z1 = fStressTable(1);
	
	for (int i=0; i<fLengths_E.Length(); i++)
		*pl++ = sqrt(e0* (*zo++) + e1* (*z1++));
}

void RGVIB2D::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Allocate(numpoints);
	
	/* fetch jacobians */
	fjacobian = fCircle->Jacobians();
	
	/* set pointers */
	double *s0 = fStressTable(0);
	double *s1 = fStressTable(1);

	double *c00 = fModuliTable(0);
	double *c11 = fModuliTable(1);
	double *c01 = fModuliTable(2);

	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s0[i] = cosi*cosi;
		s1[i] = sini*sini;
	
		/* moduli angle tables */
		c00[i] = s0[i]*s0[i];
		c11[i] = s1[i]*s1[i];
		c01[i] = s0[i]*s1[i];
	}
}
