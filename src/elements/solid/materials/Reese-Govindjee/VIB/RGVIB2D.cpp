/* $Id: RGVIB2D.cpp,v 1.1 2003-03-19 19:00:57 thao Exp $ */
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
	Material2DT(in, kPlaneStress),
	ViscVIB(in, 2, 2, 3),
	fCircle(NULL),
        fb(2), 
        fEigs(2), 
        fEigs_e(2), 
        fsigA_E(2), 
        fdtauAdepB_E(2), 
        fsigA_I(2), 
        fdtauAdepB_I(2), 
        fCalg(2), 
        fModMat(3), 
        fModulus(3), 
        fStress(2), 
	fiKAB(2,2)
{
        if (PurePlaneStress())
                fconst =0.5;
        else
	        fconst = 1.0/3.0;

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
	Material2DT::Print(out);
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
        fJ = 1.0;
        dtauAdepB(fEigs, fsigA_E, fdtauAdepB_E, Elastic);
        dtauAdepB(fEigs, fsigA_I, fdtauAdepB_I, Inelastic);

        double lambda = fdtauAdepB_E(0,1)+ fdtauAdepB_I(0,1);
        double mu = 0.5*((fdtauAdepB_E(0,0) - fdtauAdepB_E(0,1))+
                        (fdtauAdepB_I(0,0) - fdtauAdepB_I(0,1)));

        if (fConstraintOption == Material2DT::kPlaneStress)
	  {
                if (PurePlaneStress())
			IsotropicT::Set_PurePlaneStress_mu_lambda(mu, lambda);
                else
		  {
                        lambda *= 2.0*mu/(lambda + 2.0*mu);
                        double kappa = lambda + 2.0/3.0*mu;
			IsotropicT::Set_mu_kappa(mu, kappa);
		  }
	  }
        else
	  {
                        double kappa = lambda + 2.0/3.0*mu;
			IsotropicT::Set_mu_kappa(mu, kappa);
	  }
}

double RGVIB2D::StrainEnergyDensity(void)
{
	/* Calculates eigenvalues of total stretch tensor */
	Compute_b(fb);
	fb.PrincipalValues(fEigs);
	fJ = fEigs.Product();

	/*eigenvalues of elastic stretch tensor*/
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	ComputeEigs_e(fEigs, fEigs_e, fsigA_I, fdtauAdepB_I);
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
        /*Assumes s_ij() has been calculated*/ 
  
        dtauAdepB(fEigs, fsigA_E, fdtauAdepB_E, Elastic); 
        dSymMatrixT gamAB_E = fdtauAdepB_E; 
        Calgorithm(fEigs_e, fsigA_I, fdtauAdepB_I, fCalg);
	dMatrixT gamAB_I =fCalg;
        double iJ = 1.0/fJ;
	gamAB_E *= iJ;
	gamAB_I *= iJ;

	gamAB_E(0,0) -= 2.0*fsigA_E[0];
	gamAB_E(1,1) -= 2.0*fsigA_E[1];
	gamAB_I(0,0) -= 2.0*fsigA_I[0];
	gamAB_I(1,1) -= 2.0*fsigA_I[1];

        /*axial*/ 
 
        /*add the elastic part*/ 
        fModulus = fSpectralDecompSpat.EigsToRank4(gamAB_E);       
        /*add the inelastic part*/ 
        fModulus += fSpectralDecompSpat.NonSymEigsToRank4(gamAB_I); 
 
        /* shear terms */ 
        const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors(); 
        const ArrayT<dArrayT>& eigenvectors2=fSpectralDecompRef.Eigenvectors(); 
        double dlamb, coeff; 
	double sig1 = fsigA_E[0]+fsigA_I[0]; 
	double sig2 = fsigA_E[1]+fsigA_I[1]; 
	
	double& lamb1 = fEigs[0]; 
	double& lamb2 = fEigs[1]; 
	
	dlamb = lamb1 - lamb2; 
	/* modulus coefficient */ 
	if (fabs(dlamb) > kSmall) 
	  coeff = (sig1*lamb2 - sig2*lamb1)/dlamb; 
	else 
	  coeff = 0.5*(gamAB_E(0,0) - gamAB_E(0,1) +  
		       gamAB_I(0,0) - gamAB_I(0,1));           
	MixedRank4_2D(eigenvectors[0], eigenvectors[1], fModMat); 
	fModulus.AddScaled(2.0*coeff, fModMat); 
	
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
        fJ = sqrt(fEigs.Product()) ; 
         
        /*load viscoelastic principal stretches from state variable arrays*/ 
        ElementCardT& element = CurrentElement(); 
        Load(element, CurrIP()); 
	if (fFSMatSupport.RunState() == GlobalT::kFormRHS) 
	  { 
                double Jvn = sqrt(fC_vn.Det());

		/*A value of 1.2 is currently hardwired for the cut-off of Jv 
		  for finite viscosity. The cut is abrupt and no "step down"
                  procedure is implemented.  Eventually the cut-off should be
                  moved to the viscosity function*/

                if (Jvn < 1.2) 
		  { 
                    const dMatrixT& F = F_mechanical(); 
                    dSymMatrixT b_tr(NumSD()); 
                    dSymMatrixT iCvn = fC_vn; 
                    iCvn.Inverse(); 
                    /*calculate trial state;*/ 
                    b_tr.MultQBQT(F,iCvn); 
                    fSpectralDecompSpat.SpectralDecomp_Jacobi(b_tr, false);      
                    /*set initial value of elastic principal stretches  
                     * to trial values*/ 
                    fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
                    ComputeEigs_e(fEigs, fEigs_e, fsigA_I, fdtauAdepB_I); 
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
        sigA(fEigs_e, fsigA_I, Inelastic); 
        sigA(fEigs, fsigA_E, Elastic); 
 
        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);    
        fStress = fSpectralDecompSpat.EigsToRank2(fsigA_E); 
        fStress += fSpectralDecompSpat.EigsToRank2(fsigA_I); 
 
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
         
        /*load the viscoelastic principal stretches from state variable arrays*/ 
        ElementCardT& element = CurrentElement(); 
        Load(element, CurrIP()); 
 
        fSpectralDecompRef.SpectralDecomp_Jacobi(fC_v, false);   
        dArrayT Eigs_v = fSpectralDecompRef.Eigenvalues(); 
         
        double Jv = sqrt(Eigs_v.Product()); 
        double Je = sqrt(fEigs.Product())/Jv; 
 
        output[0] = Jv; 
        output[1] = Je; 
         
        sigA(fEigs_e, fsigA_I, Inelastic); 
        sigA(fEigs, fsigA_E, Elastic); 

        fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false); 
        dSymMatrixT dev_stress(NumSD()); 
        dev_stress = fSpectralDecompSpat.EigsToRank2(fsigA_I); 
        double mean_stress = fconst*dev_stress.Sum();  
        dev_stress.PlusIdentity(-mean_stress); 
 
	fietaS = fShearVisc->Function(Jv,Je);
	fietaB = fBulkVisc->Function(Jv,Je);

        double rate_visc_disp = dev_stress.ScalarProduct()*0.5*fietaS+ 
          fconst*mean_stress*fconst*fietaB; 
         
        output[2] = rate_visc_disp*fFSMatSupport.TimeStep(); 
         
}  
/***********************************************************************
 * Protected
 ***********************************************************************/
void RGVIB2D::sigA(const dArrayT& eigenstretch, dArrayT& eigenstress, 
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
	double iJ = 1/fJ;
	s0 *= l0*iJ*fThickness;
	s1 *= l1*iJ*fThickness;
}

void RGVIB2D::dtauAdepB(const dArrayT& eigenstretch, dArrayT& eigenstress, 
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
	s0 *= l0*fThickness;
	s1 *= l1*fThickness;
	
	c00 *= l0*l0*fThickness;
	c11 *= l1*l1*fThickness;
	c01 *= l0*l1*fThickness;	
	
	c00 += 2.0*s0;
	c11 += 2.0*s1;

	double iJ = 1/fJ;
	s0 *= iJ;
	s1 *= iJ;
}

void RGVIB2D::Calgorithm(const dArrayT& eigenstretch, dArrayT& eigenstress,
	dSymMatrixT& eigenmodulus, dMatrixT& Calg)
{	
        /*get sigma and dsigma/dep_e*/
	dtauAdepB(eigenstretch, eigenstress, eigenmodulus, Inelastic);
	
	double& s0 = eigenstress[0];
	double& s1 = eigenstress[1];

	/*caculate means*/
	double sm = fconst*(s0+s1);

	/*calculate elastic and viscous parts of the jacobean*/
	double Je = sqrt(eigenstretch[0]*eigenstretch[1]);
	double Jv = fJ/Je;

	fietaS = fShearVisc->Function(Jv,Je);
	fietaB = fBulkVisc->Function(Jv,Je);

	fietaS = 1/fietaS;
	fietaB = 1/fietaB;
        /*Calculate derivatives of viscosities*/
        double DietaSDep = fietaS*fShearVisc->DFuncDJv(Jv, Je)*Jv;
        double DietaBDep = fietaB*fBulkVisc->DFuncDJv(Jv, Je)*Jv;

	/*evaluate KAB^-1 where 
	 *KAB = 1+dt D/Dep_e(sigA_Idev/nD+isostress/nV)*/
	ComputeiKAB(Je, Jv, eigenstress, eigenmodulus);

	/*GAB= 1 + dt D/Dep(sigA_Idev/nD+isostress/nV+deta_de*stress)*/
	dMatrixT GAB(2,2);

	double dt = fFSMatSupport.TimeStep();
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
	Calg.MultSymAB(eigenmodulus,DAB);
}

void RGVIB2D::ComputeEigs_e(const dArrayT& eigenstretch, 
			    dArrayT& eigenstretch_e, dArrayT& eigenstress, 
			    dSymMatrixT& eigenmodulus) 
{		
	const double ctol = 1.00e-10;
		
	/*set references to principle stretches*/
	double& l0 = eigenstretch[0];
	double& l1 = eigenstretch[1];
	
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	
	double ltr0 = le0;
	double ltr1 = le1;

	double lvn0 = l0/le0;
	double lvn1 = l1/le1;
	double Jv=sqrt(lvn0*lvn1);
	
	double tol;
	
	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(ltr0);
	double ep_tr1 = 0.5*log(ltr1);
	
	/*initial guess for ep_e*/
	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;
	
	double Je;
	int counter = 0;
	
	/*initializes principle viscous stretch*/
	do 
	{
	        counter ++;
		/*calculate stresses and moduli*/
		dtauAdepB(eigenstretch_e, eigenstress, eigenmodulus, Inelastic);
	    
		/*initialize references*/
		double& s0 = eigenstress[0];
		double& s1 = eigenstress[1];
	    
		/*caculate means*/
		double sm = fconst*(s0+s1);
	    
		/*caluclate inverse viscosities*/
		Je = sqrt(le0*le1);
		Jv = fJ/Je;
	    
		fietaS = 1.0/fShearVisc->Function(Jv, Je);
		fietaB = 1.0/fBulkVisc->Function(Jv, Je);
	    
		ComputeiKAB(Jv,Je,eigenstress,eigenmodulus);
	    
		/*calculate the residual*/
		double dt = fFSMatSupport.TimeStep();
		double res0 = ep_e0 + dt*(0.5*fietaS*(s0-sm) +
					   fconst*fietaB*sm) - ep_tr0;
		double res1 = ep_e1 + dt*(0.5*fietaS*(s1-sm) +
					   fconst*fietaB*sm) - ep_tr1;
		
		/*solve for the principal strain increments*/
		double dep_e0 = -fiKAB(0,0)*res0 - fiKAB(0,1)*res1;
		double dep_e1 = -fiKAB(1,0)*res0 - fiKAB(1,1)*res1;
	    
		/*updates principal elastic stretches*/ 
		ep_e0 += dep_e0;
		ep_e1 += dep_e1;

		le0 = exp(2*ep_e0);
		le1 = exp(2*ep_e1);
		//	cout << "\n depsilon1 "<< dep_e0;
		
		if (counter > 100)
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

void RGVIB2D::ComputeiKAB(double& Jv, double& Je, 
			       dArrayT& eigenstress, dSymMatrixT& eigenmodulus)
{	
        /*Calculate derivatives of viscosities*/
        double DietaSDep_e = fietaS*(fShearVisc->DFuncDJe(Jv, Je)*Je - 
			     fShearVisc->DFuncDJv(Jv, Je)*Jv);
        double DietaBDep_e = fietaB*(fBulkVisc->DFuncDJe(Jv, Je)*Je - 
			     fBulkVisc->DFuncDJv(Jv, Je)*Jv);
	
	//	       	cout <<"\n deta: "<< DietaBDep_e;
	/*Calculate d_sigA_I/d_eps_e*/
	double& s0 = eigenstress[0];
	double& s1 = eigenstress[1];
	
	double iJ = 1.0/fJ;
	
	double c0 = eigenmodulus[0]*iJ;
	double c1 = eigenmodulus[1]*iJ;
	double c01 = eigenmodulus[2]*iJ;

	/*mean of dtauAdepB_I*/
	double sm = fconst*(s0+s1);
	double cm0 = fconst*(c0 + c01);
	double cm1 = fconst*(c01 + c1);

	dMatrixT& KAB = fiKAB;
		
	/*calculates  KAB = 1+dt*D(sigA_Idev/nD+isostress/nV)/Dep_e*/

	double dt = fFSMatSupport.TimeStep();
	KAB(0,0) = 1+0.5*fietaS*dt*(c0-cm0-DietaSDep_e*(s0-sm))+
	             fconst*fietaB*dt*(cm0 - DietaBDep_e*sm);

	KAB(1,1) = 1+0.5*fietaS*dt*(c1-cm1-DietaSDep_e*(s1-sm))+
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
