/* $Id: BoyceViscoPlasticity.cpp,v 1.5 2007-03-10 18:32:48 tdnguye Exp $ */
/* created: TDN (01/22/2001) */

#include "BoyceViscoPlasticity.h"

#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
#include <iostream.h>
#include <stdlib.h>

using namespace Tahoe;

const double third = 1.0/3.0; 
const int kNumOutputVar =4; 
static const char* Labels[kNumOutputVar] = {"VM", "Sback","lp","Lang"}; 

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
BoyceViscoPlasticity::BoyceViscoPlasticity(void):
  ParameterInterfaceT("boyce_viscoplasticity"),
  fSpectralDecompSpat(3)
{}

/*initializes history variable */
void  BoyceViscoPlasticity::PointInitialize(void)
{
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
		      
			  fFv.Identity();
			  fFv_n.Identity();
			  *fs = fs0;
			  *fs_n = fs0;

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void BoyceViscoPlasticity::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fFv_n = fFv;
		*fs_n = *fs;
		
		/* write to storage */
		Store(element, ip);
	}
}

void BoyceViscoPlasticity::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fFv = fFv_n;
		*fs = *fs_n;
		
		/* write to storage */
		Store(element, ip);
	}
}

/* describe the parameters needed by the interface */
void BoyceViscoPlasticity::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  BoyceBaseT::DefineParameters(list);

  /* common limit */
  LimitT positive(0.0, LimitT::Lower);
  LimitT lower(1.0, LimitT::Lower);

  /* hardening law parameters */
  ParameterT ref_shear_rate(ParameterT::Double, "gamma_0_dot");
  ParameterT yield_stress(ParameterT::Double, "yield_stress");
  ParameterT sat_stress(ParameterT::Double, "saturation_stress");
  ParameterT hard_mod(ParameterT::Double, "hardening_modulus");
  ParameterT stress_stiffening(ParameterT::Double, "stress_stiffening");
  ParameterT pressure_dep(ParameterT::Double, "pressure_dependence_param");
  ParameterT temperature(ParameterT::Double, "ref_temperature");
  
  ref_shear_rate.AddLimit(positive);
  yield_stress.AddLimit(positive);
  sat_stress.AddLimit(positive);
  hard_mod.AddLimit(positive);
  stress_stiffening.AddLimit(positive);
  pressure_dep.AddLimit(positive);
  temperature.AddLimit(positive);

  list.AddParameter(ref_shear_rate);
  list.AddParameter(yield_stress);
  list.AddParameter(sat_stress);
  list.AddParameter(hard_mod);
  list.AddParameter(stress_stiffening);
  list.AddParameter(pressure_dep);
  list.AddParameter(temperature);

  /* back stress params */
  ParameterT back_stress_mod(ParameterT::Double, "back_stress_modulus");
  ParameterT limiting_stretch(ParameterT::Double, "limiting_stretch");
  
  back_stress_mod.AddLimit(positive);
  limiting_stretch.AddLimit(lower);

  list.AddParameter(back_stress_mod);
  list.AddParameter(limiting_stretch);
  
  /* elastic params - could make this a choice but just neo-Hookean for now */
  ParameterT mu(ParameterT::Double, "shear_modulus");
  ParameterT kappa(ParameterT::Double, "bulk_modulus");
  mu.AddLimit(positive);
  kappa.AddLimit(positive);
  list.AddParameter(mu);
  list.AddParameter(kappa);
}

void BoyceViscoPlasticity::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  /*allows one neq process: */
  BoyceBaseT::TakeParameterList(list);

  fgammadot0 = list.GetParameter("gamma_0_dot");
  fs0 = list.GetParameter("yield_stress");
  fs_ss = list.GetParameter("saturation_stress");
  fh = list.GetParameter("hardening_modulus");
  fA = list.GetParameter("stress_stiffening");
  falpha = list.GetParameter("pressure_dependence_param");
  fT = list.GetParameter("ref_temperature");

  fmu_R = list.GetParameter("back_stress_modulus");
  flambda_L = list.GetParameter("limiting_stretch");
  
  /* elastic moduli */
  fmu = list.GetParameter("shear_modulus");
  fkappa = list.GetParameter("bulk_modulus");
  
  /*Dimension work space*/
  Initialize();
}



double BoyceViscoPlasticity::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	const dMatrixT& F = DefGrad();
	
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
	/*elastic stretch*/
	fInverse = fFv;
	fInverse.Inverse();
	fFe.MultAB(F,fInverse);

	/*Be*/
	fStretch_e.MultAAT(fFe);

	double I1 = (fStretch_e.Trace());
	double I3 = fStretch_e.Det();
	/*deviatoric part*/
		
	double energy = 0.5*fmu*(I1-1.0) * 0.25*fkappa*(I3 - log(I3) - 1.0);
	return(energy);
}

/* modulus */
const dMatrixT& BoyceViscoPlasticity::c_ijkl(void)
{
    
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
	const dMatrixT& F = DefGrad();

	/*calculate eigenvalues of total stretch tensor*/
	fStretch.MultAAT(F);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();

	const double l0 = fEigs[0];
	const double l1 = fEigs[1];
	const double l2 = fEigs[2];

	/*calculate trial state*/
	fInverse = fFv_n;
	fInverse.Inverse();
	fFe.MultAB(F,fInverse); /*trial elastic deformation gradient*/
/*
		cout << "\nIP: "<<CurrIP();
		cout << "\nF: "<<F;
		cout << "\nfFe: "<<fFe;
		cout << "\nFv_n: "<<fFv_n;
*/
	/*eigenvalues of trial elastic stretch*/
	fStretch_e.MultAAT(fFe); /*Be_tr*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch_e, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues();

	const double l0_tr = fEigs_e[0];
	const double l1_tr = fEigs_e[1];
	const double l2_tr = fEigs_e[2];

	/*calculate elastic state*/
	fInverse = fFv;
	fInverse.Inverse();
	fFe.MultAB(F,fInverse); /*trial elastic deformation gradient*/

	/*eigenvalues of  elastic stretch*/
	fStretch_e.MultAAT(fFe); /*Be*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch_e, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues();
	const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();

	const double le0 = fEigs_e[0];
	const double le1 = fEigs_e[1];
	const double le2 = fEigs_e[2];
	
	/*calculate mean stretches*/
	double leb = third*(le0+le1+le2);
	double Je23 = pow(le0*le1*le2, -third);
	double J = sqrt(l0*l1*l2);

	/*Kirchoff stress*/	
	double T0 = fmu*Je23*(le0-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
	double T1 = fmu*Je23*(le1-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
	double T2 = fmu*Je23*(le2-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);	

	/*Calg_AB*/
	Compute_Calg(fEigs, fEigs_e);

	fCalg(0,0) -= 2.0*T0;
	fCalg(1,1) -= 2.0*T1;
	fCalg(2,2) -= 2.0*T2;
	
	fModulus3D = fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
		
	double dl_tr;
	double coeff;
	
	dl_tr = l0_tr - l1_tr;
	if (fabs(dl_tr) > kSmall)
		coeff = (T0*l1_tr - T1*l0_tr)/dl_tr;
	else 
		coeff = 0.5*(fCalg(0,0)-fCalg(0,1))-T0;
	MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[1], fModMat);
	fModulus3D.AddScaled(2.0*coeff, fModMat);
    
	dl_tr = l0_tr - l2_tr;
	if (fabs(dl_tr) > kSmall)
		coeff =(T0*l2_tr - T2*l0_tr)/dl_tr;
	else 
		coeff = 0.5*(fCalg(0,0)-fCalg(0,2))-T2;	
	MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[2], fModMat);
	fModulus3D.AddScaled(2.0*coeff, fModMat);
    
	dl_tr = l1_tr - l2_tr;
	if (fabs(dl_tr) > kSmall)
		coeff  = (T1*l2_tr - T2*l1_tr)/dl_tr;
	else
		coeff = 0.5*(fCalg(1,1)-fCalg(1,2))-T1;	
	MixedRank4_3D(eigenvectors_e[1], eigenvectors_e[2], fModMat);
	fModulus3D.AddScaled(2.0*coeff, fModMat);

	if (NumSD() == 2)
	{
		fModulus[0] = fModulus3D[0];
		fModulus[1] = fModulus3D[1];
		fModulus[2] = fModulus3D[5];

		fModulus[3] = fModulus3D[6];
		fModulus[4] = fModulus3D[7];
		fModulus[5] = fModulus3D[11];
		fModulus[6] = fModulus3D[30];
		fModulus[7] = fModulus3D[31];
		fModulus[8] = fModulus3D[35];
	}
	else fModulus = fModulus3D;

	const dMatrixT& Ftotal = F_total();	
	fModulus *= 1.0/Ftotal.Det();

    return fModulus;
}

/* stresses */
const dSymMatrixT& BoyceViscoPlasticity::s_ij(void)
{
	const dMatrixT& F = DefGrad();
	
    /*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {		
		/*calculate eigenvalues of total stretch tensor*/
		fStretch.MultAAT(F); /*B*/
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch, false);	
		fEigs = fSpectralDecompSpat.Eigenvalues();

		/*calculate trial state*/
		fInverse = fFv_n;
		fInverse.Inverse();
		fFe.MultAB(F,fInverse); /*trial elastic deformation gradient*/

		/*compute polar decomposition to obtain elatic rotation tensor*/
		fSpectralDecompSpat.PolarDecomp(fFe, fRe, fStretch_e,false);
		/*eigenvalues of trial elastic stretch*/
		fStretch_e.MultAAT(fFe);  /*Be*/
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch_e, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues();
		const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();

//		cout << "\nF: "<<F;
		/*calc elastic stretch tensor*/
		ComputeEigs_e(fEigs, fEigs_e);

		const double le0 = fEigs_e[0];
		const double le1 = fEigs_e[1];
		const double le2 = fEigs_e[2];
	
		/*calculate mean stretches*/
		double leb = third*(le0+le1+le2);
		double Je23 = pow(le0*le1*le2, -third);

		/*Kirchoff stress*/	
		fEigs_Stress[0] = fmu*Je23*(le0-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
		fEigs_Stress[1] = fmu*Je23*(le1-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
		fEigs_Stress[2] = fmu*Je23*(le2-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
		
		/*calculate Fv_n+1*/
		fEigs_e[0] = sqrt(fEigs_e[0]);
		fEigs_e[1] = sqrt(fEigs_e[1]);
		fEigs_e[2] = sqrt(fEigs_e[2]);
		fStretch_e = fSpectralDecompSpat.EigsToRank2(fEigs_e);

		/*calculate the elastic deformation gradient as Ff = Re Ue*/
		fStretch_e.ToMatrix(fInverse);
		fFe.MultAB(fRe, fInverse);		
	
		/*calculate plastic deformation gradient*/
		fInverse.Inverse(fFe);
		fFv.MultAB(fInverse, F);
		
		/*store updated solution*/	
		Store(element, CurrIP());
	}	
    else 
    {
		/*elastic stretch*/
		fInverse.Inverse(fFv);
		fFe.MultAB(F,fInverse);
		fStretch_e.MultAAT(fFe); /*Be*/
		fEigs_e = fSpectralDecompSpat.Eigenvalues();		
		const double le0 = fEigs_e[0];
		const double le1 = fEigs_e[1];
		const double le2 = fEigs_e[2];
	
		/*calculate mean stretches*/
		double leb = third*(le0+le1+le2);
		double Je23 = pow(le0*le1*le2, -third);

		/*Kirchoff stress*/	
		fEigs_Stress[0] = fmu*Je23*(le0-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
		fEigs_Stress[1] = fmu*Je23*(le1-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
		fEigs_Stress[2] = fmu*Je23*(le2-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
	}

	/*Assemble*/
	fStress3D = fSpectralDecompSpat.EigsToRank2(fEigs_Stress);

    if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];
        fStress[1] = fStress3D[1];
        fStress[2] = fStress3D[5];
    }
    else fStress = fStress3D;
	
	const dMatrixT& Ftotal = F_total();	
    fStress *= 1.0/Ftotal.Det();
	return fStress;
}

/* material description */
const dMatrixT& BoyceViscoPlasticity::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& BoyceViscoPlasticity::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
  
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

int BoyceViscoPlasticity::NumOutputVariables() const {return kNumOutputVar;} 

void BoyceViscoPlasticity::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 

void BoyceViscoPlasticity::ComputeOutput(dArrayT& output)
{
	const dMatrixT& F = DefGrad();
	
    /*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	/*total stretch*/
	fStretch.MultAAT(F); /*B*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues(); /*intermediate config*/		

	/*elastic stretch*/
	fInverse.Inverse(fFv);
	fFe.MultAB(F,fInverse);
	fStretch_e.MultAAT(fFe); /*Be*/
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fStretch_e, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues(); /*intermediate config*/		

	const double l0 = fEigs[0];		
	const double l1 = fEigs[1];		
	const double l2 = fEigs[2];		
	
	/*set references to principle stretches*/
	const double le0 = fEigs_e[0];
	const double le1 = fEigs_e[1];
	const double le2 = fEigs_e[2];

	/*plastic stretch*/
	double lv0 = l0/le0;
	double lv1 = l1/le1;
	double lv2 = l2/le2;

	/*calculate mean stretches*/
	double leb = third*(le0+le1+le2);
	double lvb = third*(lv0+lv1+lv2);
	double Je23 = pow(le0*le1*le2, -third);
		
	/*calculate back stress*/
	double r = sqrt(lvb)/flambda_L;
	double coeff = fmu_R/r*third;
		
	double Sback0 = coeff*fInvL.Function(r)*(l0/le0 - lvb);
	double Sback1 = coeff*fInvL.Function(r)*(l1/le1 - lvb);
	double Sback2 = coeff*fInvL.Function(r)*(l2/le2 - lvb);
	/*calculate driving stress*/
	double Te0 = fmu*Je23*(le0-leb)/le0 - Sback0;
	double Te1 = fmu*Je23*(le1-leb)/le1 - Sback1;
	double Te2 = fmu*Je23*(le2-leb)/le2 - Sback2;
	
	output[0] = sqrt(Te0*Te0 + Te1*Te1 + Te2*Te2);
	output[1] = sqrt(Sback0*Sback0 + Sback1*Sback1 + Sback2*Sback2);
	output[2] = sqrt(third*(lv0+lv1+lv2));
	output[3] = fInvL.Function(r);
/*	output[1] = Sback0;
	output[2] = Sback1;
	output[3] = Sback2;
	output[4] = lv0;
	output[5] = lv1;
	output[6] = lv2;*/
}

/***********************************************************************
 * Protected
 ***********************************************************************/
void BoyceViscoPlasticity::Initialize(void)
{
  /*initialize internal and history variables*/
  int ndof = 3;
  int numval = ndof*ndof;

  fnstatev = 0;
  fnstatev += numval;   /*current Fv*/
  fnstatev += numval;   /*last Fv_n*/
  fnstatev ++;			/*current shear strength*/
  fnstatev ++;          /*previous shear strength*/
  	
  fstatev.Dimension(fnstatev);
  double* pstatev = fstatev.Pointer();
		
  /* assign pointers to current and last blocks of state variable array */
  fFv.Set(ndof,ndof, pstatev);
  pstatev += numval;
  fFv_n.Set(ndof,ndof, pstatev);
  pstatev += numval;
  fs = pstatev++;
  fs_n = pstatev;

 /* dimension work space */
  fInverse.Dimension(3);
  fF3D.Dimension(3);
  fStretch.Dimension(3);

  fFe.Dimension(3);
  fStretch_e.Dimension(3);
  fRe.Dimension(3);

  fEigs.Dimension(3);
  fEigs_e.Dimension(3);
  fEigs_Stress.Dimension(3);

  fRes.Dimension(4);
  fdel.Dimension(4);
  
  fK.Dimension(4);
  fG.Dimension(3);
  fM.Dimension(3);
  fH.Dimension(3);
  fCalg.Dimension(3);
  fModulus3D.Dimension(6);
  fModMat.Dimension(6);
  fStress3D.Dimension(3);
  
  fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
  fStress.Dimension(NumSD());
  
  fModulus3D = 0.0;
}

/***********************************************************************
 * Private
 ***********************************************************************/
 void BoyceViscoPlasticity::Compute_Calg(const dArrayT& eigenstretch, const dArrayT& eigenstretch_e) 
 {
	const double& l0 = eigenstretch[0];
	const double& l1 = eigenstretch[1];
	const double& l2 = eigenstretch[2];

	const double& le0 = eigenstretch_e[0];
	const double& le1 = eigenstretch_e[1];
	const double& le2 = eigenstretch_e[2];
	

	/*Kirchoff stress*/	
	double leb = third*(le0+le1+le2);
	double Je23 = pow(le0*le1*le2, -third);

	double T0 = fmu*Je23*(le0-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
	double T1 = fmu*Je23*(le1-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
	double T2 = fmu*Je23*(le2-leb) + 0.5*fkappa*(le0*le1*le2 - 1.0);
	
	/*plastic stretch*/
	double lv0 = l0/le0;
	double lv1 = l1/le1;
	double lv2 = l2/le2;

	/*calculate mean stretches*/
	double lvb = third*(lv0+lv1+lv2);
	double J = sqrt(l0*l1*l2);
	
	/*calculate back stress*/
	double r = sqrt(lvb)/flambda_L;
	double coeff = fmu_R/r*third;
		
	double Sback0 = coeff*fInvL.Function(r)*(l0/le0 - lvb);
	double Sback1 = coeff*fInvL.Function(r)*(l1/le1 - lvb);
	double Sback2 = coeff*fInvL.Function(r)*(l2/le2 - lvb);
		
	/*calculate driving stress*/
	double Te0 = fmu*Je23*(le0-leb)/le0 - Sback0;
	double Te1 = fmu*Je23*(le1-leb)/le1 - Sback1;
	double Te2 = fmu*Je23*(le2-leb)/le2 - Sback2;
	/*norm of driving stress*/
	double tau = sqrt(0.5*(Te0*Te0 + Te1*Te1 + Te2*Te2));
	
	/*pressure*/
	double p = -0.5*fkappa/J*(le0*le1*le2 - 1.0);
		
	/*plastic stretch rate/driving stress norm*/
	double s_bar = *fs + falpha*p;
	double r56 = pow(tau/s_bar, 2.5*third);
	double r16 = pow(tau/s_bar, 0.5*third);
	double gammadot = 0;
	double f = 0;
	double g = 0;
	if (tau > kSmall) 
	{ 
		gammadot = fgammadot0*exp(-fA/fT*s_bar*(1.0-r56));
		g = fh*(1.0 - *fs/fs_ss)*gammadot;
		f = gammadot/(sqrt(2)*tau);
	} 
			
	/*calculate the residual*/
	double dt = fFSMatSupport->TimeStep();
		
	/*calculate stiffness matrix*/		
	double dgamdtau = 0.0;
	double dgamdp = 0.0;
	double dgamds = 0.0;
		
	double dfdtau = 0.0;
	double dfdp = 0.0;
	double dfds = 0.0;
				
	double dgdtau = 0.0;
	double dgdp = 0.0;
	double dgds = 0.0;
		
	if (tau > kSmall) 
	{ 
		dgamdtau = gammadot*(5.0*fA/(6.0*fT*r16));
		dgamdp = -gammadot*(5.0*fA*falpha*tau/(6.0*s_bar*fT*r16) + fA*falpha/fT*(1.0-r56));
		dgamds = -gammadot*(5.0*fA*tau/(6.0*s_bar*fT*r16) + fA/fT*(1.0-r56));

		dfdtau = (dgamdtau - gammadot/tau)/(sqrt(2.0)*tau);
		dfdp = dgamdp/(sqrt(2.0)*tau);
		dfds = dgamds/(sqrt(2.0)*tau);
				
		dgdtau = fh*(1.0-*fs/fs_ss)*dgamdtau;
		dgdp = fh*(1.0-*fs/fs_ss)*dgamdp;
		dgds = -fh/fs_ss*gammadot + fh*(1.0-*fs/fs_ss)*dgamds;
	} 
	/*********************partials with respect to lambda^e_B*******************/	
	/*dp = dp/dlambda^e_B lambda^e_B*/
	double dp  = -(fkappa/J) * (le0*le1*le2);
			
	/*dSe_AB = dSe_A/dlambda^e_B lambda^e_B*/
	double dSe00 = -4.0*third*third*fmu*Je23/le0 *(le0 - 2.0*le1 - 2.0*le2);
	double dSe11 = -4.0*third*third*fmu*Je23/le1 *(le1 - 2.0*le2 - 2.0*le0);
	double dSe22 = -4.0*third*third*fmu*Je23/le2 *(le2 - 2.0*le0 - 2.0*le1);
		
	double dSe01 = -2.0*third*third*fmu*Je23/le0 *(2.0*le1 + 2.0*le0 - le2);
	double dSe02 = -2.0*third*third*fmu*Je23/le0 *(2.0*le2 + 2.0*le0 - le1);
	
	double dSe10 = -2.0*third*third*fmu*Je23/le1 *(2.0*le0 + 2.0*le1 - le2);
	double dSe12 = -2.0*third*third*fmu*Je23/le1 *(2.0*le2 + 2.0*le1 - le0);
		
	double dSe20 = -2.0*third*third*fmu*Je23/le2 *(2.0*le0 + 2.0*le2 - le1);
	double dSe21 = -2.0*third*third*fmu*Je23/le2 *(2.0*le1 + 2.0*le2 - le0);
		
	/*dSb_AB = dSb_A/dlambda^e_B lambda^e_B*/
	double dSb00 = coeff*( (lv0-lvb)*(third*lv0/lvb * fInvL.Function(r)
				- third*lv0/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) - 4.0*third*lv0*fInvL.Function(r) );
	double dSb11 = coeff*( (lv1-lvb)*(third*lv1/lvb * fInvL.Function(r)
				- third*lv1/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) - 4.0*third*lv1*fInvL.Function(r) );
	double dSb22 = coeff*( (lv2-lvb)*(third*lv2/lvb * fInvL.Function(r)
				- third*lv2/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) - 4.0*third*lv2*fInvL.Function(r) );
					
	double dSb01 = coeff*( (lv0-lvb)*(third*lv1/lvb * fInvL.Function(r)
				- third*lv1/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv1*fInvL.Function(r) );
	double dSb02 = coeff*( (lv0-lvb)*(third*lv2/lvb * fInvL.Function(r)
				- third*lv2/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv2*fInvL.Function(r) );
									
	double dSb10 = coeff*( (lv1-lvb)*(third*lv0/lvb * fInvL.Function(r)
				- third*lv0/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv0*fInvL.Function(r) );
	double dSb12 = coeff*( (lv1-lvb)*(third*lv2/lvb * fInvL.Function(r)
				- third*lv2/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv2*fInvL.Function(r) );

	double dSb20 = coeff*( (lv2-lvb)*(third*lv0/lvb * fInvL.Function(r)
				- third*lv0/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv0*fInvL.Function(r) );
	double dSb21 = coeff*( (lv2-lvb)*(third*lv1/lvb * fInvL.Function(r)
				- third*lv1/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv1*fInvL.Function(r) );

	/*dtau_B = dtau/dlambda^e_B lambda^e_B*/
	double dtau0 = 0.0;					
	double dtau1 = 0.0;					
	double dtau2 = 0.0;	
	double ct = 0.0;
	double cp = 0.0;
	if (tau > kSmall) 
	{ 
		dtau0 = (Te0*(dSe00-dSb00) + Te1*(dSe10-dSb10) + Te2*(dSe20-dSb20))/(2*tau);					
		dtau1 = (Te0*(dSe01-dSb01) + Te1*(dSe11-dSb11) + Te2*(dSe21-dSb21))/(2*tau);					
		dtau2 = (Te0*(dSe02-dSb02) + Te1*(dSe12-dSb12) + Te2*(dSe22-dSb22))/(2*tau);	

		ct = dfdtau + dfds*(dt/(1.0 - dt*dgds))*dgdtau;
		cp = dfdp + dfds*(dt/(1.0 - dt*dgds))*dgdp;
	} 
	

	/*\sum_B H_AB delta \epsilon^e_B*/
	fH(0,0) = 1.0 + dt*(f*(dSe00-dSb00) + Te0*(ct*dtau0 + cp*dp));
	fH(1,1) = 1.0 + dt*(f*(dSe11-dSb11) + Te1*(ct*dtau1 + cp*dp));
	fH(2,2) = 1.0 + dt*(f*(dSe22-dSb22) + Te2*(ct*dtau2 + cp*dp));
	
	fH(0,1) = dt*(f*(dSe01-dSb01) + Te0*(ct*dtau1 + cp*dp));
	fH(0,2) = dt*(f*(dSe02-dSb02) + Te0*(ct*dtau2 + cp*dp));
	
	fH(1,0) = dt*(f*(dSe10-dSb10) + Te1*(ct*dtau0 + cp*dp));
	fH(1,2) = dt*(f*(dSe12-dSb12) + Te1*(ct*dtau2 + cp*dp));

	fH(2,0) = dt*(f*(dSe20-dSb20) + Te2*(ct*dtau0 + cp*dp));
	fH(2,1) = dt*(f*(dSe21-dSb21) + Te2*(ct*dtau1 + cp*dp));

/*	cout << "\neigs: "<<eigenstretch;
	cout << "\neigs_e: "<<eigenstretch_e;
	cout << "\nSe2: "<< Te2-Sback2;
	cout << "\nSback2: "<< Sback2;
	cout << "\ndtau2: "<<dtau2;
	cout << "\ndSe22: "<<dSe22;
	cout << "\ndSb22: "<<dSb22;
	cout << "\nfH: "<<fH;*/
	fH.Inverse();

	/*********************partials with respect to lambda_B*******************/	
	/*dp/dlambda_B*/
	dp  = 0.5*fkappa/J * (le0*le1*le2 - 1.0);
	
	/*dSb_AB = dSb_A/dlambda_B lambda_B = -dSb_A/dlambda^e_B lambda^e_B*/
	dSb00 *= -1.0;
	dSb11 *= -1.0;
	dSb22 *= -1.0;
						
	dSb01 *= -1.0;
	dSb02 *= -1.0;
									
	dSb10 *= -1.0;
	dSb12 *= -1.0;

	dSb20 *= -1.0;
	dSb21 *= -1.0;

	/*dtau_B = dtau/dlambda^e_B lambda^e_B*/
	if (tau > kSmall) 
	{ 
		dtau0 = -(Te0*dSb00 + Te1*dSb10 + Te2*dSb20)/(2*tau);					
		dtau1 = -(Te0*dSb01 + Te1*dSb11 + Te2*dSb21)/(2*tau);					
		dtau2 = -(Te0*dSb02 + Te1*dSb12 + Te2*dSb22)/(2*tau);					
	} 

	/*\sum_B G_AB delta delta\epsilon_B; delta \epsilon_B = \delta\epsilon^e_tr_B */
	fG(0,0) = 1.0 - dt*(f*(-dSb00) + Te0*(ct*dtau0 + cp*dp));
	fG(1,1) = 1.0 - dt*(f*(-dSb11) + Te1*(ct*dtau1 + cp*dp));
	fG(2,2) = 1.0 - dt*(f*(-dSb22) + Te2*(ct*dtau2 + cp*dp));
	
	fG(0,1) = -dt*(f*(-dSb01) + Te0*(ct*dtau1 + cp*dp));
	fG(0,2) = -dt*(f*(-dSb02) + Te0*(ct*dtau2 + cp*dp));
	
	fG(1,0) = -dt*(f*(-dSb10) + Te1*(ct*dtau0 + cp*dp));
	fG(1,2) = -dt*(f*(-dSb12) + Te1*(ct*dtau2 + cp*dp));

	fG(2,0) = -dt*(f*(-dSb20) + Te2*(ct*dtau0 + cp*dp));
	fG(2,1) = -dt*(f*(-dSb21) + Te2*(ct*dtau1 + cp*dp));

//	cout << "\nfG: "<<fG;
	/*dT_A/depsilon^e_B*/
	fM(0,0) = 2.0*fmu*third*third*Je23*(4.0*le0 + le1 + le2) + fkappa*(le0*le1*le2);
	fM(1,1) = 2.0*fmu*third*third*Je23*(4.0*le1 + le2 + le0) + fkappa*(le0*le1*le2);
	fM(2,2) = 2.0*fmu*third*third*Je23*(4.0*le2 + le0 + le1) + fkappa*(le0*le1*le2);
	
	fM(0,1) = 2.0*fmu*third*third*Je23*(-2.0*le0 - 2.0*le1 + le2) + fkappa*(le0*le1*le2);
	fM(0,2) = 2.0*fmu*third*third*Je23*(-2.0*le0 - 2.0*le2 + le1) + fkappa*(le0*le1*le2);
	fM(1,2) = 2.0*fmu*third*third*Je23*(-2.0*le1 - 2.0*le2 + le0) + fkappa*(le0*le1*le2);
	
	fM(1,0) = fM(0,1);
	fM(2,0) = fM(0,2);
	fM(2,1) = fM(1,2);
//	cout << "\nfM: "<<fM;
	/*Calg_AB*/
	fCalg.MultABC(fM, fH, fG);
}

void BoyceViscoPlasticity::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e) 
{		
	const double ctol = 1.00e-12;
	
	const double l0 = eigenstretch[0];		
	const double l1 = eigenstretch[1];		
	const double l2 = eigenstretch[2];		

	/*set references to principle stretches*/
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
  
	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(le0);
	double ep_tr1 = 0.5*log(le1);
	double ep_tr2 = 0.5*log(le2);

	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;
	
	double& s = *fs;
	const double& s_n = *fs_n;
	
	/*plastic stretch*/
	double lv0 = l0/le0;
	double lv1 = l1/le1;
	double lv2 = l2/le2;

	/*calculate mean stretches*/
	double leb = third*(le0+le1+le2);
	double lvb = third*(lv0+lv1+lv2);
	double Je23 = pow(le0*le1*le2, -third);
	double J = sqrt(l0*l1*l2);
	
	/*calculate back stress*/
	double r = sqrt(lvb)/flambda_L;
	double coeff = fmu_R/r*third;
		
	double Sback0 = coeff*fInvL.Function(r)*(l0/le0 - lvb);
	double Sback1 = coeff*fInvL.Function(r)*(l1/le1 - lvb);
	double Sback2 = coeff*fInvL.Function(r)*(l2/le2 - lvb);
		
	/*calculate driving stress*/
	double Te0 = fmu*Je23*(le0-leb)/le0 - Sback0;
	double Te1 = fmu*Je23*(le1-leb)/le1 - Sback1;
	double Te2 = fmu*Je23*(le2-leb)/le2 - Sback2;
	/*norm of driving stress*/
	double tau = sqrt(0.5*(Te0*Te0 + Te1*Te1 + Te2*Te2));
	
	/*pressure*/
	double p = -0.5*fkappa/J*(le0*le1*le2 - 1.0);
		
	/*plastic stretch rate/driving stress norm*/
	double s_bar = s + falpha*p;
	double r56 = pow(tau/s_bar, 2.5*third);
	double r16 = pow(tau/s_bar, 0.5*third);
	double gammadot = 0;
	double f = 0;
	double g = 0;
	if (tau > kSmall) 
	{ 
		gammadot = fgammadot0*exp(-fA/fT*s_bar*(1.0-r56));
		g = fh*(1.0 - s/fs_ss)*gammadot;
		f = gammadot/(sqrt(2)*tau);
	} 

	/*calculate the residual*/
	double dt = fFSMatSupport->TimeStep();

	fRes[0] = ep_e0 + dt*f*Te0 - ep_tr0;
	fRes[1] = ep_e1 + dt*f*Te1 - ep_tr1;
	fRes[2] = ep_e2 + dt*f*Te2 - ep_tr2;
	fRes[3] = s - dt*g - s_n;
	double tol = sqrt(fRes[0]*fRes[0] + fRes[1]*fRes[1] + fRes[2]*fRes[2] + fRes[3]*fRes[3]);

	int iteration  = 0;	
	int max_iteration = 20;

/*	cout << "\niteration 0:";
	cout << setprecision(16)<<"\ntotal stretch: "<<l0<<"\t"<<l1<<"\t"<<l2;
	cout <<"\ntrial stretch: "<<exp(2.0*ep_tr0)<<"\t"<<exp(2.0*ep_tr1)<<"\t"<<exp(2.0*ep_tr2);
	cout <<"\nelastic stress: "<<Te0+Sback0<<"\t"<<Te1+Sback1<<"\t"<<Te2+Sback2;
	cout <<"\nsn: "<<s_n;
	cout << "\nback stress: "<<Sback0<<"\t"<<Sback1<<"\t"<<Sback2;
	cout << "\ntau: "<<tau;
	cout << "\np: "<<p;
	cout << "\ns: "<<s;
	cout << "\nsbar: "<<s_bar;
	cout << "\ngammadot: "<<gammadot;
	cout << "\nf: "<<f;
	cout << "\nfRes: "<<fRes;
	cout << "\nerror: "<<tol;
*/
	/*initializes principle viscous stretch*/
	while (tol>ctol && iteration < max_iteration) 
	{
		iteration ++;

		/*calculate stiffness matrix*/		
		double dgamdtau = gammadot*(5.0*fA/(6.0*fT*r16));
		double dgamdp = -gammadot*(5.0*fA*falpha*tau/(6.0*s_bar*fT*r16) + fA*falpha/fT*(1.0-r56));
		double dgamds = -gammadot*(5.0*fA*tau/(6.0*s_bar*fT*r16) + fA/fT*(1.0-r56));
		
		double dfdtau = (dgamdtau - gammadot/tau)/(sqrt(2.0)*tau);
		double dfdp = dgamdp/(sqrt(2.0)*tau);
		double dfds = dgamds/(sqrt(2.0)*tau);
				
		double dgdtau = fh*(1.0-s/fs_ss)*dgamdtau;
		double dgdp = fh*(1.0-s/fs_ss)*dgamdp;
		double dgds = -fh/fs_ss*gammadot + fh*(1.0-s/fs_ss)*dgamds;
		
		/*dp = dp/dlambda^e_B lambda^e_B*/
		double dp  = -(fkappa/J) * (le0*le1*le2);
			
		/*dSe_AB = dSe_A/dlambda^e_B lambda^e_B*/
		double dSe00 = -4.0*third*third*fmu*Je23/le0 *(le0 - 2.0*le1 - 2.0*le2);
		double dSe11 = -4.0*third*third*fmu*Je23/le1 *(le1 - 2.0*le2 - 2.0*le0);
		double dSe22 = -4.0*third*third*fmu*Je23/le2 *(le2 - 2.0*le0 - 2.0*le1);
		
		double dSe01 = -2.0*third*third*fmu*Je23/le0 *(2.0*le1 + 2.0*le0 - le2);
		double dSe02 = -2.0*third*third*fmu*Je23/le0 *(2.0*le2 + 2.0*le0 - le1);
	
		double dSe10 = -2.0*third*third*fmu*Je23/le1 *(2.0*le0 + 2.0*le1 - le2);
		double dSe12 = -2.0*third*third*fmu*Je23/le1 *(2.0*le2 + 2.0*le1 - le0);
		
		double dSe20 = -2.0*third*third*fmu*Je23/le2 *(2.0*le0 + 2.0*le2 - le1);
		double dSe21 = -2.0*third*third*fmu*Je23/le2 *(2.0*le1 + 2.0*le2 - le0);
		
		/*dSb_AB = dSb_A/dlambda^e_B lambda^e_B*/
		double dSb00 = coeff*( (lv0-lvb)*(third*lv0/lvb * fInvL.Function(r)
					- third*lv0/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) - 4.0*third*lv0*fInvL.Function(r) );
		double dSb11 = coeff*( (lv1-lvb)*(third*lv1/lvb * fInvL.Function(r)
					- third*lv1/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) - 4.0*third*lv1*fInvL.Function(r) );
		double dSb22 = coeff*( (lv2-lvb)*(third*lv2/lvb * fInvL.Function(r)
					- third*lv2/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) - 4.0*third*lv2*fInvL.Function(r) );
					
		double dSb01 = coeff*( (lv0-lvb)*(third*lv1/lvb * fInvL.Function(r)
					- third*lv1/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv1*fInvL.Function(r) );
		double dSb02 = coeff*( (lv0-lvb)*(third*lv2/lvb * fInvL.Function(r)
					- third*lv2/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv2*fInvL.Function(r) );
									
		double dSb10 = coeff*( (lv1-lvb)*(third*lv0/lvb * fInvL.Function(r)
					- third*lv0/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv0*fInvL.Function(r) );
		double dSb12 = coeff*( (lv1-lvb)*(third*lv2/lvb * fInvL.Function(r)
					- third*lv2/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv2*fInvL.Function(r) );

		double dSb20 = coeff*( (lv2-lvb)*(third*lv0/lvb * fInvL.Function(r)
					- third*lv0/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv0*fInvL.Function(r) );
		double dSb21 = coeff*( (lv2-lvb)*(third*lv1/lvb * fInvL.Function(r)
					- third*lv1/(flambda_L*sqrt(lvb))*fInvL.DFunction(r)) + 2.0*third*lv1*fInvL.Function(r) );

		/*dtau_B = dtau/dlambda^e_B lambda^e_B*/
		double dtau0 = (Te0*(dSe00-dSb00) + Te1*(dSe10-dSb10) + Te2*(dSe20-dSb20))/(2*tau);					
		double dtau1 = (Te0*(dSe01-dSb01) + Te1*(dSe11-dSb11) + Te2*(dSe21-dSb21))/(2*tau);					
		double dtau2 = (Te0*(dSe02-dSb02) + Te1*(dSe12-dSb12) + Te2*(dSe22-dSb22))/(2*tau);					
		
		/*Stiffness matrix*/
		fK(0,0) = 1.0 + dt*(f*(dSe00-dSb00) + (dfdtau*dtau0 + dfdp*dp)*Te0);
		fK(1,1) = 1.0 + dt*(f*(dSe11-dSb11) + (dfdtau*dtau1 + dfdp*dp)*Te1);
		fK(2,2) = 1.0 + dt*(f*(dSe22-dSb22) + (dfdtau*dtau2 + dfdp*dp)*Te2);
		fK(3,3) = 1.0 - dt*dgds;

		fK(0,1) =  dt*(f*(dSe01-dSb01) + (dfdtau*dtau1 + dfdp*dp)*Te0);
		fK(0,2) =  dt*(f*(dSe02-dSb02) + (dfdtau*dtau2 + dfdp*dp)*Te0);

		fK(1,0) =  dt*(f*(dSe10-dSb10) + (dfdtau*dtau0 + dfdp*dp)*Te1);
		fK(1,2) =  dt*(f*(dSe12-dSb12) + (dfdtau*dtau2 + dfdp*dp)*Te1);

		fK(2,0) =  dt*(f*(dSe20-dSb20) + (dfdtau*dtau0 + dfdp*dp)*Te2);
		fK(2,1) =  dt*(f*(dSe21-dSb21) + (dfdtau*dtau1 + dfdp*dp)*Te2);

		fK(3,0) = -dt*(dgdtau*dtau0 + dgdp*dp);
		fK(3,1) = -dt*(dgdtau*dtau1 + dgdp*dp);
		fK(3,2) = -dt*(dgdtau*dtau2 + dgdp*dp);

		fK(0,3) = dt*dfds*Te0;
		fK(1,3) = dt*dfds*Te1;
		fK(2,3) = dt*dfds*Te2;
	
/*		cout << "\ndp: "<<dp;
		cout << "\ndtau: "<<dtau0<<"\t"<<dtau1<<"\t"<<dtau2;
		cout << "\ndSe: "<<dSe00<<"\t"<<dSe01<< "\t"<<dSe02<<"\n"<<dSe10<<"\t"<<dSe11<<"\t"<<dSe12
			<<"\n"<<dSe20<<"\t"<<dSe21<<"\t"<<dSe22;
	cout << "\nlang: "<<r<<"\t"<<fInvL.Function(r)<<"\t"<<fInvL.DFunction(r);
		cout << "\ndSb: "<<dSb00<<"\t"<<dSb01<< "\t"<<dSb02<<"\n"<<dSb10<<"\t"<<dSb11<<"\t"<<dSb12
			<<"\n"<<dSb20<<"\t"<<dSb21<<"\t"<<dSb22;
		cout << "\nfK: "<<fK;
*/

	    /*solve for the principal strain increments*/
		fK.Inverse();
		fK.Multx(fRes,fdel,-1.0);
		
	    /*updates principal elastic stretches*/ 
	    ep_e0 += fdel[0];
	    ep_e1 += fdel[1];
	    ep_e2 += fdel[2];
		s += fdel[3];
		
		
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
	    
/*		cout << "\nfdel: "<<fdel;
		cout << "\neig: "<<le0<<"\t"<<le1<<"\t"<<le2;
		cout <<"\ns: "<<s;
*/
		/*update plastic stretch*/
		lv0 = l0/le0;
		lv1 = l1/le1;
		lv2 = l2/le2;

		/*calculate mean stretches*/
		leb = third*(le0+le1+le2);
		lvb = third*(lv0+lv1+lv2);
		Je23 = pow(le0*le1*le2, -third);
		J = sqrt(l0*l1*l2);
		
		/*calculate back stress*/
		r = sqrt(lvb)/flambda_L;
		coeff = fmu_R/r*third;
		
		Sback0 = coeff*fInvL.Function(r)*(l0/le0 - lvb);
		Sback1 = coeff*fInvL.Function(r)*(l1/le1 - lvb);
		Sback2 = coeff*fInvL.Function(r)*(l2/le2 - lvb);
		
		/*calculate driving stress*/
		Te0 = fmu*Je23*(le0-leb)/le0 - Sback0;
		Te1 = fmu*Je23*(le1-leb)/le1 - Sback1;
		Te2 = fmu*Je23*(le2-leb)/le2 - Sback2;
		
		/*norm of driving stress*/
		tau = sqrt(0.5*(Te0*Te0 + Te1*Te1 + Te2*Te2));
	
		/*pressure*/
		p = -0.5*fkappa/J*(le0*le1*le2 - 1.0);
		
		/*plastic stretch rate/driving stress norm*/
		s_bar = s + falpha*p;
		r56 = pow(tau/s_bar, 2.5*third);
		r16 = pow(tau/s_bar, 0.5*third);
		gammadot = fgammadot0*exp(-fA/fT*s_bar*(1.0-r56));
		g = fh*(1.0 - s/fs_ss)*gammadot;
		f;
		(tau > kSmall) ? f = gammadot/(sqrt(2)*tau):f=0;

	    /*update the residual*/
	    fRes[0] = ep_e0 + dt*f*Te0 - ep_tr0;
	    fRes[1] = ep_e1 + dt*f*Te1 - ep_tr1;
	    fRes[2] = ep_e2 + dt*f*Te2 - ep_tr2;
		fRes[3] = s - dt*g - s_n;

	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(fRes[0]*fRes[0] + fRes[1]*fRes[1] + fRes[2]*fRes[2] + fRes[3]*fRes[3]);

/*	cout <<"\nelastic stretch: "<<le0<<"\t"<<le1<<"\t"<<le2;
	cout <<"\nelastic stress: "<<Te0+Sback0<<"\t"<<Te1+Sback1<<"\t"<<Te2+Sback2;
	cout << "\nlang: "<<r<<"\t"<<fInvL.Function(r);
	cout << "\nback stress: "<<Sback0<<"\t"<<Sback1<<"\t"<<Sback2;
	cout << "\ntau: "<<tau;
	cout << "\np: "<<p;
	cout << "\ns: "<<s;
	cout << "\nsbar: "<<s_bar;
	cout << "\ngammadot: "<<gammadot;
	cout << "\nf: "<<f;
	cout << "\nfRes: "<<fRes;
	cout <<"\n\niteration: "<<iteration<<"\terror: "<<tol;
*/
	} 
	if (iteration >= max_iteration) 
		ExceptionT::GeneralFail("BoyceViscoPlasticity::ComputeEigs_e", 
			"number of iteration exceeds maximum of 10");
}

dMatrixT BoyceViscoPlasticity::DefGrad(void)
{
	const dMatrixT& F = F_mechanical();
	if (NumSD() == 2)
	{
		fF3D[0] = F[0];
		fF3D[1] = F[1];
		fF3D[2] = 0.0;
	    
		fF3D[3] = F[2];
		fF3D[4] = F[3];
		fF3D[5] = 0.0;
	    
		fF3D[6] = 0.0;
		fF3D[7] = 0.0;
		fF3D[8] = 1.0;
	}
	else fF3D = F;
	
	return(fF3D);
}
