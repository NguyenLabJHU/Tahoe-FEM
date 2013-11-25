/* $Id: ElasticHydrogelSuoT.cpp,v 1.1 2013-11-25 22:24:31 tahoe.xiaorui Exp $ */
/* created : RX (1/5/2012) */
#include "ElasticHydrogelSuoT.h"
#include "ParameterContainerT.h"

#include "ExceptionT.h"

#include "NeoHookean.h"
#include "MooneyRivlin.h"
#include "ArrudaBoyce.h"
#include <cmath>

using namespace Tahoe;

const double third = 1.0/3.0; 
/*const int kNumOutputVar =1;*/
const double R=8.31;
const double upsilon=1.80e-5;
const int kNumOutputVar =2;
static const char* Labels[kNumOutputVar] = {"phi","T"};


/***********************************************************************
 * Public
 ***********************************************************************/

/* constructor */
ElasticHydrogelSuoT::ElasticHydrogelSuoT(void):
	FSSolidMatT(),
	ParameterInterfaceT("ElasticHydrogelSuoT"),
	fSpectralDecompSpat(3)
{
}

int ElasticHydrogelSuoT::NumOutputVariables() const {return kNumOutputVar;} 

void ElasticHydrogelSuoT::OutputLabels(ArrayT<StringT>& labels) const 
{ 
	/*allocates space for labels*/
	labels.Dimension(kNumOutputVar); 
	
	/*copy labels*/
	for (int i = 0; i< kNumOutputVar; i++) 
		labels[i] = Labels[i]; 
} 


double ElasticHydrogelSuoT::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	const dMatrixT& F = F_total();
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
    double alpha=pow(solidfraction0,-1.0*third);
	
	/*principal stretch*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	double fmu=fPot->GetMu();
	double J=sqrt(fEigs.Product());
	double energy = 0.0;
	double I1 = fEigs[0]+fEigs[1]+fEigs[2];
	energy=0.5*fmu*(alpha*alpha*I1-3-2.0*log(alpha*alpha*alpha*J));
	
	ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
    fEigs_e = fEigs;
	fEigs_e *= alpha*alpha*pow(*fSolidFraction,2.0*third);
	
	/*deviatoric part*/
	double Je=sqrt(fEigs_e.Product());
    energy +=fPot->MeanEnergy(Je);
    
    double fT=Compute_Temperature();
    double fchi=(xL+xH)/2.0+(xH-xL)/2.0*(1.0-exp(-2.0*(fT-Ttran)/DeltaT))/(1.0+exp(-2.0*(fT-Ttran)/DeltaT));
    double phi=*fSolidFraction;

    double energymix= R*fT/(upsilon*(phi))*((1.0-phi)*log(1.0-phi)+fchi*phi*(1-phi));
    energy +=energymix;
	return(energy);
}

/* modulus */
const dMatrixT& ElasticHydrogelSuoT::c_ijkl(void)
{
	const dMatrixT& F =F_total();
	
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
    
   /* if(CurrElementNumber()==0&&CurrIP()==0)
        cout <<"\n the deformaiton gradient is:  "<<F;*/
	/*principal stretch*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	//the following is RX added
	const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();
	
	/*elastic part*/
	//TO DO:: call function that calculates solid fraction
	//the following is RX add 
	double J=sqrt(fEigs.Product());
	
   /*load the history variables*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    /*if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
	{
    	CalculatePhi(*fSolidFraction, *fSolidFraction_n, J);
		Store(element, CurrIP());
	}*/
    CalculatePhi(*fSolidFraction, *fSolidFraction_n, J);
    Store(element, CurrIP());


	fEigs_e = fEigs;
    double alpha=pow(solidfraction0,-1.0*third);
	fEigs_e *= alpha*alpha*pow(*fSolidFraction,2.0*third);
	
	/*deviatoric part*/
	double Je = sqrt(fEigs_e.Product());
	fEigs_dev = fEigs_e;
	fEigs_dev *= pow(Je, -2.0*third);
	
    /*fPot->DevStress(fEigs_dev, ftau);
    ftau += fPot->MeanStress(Je);  */  
     double fmu=fPot->GetMu();
	 int nsd = ftau.Length();
     ftau[0] = fmu*(alpha*alpha*fEigs[0]-1.0);
     ftau[1] = fmu*(alpha*alpha*fEigs[1]-1.0);
  
    if (nsd == 3)
    ftau[2] = fmu*(alpha*alpha*fEigs[2]-1.0);
	ftau += fPot->MeanStress(Je);
	
    //cout<<"\n the stress is: "<< ftau[0];
	//TO DO:: calculate Gamma
    //fPot->DevMod(fEigs_dev,fDtauDe);
	    double fT=Compute_Temperature();
    double fchi=(xL+xH)/2.0+(xH-xL)/2.0*(1.0-exp(-2.0*(fT-Ttran)/DeltaT))/(1.0+exp(-2.0*(fT-Ttran)/DeltaT));
	double fkappa=fPot->GetKappa();
	double factor=pow(alpha,6.0)*fkappa*upsilon*(*fSolidFraction)*(*fSolidFraction)*J*J/(R*fT*(1/(*fSolidFraction-1)+1+2*fchi*(*fSolidFraction))-fkappa*upsilon/2*(3*pow(alpha,6.0)*(*fSolidFraction)*(*fSolidFraction)*J*J-1));
	fDtauDe[0]=2.0*fmu*fEigs[0]*alpha*alpha;
	fDtauDe[1]=2.0*fmu*fEigs[1]*alpha*alpha;
	if (nsd == 2)
  {
    fDtauDe[2]=0;
  }
  else 
  {
   fDtauDe[2]=2.0*fmu*fEigs[2]*alpha*alpha;
    fDtauDe[3]=0;
	 fDtauDe[4]=0;
	  fDtauDe[5]=0;
  }

	/*double factor=fkappa*mu*(*fSolidFraction)*J/(R*fT*(1/(*fSolidFraction-1)+1+2*fchi*(*fSolidFraction))-fkappa*mu*(2*(*fSolidFraction)*J-1));*/
	
    fDtauDe += fPot->MeanMod(Je)*(1.0+factor);
	
    dSymMatrixT& Gamma = fDtauDe;
    Gamma(0,0) -= 2.0*ftau[0];
    Gamma(1,1) -= 2.0*ftau[1];
    Gamma(2,2) -= 2.0*ftau[2];
	
	
	fModulus3D = fSpectralDecompSpat.EigsToRank4(Gamma);	
	double dl, coeff;
	//	cout << "\nfModulus3D: "<<fModulus3D;
	
    double& l0 = fEigs[0];
    double& l1 = fEigs[1];
    double& l2 = fEigs[2];
	
	dl = l0 - l1;
    if (fabs(dl) > kSmall)
		coeff = (ftau[0]*l1 - ftau[1]*l0)/dl;
    else
		coeff = 0.5*(Gamma(0,0)-Gamma(0,1))-ftau[0];
    MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    dl = l0 - l2;
    if (fabs(dl) > kSmall)
		coeff = (ftau[0]*l2 - ftau[2]*l0)/dl;
    else 
		coeff = 0.5*(Gamma(0,0)-Gamma(0,2))-ftau[2];	
    MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
    fModulus3D.AddScaled(2.0*coeff, fModMat);
    
    dl = l1 - l2;
	if (fabs(dl) > kSmall)
		coeff  = (ftau[1]*l2 - ftau[2]*l1)/dl;
    else
		coeff = 0.5*(Gamma(1,1)-Gamma(1,2))-ftau[1];	
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
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
const dSymMatrixT& ElasticHydrogelSuoT::s_ij(void)
{
	const dMatrixT& F = F_total();
	
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

	/*principal stretch*/
	fb.MultAAT(fF3D);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	
	/*elastic part*/
	//TO DO:: call function that calculates solid fraction
	//the following is RX add 
	double J=sqrt(fEigs.Product());
	
    /*load the history variables*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
   /* if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
	{
    	CalculatePhi(*fSolidFraction, *fSolidFraction_n, J);
		Store(element, CurrIP());
	}*/
	  	CalculatePhi(*fSolidFraction, *fSolidFraction_n, J);
		Store(element, CurrIP());
	
	
	fEigs_e = fEigs;
    double alpha=pow(solidfraction0,-1.0*third);
	fEigs_e *= alpha*alpha*pow(*fSolidFraction,2.0*third);
	
	/*deviatoric part*/
	double Je = sqrt(fEigs_e.Product());
	fEigs_dev = fEigs_e;
	fEigs_dev *= pow(Je, -2.0*third);
	
	/*fPot->DevStress(fEigs_dev, ftau);*/
	 double fmu=fPot->GetMu();
	 int nsd = ftau.Length();
     ftau[0] = fmu*(alpha*alpha*fEigs[0]-1.0);
     ftau[1] = fmu*(alpha*alpha*fEigs[1]-1.0);
  
    if (nsd == 3)
    ftau[2] = fmu*(alpha*alpha*fEigs[2]-1.0);

	ftau += fPot->MeanStress(Je);
	
	fStress3D = fSpectralDecompSpat.EigsToRank2(ftau);

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
const dMatrixT& ElasticHydrogelSuoT::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
	
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& ElasticHydrogelSuoT::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F_total();
	
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

void ElasticHydrogelSuoT::ComputeOutput(dArrayT& output)
{
   /*load the history variables*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
	output[0] = *fSolidFraction;
    double Temp=Compute_Temperature();
    output[1]=Temp-273.0;
}

double ElasticHydrogelSuoT::Pressure(void) const
{
    return fStress.Trace()/3.0;
}


/***********************************************************************
 * Protected
 ***********************************************************************/
//TO DO: Redo interface
/* information about subordinate parameter lists */
void ElasticHydrogelSuoT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);
	/*material parameters for matrix*/
	sub_list.AddSub("gel_potential", ParameterListT::Once);	
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ElasticHydrogelSuoT::NewSub(const StringT& name) const
{
	PotentialT* pot = NULL;
	if (name == "neo-hookean")
		pot = new NeoHookean;
	else if (name == "mooney-rivlin")
		pot = new MooneyRivlin;
	else if (name == "arruda-boyce")
		pot = new ArrudaBoyce;
	if (pot)
		return pot;
		
	if (name == "gel_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);
		
		/* choice of parameters */
		choice->AddSub("neo-hookean");
		choice->AddSub("mooney-rivlin");
		choice->AddSub("arruda-boyce");
		return(choice);
	}
}

void ElasticHydrogelSuoT::DefineParameters(ParameterListT& list) const
{
	SolidMaterialT::DefineParameters(list);
	/*list.AddParameter(ParameterT::Double, "density"); this is added*/
	/*list.AddParameter(ParameterT::Double, "reference_temperature");*/
	/*list.AddParameter(ParameterT::Double, "FH_parameter");*/
    list.AddParameter(ParameterT::Double, "FH_parameter1");
    list.AddParameter(ParameterT::Double, "FH_parameter2");
    list.AddParameter(ParameterT::Double, "FH_parameter3");
    list.AddParameter(ParameterT::Double, "FH_parameter4");
    list.AddParameter(ParameterT::Double, "initial_solidfraction");
    
}

/* accept parameter list */
void ElasticHydrogelSuoT::TakeParameterList(const ParameterListT& list)
{
	FSSolidMatT::TakeParameterList(list);
	StringT caller = "ElasticHydrogelSuoT::TakeParameterList";
	
	const ParameterListT& gel_pot = list.GetListChoice(*this, "gel_potential");
	if(gel_pot.Name() == "neo-hookean")
		fPot = new NeoHookean;
	else if(gel_pot.Name() == "mooney-rivlin")
		fPot = new MooneyRivlin;
	else if(gel_pot.Name() == "arruda-boyce")
	     fPot = new ArrudaBoyce;
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	if (!fPot) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", gel_pot.Name().Pointer());			

	fPot->TakeParameterList(gel_pot);
	
	/*fT = list.GetParameter("reference_temperature");
	fchi = list.GetParameter("FH_parameter");*/
    xL=list.GetParameter("FH_parameter1");
    xH=list.GetParameter("FH_parameter2");
    Ttran=list.GetParameter("FH_parameter3");
    DeltaT=list.GetParameter("FH_parameter4");
    solidfraction0=list.GetParameter("initial_solidfraction");
	
	/*dimension workspace*/

	fF3D.Dimension(3);
	fb.Dimension(3);
	fbe.Dimension(3);
	fEigs_dev.Dimension(3);
	fEigs.Dimension(3);
	fEigs_e.Dimension(3);

	ftau.Dimension(3);
	fStress.Dimension(NumSD());
	fStress3D.Dimension(3);

	fDtauDe.Dimension(3);
	fModulus3D.Dimension(6);
	fModMat.Dimension(6);
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));	
	
	/*dimension state variables*/
	fnstatev = 0;
	fnstatev ++; /*fSolidFraction*/
	fnstatev ++; /*fSolidFraction_n*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	fSolidFraction = pstatev;
	pstatev ++;
	fSolidFraction_n = pstatev;
}


/*initializes history variable */
void  ElasticHydrogelSuoT::PointInitialize(void)
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
		      
			  *fSolidFraction = solidfraction0;
			  *fSolidFraction_n = solidfraction0;

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void ElasticHydrogelSuoT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP();ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		*fSolidFraction_n = *fSolidFraction;

		/* write to storage */
		Store(element, ip);
	}
}

void ElasticHydrogelSuoT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		*fSolidFraction = *fSolidFraction_n;
		
		/* write to storage */
		Store(element, ip);
	}
}

void ElasticHydrogelSuoT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}
void ElasticHydrogelSuoT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;

}

/*************************************************************************
 * Private
 *************************************************************************/
 
/*void ElasticHydrogelSuoT::CalculatePhi(double& phi,  const double& phi_n, const double& J)*/
void ElasticHydrogelSuoT::CalculatePhi(double& phi,  const double& phi_n, const double& J)
{ 

	const double ctol=1.00e-12;
	double maxiteration=50;

	/*double chempo=0;*/
    double fT=Compute_Temperature();
    double fchi=(xL+xH)/2.0+(xH-xL)/2.0*(1.0-exp(-2.0*(fT-Ttran)/DeltaT))/(1.0+exp(-2.0*(fT-Ttran)/DeltaT));
	phi = phi_n;
	
	//This is something I added to give an initialize J, which is quite not goood//
	

	//see the following right or wrong
	//bulid says there is no GetParamter funciton in PontentialT.cpp
	double fkappa=fPot->GetKappa();
    double alpha=pow(solidfraction0,-1.0*third);
	double res=R*fT*(log(1-phi)+phi+fchi*phi*phi)-fkappa*upsilon*phi/2*(pow(alpha,6.0)*phi*J*phi*J-1);
	double error=sqrt(res*res);
  /*  cout <<"\n temperature is:  "<<fT;
	cout <<"\n flory huggins parameter is:  "<<fchi; 
	cout <<"\n J is:  "<<J;
	cout <<"\n residule is: "<<res; */
	int iteration=0;
	while (error>ctol && iteration <  maxiteration)
	{ 
		iteration ++;
		
		double K=R*fT*(1/(phi-1)+1+2*fchi*phi)-fkappa*upsilon/2*(3*phi*phi*J*J*pow(alpha,6.0)-1);
		double dphi = -res/K;
		phi = phi + dphi;
		
		res=R*fT*(log(1-phi)+phi+fchi*phi*phi)-fkappa*upsilon*phi/2*(pow(alpha,6.0)*phi*J*phi*J-1);
		error=sqrt(res*res);
	}
/*    if(CurrElementNumber()==0&&CurrIP()==0)
    {cout <<"\n polymer concentration is:  "<<phi;
        cout <<"\n residule is: "<<res;} */
	if (iteration >= maxiteration) 
	{
		ExceptionT::GeneralFail("ElasticHydrogelSuoT::ComputePhi", 
			"number of iteration exceeds maximum");
	}
}

/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void ElasticHydrogelSuoT::MixedRank4_2D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 2 ||
	    b.Length() != 2 ||
	    rank4_ab.Rows() != 3 ||
	    rank4_ab.Cols() != 3) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = B(1.);
//	z4 = B(2.);

	z1 = a[0];
	z2 = a[1];
	z3 = b[0];
	z4 = b[1];

	z5 = z1*z1;
	z6 = z2*z2;
	z7 = z3*z3;
	z8 = 2.*z1*z2*z3*z4;
	z9 = z4*z4;
	z3 = 2.*z3*z4;
	z4 = z3*z5;
	z3 = z3*z6;
	z10 = 2.*z1*z2*z7;
	z11 = 2.*z5*z7;
	z7 = z6*z7;
	z1 = 2.*z1*z2*z9;
	z2 = z5*z9;
	z5 = 2.*z6*z9;
	z4 = z10 + z4;
	z1 = z1 + z3;
	z2 = z2 + z7 + z8;
	z3 = 0.5*z4;
	z1 = 0.5*z1;
	z2 = 0.5*z2;

	//{{z11, z8, z3}, 
	// {z8, z5, z1}, 
	// {z3, z1, z2}}

	double* p = rank4_ab.Pointer();
	*p++ = z11;
    *p++ = z8;
    *p++ = z3;
    *p++ = z8;
    *p++ = z5;
    *p++ = z1;
    *p++ = z3;
    *p++ = z1;
    *p   = z2;
}

void ElasticHydrogelSuoT::MixedRank4_3D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 3 ||
	    b.Length() != 3 ||
	    rank4_ab.Rows() != 6 ||
	    rank4_ab.Cols() != 6) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = A(3.);
//	z4 = B(1.);
//	z5 = B(2.);
//	z6 = B(3.);

	z1 = a[0];	
	z2 = a[1];
	z3 = a[2];
	z4 = b[0];
	z5 = b[1];
	z6 = b[2];
	z7 = z1*z1;
	z8 = z2*z2;
	z9 = z3*z3;
	z10 = 2.*z1*z4;
	z11 = z10*z2;
	z12 = z2*z4;
	z13 = z1*z12;
	z14 = z4*z4;
	z15 = z10*z5;
	z15 = z15*z3;
	z16 = z11*z5;
	z17 = z12*z5;
	z18 = 2.*z17;
	z17 = z17*z3;
	z18 = z18*z3;
	z19 = z1*z4*z5;
	z19 = z19*z3;
	z20 = z5*z5;
	z11 = z11*z6;
	z13 = z13*z6;
	z10 = z10*z3*z6;
	z12 = z12*z3*z6;
	z21 = 2.*z12;
	z22 = z1*z2*z5*z6;
	z23 = 2.*z22;
	z24 = z1*z3*z5*z6;
	z25 = 2.*z24;
	z26 = 2.*z2*z3*z5*z6;
	z27 = z6*z6;
	z28 = 2.*z1*z14;
	z29 = 2.*z1*z2;
	z30 = z14*z2;
	z11 = z11 + z15;
	z15 = z1*z20*z3;
	z31 = 2.*z2*z20*z3;
	z18 = z18 + z23;
	z21 = z21 + z25;
	z23 = z1*z2*z27;
	z1 = 2.*z1*z27*z3;
	z25 = 2.*z2*z27*z3;
	z2 = z2*z28;
	z28 = z28*z3;
	z29 = z20*z29;
	z3 = z3*z30;
	z30 = 2.*z14*z7;
	z32 = z20*z7;
	z33 = z27*z7;
	z34 = 2.*z4*z5*z7;
	z35 = 2.*z4*z6*z7;
	z7 = z5*z6*z7;
	z36 = z14*z8;
	z37 = 2.*z20*z8;
	z38 = z27*z8;
	z39 = 2.*z4*z5*z8;
	z40 = z4*z6*z8;
	z8 = 2.*z5*z6*z8;
	z14 = z14*z9;
	z20 = z20*z9;
	z27 = 2.*z27*z9;
	z41 = z4*z5*z9;
	z4 = 2.*z4*z6*z9;
	z5 = 2.*z5*z6*z9;
	z6 = 0.5*z11;
	z9 = 0.5*z18;
	z11 = 0.5*z21;
	z2 = z2 + z34;
	z18 = z28 + z35;
	z3 = z13 + z19 + z3 + z7;
	z7 = z16 + z32 + z36;
	z13 = z29 + z39;
	z15 = z15 + z17 + z22 + z40;
	z8 = z31 + z8;
	z14 = z10 + z14 + z33;
	z17 = z20 + z26 + z38;
	z12 = z12 + z23 + z24 + z41;
	z1 = z1 + z4;
	z4 = z25 + z5;
	z2 = 0.5*z2;
	z5 = 0.5*z18;
	z3 = 0.5*z3;
	z7 = 0.5*z7;
	z13 = 0.5*z13;
	z15 = 0.5*z15;
	z8 = 0.5*z8;
	z14 = 0.5*z14;
	z17 = 0.5*z17;
	z12 = 0.5*z12;
	z1 = 0.5*z1;
	z4 = 0.5*z4;
	
	//{{z30, z16, z10,  z6,  z5,  z2}, 
	// {z16, z37, z26,  z8,  z9, z13}, 
	// {z10, z26, z27,  z4,  z1, z11}, 
	// { z6,  z8,  z4, z17, z12, z15}, 
	// { z5,  z9,  z1, z12, z14,  z3},
	// { z2, z13, z11, z15,  z3,  z7}}
	
	double* p = rank4_ab.Pointer();
    *p++ = z30;
    *p++ = z16;
    *p++ = z10;
    *p++ = z6;
    *p++ = z5;
    *p++ = z2;
    *p++ = z16;
    *p++ = z37;
    *p++ = z26;
    *p++ = z8;
    *p++ = z9;
    *p++ = z13;
    *p++ = z10;
    *p++ = z26;
    *p++ = z27;
    *p++ = z4;
    *p++ = z1;
    *p++ = z11;
    *p++ = z6;
    *p++ = z8;
    *p++ = z4;
    *p++ = z17;
    *p++ = z12;
    *p++ = z15;
    *p++ = z5;
    *p++ = z9;
    *p++ = z1;
    *p++ = z12;
    *p++ = z14;
    *p++ = z3;
    *p++ = z2;
    *p++ = z13;
    *p++ = z11;
    *p++ = z15;
    *p++ = z3;
    *p  = z7;
}

