/* $Id: IsoVECorneaModel.cpp,v 1.3 2005-04-20 23:45:17 thao Exp $ */
/* created: paklein (11/08/1997) */
#include "IsoVECorneaModel.h"

#include <math.h>
#include <iostream.h>
#include "ParameterContainerT.h"
#include "toolboxConstants.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ifstreamT.h"

/*potential*/
#include "WormLikeChain.h"
#include "WLCwRep.h"
#include "FungType.h"
#include "FungwRep.h"

/* point generators */
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"
#include "FCCPtsT.h"

using namespace Tahoe;

const int kNumOutputVar = 1; 
static const char* Labels[kNumOutputVar] = {"Dvisc"}; 
static const double third = 1.0/3.0;
static const double pi = 3.14159274641021;

/* constructors */
IsoVECorneaModel::IsoVECorneaModel(void):
	RGViscoelasticityT(),
	ParameterInterfaceT("Isotropic_Viscoelastic_Cornea_Model"),
	fSphere(NULL),
	fPot_EQ(NULL),
	fPot_NEQ(NULL),
	fNumSD(3),
	fNumStress(3),
	fNumModuli(6),
	fModulus(6),
	fStress(3),
	fSpectralDecompSpat(3),
	fb(3),
	fbe(3),
	fiC(3),
	fEigs(3),
	fEigs_e(3),
	fdWdE_EQ(3),
	fdWdE_NEQ(3),
	fddWddE_EQ(3),
	fddWddE_NEQ(3),
	fCalg(3),
	fModMat(6),
	fiKAB(3)	
{
}

/* destructor */
IsoVECorneaModel::~IsoVECorneaModel(void) 
{ 
  delete fSphere; 
  delete fPot_NEQ;
  delete fPot_EQ;
}

/* information about subordinate parameter lists */
void IsoVECorneaModel::DefineSubs(SubListT& sub_list) const
{
  /* inherited */
  RGViscoelasticityT::DefineSubs(sub_list);

  /* choice of integration schemes */
  sub_list.AddSub("sphere_integration_choice", ParameterListT::Once, true);
  sub_list.AddSub("equilibrium_potential_choice", ParameterListT::Once, true);
  sub_list.AddSub("nonequilibrium_potential_choice", ParameterListT::Once, true);
  sub_list.AddSub("viscosity_choice", ParameterListT::Once, true);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* IsoVECorneaModel::NewSub(const StringT& name) const
{
  /* inherited */
  ParameterInterfaceT* sub = RGViscoelasticityT::NewSub(name);
  if (sub) 
  {
    return sub;
  }
  else if (name == "sphere_integration_choice")
  {
    ParameterContainerT* choice = new ParameterContainerT(name);
    choice->SetListOrder(ParameterListT::Choice);

    /* bound */
    LimitT lower(1, LimitT::LowerInclusive);
  
    /* like the grid on a sphere */
    ParameterContainerT lat_long("latitude_longitude");
    ParameterT n_phi(ParameterT::Integer, "n_phi");
    n_phi.AddLimit(lower);
    ParameterT n_theta(ParameterT::Integer, "n_theta");
    n_theta.AddLimit(lower);
    lat_long.AddParameter(n_phi);
    lat_long.AddParameter(n_theta);
    choice->AddSub(lat_long);
    
    /* icosahedral points */
    ParameterContainerT ico("icosahedral");
    ParameterT ico_points(ParameterT::Integer, "points");
    ico_points.AddLimit(6, LimitT::Only);
    ico_points.AddLimit(10, LimitT::Only);
    ico_points.AddLimit(40, LimitT::Only);
    ico_points.AddLimit(160, LimitT::Only);
    ico.AddParameter(ico_points);
    choice->AddSub(ico);
    
    /* FCC point arrangement */
    ParameterContainerT fcc("fcc_points");
    ParameterT n_shells(ParameterT::Integer, "shells");
    n_shells.AddLimit(lower);
    fcc.AddParameter(n_shells);
    fcc.AddParameter(ParameterT::Double, "nearest_neighbor_distance");
    choice->AddSub(fcc);
    return choice;
  }
  else if (name == "equilibrium_potential_choice")
  {
    ParameterContainerT* choice = new ParameterContainerT(name);
     choice->SetListOrder(ParameterListT::Choice);

   /* bound */
    LimitT lower(0.0, LimitT::Lower);

    /* worm like chain statistics */
    ParameterContainerT wlc("worm_like_chain");
    ParameterT N(ParameterT::Double, "chain_density");
    N.AddLimit(lower);
    ParameterT T(ParameterT::Double, "temperature");
    T.AddLimit(lower);
    ParameterT A(ParameterT::Double, "persistence_length");
    A.AddLimit(lower);
    ParameterT R0(ParameterT::Double, "initial_coil_length");
    R0.AddLimit(lower);
    wlc.AddParameter(N);
    wlc.AddParameter(T);
    wlc.AddParameter(A);
    wlc.AddParameter(R0);
    choice->AddSub(wlc);
 
    /* worm like chain statistics with power law repulsion*/
    ParameterContainerT wlcrep("worm_like_chain_power_repulsion");
    ParameterT N2(ParameterT::Double, "chain_density");
    N2.AddLimit(lower);
    ParameterT T2(ParameterT::Double, "temperature");
    T2.AddLimit(lower);
    ParameterT A2(ParameterT::Double, "persistence_length");
    A2.AddLimit(lower);
    ParameterT R02(ParameterT::Double, "initial_coil_length");
    R02.AddLimit(lower);
    ParameterT n(ParameterT::Double, "repulsion_exponent");
    n.AddLimit(lower);
    wlcrep.AddParameter(N2);
    wlcrep.AddParameter(T2);
    wlcrep.AddParameter(A2);
    wlcrep.AddParameter(R02);
    wlcrep.AddParameter(n);
    choice->AddSub(wlcrep);

    /* fung-type chains */
    ParameterContainerT fungtype("fung_type_chain");
    ParameterT C1(ParameterT::Double, "multiplier_C1");
    C1.AddLimit(lower);
    ParameterT beta(ParameterT::Double, "exponent_beta");
    beta.AddLimit(lower);
    fungtype.AddParameter(C1);
    fungtype.AddParameter(beta);
    choice->AddSub(fungtype);

    /* fung-type chains with repulsion */
    ParameterContainerT fungrep("fung_type_power_repulsion");
    ParameterT C12(ParameterT::Double, "multiplier_C1");
    C12.AddLimit(lower);
    ParameterT beta2(ParameterT::Double, "exponent_beta");
    beta2.AddLimit(lower);
    ParameterT n2(ParameterT::Double, "repulsion_exponent");
    n2.AddLimit(lower);
    fungrep.AddParameter(C12);
    fungrep.AddParameter(beta2);
    fungrep.AddParameter(n2);
    choice->AddSub(fungrep);
    return choice;
  }
  else if (name == "nonequilibrium_potential_choice")
  {
    ParameterContainerT* choice = new ParameterContainerT(name);
     choice->SetListOrder(ParameterListT::Choice);

   /* bound */
    LimitT lower(0.0, LimitT::Lower);

    /* worm like chain statistics */
    ParameterContainerT wlc("worm_like_chain");
    ParameterT N(ParameterT::Double, "chain_density");
    N.AddLimit(lower);
    ParameterT T(ParameterT::Double, "temperature");
    T.AddLimit(lower);
    ParameterT A(ParameterT::Double, "persistence_length");
    A.AddLimit(lower);
    ParameterT R0(ParameterT::Double, "initial_coil_length");
    R0.AddLimit(lower);
    wlc.AddParameter(N);
    wlc.AddParameter(T);
    wlc.AddParameter(A);
    wlc.AddParameter(R0);
    choice->AddSub(wlc);
 
    /* worm like chain statistics with power law repulsion*/
    ParameterContainerT wlcrep("worm_like_chain_power_repulsion");
    ParameterT N2(ParameterT::Double, "chain_density");
    N2.AddLimit(lower);
    ParameterT T2(ParameterT::Double, "temperature");
    T2.AddLimit(lower);
    ParameterT A2(ParameterT::Double, "persistence_length");
    A2.AddLimit(lower);
    ParameterT R02(ParameterT::Double, "initial_coil_length");
    R02.AddLimit(lower);
    ParameterT n(ParameterT::Double, "repulsion_exponent");
    n.AddLimit(lower);
    wlcrep.AddParameter(N2);
    wlcrep.AddParameter(T2);
    wlcrep.AddParameter(A2);
    wlcrep.AddParameter(R02);
    wlcrep.AddParameter(n);
    choice->AddSub(wlcrep);

    /* fung-type chains */
    ParameterContainerT fungtype("fung_type_chain");
    ParameterT C1(ParameterT::Double, "multiplier_C1");
    C1.AddLimit(lower);
    ParameterT beta(ParameterT::Double, "exponent_beta");
    beta.AddLimit(lower);
    fungtype.AddParameter(C1);
    fungtype.AddParameter(beta);
    choice->AddSub(fungtype);

    /* fung-type chains with repulsion */
    ParameterContainerT fungrep("fung_type_power_repulsion");
    ParameterT C12(ParameterT::Double, "multiplier_C1");
    C12.AddLimit(lower);
    ParameterT beta2(ParameterT::Double, "exponent_beta");
    beta2.AddLimit(lower);
    ParameterT n2(ParameterT::Double, "repulsion_exponent");
    n2.AddLimit(lower);
    fungrep.AddParameter(C12);
    fungrep.AddParameter(beta2);
    fungrep.AddParameter(n2);
    choice->AddSub(fungrep);
    return choice;
  }
  else if (name == "viscosity_choice")
  {
    ParameterContainerT* choice = new ParameterContainerT(name);
     choice->SetListOrder(ParameterListT::Choice);

    /* bound */
    LimitT lower(0.0, LimitT::Lower);

    /* worm like chain statistics */
    ParameterContainerT linearvisc("linear_viscosity");
    ParameterT tauS(ParameterT::Double, "initial_shear_relaxation_time");
    tauS.AddLimit(lower);
    ParameterT tauB(ParameterT::Double, "initial_bulk_relaxation_time");
    tauB.AddLimit(lower);
    linearvisc.AddParameter(tauS);
    linearvisc.AddParameter(tauB);
    choice->AddSub(linearvisc);
    return choice;
  }
}

/* accept parameter list */
void IsoVECorneaModel::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
  RGViscoelasticityT::TakeParameterList(list);
  const ParameterListT& eqpotential = list.GetListChoice(*this, "equilibrium_potential_choice");
  const ParameterListT& neqpotential = list.GetListChoice(*this, "nonequilibrium_potential_choice");
  const ParameterListT& points = list.GetListChoice(*this, "sphere_integration_choice");
  const ParameterListT& viscosities = list.GetListChoice(*this, "viscosity_choice");
  
  if (eqpotential.Name() == "worm_like_chain")
  {
      double N = eqpotential.GetParameter("chain_density");
      double T = eqpotential.GetParameter("temperature");
      double A = eqpotential.GetParameter("persistence_length");
      fR0_EQ = eqpotential.GetParameter("initial_coil_length");
      double L = fR0_EQ*fR0_EQ/(2*A);
      double k = 1.3806503e-23; 
		
      /*set wlc potential function*/
      double K = 0.75/pi*(N*k*T/A);
      fPot_EQ = new WormLikeChain(K,L);
      if (!fPot_EQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      double p = sqrt(2*A/L);
      fC0_EQ = (N*k*T/A) * ( p - 0.25 + 0.25/((1-p)*(1-p)) );
  }
  else if (eqpotential.Name() == "worm_like_chain_power_repulsion")
  {
      double N = eqpotential.GetParameter("chain_density");
      double T = eqpotential.GetParameter("temperature");
      double A = eqpotential.GetParameter("persistence_length");
      double R0 = eqpotential.GetParameter("initial_coil_length");
      double n = eqpotential.GetParameter("repulsion_exponent");
      double L = R0*R0/(2*A);
      double k = 1.3806503e-23; 
		
      /*set wlc potential function*/
      double K = 0.75/pi*(N*k*T/A);
      double p = sqrt(2*A/L);
      double C = 0.75/pi*(N*k*T/A) * ( p - 0.25 + 0.25/((1-p)*(1-p)) )*pow(R0,(n+1.0))/n;
      fPot_EQ = new WLCwRep(K,L,C,n);
      if (!fPot_EQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0_EQ = 0.0;
      fR0_EQ = R0;
  }
  else if (eqpotential.Name() == "fung_type_chain")
    {
      double C1 = eqpotential.GetParameter("multiplier_C1");
      double beta = eqpotential.GetParameter("exponent_beta");

      /*set fung type potential*/
      double K = 0.75/pi*C1;
      fPot_EQ = new FungType(K,beta);
      if (!fPot_EQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0_EQ = C1*beta;
      fR0_EQ = 1.0;      
    }
  else if (eqpotential.Name() == "fung_type_power_repulsion")
    {
      double C1 = eqpotential.GetParameter("multiplier_C1");
      double beta = eqpotential.GetParameter("exponent_beta");
      double n = eqpotential.GetParameter("repulsion_exponent");

      /*set fung type potential*/
      double K = 0.75/pi*C1;
      double C = 0.75/pi*C1*beta/n;
	  fPot_EQ = new FungwRep(K,beta,C,n);
      if (!fPot_EQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0_EQ = 0.0;
      fR0_EQ = 1.0;      
    }
  else
    ExceptionT::GeneralFail("IsoVECorneaModel::TakeParameterList", "unrecognized potential \"%s\"", eqpotential.Name().Pointer());

  if (neqpotential.Name() == "worm_like_chain")
  {
      double N = neqpotential.GetParameter("chain_density");
      double T = neqpotential.GetParameter("temperature");
      double A = neqpotential.GetParameter("persistence_length");
      fR0_NEQ = neqpotential.GetParameter("initial_coil_length");
      double L = fR0_NEQ*fR0_NEQ/(2*A);
      double k = 1.3806503e-23; 
		
      /*set wlc potential function*/
      double K = 0.75/pi*(N*k*T/A);
      fPot_NEQ = new WormLikeChain(K,L);
      if (!fPot_NEQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      double p = sqrt(2*A/L);
      fC0_NEQ = (N*k*T/A) * ( p - 0.25 + 0.25/((1-p)*(1-p)) );
  }
  else if (neqpotential.Name() == "worm_like_chain_power_repulsion")
  {
      double N = neqpotential.GetParameter("chain_density");
      double T = neqpotential.GetParameter("temperature");
      double A = neqpotential.GetParameter("persistence_length");
      double R0 = neqpotential.GetParameter("initial_coil_length");
      double n = neqpotential.GetParameter("repulsion_exponent");
      double L = R0*R0/(2*A);
      double k = 1.3806503e-23; 
		
      /*set wlc potential function*/
      double K = 0.75/pi*(N*k*T/A);
      double p = sqrt(2*A/L);
      double C = 0.75/pi*(N*k*T/A) * ( p - 0.25 + 0.25/((1-p)*(1-p)) )*pow(R0,(n+1.0))/n;
      fPot_NEQ = new WLCwRep(K,L,C,n);
      if (!fPot_NEQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0_NEQ = 0.0;
      fR0_NEQ = R0;
  }
  else if (neqpotential.Name() == "fung_type_chain")
  {
      double C1 = neqpotential.GetParameter("multiplier_C1");
      double beta = neqpotential.GetParameter("exponent_beta");

      /*set fung type potential*/
      double K = 0.75/pi*C1;
      fPot_NEQ = new FungType(K,beta);
      if (!fPot_NEQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0_NEQ = C1*beta;
      fR0_NEQ = 1.0;      
  }
  else if (neqpotential.Name() == "fung_type_power_repulsion")
    {
      double C1 = neqpotential.GetParameter("multiplier_C1");
      double beta = neqpotential.GetParameter("exponent_beta");
      double n = neqpotential.GetParameter("repulsion_exponent");

      /*set fung type potential*/
      double K = 0.75/pi*C1;
      double C = 0.75/pi*C1*beta/n;
	  fPot_NEQ = new FungwRep(K,beta,C,n);
      if (!fPot_NEQ) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0_NEQ = 0.0;
      fR0_NEQ = 1.0;      
    }
  else
    ExceptionT::GeneralFail("IsoVECorneaModel::TakeParameterList", "unrecognized potential \"%s\"", neqpotential.Name().Pointer());

  if (points.Name() == "latitude_longitude")
  {
      int n_phi = points.GetParameter("n_phi");
      int n_theta = points.GetParameter("n_theta");
      fSphere = new LatLongPtsT(n_phi, n_theta);
  }
  else if (points.Name() == "icosahedral")
  {
      int np = points.GetParameter("points");
      fSphere = new IcosahedralPtsT(np);
  }
  else if (points.Name() == "fcc_points")
  {
      int num_shells = points.GetParameter("shells");
      double bond_length = points.GetParameter("nearest_neighbor_distanc\
e");
      fSphere = new FCCPtsT(num_shells, bond_length);
  }
  else
    ExceptionT::GeneralFail("IsoVECorneaModel::TakeParameterList", "unrecognized point scheme \"%s\"", points.Name().Pointer());

   /* set tables */
  Construct();
  
  if (viscosities.Name() == "linear_viscosity")
  {
      double tauS = viscosities.GetParameter("initial_shear_relaxation_time");
      double tauB = viscosities.GetParameter("initial_bulk_relaxation_time");
	/* initial nonequilibrium moduli */
	fEigs_e = 1.0;
	ddWddE(fEigs_e, fdWdE_NEQ, fddWddE_NEQ, fPot_NEQ, fR0_NEQ, fC0_NEQ);
	double lambda = fddWddE_NEQ(0,1);
	double mu = 0.5*(fddWddE_NEQ(0,0) - fddWddE_NEQ(0,1));
	double kappa = lambda + 2.0/3.0*mu;
	
	fietaS = 1.0/(mu*tauS);
	fietaB = 1.0/(kappa*tauB);
  }
  else
    ExceptionT::GeneralFail("IsoVECorneaModel::TakeParameterList", "unrecognized potential \"%s\"", viscosities.Name().Pointer());

}

/* strain energy density */
double IsoVECorneaModel::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	Compute_b(fb);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
	fEigs = fSpectralDecompSpat.Eigenvalues();
	
	/* stretched bonds */
	ComputeLengths(fEigs, fR0_EQ);
	/* update potential table */
	fPot_EQ->MapFunction(fLengths,fU);
	/* sum contributions to equilibrium entropic component*/
	double  energy = 0.0;
	double* pU = fU.Pointer();
	double* pj = fjacobian.Pointer();
	double iR0 = 1.0/fR0_EQ;	
	for (int j = 0; j < fLengths.Length(); j++)
		energy += (*pU++)*(*pj++)*iR0;
	/*add contribution from repulsion component*/
	double I3 = fEigs.Product();
	energy -= 0.5*fC0_EQ*log(I3);
	
	/*adds nonequilibrium part */
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
	const dMatrixT& F = F_mechanical();	
	fiC = fC_v.Inverse();
	fbe.MultQBQT(F,fiC);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues();
	I3 = fEigs_e.Product();
	/* stretched bonds */
	ComputeLengths(fEigs_e, fR0_NEQ);
	/* update potential table */
	fPot_NEQ->MapFunction(fLengths,fU);
	/* sum contributions to equilibrium entropic component*/
	pU = fU.Pointer();
	pj = fjacobian.Pointer();
	iR0 = 1.0/fR0_NEQ;	
	for (int j = 0; j < fLengths.Length(); j++)
		energy += (*pU++)*(*pj++)*iR0;
	/*add contribution from repulsion component*/
	energy -= 0.5*fC0_NEQ*log(I3);
	
	return(energy);
}

/* modulus */
const dMatrixT& IsoVECorneaModel::c_ijkl(void)
{
    
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
    
	const dMatrixT& F = F_mechanical();
    Compute_b(fb);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
    const ArrayT<dArrayT>& eigenvectors=fSpectralDecompSpat.Eigenvectors();

	/*add equilibrium component*/
	ddWddE(fEigs, fdWdE_EQ, fddWddE_EQ, fPot_EQ, fR0_EQ, fC0_EQ); 
    fddWddE_EQ(0,0) -= 2.0*fdWdE_EQ[0];
    fddWddE_EQ(1,1) -= 2.0*fdWdE_EQ[1];
    fddWddE_EQ(2,2) -= 2.0*fdWdE_EQ[2];
   
	fModulus = fSpectralDecompSpat.EigsToRank4(fddWddE_EQ);	
  
	double dl, coeff;

    double& l0 = fEigs[0];
    double& l1 = fEigs[1];
    double& l2 = fEigs[2];
	
	dl = l0 - l1;
    if (fabs(dl) > kSmall)
		coeff = (fdWdE_EQ[0]*l1 - fdWdE_EQ[1]*l0)/dl;
    else 
		coeff = (fddWddE_EQ(0,0)-fddWddE_EQ(0,1));
    MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
    fModulus.AddScaled(2.0*coeff, fModMat);
    
    dl = l0 - l2;
    if (fabs(dl) > kSmall)
      coeff = (fdWdE_EQ[0]*l2 - fdWdE_EQ[2]*l0)/dl;
    else 
      coeff = (fddWddE_EQ(0,0)-fddWddE_EQ(0,2));	
    MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
    fModulus.AddScaled(2.0*coeff, fModMat);
    
    dl = l1 - l2;
   if (fabs(dl) > kSmall)
		coeff  = (fdWdE_EQ[1]*l2 - fdWdE_EQ[2]*l1)/dl;
    else
      coeff = (fddWddE_EQ(1,1)-fddWddE_EQ(1,2));	
    MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
    fModulus.AddScaled(2.0*coeff, fModMat);
		
	/*calculates trial stretches*/
	fiC = fC_vn;
	fbe.MultQBQT(F, fiC.Inverse());
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues();  /*trial elastic stretch*/

	double l0_tr = fEigs_e[0];
	double l1_tr = fEigs_e[1];
	double l2_tr = fEigs_e[2];
	const ArrayT<dArrayT>& eigenvectors_e=fSpectralDecompSpat.Eigenvectors();

	/*calculates elastic stretches*/
	fiC = fC_v;
	fbe.MultQBQT(F, fiC.Inverse());
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
	ddWddE(fEigs_e, fdWdE_NEQ, fddWddE_NEQ, fPot_NEQ, fR0_NEQ, fC0_NEQ); 

	ComputeiKAB(fddWddE_NEQ, fietaS, fietaB);
	
	fCalg(0,0) = fddWddE_NEQ(0,0)*fiKAB(0,0) + fddWddE_NEQ(0,1)*fiKAB(1,0) + fddWddE_NEQ(0,2)*fiKAB(2,0) 
				- 2.0*fdWdE_NEQ[0];
	fCalg(1,0) = fddWddE_NEQ(1,0)*fiKAB(0,0) + fddWddE_NEQ(1,1)*fiKAB(1,0) + fddWddE_NEQ(1,2)*fiKAB(2,0);
	fCalg(2,0) = fddWddE_NEQ(2,0)*fiKAB(0,0) + fddWddE_NEQ(2,1)*fiKAB(1,0) + fddWddE_NEQ(2,2)*fiKAB(2,0);
	fCalg(0,1) = fddWddE_NEQ(0,0)*fiKAB(0,1) + fddWddE_NEQ(0,1)*fiKAB(1,1) + fddWddE_NEQ(0,2)*fiKAB(2,1);
	fCalg(1,1) = fddWddE_NEQ(1,0)*fiKAB(0,1) + fddWddE_NEQ(1,1)*fiKAB(1,1) + fddWddE_NEQ(1,2)*fiKAB(2,1) 
				- 2.0*fdWdE_NEQ[1];
	fCalg(2,1) = fddWddE_NEQ(2,0)*fiKAB(0,1) + fddWddE_NEQ(2,1)*fiKAB(1,1) + fddWddE_NEQ(2,2)*fiKAB(2,1);
	fCalg(0,2) = fddWddE_NEQ(0,0)*fiKAB(0,2) + fddWddE_NEQ(0,1)*fiKAB(1,2) + fddWddE_NEQ(0,2)*fiKAB(2,2);
	fCalg(1,2) = fddWddE_NEQ(1,0)*fiKAB(0,2) + fddWddE_NEQ(1,1)*fiKAB(1,2) + fddWddE_NEQ(1,2)*fiKAB(2,2);
	fCalg(2,2) = fddWddE_NEQ(2,0)*fiKAB(0,2) + fddWddE_NEQ(2,1)*fiKAB(1,2) + fddWddE_NEQ(2,2)*fiKAB(2,2) 
				- 2.0*fdWdE_NEQ[2];
					
	fModulus += fSpectralDecompSpat.NonSymEigsToRank4(fCalg);
    
	double dl_tr;
	
	dl_tr = l0_tr - l1_tr;
	if (fabs(dl_tr) > kSmall)
		coeff = (fdWdE_NEQ[0]*l1_tr - fdWdE_NEQ[1]*l0_tr)/dl_tr;
	else 
		coeff = (fCalg(0,0)-fCalg(0,1));
	MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[1], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);
    
	dl_tr = l0_tr - l2_tr;
	if (fabs(dl_tr) > kSmall)
		coeff =(fdWdE_NEQ[0]*l2_tr - fdWdE_NEQ[2]*l0_tr)/dl_tr;
	else 
		coeff = (fCalg(0,0)-fCalg(0,2));	
	MixedRank4_3D(eigenvectors_e[0], eigenvectors_e[2], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);

	dl_tr = l1_tr - l2_tr;
	if (fabs(dl_tr) > kSmall)
		coeff  = (fdWdE_NEQ[1]*l2_tr - fdWdE_NEQ[2]*l1_tr)/dl_tr;
	else
		coeff = (fCalg(1,1)-fCalg(1,2));	
	MixedRank4_3D(eigenvectors_e[1], eigenvectors_e[2], fModMat);
	fModulus.AddScaled(2.0*coeff, fModMat);
	
	fModulus *= 1.0/sqrt(fEigs.Product());
//	cout<<"\nmodulus: "<< fModulus;
    return fModulus;
}

/* stresses */
const dSymMatrixT& IsoVECorneaModel::s_ij(void)
{
    /* stretch tensor */
    Compute_b(fb);
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();
	/*sum equilibrium components*/
	dWdE(fEigs, fdWdE_EQ, fPot_EQ, fR0_EQ, fC0_EQ); 
	fStress = fSpectralDecompSpat.EigsToRank2(fdWdE_EQ);

	/*load inelastic stretch*/
	const dMatrixT& F = F_mechanical();
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
	{		
		/*calc trial state*/
		fiC = fC_vn;
		fbe.MultQBQT(F, fiC.Inverse());
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues();  /*trial strains*/
		ComputeEigs_e(fEigs, fEigs_e, fdWdE_NEQ, fddWddE_NEQ, 
			fPot_NEQ, fR0_NEQ, fC0_NEQ, fietaS, fietaB);
		/*sum nonequilibrium component*/
		dWdE(fEigs_e, fdWdE_NEQ, fPot_NEQ, fR0_NEQ, fC0_NEQ); 
		fStress += fSpectralDecompSpat.EigsToRank2(fdWdE_NEQ);

		/*Calculate Cv*/
		fbe = fSpectralDecompSpat.EigsToRank2(fEigs_e); /*be which is colinear with btr*/
		fC_v.MultQTBQ(F, fbe.Inverse()); 
		Store(element, CurrIP());
	}	
	else 
	{
		/*calc elastic stretch*/
		fiC = fC_v;
		fbe.MultQBQT(F, fiC.Inverse());
		fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
		fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
		
		dWdE(fEigs_e, fdWdE_NEQ, fPot_NEQ, fR0_NEQ, fC0_NEQ); 
		fStress += fSpectralDecompSpat.EigsToRank2(fdWdE_NEQ);
	}
	
    fStress *= 1.0/sqrt(fEigs.Product());
//	cout << "\nstress: "<<fStress;
	return fStress;
}

/* material description */
const dMatrixT& IsoVECorneaModel::C_IJKL(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
    return fModulus;	
}

const dSymMatrixT& IsoVECorneaModel::S_IJ(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
    return fStress;
}

int IsoVECorneaModel::NumOutputVariables() const {return kNumOutputVar;} 

void IsoVECorneaModel::OutputLabels(ArrayT<StringT>& labels) const 
{ 
     /*allocates space for labels*/
     labels.Dimension(kNumOutputVar); 
  
     /*copy labels*/
     for (int i = 0; i< kNumOutputVar; i++) 
       labels[i] = Labels[i]; 
} 

void IsoVECorneaModel::ComputeOutput(dArrayT& output)
{
    /* spectral decomposition */
    Compute_b(fb);
    fSpectralDecompSpat.SpectralDecomp_Jacobi(fb, false);	
    fEigs = fSpectralDecompSpat.Eigenvalues();

	const dMatrixT& F = F_mechanical();
	/*load the viscoelastic principal stretches from state variable arrays*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());
   
	output[0] = 0.0;
	
	/*calc elastic stretch*/
	fiC = fC_v;
	fbe.MultQBQT(F, fiC.Inverse());
	fSpectralDecompSpat.SpectralDecomp_Jacobi(fbe, false);	
	fEigs_e = fSpectralDecompSpat.Eigenvalues(); 
	dWdE(fEigs_e, fdWdE_NEQ, fPot_NEQ, fR0_NEQ, fC0_NEQ);
	fStress = fSpectralDecompSpat.EigsToRank2(fdWdE_NEQ);
	double p = third*(fStress[0]+fStress[1]+fStress[2]);
	fStress[0] -= p;
	fStress[1] -= p;
	fStress[2] -= p;
	output[0] += 0.5*(0.5*fietaS*fStress.ScalarProduct()+fietaB*p*p);
}

/***********************************************************************
* Protected
***********************************************************************/

void IsoVECorneaModel::dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress, 
					const C1FunctionT* potential, const double R0, const double C0)
{
	ComputeLengths(eigenstretch2, R0);
	potential->MapDFunction(fLengths, fdU);
	
	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pr  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();

	double *p0  = fStressTable(0);
	double *p1  = fStressTable(1);
	double *p2  = fStressTable(2);
	
	/* PK2 principal values */	
	/*contribution from entropic part*/
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	double& s2 = eigenstress[2] = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pdU++)*R0/(*pr++)*(*pj++);
		s0 += factor*(*p0++);
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
	}
	s0 *= eigenstretch2[0];
	s1 *= eigenstretch2[1];
	s2 *= eigenstretch2[2];
	
	/*contribution from repulsion part*/
	s0 -= C0;
	s1 -= C0;
	s2 -= C0;	
}

/*calculates Dtau_A/Depsilon_B*/
void IsoVECorneaModel::ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
	dSymMatrixT& eigenmod, const C1FunctionT* potential, const double R0, const double C0)
{
	ComputeLengths(eigenstretch2, R0);
	potential->MapDFunction(fLengths, fdU);
	potential->MapDDFunction(fLengths, fddU);

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pr   = fLengths.Pointer();
	double* pj   = fjacobian.Pointer();

	/* stress */
	double* ps0  = fStressTable(0);
	double* ps1  = fStressTable(1);
	double* ps2  = fStressTable(2);

	/* modulus */
	double* pc00  = fModuliTable(0);
	double* pc11  = fModuliTable(1);
	double* pc22  = fModuliTable(2);

	double* pc12  = fModuliTable(3);
	double* pc02  = fModuliTable(4);
	double* pc01  = fModuliTable(5);
	
	/* PK2 principal values */	
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	double& s2 = eigenstress[2] = 0.0;

	double& c00 = eigenmod(0,0) = 0.0;
	double& c11 = eigenmod(1,1) = 0.0;
	double& c22 = eigenmod(2,2) = 0.0;

	double& c12 = eigenmod(1,2) = 0.0;
	double& c02 = eigenmod(0,2) = 0.0;
	double& c01 = eigenmod(0,1) = 0.0;

	/*contribution from entropic part*/
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pdU)*R0/(*pr)*(*pj);
		double cfactor = ((*pddU++) - (*pdU++)/(*pr))*(R0*R0*R0/((*pr)*(*pr)))*(*pj++);
		pr++;

		s0 += sfactor*(*ps0++);
		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);

		c00 += cfactor*(*pc00++);
		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);

		c12 += cfactor*(*pc12++);
		c02 += cfactor*(*pc02++);
		c01 += cfactor*(*pc01++);
	}
	
	/*contribution from repulsion part*/
	const double& l0 = eigenstretch2[0];
	const double& l1 = eigenstretch2[1];
	const double& l2 = eigenstretch2[2];
	
	s0 *= l0;
	s1 *= l1;
	s2 *= l2;

	s0 -= C0;
	s1 -= C0;
	s2 -= C0;	

	c00 *= l0*l0;
	c11 *= l1*l1;
	c22 *= l2*l2;
	c01 *= l1*l0;
	c02 *= l0*l2;
	c12 *= l1*l2;
	c00 += 2.0*C0 + 2.0*s0;	
	c11 += 2.0*C0 + 2.0*s1;	
	c22 += 2.0*C0 + 2.0*s2;	
}
/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void IsoVECorneaModel::ComputeLengths(const dArrayT& eigs, const double R0)
{
	double l0 = eigs[0];
	double l1 = eigs[1];
	double l2 = eigs[2];

	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s0 = fStressTable(0);
	double* s1 = fStressTable(1);
	double* s2 = fStressTable(2);
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = R0*sqrt(l0*(*s0++) + l1*(*s1++) + l2*(*s2++));
}

/* allocate memory for all the tables */
void IsoVECorneaModel::Dimension(int numbonds)
{
	/* length table */
	fLengths.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* jacobian table */
	fjacobian.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumModuli, numbonds);	
}


/***********************************************************************
* Private
***********************************************************************/
void IsoVECorneaModel::ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
			     dArrayT& eigenstress, dSymMatrixT& eigenmodulus, 
				 const C1FunctionT* potential, const double R0, const double C0,
				 const double ietaS, const double ietaB) 
{		
	const double ctol = 1.00e-14;
		
	/*set references to principle stretches*/
//    cout << "\neigenstretch: "<<eigenstretch; 
//    cout << "\neigenstretch_e: "<<eigenstretch_e; 
	 
	double& le0 = eigenstretch_e[0];
	double& le1 = eigenstretch_e[1];
	double& le2 = eigenstretch_e[2];
  
	double tol;

	/*initialize principle elastic and trial elastic log strains */
	double ep_tr0 = 0.5*log(le0);
	double ep_tr1 = 0.5*log(le1);
	double ep_tr2 = 0.5*log(le2);

	double ep_e0 = ep_tr0;		
	double ep_e1 = ep_tr1;	
	double ep_e2 = ep_tr2;
	

	/*initializes principle viscous stretch*/
	do 
	{
		ddWddE(eigenstretch_e, eigenstress, eigenmodulus, potential, R0, C0);
		
	    double s0 = eigenstress[0];
	    double s1 = eigenstress[1];
	    double s2 = eigenstress[2];
	    double sm = third*(s0+s1+s2);
		
		ComputeiKAB(eigenmodulus, ietaS, ietaB);
	    
	    /*calculate the residual*/
	    double dt = fFSMatSupport->TimeStep();
	    double res0 = ep_e0 + dt*(0.5*ietaS*(s0-sm) + third*ietaB*sm) - ep_tr0;
	    double res1 = ep_e1 + dt*(0.5*ietaS*(s1-sm) + third*ietaB*sm) - ep_tr1;
	    double res2 = ep_e2 + dt*(0.5*ietaS*(s2-sm) + third*ietaB*sm) - ep_tr2;
		

	    /*solve for the principal strain increments*/
	    double dep_e0=-fiKAB(0,0)*res0-fiKAB(0,1)*res1-fiKAB(0,2)*res2;
	    double dep_e1=-fiKAB(1,0)*res0-fiKAB(1,1)*res1-fiKAB(1,2)*res2;
	    double dep_e2=-fiKAB(2,0)*res0-fiKAB(2,1)*res1-fiKAB(2,2)*res2;
	    
	    /*updates principal elastic stretches*/ 
	    ep_e0 += dep_e0;
	    ep_e1 += dep_e1;
	    ep_e2 += dep_e2;
	    
	    le0 = exp(2.0*ep_e0);
	    le1 = exp(2.0*ep_e1);
	    le2 = exp(2.0*ep_e2);
		
	    
	    /*Check that the L2 norm of the residual is less than tolerance*/
	    tol = sqrt(res0*res0 + res1*res1+res2*res2);
	}while (tol>ctol); 
}

void IsoVECorneaModel::ComputeiKAB(const dSymMatrixT& eigenmodulus, const double ietaS, const double ietaB)
{	
	/*deviatoric values*/
	double c0 = eigenmodulus(0,0);
	double c1 = eigenmodulus(1,1);
	double c2 = eigenmodulus(2,2);

	double c12 = eigenmodulus(1,2);
	double c02 = eigenmodulus(0,2);
	double c01 = eigenmodulus(0,1);
	
	/*mean value*/
	double cm0 = third*(c0+c01+c02);
	double cm1 = third*(c1+c01+c12);
	double cm2 = third*(c2+c02+c12);

	dMatrixT& KAB = fiKAB;
		
	/*calculates  KAB = 1+dt*D(dWdE_Idev/nD+isostress/nV)/Dep_e*/

	double dt = fFSMatSupport->TimeStep();
	KAB(0,0) = 1+0.5*ietaS*dt*(c0-cm0)+third*ietaB*dt*cm0;
	KAB(1,1) = 1+0.5*ietaS*dt*(c1-cm1)+third*ietaB*dt*cm1;
	KAB(2,2) = 1+0.5*ietaS*dt*(c2-cm2)+third*ietaB*dt*cm2;

	KAB(1,2) = 0.5*ietaS*dt*(c12-cm2)+third*ietaB*dt*cm2;
	KAB(0,2) = 0.5*ietaS*dt*(c02-cm2)+third*ietaB*dt*cm2;
	KAB(0,1) = 0.5*ietaS*dt*(c01-cm1)+third*ietaB*dt*cm1;
       
	KAB(2,1) = 0.5*ietaS*dt*(c12-cm1)+third*ietaB*dt*cm1;
	KAB(2,0) = 0.5*ietaS*dt*(c02-cm0)+third*ietaB*dt*cm0;
	KAB(1,0) = 0.5*ietaS*dt*(c01-cm0)+third*ietaB*dt*cm0;
	
//	cout << "\nKAB: "<<KAB;
	/*inverts KAB*/
	fiKAB.Inverse(KAB);
}


/* Initialize angle tables */
void IsoVECorneaModel::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fSphere->SpherePoints(0.0, 0.0);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Dimension(numpoints);
	
	/* fetch jacobians */
	fjacobian = fSphere->Jacobians();
	
	/* set pointers */
	double *s0 = fStressTable(0);
	double *s1 = fStressTable(1);
	double *s2 = fStressTable(2);

	double *c00 = fModuliTable(0);
	double *c11 = fModuliTable(1);
	double *c22 = fModuliTable(2);

	double *c12 = fModuliTable(3);
	double *c02 = fModuliTable(4);
	double *c01 = fModuliTable(5);

	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);

		double xsi0 = xsi[0];
		double xsi1 = xsi[1];
		double xsi2 = xsi[2];
		
		/* stress angle tables */
		s0[i] = xsi0*xsi0;
		s1[i] = xsi1*xsi1;
		s2[i] = xsi2*xsi2;
	
		/* moduli angle tables */
		c00[i] = s0[i]*s0[i];
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];

		c12[i] = s1[i]*s2[i];
		c02[i] = s0[i]*s2[i];
		c01[i] = s0[i]*s1[i];
	}
}
