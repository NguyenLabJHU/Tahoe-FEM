/* $Id: IsoCorneaModel.cpp,v 1.2 2006-11-12 18:26:54 thao Exp $ */
/* created: paklein (11/08/1997) */

#include "IsoCorneaModel.h"
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

#ifdef VIB_MATERIAL
/* point generators */
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"
#include "FCCPtsT.h"

using namespace Tahoe;
static const double pi = 3.14159274641021;
/* constructors */
IsoCorneaModel::IsoCorneaModel(void):
	OgdenIsotropicT(),
	ParameterInterfaceT("Isotropic_Cornea_Model"),
	fNumSD(3),
	fNumStress(3),
	fNumModuli(6),
	fSphere(NULL),
	fPotential(NULL)
{
#ifndef VIB_MATERIAL
	ExceptionT::BadInputValue("IsoCorneaModel::IsoCorneaModel", 
		"VIB_MATERIAL must be enabled");
#endif
}

/* destructor */
IsoCorneaModel::~IsoCorneaModel(void) 
{ 
  delete fSphere; 
  delete fPotential; 
}

/* information about subordinate parameter lists */
void IsoCorneaModel::DefineSubs(SubListT& sub_list) const
{
  /* inherited */
  OgdenIsotropicT::DefineSubs(sub_list);

  /* choice of integration schemes */
  sub_list.AddSub("sphere_integration_choice", ParameterListT::Once, true);
  sub_list.AddSub("potential_choice", ParameterListT::Once, true);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* IsoCorneaModel::NewSub(const StringT& name) const
{
  /* inherited */
  ParameterInterfaceT* sub = OgdenIsotropicT::NewSub(name);
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
  else if (name == "potential_choice")
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
    ParameterContainerT fungrep("fung_type_chain_power_repulsion");
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
}

/* accept parameter list */
void IsoCorneaModel::TakeParameterList(const ParameterListT& list)
{
  /* inherited */
    OgdenIsotropicT::TakeParameterList(list);
  const ParameterListT& potential = list.GetListChoice(*this, "potential_choice");
  //  cout << "\npotential name: "<<potential.Name()<<endl;

  if (potential.Name() == "worm_like_chain")
  {
      double N = potential.GetParameter("chain_density");
      double T = potential.GetParameter("temperature");
      double A = potential.GetParameter("persistence_length");
      double R0 = potential.GetParameter("initial_coil_length");
      double L = R0*R0/(2*A);
      double k = 1.3806503e-23; 
		
      /*set wlc potential function*/
      double K = 0.75/pi*(N*k*T/A);
      fPotential = new WormLikeChain(K,L);
      if (!fPotential) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      double p = sqrt(2*A/L);
      fC0 = (N*k*T/A) * ( p - 0.25 + 0.25/((1-p)*(1-p)) );
      fR0 = R0;
  }
  else if (potential.Name() == "worm_like_chain_power_repulsion")
  {
      double N = potential.GetParameter("chain_density");
      double T = potential.GetParameter("temperature");
      double A = potential.GetParameter("persistence_length");
      double R0 = potential.GetParameter("initial_coil_length");
      double n = potential.GetParameter("repulsion_exponent");
      double L = R0*R0/(2*A);
      double k = 1.3806503e-23; 
		
      /*set wlc potential function*/
      double K = 0.75/pi*(N*k*T/A);
      double p = sqrt(2*A/L);
      double C = 0.75/pi*(N*k*T/A) * ( p - 0.25 + 0.25/((1-p)*(1-p)) )*pow(R0,(n+1.0))/n;
      fPotential = new WLCwRep(K,L,C,n);
      if (!fPotential) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0 = 0.0;
      fR0 = R0;
  }
  else if (potential.Name() == "fung_type_chain")
    {
      double C1 = potential.GetParameter("multiplier_C1");
      double beta = potential.GetParameter("exponent_beta");

      /*set fung type potential*/
      double K = 0.75/pi*C1;
      fPotential = new FungType(K,beta);
      if (!fPotential) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0 =2.0* C1*beta;
      fR0 = 1.0;      
    }
  else if (potential.Name() == "fung_type_chain_power_repulsion")
    {
      double C1 = potential.GetParameter("multiplier_C1");
      double beta = potential.GetParameter("exponent_beta");
      double n = potential.GetParameter("repulsion_exponent");

      /*set fung type potential*/
      double K = 0.75/pi*C1;
      double C = 0.75/pi*2.0*C1*beta/n;
	  fPotential = new FungwRep(K,beta,C,n);
      if (!fPotential) throw ExceptionT::kOutOfMemory;

      /*set parameters for repulsive potential*/
      fC0 = 0.0;
      fR0 = 1.0;      
    }
  else
    ExceptionT::GeneralFail("IsoCorneaModel::TakeParameterList", "unrecognized potential \"%s\"", potential.Name().Pointer());

  const ParameterListT& points = list.GetListChoice(*this, "sphere_integration_choice");
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
    ExceptionT::GeneralFail("IsoCorneaModel::TakeParameterList", "unrecognized point scheme \"%s\"", points.Name().Pointer());

  /* set tables */
  Construct();
}

/* strain energy density */
double IsoCorneaModel::StrainEnergyDensity(void)
{
	/* stretch */
	Compute_C(fC);
	
	/* principal stretches */
	fC.PrincipalValues(fEigs);
	
	/* stretched bonds */
	ComputeLengths(fEigs);

	/* update potential table */
	fPotential->MapFunction(fLengths,fU);

	/* sum contributions to entropic component*/
	double  energy = 0.0;
	double* pU = fU.Pointer();
	double* pj = fjacobian.Pointer();
	double iR0 = 1.0/fR0;	
	for (int i = 0; i < fLengths.Length(); i++)
		energy += (*pU++)*(*pj++)*iR0;
		
	/*add contribution from repulsion component*/
	double J = fEigs.Product();
	energy -= 0.5*fC0*log(J);
	
	return energy;
}

/***********************************************************************
* Protected
***********************************************************************/

/* principal values given principal values of the stretch tensors,
 * i.e., the principal stretches squared */
void IsoCorneaModel::dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch2);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

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
		double factor = (*pdU++)*fR0/(*pr++)*(*pj++);
		s0 += factor*(*p0++);
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
	}
	/*contribution from repulsion part*/
	s0 -= fC0/eigenstretch2[0];
	s1 -= fC0/eigenstretch2[1];
	s2 -= fC0/eigenstretch2[2];	
}

void IsoCorneaModel::ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
	dSymMatrixT& eigenmod)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch2);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
	fPotential->MapDDFunction(fLengths, fddU);	

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
		double sfactor = (*pdU)*fR0/(*pr)*(*pj);
		double cfactor = ((*pddU++) - (*pdU++)/(*pr))*(fR0*fR0*fR0/((*pr)*(*pr)))*(*pj++);
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
	
	s0 -= fC0/l0;
	s1 -= fC0/l1;
	s2 -= fC0/l2;	

	c00 += 2.0*fC0/(l0*l0);	
	c11 += 2.0*fC0/(l1*l1);	
	c22 += 2.0*fC0/(l2*l2);	
	
}

/* allocate memory for all the tables */
void IsoCorneaModel::Dimension(int numbonds)
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

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void IsoCorneaModel::ComputeLengths(const dArrayT& eigs)
{
	double C0 = eigs[0];
	double C1 = eigs[1];
	double C2 = eigs[2];

	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s0 = fStressTable(0);
	double* s1 = fStressTable(1);
	double* s2 = fStressTable(2);
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = fR0*sqrt(C0*(*s0++) + C1*(*s1++) + C2*(*s2++));
}

/***********************************************************************
* Private
***********************************************************************/

/* Initialize angle tables */
void IsoCorneaModel::Construct(void)
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
#endif /*VIB_MATERIAL*/
