/* $Id: OgdenMaterialT.cpp,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: tdn (3/17/2003) */
#include "OgdenMaterialT.h"
#include "PotentialT.h"
#include "NeoHookean.h"

#include "ifstreamT.h"

#include <iostream.h>
#include <math.h>
#include <iomanip.h>

using namespace Tahoe;

/* constructor */
OgdenMaterialT::OgdenMaterialT(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("Ogden_material"),
	OgdenBaseT(in, support),
	fthird(1.0/3.0)
{
ExceptionT::GeneralFail("OgdenMaterialT::OgdenMaterialT", "out of date");
#if 0
  /*read in potential code*/
  int code;
  in >> code;
  switch(code)
  {
     case PotentialT::kNeoHookean: 
     {
       fPot = new NeoHookean(in);
       break;
     }
     default:
     {
       throw ExceptionT::kBadInputValue;
     }
  }
  cout.precision(12);
#endif
}

OgdenMaterialT::~OgdenMaterialT(void)
{
  delete fPot;
}

double OgdenMaterialT::StrainEnergyDensity(void)
{
    /*calculates deviatoric and volumetric part of the total stretch */
  Compute_b(fb);
  fb.PrincipalValues(fEigs);
  double J = sqrt(fEigs.Product());

  dArrayT eigenstretch_bar(3);
  if (NumSD() == 2)
  {
    eigenstretch_bar[0]=fEigs[0];
    eigenstretch_bar[1]=fEigs[1];
    eigenstretch_bar[2]=1.0;
  }
  else eigenstretch_bar = fEigs;
  eigenstretch_bar *= pow(J, -2.0*fthird);
  
  double energy =0.0;
  energy = fPot->Energy(eigenstretch_bar, J);

  return(energy);
}

/* principal values of the PK2 stress given principal values of the stretch 
 * tensors, i.e., the principal stretches squared */

void OgdenMaterialT::dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress)
{
  double J = sqrt(eigenstretch.Product());

  dArrayT eigenstretch_bar(3);
  if (NumSD() == 2)
  {
    eigenstretch_bar[0]=eigenstretch[0];
    eigenstretch_bar[1]=eigenstretch[1];
    eigenstretch_bar[2]=1.0;
  }
  else eigenstretch_bar = eigenstretch;
  eigenstretch_bar *= pow(J, -2.0*fthird);

  /*evaluates Kirchoff stress*/
  fPot->DevStress(eigenstretch_bar, eigenstress);
  eigenstress += fPot->MeanStress(J);
  eigenstress /= J;
}

void OgdenMaterialT::ddWddE(const dArrayT& eigenstretch, dArrayT& eigenstress,
			  dSymMatrixT& eigenmod)
{
  double J = sqrt(eigenstretch.Product());

  dArrayT eigenstretch_bar(3);
  if (NumSD() == 2)
  {
    eigenstretch_bar[0]=eigenstretch[0];
    eigenstretch_bar[1]=eigenstretch[1];
    eigenstretch_bar[2]=1.0;
  }
  else eigenstretch_bar = eigenstretch;
  eigenstretch_bar *= pow(J, -2.0*fthird);
  
  /*evaluates Kirchoff stress*/
  fPot->DevStress(eigenstretch_bar, eigenstress);
  eigenstress += fPot->MeanStress(J);

  /*evaluates dtau_de*/
  fPot->DevMod(eigenstretch_bar,eigenmod);
  eigenmod += fPot->MeanMod(J);

  /*  cout << "\neigenmod: "<<eigenmod;
      cout << "\neigenstress: "<<eigenstress;*/

  if (NumSD() == 2)
  {
    eigenmod[0] -= 2.0*eigenstress[0];
    eigenmod[1] -= 2.0*eigenstress[1];
  }
  else
  {
    /*transform moduli to 1/lam dS_dlam*/
    eigenmod[0] -= 2.0*eigenstress[0];
    eigenmod[1] -= 2.0*eigenstress[1];
    eigenmod[2] -= 2.0*eigenstress[2];
  }

  eigenmod /= J;
  eigenstress /= J;
}

