/* created: TDN (01/22/2001) */

#ifndef _PotentialT_
#define _PotentialT_

#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "fstreamT.h"

namespace Tahoe {

class fstreamT;

class PotentialT
{
  public:
  enum TypesT {kNeoHookean = 1,
		kOgden = 2};

  /* print parameters */
  virtual void Print(ostream& out) const = 0;
  virtual void PrintName(ostream& out) const = 0;
  virtual void Initialize(void) = 0;
  
  /*free energy density*/
  virtual double Energy(const dArrayT& lambda_bar,const double& J) = 0;
  
  /*Kirchoff stress measures*/
  virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau) = 0;
  virtual double MeanStress(const double& J) = 0;
  
  /*derivative of Kirchoff stress with log strain*/
  virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus) = 0;
  virtual double MeanMod(const double& J) = 0;

  virtual const double Kappa() = 0;
  virtual const double Mu() = 0;
};
}
#endif /* _PotentialT_ */
