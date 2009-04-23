/* created: TDN (01/22/2001) */

#ifndef _PotentialT_
#define _PotentialT_

#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "fstreamT.h"
/* base class */
#include "ParameterInterfaceT.h"

namespace Tahoe {

class fstreamT;

class PotentialT:public ParameterInterfaceT
{
  public:
  PotentialT();
  ~PotentialT();
  
  enum TypesT {kNeoHookean = 1, kMooneyRivlin = 2, kVerondaWestmann = 3};
  
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
  
  /*free energy density*/
  virtual double Energy(const dArrayT& lambda_bar, const double& J, double temperature = 1.0) = 0;
  
  /*Kirchoff stress measures*/
  virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau, double temperature = 1.0) = 0;
  virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus, double temperature = 1.0) = 0;
  
  /*derivative of Kirchoff stress with log strain*/
  virtual double MeanMod(const double& J);
  virtual double MeanEnergy(const double& J);
  virtual double MeanStress(const double& J);

  void SetKappa(double kappa);
  const double GetKappa(void) const; 
  const double GetMu(void) const; 
  
  protected:
  double fKappa;
  double fMu;

};
inline const double PotentialT::GetKappa(void) const {return fKappa; } 
inline const double PotentialT::GetMu(void) const {return fMu; } 

}
#endif /* _PotentialT_ */
