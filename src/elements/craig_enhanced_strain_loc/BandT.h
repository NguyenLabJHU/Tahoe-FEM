#include "dArrayT.h"
#include "ArrayT.h"
//#include "AutoArrayT.h"
#include "SSEnhLocCraigT.h"
#include "iAutoArrayT.h"

#include "SmallStrainT.h"

#ifndef _BAND_T_H_
#define _BAND_T_H_

namespace Tahoe {

  /* forward delaration */
  class SSEnhLocCraigT;

class BandT
  {
  public:

    BandT(const dArrayT normal, const dArrayT slipDir, const dArrayT
    perpSlipDir, dArrayT &coord, double h_delta, double residCohesion, ArrayT<dSymMatrixT> stressList, SSEnhLocCraigT *element); 
    
    virtual const iAutoArrayT& ActiveNodes() const;
    virtual const dArrayT& Normal() const;
    virtual const dArrayT& SlipDir() const;
    virtual const dArrayT& PerpSlipDir() const;
    virtual double H_delta() const;
    virtual double ResidualCohesion() const;
    virtual double Jump() const;
    virtual double JumpIncrement() const;
    //virtual void IncrementJump ();
    virtual void StoreJumpIncrement(double increment);
    //void CloseStep();
    virtual dSymMatrixT Stress_List(int ip);
    virtual void IncrementStress(dSymMatrixT stressIncr, int ip);
    //virtual void UpdateCohesion();
    virtual void SetEffectiveSoftening(double effectiveSoftening);
    virtual double EffectiveSoftening();
    virtual void SetActive(bool active);
    virtual bool IsActive();
    virtual void FlipSlipDir();
    virtual void CloseStep();
    virtual dArrayT& Coords();
 

  private:

  virtual void SetEndPoints(dArrayT& coord);
  virtual void ActivateNodes(dArrayT& coord);

    int kNSD;

    dArrayT fNormal;
    dArrayT fSlipDir;
    dArrayT fPerpSlipDir;
    dArrayT fCoords; //coords of one point on band
    double fLength; //not used ?
    double fJump;
    double fJumpIncrement;
    ArrayT <dArrayT> fEndPoints; //not used yet
    iAutoArrayT fActiveNodes;
    ArrayT<dSymMatrixT> fStress_List;
    double fResidualCohesion;
    double fH_delta;
    double fEffectiveSoftening;
    double fIsBandActive;
    SSEnhLocCraigT *currentElement;

  };
  
}

#endif //_BAND_T_H_
