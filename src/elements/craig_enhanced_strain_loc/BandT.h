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
    
    const iAutoArrayT& ActiveNodes() const;
    const dArrayT& Normal() const;
    const dArrayT& SlipDir() const;
    const dArrayT& PerpSlipDir() const;
    double H_delta() const;
    double ResidualCohesion() const;
    double Jump() const;
    double JumpIncrement() const;
    void IncrementJump ();
    void StoreJumpIncrement(double increment);
    //void CloseStep();
    dSymMatrixT Stress_List(int ip);
    void IncrementStress(dSymMatrixT stressIncr, int ip);
    void UpdateCohesion();
    void SetEffectiveSoftening(double effectiveSoftening);
    double EffectiveSoftening();
    void SetActive(bool active);
    bool IsActive();
    void FlipSlipDir();


  private:

  void SetEndPoints(dArrayT& coord);
  void ActivateNodes(dArrayT& coord);

    int kNSD;

    dArrayT fNormal;
    dArrayT fSlipDir;
    dArrayT fPerpSlipDir;
    double fLength; //not used ?
    double fJump;
    double fJumpIncrement;
    ArrayT <dArrayT> fEndPoints; //not used yet
    iAutoArrayT fActiveNodes;
    SSEnhLocCraigT *currentElement;
    ArrayT<dSymMatrixT> fStress_List;
    double fResidualCohesion;
    double fH_delta;
    double fEffectiveSoftening;
    double fIsBandActive;

  };
  
}

#endif //_BAND_T_H_
