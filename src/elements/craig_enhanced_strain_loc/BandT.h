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

    BandT(const dArrayT normal, const dArrayT slipDir, const dArrayT perpSlipDir, dArrayT &coord, SSEnhLocCraigT *element); 
    
    const iAutoArrayT& ActiveNodes() const;
    const dArrayT& Normal() const;
    const dArrayT& SlipDir() const;
    const dArrayT& PerpSlipDir() const;
    double Jump() const;
    void IncrementJump (const double increment);



  private:

  void SetEndPoints(dArrayT& coord);
  void ActivateNodes(dArrayT& coord);

    int kNSD;

    dArrayT fNormal;
    dArrayT fSlipDir;
    dArrayT fPerpSlipDir;
    double fLength; //not used ?
    double fJump;
    double fJumpIncrement; //not used ?
    ArrayT <dArrayT> fEndPoints; //not used yet
    iAutoArrayT fActiveNodes;
    SSEnhLocCraigT *currentElement;
  };
  
}

#endif //_BAND_T_H_
