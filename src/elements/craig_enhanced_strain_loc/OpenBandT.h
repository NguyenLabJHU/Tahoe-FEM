#include "dArrayT.h"
#include "ArrayT.h"
//#include "AutoArrayT.h"
#include "SSEnhLocOpenT.h"
#include "iAutoArrayT.h"

#include "SmallStrainT.h"

#ifndef _OPENBAND_T_H_
#define _OPENBAND_T_H_

namespace Tahoe {

  /* forward delaration */
  class SSEnhLocOpenT;

class OpenBandT
  {
  public:
  
	//enum BandStateT { kInactive = 0, kElastic, kDamage};
  
    OpenBandT(const dArrayT normal, dArrayT shearDir, dArrayT &coord, double residCohesion,
	 ArrayT<dSymMatrixT> stressList, SSEnhLocOpenT *element); 
    
    virtual const iAutoArrayT& ActiveNodes() const;
    virtual const dArrayT& Normal() const;
	virtual const dArrayT& ShearDir() const;
    //virtual const dArrayT& SlipDir() const;
    //virtual const dArrayT& PerpSlipDir() const;
    //virtual double H_delta() const;
	virtual double InitialCohesion() const;
	virtual double Cohesion() const;
	virtual void UpdateCohesion(double cohesion);
    //virtual double ResidualCohesion() const;
    virtual dArrayT Jump() const;
    virtual dArrayT LastJump() const;
	//virtual double JumpIncrement() const;
    virtual void UpdateJump ();
    virtual void StoreJump(dArrayT jump);
    //void CloseStep();
    virtual dSymMatrixT Stress_List(int ip);
    virtual void IncrementStress(dSymMatrixT stressIncr, int ip);
    //virtual void UpdateCohesion();
    //virtual void SetEffectiveSoftening(double effectiveSoftening);
    //virtual double EffectiveSoftening();
    virtual void SetBandState(int state);
    virtual int BandState();
    virtual void FlipSlipDir();
    virtual void CloseStep();
    virtual dArrayT& Coords();
 

  private:

  virtual void SetEndPoints(dArrayT& coord);
  virtual void ActivateNodes(dArrayT& coord);

    int kNSD;

    dArrayT fNormal;
	dArrayT fShearDir; //new
    //dArrayT fSlipDir;
    //dArrayT fPerpSlipDir;
    dArrayT fCoords; //coords of one point on band
    dArrayT fJump;
    dArrayT fLastJump;
    //ArrayT <dArrayT> fEndPoints; //not used yet
    iAutoArrayT fActiveNodes;
    ArrayT<dSymMatrixT> fStress_List;
	double fCohesion;
	double fLastCohesion;
	double fInitialCohesion;
    //double fResidualCohesion;
    //double fH_delta;
    //double fEffectiveSoftening;
    int fBandState;
    SSEnhLocOpenT *currentElement;

  };
  
}

#endif //_OPENBAND_T_H_
