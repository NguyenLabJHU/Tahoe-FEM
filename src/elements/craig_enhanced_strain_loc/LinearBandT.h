#include "dArrayT.h"
#include "ArrayT.h"
#include "SSEnhLocLinearT.h"
#include "iAutoArrayT.h"

#include "SmallStrainT.h"

#ifndef _LINEAR_BAND_T_H_
#define _LINEAR_BAND_T_H_

namespace Tahoe {

  /* forward delaration */
  //class SSEnhLocLinearT;

class LinearBandT
  {
  public:
  
	friend class SSEnhLocLinearT;
  
  public:

    LinearBandT(const dArrayT normal, const dArrayT slipDir, const dArrayT
    perpSlipDir, dArrayT &coord, double h_delta, 
	ArrayT<dSymMatrixT> stressList, SSEnhLocLinearT *element); 
    
    virtual const iAutoArrayT& ActiveNodes() const;
    virtual const dArrayT& Normal() const;
    virtual const dArrayT& SlipDir(int bandIP) const;
    virtual const dArrayT& PerpSlipDir(int bandIP) const;
    virtual double H_delta(int bandIP) const;
    virtual double ResidualCohesion(int ip) const;
    virtual double Jump(int ip) const;
    virtual double JumpIncrement(int ip) const;
    virtual void IncrementJump (int ip);
    virtual void StoreJumpIncrement(double increment, int bandIP);
    //void CloseStep();
    virtual dSymMatrixT Stress_List(int ip);
	virtual double IPBandCoord(int ip);
	virtual dArrayT SurfaceIPTraction(int ip) {return fSurfaceIPTractions[ip];};
	virtual dArrayT SurfaceIPCoords(int ip) {return fSurfaceIPCoords[ip];};
	virtual double DistanceBetweenIPs() {return fDistanceBetweenSurfaceIPs;};
    virtual void IncrementStress(dSymMatrixT stressIncr, int ip);
	virtual void IncrementTractionAtBandIP(dArrayT increment, int bandIP); 
    //virtual void UpdateCohesion();
    virtual void SetEffectiveSoftening(double effectiveSoftening, int bandIP);
    virtual double EffectiveSoftening(int bandIP);
    virtual void SetActive(int ip, bool active);
    virtual bool IsActive(int ip);
    virtual void FlipSlipDir(int ip);
    virtual void CloseStep();
    virtual dArrayT& Coords();
	virtual bool IsBulkIPActive(int ip) {return fIsBulkIPActive[ip];};
	virtual int NumSurfaceIPs() {return kNumSurfaceIPs;};
	
  protected:
  
  ArrayT<dSymMatrixT> NodalStressList();

  private:
  
  virtual void SetEndPoints(dArrayT& coord);
  virtual void ActivateNodes(dArrayT& coord);
  virtual void ActivateBulkIPs(dArrayT& coord);


    int kNSD;
	static const int kNumSurfaceIPs = 2;

    dArrayT fNormal;
    ArrayT<dArrayT> fSlipDir;
    ArrayT<dArrayT> fPerpSlipDir;
    dArrayT fCoords; //coords of one point on band
    dArrayT fJump;
    dArrayT fJumpIncrement;
    iAutoArrayT fActiveNodes;
    ArrayT<dSymMatrixT> fStress_List;
	dArrayT fIPBandCoords;
	ArrayT<dArrayT> fSurfaceIPTractions;
	ArrayT<dArrayT> fSurfaceIPCoords;
	double fDistanceBetweenSurfaceIPs;
    dArrayT fResidualCohesion;
    dArrayT fH_delta;
    dArrayT fEffectiveSoftening;
    ArrayT<bool> fIsBandActive;
    SSEnhLocLinearT *fCurrentElement;
	ArrayT<bool> fIsBulkIPActive;

  };
  
}

#endif //_LINEAR_BAND_T_H_
