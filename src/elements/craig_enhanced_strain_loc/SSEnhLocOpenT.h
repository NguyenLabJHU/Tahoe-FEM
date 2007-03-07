#ifndef _SS_ENH_LOC_OPEN_T_H
#define _SS_ENH_LOC_OPEN_T_H

/* base class */
#include "SmallStrainT.h"

#include "SolidElementT.h"
#include "OpenBandT.h"

#include "HookeanMatT.h"
#include "MapT.h"

#include "ofstreamT.h"

namespace Tahoe {

  /* forward declarations */
  class SSMatSupportT;
  
  class SSEnhLocOpenT: public SmallStrainT //: public SSEnhLocCraigT
    {
    public:
	
	  enum BandStateT { kInactive = 0, kElastic, kDamage};	

      /* constructor */
      SSEnhLocOpenT(const ElementSupportT& support);

      /** describe the parameters needed by the interface */
      virtual void DefineParameters(ParameterListT& list) const;
	  
	  /* information about subordinate parameter lists */
	  virtual void DefineSubs(SubListT& sub_list) const;

      /** accept parameter list */
      virtual void TakeParameterList(const ParameterListT& list);
	  
	  enum BVPTypeT {kNonhomogeneous = 0,
			kHomogeneous,
			kPreFailed,
			kNumBVPTypes};

    protected:
	  void FormKd(double constK);
      virtual void FormStiffness(double constK);
      virtual dArrayT CalculateJump();
      //virtual bool IsBandActive();
	  

	  virtual dArrayT CalculateElasticJump();
	  virtual dArrayT CalculatePlasticJump();
	  virtual double EffectiveTraction(dArrayT jump);
	  
	  virtual dArrayT ElasticTractionBalance(dArrayT jump);
	  virtual dArrayT PlasticTractionBalance(dArrayT jump);
	  virtual dMatrixT ElasticBalanceGradient(dArrayT jump);
	  virtual dMatrixT PlasticBalanceGradient(dArrayT jump);
	  
	  virtual double CurrentCohesion(dArrayT jump);
	  virtual double EffectiveJump(dArrayT jump);

	  virtual dArrayT CohesionlessTractionBalance(dArrayT jump);
	  virtual dMatrixT CohesionlessBalanceGradient(dArrayT jump);
	  
	  virtual double DNormalStressDNormalSlip();
	  virtual double DNormalStressDShearSlip();
	  virtual double DShearStressDNormalSlip();
	  virtual double DShearStressDShearSlip();
	  virtual double QuadContraction(dMatrixT c_ijkl, 
		dArrayT vector1, dArrayT vector2, dArrayT vector3, dArrayT vector4);
	  virtual dArrayT AvgGradActive(int ndof);

      virtual OpenBandT* FormNewBand(dArrayT normal, dArrayT slipDir,
				 dArrayT perpSlipDir, dArrayT coords, double area);
      
      virtual void CloseStep(void);
	  bool IsElementTraced();
	  bool IsElementTraced(int elementNumber);
	  bool IsElementLocalized();
	  bool TraceElement();
      virtual void ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs);
      dArrayT InterceptCoords(dArrayT& localizedEleCoord,
		dArrayT& nodalCoord1, dArrayT& nodalCoord2);
	  virtual dArrayT Centroid();
	  virtual void AddNewEdgeElements(int elementNumber);
	
	  virtual void GetElement(int elementNumber);
	
      virtual void LoadBand(int elementNumber);
      
	  /*Functions for nonstandard bvp types */
	  virtual void PreFailElements();
	  
	  
	  virtual dSymMatrixT StressIncrOnBand(dArrayT jumpIncrement);
      virtual dSymMatrixT LastStressOnBand();
      virtual double NormalStress(dArrayT jump);
	  virtual double ShearStress(dArrayT jump);

      virtual dSymMatrixT FormdGdSigma(int ndof);
      virtual dMatrixT FormdGsdSigma(int ndof);

	  virtual dSymMatrixT FormGradActiveTensorJumpIncr(dArrayT jumpIncr, int ndof, int ip);
	  virtual dSymMatrixT FormGradActiveTensorShearDir(int ndof, int ip);
	  virtual dMatrixT FormGradActiveTensorBandDirs(int ndof, int ip);
			
	protected:
    

	MapT<int, OpenBandT*>* TracedElements()
	  {return &fTracedElements;}  
	
	
    private:

      OpenBandT* fOpenBand;

      /*band parameters*/

	  double fFrictionCoeff;
	  double fAlpha_sigma;
	  double fAlpha_zeta;
	  double fZeta_star;
    
	/* static variables creating linking error - investigate */
    static bool fLocalizationHasBegun;
	static bool fSeedElementsSet;
	static double fDetAMin;
	static int fLeastDetEle;

	//bool fLocalizationHasBegun;	
	//bool fSeedElementsSet;
	//double fDetAMin;
	//int fLeastDetEle;
	  
	MapT<int, OpenBandT*> fTracedElements;
	  

	int fFirstElementToLocalize; // -1 if loading is nonhomogeneous, element number otherwise
	
	int fBVPType;
	dArrayT fPreFailedNormal;
	dArrayT fPreFailedSlipDir;

	bool fMultiBand;

	
	
	/*elements with one edge on band of traced elements
	 * with that coordinate of the end of band*/
	iAutoArrayT fEdgeOfBandElements;
	AutoArrayT<dArrayT> fEdgeOfBandCoords;
	  
    }; //end class declaration

}//end namespace Tahoe

#endif // _SS_ENH_LOC_OPEN_T_H_


