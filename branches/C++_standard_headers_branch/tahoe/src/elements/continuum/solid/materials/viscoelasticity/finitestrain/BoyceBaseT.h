/* $Id: BoyceBaseT.h,v 1.2 2007-07-25 14:47:29 tdnguye Exp $ */
/* created : TDN (1/22/2001) */
#ifndef _BOYCE_BASE_T_H_
#define _BOYCE_BASE_T_H_

/* base classes */
#include "FSSolidMatT.h"
#include "SpectralDecompT.h"

namespace Tahoe {

/** base class for nonlinear finite deformation viscoelasticity **/

class BoyceBaseT: public FSSolidMatT
{
  public:
  
	/* constructor */
	BoyceBaseT(void);

	enum LocIntCodeT {
    kExplicit = 0,
          kImplicit = 1
	};

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const {
		ExceptionT::GeneralFail("BoyceBaseT::Pressure", "not implemented");
		return 0.0;
	};

	/** return true if the material has history variables */
	virtual bool HasHistory(void) const { return true; };
	
	/*Initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; 
	virtual void PointInitialize(void);              

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time
	/* apply pre-conditions at the current time step */
	virtual void InitStep(void){ FSSolidMatT::InitStep(); };
	
	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	
 protected:
	
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  		
  protected:
	/*internal state variables. Dimension numprocess<nsd x nsd>*/
	dMatrixT fFv;
	dMatrixT fFv_n;
	
	/*number of state variables*/
	int fnstatev;
	
	/* internal state variables array*/
	dArrayT fstatev;

   private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompRef;
	
	
};

}

#endif /* _BOYCE_BASE_T_H_ */

