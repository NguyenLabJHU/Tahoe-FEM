/* $Id: RGViscoelasticityT.h,v 1.1 2003-04-05 20:06:59 thao Exp $ */
/* created : TDN (1/22/2001) */
#ifndef _RG_VISCO_T_H_
#define _RG_VISCO_T_H_

/* base classes */
#include "FSSolidMatT.h"
#include "IsotropicT.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

/** base class for large deformation isotropic material following
 * Ogden's formulation */
class RGViscoelasticityT: public FSSolidMatT, public IsotropicT
{
  public:
  
	/* constructor */
	RGViscoelasticityT(ifstreamT& in, const FSMatSupportT& support);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const 
	{
		cout << "\n RGViscoelasticityT::Pressure: not implemented" << endl;
		throw ExceptionT::kGeneralFail;
		return 0.0;
	};

	/* print parameters */	
	virtual	void Print(ostream& out) const;	
	virtual void PrintName(ostream& out) const;
	
	/** initialization called immediately after constructor. This function
	 * dimensions and set source for viscous history variables */
	virtual void Initialize(void);
	
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

 protected:
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  		

  protected:
	/*internal state variables*/
	dSymMatrixT     fC_v;
	dSymMatrixT     fC_vn;
		
	/*number of state variables*/
	int fnstatev;
	
	/* internal state variables array*/
	dArrayT fstatev;
	
	/*dof of internal variables*/
	int fndof;
};
}

#endif /* _RG_VISCO_T_H_ */

