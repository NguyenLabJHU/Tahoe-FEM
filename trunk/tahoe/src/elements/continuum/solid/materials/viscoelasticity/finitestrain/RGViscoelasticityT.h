/* $Id: RGViscoelasticityT.h,v 1.3 2006-08-03 23:13:34 tdnguye Exp $ */
/* created : TDN (1/22/2001) */
#ifndef _RG_VISCO_T_H_
#define _RG_VISCO_T_H_

/* base classes */
#include "FSSolidMatT.h"
#include "SpectralDecompT.h"

namespace Tahoe {

/** base class for nonlinear finite deformation viscoelasticity **/

class RGViscoelasticityT: public FSSolidMatT
{
  public:
  
	/* constructor */
	RGViscoelasticityT(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const {
		ExceptionT::GeneralFail("RGViscoelasticityT::Pressure", "not implemented");
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

	/*Returns eigenvalues of viscous deformation gradient
	Assumes that current values of Cv and Cvn have been loaded using Load(ElementCardT& element, int ip)*/
	const dArrayT& Compute_Eigs_v(const int process_id);
	const dArrayT& Compute_Eigs_vn(const int process_id);
	
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/* Dimension internal state variables*/
	void SetStateVariables (const int numprocess);
	
 protected:
	
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  		
  protected:
	/*internal state variables. Dimension numprocess<nsd x nsd>*/
	ArrayT<dSymMatrixT> fC_v;
	ArrayT<dSymMatrixT> fC_vn;
	
	/* number of nonequilibrium processes*/
	/* must be set in derived classes before TakeParameterList is called*/
	/* default value is 1*/
	int fNumProcess;
	
	/*number of state variables*/
	int fnstatev;
	
	/* internal state variables array*/
	dArrayT fstatev;

   private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompRef;
	
	
};

}

#endif /* _RG_VISCO_T_H_ */

