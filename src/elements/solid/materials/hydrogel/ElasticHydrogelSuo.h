/* $Id: ElasticHydrogelSuo.h,v 1.1 2013-02-01 17:03:24 tahoe.xiaorui Exp $ */
/* created : RX (1/5/2012) */
#ifndef _ElasticHydrogelSuo
#define _ElasticHydrogelSuo

/* base classes */
#include "FSSolidMatT.h"
#include "SpectralDecompT.h"
#include "PotentialT.h"
#include "C1FunctionT.h"
#include <cmath>
#include "RGSplitT2.h"


namespace Tahoe {

/** base class for nonlinear finite deformation viscoelasticity **/

class ElasticHydrogelSuo: public FSSolidMatT 
{
  public:
  
	/* constructor */
	ElasticHydrogelSuo(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij.  */
/*	virtual double Pressure(void); */
	
	/** return true if the material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/*Initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; 
	virtual void PointInitialize(void);              

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time

	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
	virtual const dMatrixT& C_IJKL(void);
	virtual const dSymMatrixT& S_IJ(void);
	
	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*virtual void Load(ElementCardT& element, int ip)
	
	/*@}*/
	//add these two
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	virtual void DefineParameters(ParameterListT& list) const;

	
protected:
	/*virtual double ComputeSwelling(const dArrayT& eigenstretch); */   
	void CalculatePhi(double& phi,  const double& phi_n, const double& J);
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
	virtual void Load(ElementCardT& element, int ip);
	virtual void Store(ElementCardT& element, int ip);	

 	
protected:
    double fchi;
	/* return values */
	double fT;
	dMatrixT fModulus;
	dSymMatrixT fStress;
	
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;
	
	/*state variable to store solid fraction*/
	int fnstatev;
	dArrayT fstatev;    /*internal state variable array*/
	double* fSolidFraction;
	double* fSolidFraction_n;
	
	/*work space*/
	dSymMatrixT fb;
	dSymMatrixT fbe;
	dMatrixT fF3D;
	
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_dev;
	
	dArrayT	    ftau;
	dSymMatrixT fDtauDe;

	dSymMatrixT fStress3D;
	dMatrixT    fModulus3D;
	dMatrixT    fModMat;
	
	/*potential*/
	PotentialT* fPot;

	
};

}

#endif /* _RG_VISCO_T_H_ */

