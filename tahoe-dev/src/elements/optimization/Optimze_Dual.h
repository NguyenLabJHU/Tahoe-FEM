/* $Id: Optimze_Dual.h,v 1.1 2009-04-23 03:03:43 thao Exp $ */
/*Class to calculate objective function and gradient for inverse elasticity problems using the adjoint method */
/*A Oberai, NH Gakhale, GR Feijoo (2003) Inverse Problems 19:297-313*/

/*currently parameters are read in and acted on at the element level only.*/
/*assumes a homogeneous distribution*/

#ifndef _Optimize_Dual_
#define _Optimize_Dual_

/* base class */
#include "ParameterInterfaceT.h"

#include "ofstreamT.h"
#include "ifstreamT.h"

namespace Tahoe {

/** Interface for linear strain deformation and field gradients */
class Optimize_Dual: virtual public ParameterInterfaceT
{
  public:
      
	/** constructor */
	Optimize_Dual(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

  protected:
    
	virtual void InitParams(void);
	
	virtual void ReadData(void);
		
   /*write Dakota output file*/
	void Write_Dakota_Output(ofstreamT& output);

   /*returns the number of parameters*/
	const int NumParams(void) const {return fParamVals.Length();}
	
	/*return values both const and non const*/
	const dArrayT& ParamVals(void) const {return fParamVals;};


	const ArrayT<StringT>& ParamLabels(void) const {return fParamLabels;};
	
	/** calculate the internal force contribution*/
	virtual void FormKd(double constK);

	/** calculate the internal force contribution*/
	virtual void FormStiffness(double constK);

	/** calculate objective function*/
	virtual void Compute_Cost(void);
	
	/** calculate objective function*/
	virtual void Compute_Gradients(void);

	/**calculate gradients*/
	
	virtual void SetLocalArrays(void);
	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);
	
	bool NextElement(void);
  protected:
	SSOptimize_MatSupportT* fSSOptimize_MatSupport;
	
	SSOptimize_MatT* fOptimize_Mat; 

	/*The primal problem*/
	SmallStrainT* fPrimal_Element;
	
	/** \name return values */
	/*@{*/
  	ArrayT<dSymMatrixT> fDual_Strain_List;
//  	ArrayT<dSymMatrixT> fDual_Strain_last_List;
	/*@}*/

	LocalArrayT fLocPrimalDisp;      /**< nodal primal disp field */
	LocalArrayT fLocPrimalDisp_last; /**< last  primal disp field */
	LocalArrayT fLocData;      /**< measurement data field */
		
	/*regularization?*/
	double fCostFunction;
	double fconstraint;
	dArrayT fGradients;  /*gradients of the cost function*/
	dArrayT fConstraint_Grad;
	
	dArrayT fParamVals;		/*parameter values*/
	ArrayT<StringT> fParamLabels; /*parameter labels*/
	
	private:
	/*work spaces*/
	dArrayT fip_residual;
	dArrayT fip_data;
	
	StringT fOutFile;
	ofstreamT fOutput;
	
	StringT fDataFileRoot;
	StringT fDataFile;
	ifstreamT fDataInput;
	
	dArray2DT fData;
	dArray2DT fData_Coords;
	
	double fweight_cost;

};


#endif /* _Optimize_Dual_ */
