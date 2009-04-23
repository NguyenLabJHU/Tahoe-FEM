/* $Id: SS_Optimize_Dual.h,v 1.1 2009-04-23 03:03:44 thao Exp $ */
/*Class to calculate objective function and gradient for inverse elasticity problems using the adjoint method */
/*A Oberai, NH Gakhale, GR Feijoo (2003) Inverse Problems 19:297-313*/

/*currently parameters are read in and acted on at the element level only.*/
/*assumes a homogeneous distribution*/

#ifndef _SS_Optimize_Dual_
#define _SS_Optimize_Dual_

/* base class */
#include "SmallStrainT.h"
#include "SSOptimize_MatT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"


namespace Tahoe {

/** Interface for linear strain deformation and field gradients */
class SS_Optimize_Dual: public SmallStrainT
{
  public:
      
	/** constructor */
	SS_Optimize_Dual(const ElementSupportT& support);

	/*destructor*/
	~SS_Optimize_Dual(void);

	/** destructor */
//	~SS_Optimize_Dual(void);

	virtual void InitialCondition(void);
	
	virtual void InitStep(void);
	
	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearDualStrain(void) const;
	const dSymMatrixT& LinearDualStrain(int ip) const;
	/*@}*/

	virtual void WriteOutput(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

  protected:
  
   /*write Dakota output file*/
   void Write_Dakota_Output(ofstreamT& output);
  
   /*returns the number of parameters*/
	const int NumParams(void) const {return fParamVals.Length();}
	
	/*return values both const and non const*/
	const dArrayT& ParamVals(void) const {return fParamVals;};


	const ArrayT<StringT>& ParamLabels(void) const {return fParamLabels;};
	
	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param name list identifier
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

	/** calculate the internal force contribution*/
	virtual void FormKd(double constK);

	void FormStiffness(double constK);

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
	dSymMatrixT fmat;

};

inline const dSymMatrixT& SS_Optimize_Dual::LinearDualStrain(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
		ExceptionT::GeneralFail("SmallStrainT::LinearStrain", "material %d did not specify this need", 
			mat_num + 1);
#endif

	return fDual_Strain_List[CurrIP()];
}

inline const dSymMatrixT& SS_Optimize_Dual::LinearDualStrain(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
		ExceptionT::GeneralFail("SmallStrainT::LinearStrain", "material %d did not specify this need", 
			mat_num + 1);
#endif
	return fDual_Strain_List[ip];
}

} /* namespace Tahoe */

#endif /* _SMALLSTRAIN_T_H_ */
