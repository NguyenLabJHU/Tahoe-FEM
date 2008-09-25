/* $Id: NLV_Nfibers.h,v 1.1 2008-09-25 13:32:32 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _NLV_Nfibers_
#define _NLV_Nfibers_ 
 
/* base class */
#include "FSFiberMatViscT.h"
#include "C1FunctionT.h"
#include "Array2DT.h"

namespace Tahoe {

class NLV_Nfibers: public FSFiberMatViscT
{
   public:
  
	/* constructor/destructor */
	NLV_Nfibers(void);
	
/* destructor */
	~NLV_Nfibers(void);
	
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** describe the parameters */
	virtual void DefineParameters(ParameterListT& list) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** types of analysis */
	enum FiberPot {
	         kStandard = 0,      // W(I) = k(I-1)^2
		      kExponential = 1,  // W(I) = k1/k2 (exp(k2(I-1)) - k2 I) veronda westmann
		     kExpQuad = 2        // W(I) = k1/2k2 (exp(k2(I-1)^2) - 2 k2 I) holzapfel gasser
		   };
		

	/****************************** fiber ***********************************************/
	/*computes fiber stress in local frame*/
	/*redefined from base class*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress, 
				const int process_index);
	
	/*computes  moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress, dMatrixT& Mod, 
				const int process_index);

	/*compute the algorithmic moduli dSNEQ/dCv deltaCv/deltadC in local fiber coord sys*/
	virtual void ComputeCalg (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v,  dMatrixT& Calg, const int process_index); 
	
	/*local newton loop for viscous stretch tensor*/ 
	virtual void Compute_Cv(const dSymMatrixT& FiberStretch, const dSymMatrixT& FiberStretch_vn, 
			dSymMatrixT& FiberStretch_v, const int process_index);
	
protected:
	

 	/* work space needed to compute fiber viscoelastic behavior*/
	/*dimension fNumFibProcess + 1*/
//	int fPotType_f;
	Array2DT<C1FunctionT*> fPot_f;
	Array2DT<C1FunctionT*> fVisc_f;

	bool fsame;
	/*workspaces*/
	dArrayT fResidual;

	dArrayT fVec;
	dSymMatrixT fFlowStress;
	dSymMatrixT fMod1;
	dSymMatrixT fMod2;
	dMatrixT fiVisc;
	dMatrixT fiK_f;
	dMatrixT fG;
	
};
	
}
#endif /* _NLV_Nfibers_ */
