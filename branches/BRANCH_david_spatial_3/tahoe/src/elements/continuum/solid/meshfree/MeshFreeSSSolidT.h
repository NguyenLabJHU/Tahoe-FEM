/* $Id: MeshFreeSSSolidT.h,v 1.9.30.1 2005-07-25 02:36:56 paklein Exp $ */
/* created: paklein (09/11/1998) */
#ifndef _MF_SMALLSTRAIN_T_H_
#define _MF_SMALLSTRAIN_T_H_

/* base classes */
#include "SmallStrainT.h"

/* direct members */
#include "nVariMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
class MeshFreeShapeFunctionT;
class MeshFreeFractureSupportT;

/** small strain elasticity with MLS shapefunctions for the
 * field (displacement) representation
 * \note clean up code governing when crack growth algorithm
 * is used, initiation criteria, etc. (PAK 09/28/1999) */
class MeshFreeSSSolidT: public SmallStrainT
{
public:

	/** constructor */
	MeshFreeSSSolidT(const ElementSupportT& support);

	/** destructor */
	~MeshFreeSSSolidT(void);
	
	/** append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** appends group connectivities to the array */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/** write output */
	virtual void WriteOutput(void);

	/** returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* returns 1 if DOF's are interpolants of the nodal values */
	virtual int InterpolantDOFs(void) const;

	/* retrieve nodal unknowns */
	virtual void NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const;

	/* weight the computational effort of every node */
	virtual void WeightNodalCost(iArrayT& weight) const;

	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state

	/** accessors */
	MeshFreeSupportT& MeshFreeSupport(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* initialization functions */
	virtual void SetShape(void);

	/* increment current element */
	virtual bool NextElement(void);

	/* driver for calculating output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

private:

	/* write displacement field and gradients */
	 virtual void WriteField(void); //TEMP?
	
private:

	/** meshless shape functions */
	MeshFreeShapeFunctionT* fMFShapes;

	/** support for meshless calculations */
	MeshFreeFractureSupportT* fMFFractureSupport;

	/** make field at bounding nodes nodally exact */
	bool fAutoBorder;

	/* dynamic wrapper */
	nVariMatrixT<double> fB_wrap;
	
	/* field ready flag */
	bool fFieldSet;

	/** pointer to list parameters needed to construct meshless shape functions. This
	 * pointer is set during MeshFreeSSSolidT::TakeParamaterListT and used during
	 * MeshFreeSSSolidT::SetShape */
	const ParameterListT* fMeshfreeParameters;
};

} // namespace Tahoe 
#endif /* _MF_SMALLSTRAIN_T_H_ */
