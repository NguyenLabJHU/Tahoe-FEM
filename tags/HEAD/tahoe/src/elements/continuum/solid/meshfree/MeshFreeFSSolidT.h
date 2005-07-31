/* $Id: MeshFreeFSSolidT.h,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (09/16/1998)                                          */
/* large deformation elasticity with MLS shapefunctions for the           */
/* field (displacement) representation                                    */
/* NOTE: clean up code governing when crack growth algorithm              */
/* is used, initiation criteria, etc. (PAK 09/28/1999)                    */

#ifndef _EFG_FDELASTIC_T_H_
#define _EFG_FDELASTIC_T_H_

/* base classes */
#include "TotLag_FSSolidT.h"
#include "MeshFreeFractureSupportT.h"

/* direct members */
#include "nVariMatrixT.h"

class MeshFreeFSSolidT: public TotLag_FSSolidT,
	public MeshFreeFractureSupportT
{
public:

	/* constructor */
	MeshFreeFSSolidT(FEManagerT& fe_manager);
	
	/* data initialization */
	virtual void Initialize(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/* appends group connectivities to the array */
virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/* write output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/* returns true if the internal force has been changed since
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
	virtual void ResetStep(void); // restore last converged state
					
protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;

	/* initialization functions */
	virtual void SetShape(void);

	/* increment current element */
	virtual bool NextElement(void);

	/* driver for nodal value calculations */
	virtual void ComputeNodalValues(const iArrayT& codes);

private:

	/* write displacement field and gradients */
	 virtual void WriteField(void); //TEMP?
	
protected:

	/* wrappers */
	nVariMatrixT<double>  fStressStiff_wrap;
	nVariMatrixT<double>  fB_wrap;
	nVariMatrixT<double>  fGradNa_wrap;
	nVariArray2DT<double> fDNa_x_wrap;
};

#endif /* _EFG_FDELASTIC_T_H_ */
