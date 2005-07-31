/* $Id: DiffusionElementT.h,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (10/02/1999)                                          */

#ifndef _DIFFUSE_T_H_
#define _DIFFUSE_T_H_

/* base class */
#include "ContinuumElementT.h"

/* direct members */
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ShapeFunctionT;
class DiffusionMaterialT;
class StringT;

class DiffusionElementT: public ContinuumElementT
{
public:
	
	enum OutputCodeT {iNodalCoord = 0, // (reference) nodal coordinates
                       iNodalDisp = 1, // nodal "displacements"
                    iMaterialData = 2};// material model output

	/* constructor */
	DiffusionElementT(FEManagerT& fe_manager);
	
	/* data initialization */
	virtual void Initialize(void);

	/* set the controller */
	virtual void SetController(eControllerT* controller);

	/* compute nodal force */
	virtual void AddNodalForce(int node, dArrayT& force);

	/* returns the stored energy */
	virtual double InternalEnergy(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out);

	/* initialization functions */
	virtual void SetLocalArrays(void);
	virtual void SetShape(void);

	/* construct the effective mass matrix */
	virtual void LHSDriver(void);

	/* form the residual force vector */
	virtual void RHSDriver(void);

	/* increment current element */
	virtual bool NextElement(void);	
	
	/* form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/* calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/* return a pointer to a new material list */
	virtual MaterialListT* NewMaterialList(int size) const;

	/* driver for nodal value calculations */
	virtual void ComputeNodalValues(const iArrayT& codes);

private:

	/* construct output labels array */
	virtual void SetOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts);
	virtual void GenerateOutputLabels(const iArrayT& counts,
		ArrayT<StringT>& labels) const;

protected:

	/* run time */
	DiffusionMaterialT* fCurrMaterial;

	/* arrays with local ordering */
	LocalArrayT fLocVel; // "temperature rate"
	
	/* work space */
	dMatrixT fD; /* constitutive matrix          */
	dMatrixT fB; /* "strain-displacement" matrix */
	dArrayT  fq; /* heat flow = k_ij T,j         */

	/* parameters */
	static const int NumOutputCodes;
};

#endif /* _DIFFUSE_T_H_ */
