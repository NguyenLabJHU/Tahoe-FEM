/* $Id: CSEAnisoT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (11/19/1997)                                          */
/* Cohesive surface elements with vector traction potentials,             */
/* i.e., like Xu-Needleman's potential.                                   */

#ifndef _CSE_ANISO_T_H_
#define _CSE_ANISO_T_H_

/* base class */
#include "CSEBaseT.h"

/* direct members */
#include "pArrayT.h"

/* forward declarations */
class SurfacePotentialT;

class CSEAnisoT: public CSEBaseT
{
public:

	/* constructor */
	CSEAnisoT(FEManagerT& fe_manager);

	/* destructor */
	~CSEAnisoT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialize class data */
	virtual void Initialize(void);

protected:

	/* tangent matrix and force vector */
	virtual void LHSDriver(void);
	virtual void RHSDriver(void);

	/* driver for nodal value calculations */
	virtual void ComputeNodalValues(const iArrayT& codes);

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

private:

	/* nodal value calculations */
	virtual void SetOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;

	/* construct output labels array */
	virtual void GenerateOutputLabels(const iArrayT& codes,
		ArrayT<StringT>& labels) const;

	/* operations with pseudo rank 3 (list in j) matrices */
	void u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
		dMatrixT& Qu);

	void Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
		dMatrixT& Qu);

protected:

	/* shape functions (wrt current config) */
	SurfaceShapeT* fCurrShapes;

	/* cohesive surface potentials */
	pArrayT<SurfacePotentialT*> fSurfPots;
	
	/* coordinate transformation */
	dMatrixT fQ;     // t'_i = Q_ji t_j, where t' is in the local frame
	dArrayT  fdelta; // gap vector (local frame)
	dArrayT  fT;     // traction vector (global frame)
	dMatrixT fddU;	 // surface stiffness (local frame)
	
	ArrayT<dMatrixT> fdQ; // list representation of rank 3 of dQ_ij/du_k
	
	/* work space (for tangent) */
	dMatrixT fnsd_nee_1;
	dMatrixT fnsd_nee_2;
};

#endif /* _CSE_ANISO_T_H_ */
