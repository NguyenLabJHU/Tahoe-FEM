/* $Id: CSEAnisoT.h,v 1.3 2001-02-27 00:06:24 paklein Exp $ */
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
	CSEAnisoT(FEManagerT& fe_manager, bool rotate);

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

	/* nodal value calculations */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;

	/* compute output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values);

	/* construct output labels array */
	virtual void GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const;

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

private:

	/* operations with pseudo rank 3 (list in j) matrices */
	void u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
		dMatrixT& Qu);

	void Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
		dMatrixT& Qu);

protected:

	/* shape function - if (fRotate) then wrt current config */
	bool fRotate;
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
