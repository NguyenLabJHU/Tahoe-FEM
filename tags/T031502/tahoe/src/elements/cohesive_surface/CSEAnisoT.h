/* $Id: CSEAnisoT.h,v 1.9 2002-02-18 19:09:42 cjkimme Exp $ */
/* created: paklein (11/19/1997) */

#ifndef _CSE_ANISO_T_H_
#define _CSE_ANISO_T_H_

/* base class */
#include "CSEBaseT.h"

/* direct members */
#include "pArrayT.h"
#include "RaggedArray2DT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

/* forward declarations */
class SurfacePotentialT;

/* Cohesive surface elements with vector argument cohesive relations. */
class CSEAnisoT: public CSEBaseT
{
public:

	/* constructor */
	CSEAnisoT(FEManagerT& fe_manager, bool rotate);

	/* destructor */
	~CSEAnisoT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** initialize class data */
	virtual void Initialize(void);

	/** close current time increment */
	virtual void CloseStep(void);

	/** write restart data to the output stream. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. */
	virtual void ReadRestart(istream& in);

protected:

	/* tangent matrix and force vector */
	virtual void LHSDriver(void);
	virtual void RHSDriver(void);

	/* nodal value calculations */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
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
	iArrayT fNumStateVariables;
	pArrayT<SurfacePotentialT*> fSurfPots;

	/** state variable storage array. 
	 * Array has dimensions: [nel] x [nip * nvar] */
	RaggedArray2DT<double> fStateVariables;
	RaggedArray2DT<double> fStateVariables_last;
	
	/* coordinate transformation */
	dMatrixT fQ;     // t'_i = Q_ji t_j, where t' is in the local frame
	dArrayT  fdelta; // gap vector (local frame)
	dArrayT  fT;     // traction vector (global frame)
	dMatrixT fddU;	 // surface stiffness (local frame)
	
	ArrayT<dMatrixT> fdQ; // list representation of rank 3 of dQ_ij/du_k
	
	/* work space (for tangent) */
	dMatrixT fnsd_nee_1;
	dMatrixT fnsd_nee_2;

	/* variables for calculating nodal info */
	/* Added by cjkimme 11/07/01 */
	bool fCalcNodalInfo;
	int fNodalInfoCode;
        dArray2DT fNodalQuantities;
	LocalArrayT fNodalValues;
	int fBulkGroup;

	/* change the time Step */
	/*double fNewTimeStep;
	bool fNeedsNewTimeStep;
	bool NewTimeStepQ(void);
	double NewTimeStep(void);*/
};

#endif /* _CSE_ANISO_T_H_ */
