/* $Id: MultiManagerT.h,v 1.5 2004-07-15 08:31:03 paklein Exp $ */
#ifndef _MULTI_MANAGER_H_
#define _MULTI_MANAGER_H_

/* element configuration header */
#include "ElementsConfig.h"
#ifdef BRIDGING_ELEMENT

/* base class  */
#include "FEManagerT.h"

namespace Tahoe {

class FEManagerT_bridging;

/** wrapper for running multiple FEManagerT with a single solver */
class MultiManagerT: public FEManagerT
{
public:

	/** constructor */
	MultiManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
		FEManagerT_bridging* fine, FEManagerT_bridging* coarse);

	/** destructor */
	virtual ~MultiManagerT(void);
	
	/** initialize members */
	virtual void Initialize(InitCodeT init = kFull);

	/** (re-)set the equation number for the given group */
	virtual void SetEquationSystem(int group, int start_eq_shift = 0);

	/** \name solution steps */
	/*@{*/
	/** initialize the current time increment for all groups */
	virtual ExceptionT::CodeT InitStep(void);

	/** close the current time increment for all groups */
	virtual ExceptionT::CodeT CloseStep(void);

	/** called for all groups if the solution procedure for any group fails */
	virtual ExceptionT::CodeT ResetStep(void);
	/*@}*/

	/** \name solution messaging */
	/*@{*/
	/** compute LHS-side matrix and assemble to solver */
	virtual void FormLHS(int group, GlobalT::SystemTypeT sys_type) const;

	/** compute RHS-side */
	virtual void FormRHS(int group) const;

	/** send update of the solution to the NodeManagerT */
	virtual void Update(int group, const dArrayT& update);

	/** system relaxation */
	virtual GlobalT::RelaxCodeT RelaxSystem(int group) const;
	/*@}*/

	/** \name output */
	/*@{*/
	/** initiate the process of writing output from all output sets 
	 * \param time time label associated with the output data */
	virtual void WriteOutput(double time);

	/** (temporarily) direct output away from main out */
	virtual void DivertOutput(const StringT& outfile);

	/** restore outputs to their regular destinations */
	virtual void RestoreOutput(void);
	/*@}*/

	/** \name load control functions (returns true if successful) */
	/*@{*/
	virtual bool DecreaseLoadStep(void);
	virtual bool IncreaseLoadStep(void);
	/*@}*/

private:

	/** \name sub-managers */
	/*@{*/
	FEManagerT_bridging* fFine;
	FEManagerT_bridging* fCoarse;	
	/*@}*/
	
	/** \name equations for each sub-manager */
	/*@{*/
	iArrayT fEqnos1;
	iArrayT fEqnos2;
	/*@}*/

	/** work space */
	dArray2DT fFieldAtGhosts;
	
	/** \name coarse/fine output */
	/*@{*/ 
	iArray2DT fAtomConnectivities;
	int fOutputID;
	bool fDivertOutput;
	/*@}*/
	
	/** \name workspace for cross terms */
	/*@{*/
	/** transpose of data for interpolating data from the coarse scale
	 * onto the fine scale points */
	//InterpolationDataT fFollowerCellTranspose;

	dArray2DT fR_U; /**< coarse scale forces */
	iArrayT   fR_U_eqnos; /**< equations for assembly for MultiManagerT::fR_U */
	dArray2DT fR_Q; /**< fine scale forces */
	iArrayT   fR_Q_eqnos; /**< equations for assembly for MultiManagerT::fR_Q */
	
	const FieldT* fFineField;
	const FieldT* fCoarseField;
	/*@}*/
	
	/** \name keep/omit cross terms */
	/*@{*/
	bool fFineToCoarse; /**< fine scale contribution to coarse scale equations */ 
	bool fCoarseToFine; /**< coarse scale contribution to fine scale equations */ 
	int  fCorrectOverlap; /**< adjust C-B bond densities to account for overlap */ 
	double fCBTikhonov; /**< regularization used to solve C-B bond densities in the overlap */
	double fK2; /**< penalization constant for density different from 1.0 */
	/*@}*/
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _MULTI_MANAGER_H_ */
