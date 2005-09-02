/* $Id: BridgingScaleManagerT.h,v 1.7 2005-09-02 02:07:50 d-farrell2 Exp $ */
#ifndef _BRIDGING_SCALE_MANAGER_H_
#define _BRIDGING_SCALE_MANAGER_H_

/* element configuration header */
#include "ElementsConfig.h"
#include "DevelopmentElementsConfig.h"
#include "dSPMatrixT.h"
#include "nMatrixT.h"

#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)

/* base class  */
#include "MultiManagerT.h"

namespace Tahoe {

/* forward declarations */
class FEManagerT_THK;

/** manager for dynamic bridging scale calculations */
class BridgingScaleManagerT: public MultiManagerT
{
public:

	/** constructor */
	BridgingScaleManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, TaskT task);

	/** destructor */
	virtual ~BridgingScaleManagerT(void);

	/** solve all the time sequences */
	virtual void Solve(void);

	/** (re-)set system to initial conditions */
	virtual ExceptionT::CodeT InitialCondition(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:
	
	/** Determine basic solution information */
	virtual void Initialize(void);
	
	/** Determine basic BSM solution information, initial conditions */
	virtual void InitBSM(void);
	
	/** Determine basic MD/THK solution information, initial conditions */
	virtual void InitMDTHK(void);

	/** Determine BSM solution */
	virtual void SolveBSM(void);
	
	/** Determine MD/THK solution */
	virtual void SolveMDTHK(void);

	/** calculate internal+external force for the given displacement u */
	const dArray2DT& TotalForce(const StringT& field_name, const dArray2DT& field_values, 
		FEManagerT_bridging& bridging, dArray2DT& rhs_2D) const;

private:

	/** cast of MultiManagerT::fFine to FEManagerT_THK */
	FEManagerT_THK* fFine_THK;
	
	// Misc. parameters for both BSM & MD/THK
	int fNSD, fNND, fNumgatoms, fNumbatoms;
	dArray2DT fGadisp, fGavel, fGaacc;	// Ghost atom kinematic info
	dArray2DT fBadisp, fBavel, fBaacc;	// Boundary atom kinematic info
	dArray2DT fBoundghostdisp, fBoundghostvel, fBoundghostacc, fTHKforce; // boundary & ghost atom information and THK force
	iArrayT fBoundaryghostatoms, fAllatoms, fGatoms, fBatoms, fBoundatoms; // arrays of numbers of boundary and ghost atoms
	TimeManagerT* fFine_time_manager;
	
	// Misc. parameters for BSM only
	dArray2DT fRHS_2D_true, fFubig, fTotalu, fFu, fProjectedu, fNtfproduct;
	iArrayT fActiveFENodes;
	dSPMatrixT fNtF;
	CommManagerT* fFine_comm_manager;
	int fFubig_ID;
	TimeManagerT* fCoarse_time_manager;
	dArrayT fFx, fTempx;	// some arrays for the NtF calculation (perhaps don't need?)
	nMatrixT<int> fGhostonmap, fGhostoffmap;
		
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
#endif /* _BRIDGING_SCALE_MANAGER_H_ */
