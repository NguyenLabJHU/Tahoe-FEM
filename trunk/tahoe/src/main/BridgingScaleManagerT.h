/* $Id: BridgingScaleManagerT.h,v 1.3 2007-08-14 16:22:43 d-farrell2 Exp $ */
#ifndef _BRIDGING_SCALE_MANAGER_H_
#define _BRIDGING_SCALE_MANAGER_H_

/* element configuration header */
#include "ElementsConfig.h"
#include "dSPMatrixT.h"
#include "nMatrixT.h"

#if defined(BRIDGING_ELEMENT)

/* base class  */
#include "MultiManagerT.h"

/* File related to calculation using the Wagner-Karpov-Liu Bridging Scale Method
* If you make use of this code, please cite the following publications (they are also handy references for the method and implementation):
* 
*		1)	Wagner GJ, WK Liu. "Coupling of Atomistic and Continuum Simulations using a Bridging Scale Decomposition", Journal of Computational Physics, 190:249-274 (2003)
*		2)	Wagner GJ, EG Karpov, WK Liu. "Molecular Dynamics Boundary Conditions for Regular Crystal Lattices", CMAME, 193(17-20):1579-1601 (2004)
*		3)	Park HS, WK Liu. "An Introduction and Tutorial on Multiple Scale Analysis in Solids", CMAME, 193(17-20):1733-1772 (2004)
*		4)	Park HS, EG Karpov, PA Klein, WK Liu. "Three-Dimensional Bridging Scale Analysis of Dynamic Fracture", Journal of Computational Physics, 207(2):588-609 (2005)
*		5)	Park HS, EG Karpov, WK Liu, PA Klein. "The Bridging Scale for Two-Dimensional Atomistic/Continuum Coupling", Philosophical Magazine, 85(1):79-113 (2005)
*		6)	Farrell DE, HS Park, WK Liu. "Implementation Aspects of the Bridging Scale Method and Application to Intersonic Crack Propagation", IJNME 71:583-605 (2007)
*		7)	Farrell DE, EG Karpov, WK Liu. "Algorithms for Bridging Scale Method Parameters", Computational Mechanics, in print DOI: 10.1007/s00466-007-0156-z (2007)
*/

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
	
	/** Determine basic BSM solution information, initial conditions for Beta type THK */
	virtual void InitBetaBSM(void);
	
	/** Determine basic MD/THK solution information, initial conditions */
	virtual void InitMDTHK(void);

	/** Determine BSM solution, theta THK */
	virtual void SolveBSM(void);
	
	/** Determine BSM solution, beta THK */
	virtual void SolveBetaBSM(void);
	
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
	dArray2DT fRHS_2D_true, fFubig, fTotalu, fTotalv, fFu, fProjectedu, fNtfproduct;
	iArrayT fActiveFENodes;
	dSPMatrixT fNtF;
	CommManagerT* fFine_comm_manager;
	int fFubig_ID;
	TimeManagerT* fCoarse_time_manager;
	dArrayT fFx, fTempx;	// some arrays for the NtF calculation (perhaps don't need?)
	nMatrixT<int> fGhostonmap, fGhostoffmap;
	
	
	dArray2DT fFineU, fFineV;
		
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _BRIDGING_SCALE_MANAGER_H_ */
