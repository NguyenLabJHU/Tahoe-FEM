/* $Id: AtomicConveyorT.h,v 1.1.2.1 2003-09-18 21:03:44 cjkimme Exp $ */
#ifndef _ATOMIC_CONVEYOR_T_H_
#define _ATOMIC_CONVEYOR_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "ScheduleT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "BasicFieldT.h"

namespace Tahoe {

/** Atomistic conveyor belt for steady-state crack propagation.
 * Crack propagates from left to right, and if crack tip gets too
 * close to the right boundary, material is cut from left side
 * and pasted onto right side.
 *
 * A number of assumptions/restrictions are placed on the form of
 * the crystal that will work with this controller. A few of them
 * are:
 *
 *    The crystal sample is presumed to be a parallelpipedal strip.
 *    In N = 2 or 3 dimensions, the Nth cartesian coordinate axis 
 *       specifies the direction of equal and opposite tensile pull
 *       applied to the sample faces normal to this direction.
 *    The initial crack is along the midpoint of the sample and
 *       propogates along the (N - 1) cartesian coordinate axis.
 *    There is only one kind of interatomic interaction. 
 */
class AtomicConveyorT: public KBC_ControllerT
{
public:	

	/** constructor */
	AtomicConveyorT(NodeManagerT& node_manager, BasicFieldT& field);

	/** destructor */
 	~AtomicConveyorT(void);

	/** initialize data. Must be called immediately after construction */
	virtual void Initialize(ifstreamT& in);
	
	virtual void InitStep(void);

	/** generate velocities near the crack tip */
	virtual void InitialCondition(void);
	
	/** write parameters to .out file */
	virtual void WriteParameters(ostream& out) const;

protected:
	
	/** enforce boundary conditions/damping */
	void SetBCCards(void);
	
	/** locate nodes in specfied region */
	void NodesInRegion(dArrayT& xmin, dArrayT& xmax, iArrayT& nodes);
	/*@}*/

protected:

	/** the Field */
	BasicFieldT& fField;

	/** dimensionality of space */
	int fSD;
	
	/** dimensions of system */
	dArrayT fBoxSize;

	/** data schedule for the top/bottom displacement loading */
	const ScheduleT* fSchedule;
	double fScale;
	
	/** nodes on boundary providing the loading */
	iArrayT fTopNodes, fBottomNodes;
	
	/** regions containing above nodes if not specified by nodeset ID */
	dArrayT fTopMin, fTopMax, fBottomMin, fBottomMax;
	
	/** nodes on edges that are damped */
	iArrayT iDampedNodesLeft, iDampedNodesRight;
	
	/** initial and current length of crack tip */
	double fPreCrackLength, fCrackLength;
	
	/** height of crack face */
	double fCrackHeight;
	
	/** height tolerance for atoms on crack face of order of 
	 *  the interatomic potential range */
	double fTipRegionHeight;
	
	/** find the crack tip every this many timesteps */
	int nTipIncs;
	
	/** opening displacement and tolerance for locating crack tip */
	double fMaxOpening, fOpeningTol;
	
	/** paste on new material whenever tip is within this distance of end */
	double fMinTipApproach; 

	/** current boundaries of system */
	dArrayT fBoxMin, fBoxMax;
	
	/** reference boundaries of system */
	dArrayT fBoxMin_0, fBoxMax_0;
	
	/** damp atoms within this distance of right or left */
	double fDelta_x;
	
	/** maximum damping constant--damping is interpolated between 0 and this value */
	double fMaxBeta;
	
	/** Initial strain */
	double initStrain;

	/** Element group of the crystal */
	int iElementGroup;
};

} // namespace Tahoe 
#endif /* _ATOMIC_CONVEYOR_T_H_ */
