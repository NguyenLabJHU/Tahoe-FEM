/* $Id: ConveyorT.h,v 1.1 2003-08-06 06:54:15 paklein Exp $ */
#ifndef _CONVEYOR_T_H_
#define _CONVEYOR_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "AutoArrayT.h"
#include "iArray2DT.h"

namespace Tahoe {

/** forward declarations */
class FieldT;
class KBC_PrescribedT;

/** conveyor belt */
class ConveyorT: public KBC_ControllerT
{
public:

	/** constructor */
	ConveyorT(NodeManagerT& node_manager, FieldT& field);

	/** initialization */
	virtual void Initialize(ifstreamT& in);

	/** write parameters */
	virtual void WriteParameters(ostream& out) const;

	/** not implemented - there's no going back */
	virtual void Reset(void);

	/** set to initial conditions */
	virtual void InitialCondition(void);

	/** open time interva; */
	virtual void InitStep(void);

	/** computing residual force */
	virtual void FormRHS(void);

	/** apply the update to the solution. Does nothing by default. */
	virtual void Update(const dArrayT& update);

	/** signal that the solution has been found */
	virtual void CloseStep(void);

	/* returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

protected:

	/** enum to define tracking criterion */
	enum TrackingTypeT {
        kMax = 1,
        kMin = 2,
   kLeftMost = 3,
  kRightMost = 4
	};

	/** locate new tracking point */
	double TrackPoint(TrackingTypeT tracking_type, double threshold);

	/** reset system to new center
	 * \return true if the system focus has been changed */
	bool SetSystemFocus(double focus);

	/** mark elements linking left to right edge as inactive */
	void MarkElements(void);

	/** deactivate elements to create a pre-crack */
	void CreatePrecrack(void);
	
protected:

	/** the field */
	FieldT& fField;

	/** \name prescribed dimensions */
	/*@{*/
	/** repeating length of the mesh in the x-direction */
	double fMeshRepeatLength; 

	/** distance to shift the window to enforce ConveyorT::fRightMinSpacing */
	double fWindowShiftDistance;

	/** minimum distance between tracking point and the right edge of the domain.
	 * If the tracking point enters this zone, the window over the system is shifted
	 * by ConveyorT::fShiftDistance. */
	double fRightMinSpacing;
	/*@}*/

	/** \name stretching boundary
	 * Stretching as a displacement, velocity, or acceleration jump 
	 * split equally between the upper and lower surfaces. */
	/*@{*/
	KBC_CardT::CodeT fULBC_Code;
	double           fULBC_Value;
	int              fULBC_ScheduleNumber;
	const ScheduleT* fULBC_Schedule;
	/*@}*/

	/** \name boundary condition for the far right edge */
	/*@{*/
	KBC_PrescribedT* fRightEdge;
	AutoArrayT<int>  fShiftedNodes;
	/*@}*/

	/** \name crack tip element group */
	/*@{*/
	int fTipElementGroup; /**< number of the element group controlling the tip position */
	double fTipX_0; /**< initial x-coordinate of the crack tip */
	double fTipY_0; /**< cleavage plane position */
	
	int fTipOutputCode; /**< output flag to generate data to locate the tip */
	int fTipColumnNum;  /**< column of output variable to locate tip */
	double fTipThreshold; /**< threshold value which defines the tip position */
	/*@}*/

	/** \name edge damping */
	/*@{*/
	double fDampingWidth;
	double fDampingCoefficient;
	AutoArrayT<int> fDampingNodes;
	dArray2DT fDampingForce;
	dArray2DT fDampingCoeff;	
	iArray2DT fDampingEqnos;
	bool fDampingReset;
	/*@}*/

	/** \name nodes at upper and lower boundaries */
	/*@{*/
	iArrayT fBottomNodes;
	iArrayT fTopNodes;
	/*@}*/

	/** \name derived dimensions */
	/*@{*/
	double fX_Left;  /**< left-most edge of the undeformed mesh */
	double fX_Right; /**< right-most edge of the undeformed mesh */
	double fX_PeriodicLength; /** periodic length of the system */

	/** width of the dead element zone at either end of the domain. */
	double fWidthDeadZone;
	/*@}*/

	/** \name tracking point */
	/*@{*/
	int fTrackingInterval;

	int fTrackingCount;
	double fTrackingPoint;
	double fTrackingPoint_last;
	
	/** output stream for tracking point position at all fTrackingInterval */
	ofstreamT fTrackingOutput;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _CONVEYOR_T_H_ */
