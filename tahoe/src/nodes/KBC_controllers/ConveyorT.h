/* $Id: ConveyorT.h,v 1.6 2004-12-21 17:24:51 thao Exp $ */
#ifndef _CONVEYOR_T_H_
#define _CONVEYOR_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "AutoArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "iArrayT.h"
#include "ofstreamT.h"

namespace Tahoe {

/** forward declarations */
class FieldT;

/** conveyor belt */
class ConveyorT: public KBC_ControllerT
{
public:

	/** constructor */
	ConveyorT(const BasicSupportT& support, FieldT& field);

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

	/** \name restart functions */
	/*@{*/
	virtual void ReadRestart(ifstreamT& in);
	virtual void WriteRestart(ofstreamT& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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
	virtual bool SetSystemFocus(double focus);

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
	KBC_ControllerT* fRightEdge;
	AutoArrayT<int>  fShiftedNodes;
	/*@}*/

	/** \name crack tip element group */
	/*@{*/
	int fTipElementGroup; /**< number of the element group controlling the tip position */
	double fTipX_0; /**< initial x-coordinate of the crack tip */
	double fTipY_0; /**< cleavage plane position */

	TrackingTypeT fTrackingType;
	double fTipThreshold; /**< threshold value which defines the tip position */

	StringT fTipOutputVariable;
	int fTipOutputCode; /**< output flag to generate data to locate the tip */
	int fTipColumnNum;  /**< column of output variable to locate tip */
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

	double fX_Left_last;  /**< left-most edge of the undeformed mesh */
	double fX_Right_last; /**< right-most edge of the undeformed mesh */
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
