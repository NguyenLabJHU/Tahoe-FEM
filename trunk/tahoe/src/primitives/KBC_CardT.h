/* $Id: KBC_CardT.h,v 1.7 2005-08-03 07:52:22 paklein Exp $ */
/* created: paklein (05/23/1996) */
#ifndef _KBC_CARD_T_H_
#define _KBC_CARD_T_H_

namespace Tahoe {

/* forward declaration */
class ScheduleT;

/** container to hold kinematic boundary condition specifications */
class KBC_CardT
{
public:

	friend class NodeManagerT;
	
	/** codes */
	enum CodeT {kFix = 0, /**< fixed to 0 */
                kDsp = 1, /**< prescribed displacement (0th time derivative) */
                kVel = 2, /**< prescribed velocity (1st time derivative) */
                kAcc = 3, /**< prescribed acceleration (2nd time derivative) */
                kNull= 4  /**< component is marked as inactive, but integrator does not modify */
                };

	/** \name constructor */
	/*@{*/
	KBC_CardT(void);
	KBC_CardT(int node, int dof, CodeT code, const ScheduleT* schedule, double value);
	/*@}*/

	/** modifier */
	void SetValues(int node, int dof, CodeT code, const ScheduleT* schedule, double value);

	/** \name accessors */
	/*@{*/
	int Node(void) const;
	int DOF(void) const;
	CodeT Code(void) const;
	const ScheduleT* Schedule(void) const { return fSchedule; };
	/*@}*/

	/* returns the value of the BC */
	double Value(void) const;

	/* input operator for codes */
	static CodeT int2CodeT(int i);
	
protected:

	int      fnode;
	int      fdof;
	CodeT    fcode;
	double   fvalue;			
	const ScheduleT* fSchedule;
};

/* in-lines */

/* accessors */
inline int KBC_CardT::Node(void) const   { return fnode; }
inline int KBC_CardT::DOF(void) const    { return fdof;  }
inline KBC_CardT::CodeT KBC_CardT::Code(void) const   { return fcode; }

} // namespace Tahoe 
#endif /* _KBC_CARD_T_H_ */
