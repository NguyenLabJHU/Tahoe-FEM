/* $Id: KBC_CardT.h,v 1.6 2004-07-15 08:31:36 paklein Exp $ */
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
	enum CodeT {kFix = 0,
                kDsp = 1,
                kVel = 2,
                kAcc = 3,
                kNull= 4};

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
