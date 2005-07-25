/* $Id: KBC_CardT.h,v 1.6.30.1 2005-07-25 02:37:24 paklein Exp $ */
/* created: paklein (05/23/1996) */
#ifndef _KBC_CARD_T_H_
#define _KBC_CARD_T_H_

/* direct members */
#include "StringT.h"

namespace Tahoe {

/* forward declaration */
class ScheduleT;
class StringT;

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

	/** usage mode */
	enum ModeT {
		kUndefined,
		kNode,
		kSet
	};

	/** \name constructor */
	/*@{*/
	KBC_CardT(void);
//	KBC_CardT(int node, int dof, CodeT code, const ScheduleT* schedule, double value);
	/*@}*/

	/** \name modifier */
	/*@{*/
	void SetValues(int node, int dof, CodeT code, const ScheduleT* schedule, double value);
	void SetValues(const StringT& ID, int dof, CodeT code, const ScheduleT* schedule, double value);
	/*@}*/

	/** \name accessors */
	/*@{*/
	int Node(void) const;
	const StringT& ID(void) const;
	int DOF(void) const;
	CodeT Code(void) const;
	ModeT Mode(void) const;
	const ScheduleT* Schedule(void) const { return fSchedule; };
	/*@}*/

	/* returns the value of the BC */
	double Value(void) const;

	/* input operator for codes */
	static CodeT int2CodeT(int i);
	
protected:

	/** \name node or node set */
	/*@{*/
	int fnode;
	StringT fID;
	/*@}*/

	int      fdof;
	CodeT    fcode;
	ModeT    fmode;
	double   fvalue;			
	const ScheduleT* fSchedule;
};

/* in-lines */

/* accessors */
inline int KBC_CardT::Node(void) const   { return fnode; }
inline const StringT& KBC_CardT::ID(void) const   { return fID; }
inline int KBC_CardT::DOF(void) const    { return fdof;  }
inline KBC_CardT::CodeT KBC_CardT::Code(void) const   { return fcode; }
inline KBC_CardT::ModeT KBC_CardT::Mode(void) const   { return fmode; }

} // namespace Tahoe 
#endif /* _KBC_CARD_T_H_ */
