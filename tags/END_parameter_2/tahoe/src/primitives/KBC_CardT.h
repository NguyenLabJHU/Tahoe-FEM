/* $Id: KBC_CardT.h,v 1.5.24.1 2004-03-22 18:36:18 paklein Exp $ */
/* created: paklein (05/23/1996) */
#ifndef _KBC_CARD_T_H_
#define _KBC_CARD_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

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

	/** \name modifiers */
	/*@{*/
	void SetValues(istream& in);
	void SetValues(int node, int dof, CodeT code, const ScheduleT* schedule, double value);
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	int Node(void) const;
	int DOF(void) const;
	CodeT Code(void) const;
	const ScheduleT* Schedule(void) const { return fSchedule; };
	/*@}*/

	/* returns the value of the BC */
	double Value(void) const;

	/* I/O */
	static void WriteHeader(ostream& out);
	void WriteValues(ostream& out) const;

	/* input operator for codes */
	static CodeT int_to_CodeT (int i);
	friend istream& operator>>(istream& in, KBC_CardT::CodeT& code);
	
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
