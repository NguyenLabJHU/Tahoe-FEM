/* $Id: KBC_CardT.h,v 1.2.2.2 2002-04-30 01:30:22 paklein Exp $ */
/* created: paklein (05/23/1996) */

#ifndef _KBC_CARD_T_H_
#define _KBC_CARD_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ScheduleT;

class KBC_CardT
{
public:

	friend class NodeManagerT;
	
	/* codes */
	enum CodeT {kFix = 0,
                kDsp = 1,
                kVel = 2,
                kAcc = 3};

	/* constructor */
	KBC_CardT(void);
	KBC_CardT(int node, int dof, CodeT code, int nLTF, double value);

	/* modifiers */
	void SetValues(istream& in);
	void SetValues(int node, int dof, CodeT code, int nLTF, double value);
	void SetSchedule(const ScheduleT* schedule);
	
	/* accessors */
	int Node(void) const;
	int DOF(void) const;
	CodeT Code(void) const;
	int ScheduleNum(void) const;
		
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
	int      fSchedNum;
	double   fvalue;			
	const ScheduleT* fSchedule;
};

/* in-lines */

/* accessors */
inline int KBC_CardT::Node(void) const   { return fnode; }
inline int KBC_CardT::DOF(void) const    { return fdof;  }
inline KBC_CardT::CodeT KBC_CardT::Code(void) const   { return fcode; }
inline int KBC_CardT::ScheduleNum(void) const { return fSchedNum; }

#endif /* _KBC_CARD_T_H_ */
