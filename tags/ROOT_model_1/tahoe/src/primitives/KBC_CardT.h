/* $Id: KBC_CardT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/23/1996)                                          */

#ifndef _KBC_CARD_T_H_
#define _KBC_CARD_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class LoadTime;

class KBC_CardT
{
public:

	friend class NodeManagerPrimitive;
	
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
	void SetSchedule(const LoadTime* LTfPtr);
	
	/* accessors */
	int Node(void) const;
	int DOF(void) const;
	CodeT Code(void) const;
	int LTfNum(void) const;
		
	/* returns the value of the BC */
	double Value(void) const;

	/* I/O */
	void WriteHeader(ostream& out) const;
	void WriteValues(ostream& out) const;

	/* input operator for codes */
	friend istream& operator>>(istream& in, KBC_CardT::CodeT& code);
	
protected:

	int      fnode;
	int      fdof;
	CodeT    fcode;
	int      fnLTf;
	double   fvalue;			
	const LoadTime* fLTfPtr;
};

/* in-lines */

/* accessors */
inline int KBC_CardT::Node(void) const   { return fnode; }
inline int KBC_CardT::DOF(void) const    { return fdof;  }
inline KBC_CardT::CodeT KBC_CardT::Code(void) const   { return fcode; }
inline int KBC_CardT::LTfNum(void) const { return fnLTf; }

#endif /* _KBC_CARD_T_H_ */
