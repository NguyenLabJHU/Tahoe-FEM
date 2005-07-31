/* $Id: IC_CardT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (07/16/1997)                                          */
/* Container class for kinematic initial condition data.                  */
/* Handles mainly I/O and provides access to data via                     */
/* (inline) accessors.                                                    */

#ifndef _IC_CARD_T_H_
#define _IC_CARD_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;

class IC_CardT
{
public:

	/* codes */
	enum CodeT {kDsp = 1,
	            kVel = 2,
	            kAcc = 3};

	/* constructor */
	IC_CardT(void);

	/* modifiers */
	void SetValues(ifstreamT& in);
	void SetValues(int node, int dof, CodeT code, double value);
	
	/* accessors */
	int Node(void) const;
	int DOF(void) const;
	CodeT Code(void) const;
	double Value(void) const;

	/* I/O */
	void WriteHeader(ostream& out) const;
	void WriteValues(ostream& out) const;
	
	/* input operator for codes */
	friend istream& operator>>(istream& in, IC_CardT::CodeT& code);
	
private:

	int    fnode;
	int    fdof;
	CodeT  fcode;
	double fvalue;			
};

/* inline functions */

/* accessors */
inline int IC_CardT::Node(void) const     { return fnode;  }
inline int IC_CardT::DOF(void) const      { return fdof;   }
inline IC_CardT::CodeT IC_CardT::Code(void) const { return fcode;  }
inline double IC_CardT::Value(void) const { return fvalue; }

#endif /* _IC_CARD_T_H_ */
