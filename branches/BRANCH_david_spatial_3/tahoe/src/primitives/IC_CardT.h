/* $Id: IC_CardT.h,v 1.6.30.1 2005-07-25 02:37:24 paklein Exp $ */
/* created: paklein (07/16/1997) */
#ifndef _IC_CARD_T_H_
#define _IC_CARD_T_H_

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** container class for kinematic initial condition data.
 * Handles mainly I/O and provides access to data via (inline) accessors */
class IC_CardT
{
public:

	/** usage mode */
	enum ModeT {
		kUndefined,
		kNode,
		kSet
	};

	/** constructor */
	IC_CardT(void);

	/** \name modifier */
	/*@{*/
	void SetValues(int node, int dof, int order, double value);
	void SetValues(const StringT& ID, int dof, int order, double value);
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	int Node(void) const;
	const StringT& ID(void) const;
	int DOF(void) const;
	int Order(void) const;
	double Value(void) const;
	ModeT Mode(void) const { return fmode; };
	/*@}*/

private:

	/** \name node or node set */
	/*@{*/
	int fnode;
	StringT fID;
	/*@}*/

	int    fdof;
	int    forder; /**< time derivative */
	double fvalue;
	ModeT  fmode;
};

/* inline functions */

/* accessors */
inline int IC_CardT::Node(void) const     { return fnode;  }
inline const StringT& IC_CardT::ID(void) const { return fID; };
inline int IC_CardT::DOF(void) const      { return fdof;   }
inline int IC_CardT::Order(void) const    { return forder; }
inline double IC_CardT::Value(void) const { return fvalue; }

} // namespace Tahoe 
#endif /* _IC_CARD_T_H_ */
