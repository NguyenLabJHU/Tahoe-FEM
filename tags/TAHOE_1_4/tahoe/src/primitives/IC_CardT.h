/* $Id: IC_CardT.h,v 1.5 2002-07-05 22:28:33 paklein Exp $ */
/* created: paklein (07/16/1997) */

#ifndef _IC_CARD_T_H_
#define _IC_CARD_T_H_

#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** container class for kinematic initial condition data.
 * Handles mainly I/O and provides access to data via (inline) accessors */
class IC_CardT
{
public:

	/** constructor */
	IC_CardT(void);

	/** \name modifiers */
	/*@{*/
	void SetValues(ifstreamT& in);
	void SetValues(int node, int dof, int order, double value);
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	int Node(void) const;
	int DOF(void) const;
	int Order(void) const;
	double Value(void) const;
	/*@}*/

	/** \name I/O methods */
	/*@{*/
	static void WriteHeader(ostream& out);
	void WriteValues(ostream& out) const;
	/*@}*/
	
private:

	int    fnode;
	int    fdof;
	int    forder; /**< time derivative */
	double fvalue;			
};

/* inline functions */

/* accessors */
inline int IC_CardT::Node(void) const     { return fnode;  }
inline int IC_CardT::DOF(void) const      { return fdof;   }
inline int IC_CardT::Order(void) const    { return forder; }
inline double IC_CardT::Value(void) const { return fvalue; }

} // namespace Tahoe 
#endif /* _IC_CARD_T_H_ */
