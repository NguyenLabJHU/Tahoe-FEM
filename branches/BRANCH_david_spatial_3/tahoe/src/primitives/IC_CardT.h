/* $Id: IC_CardT.h,v 1.6 2004-07-15 08:31:36 paklein Exp $ */
/* created: paklein (07/16/1997) */
#ifndef _IC_CARD_T_H_
#define _IC_CARD_T_H_

namespace Tahoe {

/** container class for kinematic initial condition data.
 * Handles mainly I/O and provides access to data via (inline) accessors */
class IC_CardT
{
public:

	/** constructor */
	IC_CardT(void);

	/** modifier */
	void SetValues(int node, int dof, int order, double value);
	
	/** \name accessors */
	/*@{*/
	int Node(void) const;
	int DOF(void) const;
	int Order(void) const;
	double Value(void) const;
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
