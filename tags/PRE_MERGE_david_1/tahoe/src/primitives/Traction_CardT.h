/* $Id: Traction_CardT.h,v 1.5 2004-07-15 08:31:36 paklein Exp $ */
/* created: paklein (05/29/1996) */

#ifndef _TRACTION_T_H_
#define _TRACTION_T_H_

#include "Environment.h"

/* direct members */
#include "iArrayT.h"
#include "LocalArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ScheduleT;
class DomainIntegrationT;
class ElementSupportT;

/** natural boundary condition information */
class Traction_CardT
{
public:

	/** coordinate system for traction vector */
	enum CoordSystemT {kCartesian = 0, /**< x,y,z components of traction */
	                       kLocal = 1  /**< last component is normal component */
	                       }; 
	static CoordSystemT int2CoordSystemT(int i);

	/** constructor */
	Traction_CardT(void);

	/* modifiers */
	void EchoValues(const ElementSupportT& support, const DomainIntegrationT& domain,
		int element, int ndof, ifstreamT& in, ostream& out);

	void EchoValues(const ElementSupportT& support, int elem, int facet, int nLTf,
		 CoordSystemT coord_sys, const iArrayT& locnodenums, const dArray2DT& valuesT,
		 ostream& out);

	/** define traction values */
	void SetValues(const ElementSupportT& support, int elem, int facet, int nLTf,
		 CoordSystemT coord_sys, const iArrayT& locnodenums, const dArray2DT& valuesT);

	/* return the element and facet number specified for the force */
	void Destination(int& elem_num, int& facet_num) const;
	CoordSystemT CoordSystem(void) const;

	/* return the traction value: (ndof x nnd) */
	void CurrentValue(LocalArrayT& traction) const;
	
	/* local numbers of facet nodes */
	const iArrayT& LocalNodeNumbers(void) const;
	
	/* read/write (global) node numbers and equations */
	iArrayT& Nodes(void);
	iArrayT& Eqnos(void);

	/* read only */
	const iArrayT& Nodes(void) const;
	const iArrayT& Eqnos(void) const;
	
	/* write the standard header */
	void WriteHeader(ostream& out, int ndof) const;

	/* input operator for codes */
	friend istream& operator>>(istream& in, Traction_CardT::CoordSystemT& code);

private:
	
	/* parameters */
	int fElemNum;
	int fFacetNum;
	CoordSystemT fCoordSystem;
	const ScheduleT* fLTfPtr;
	LocalArrayT fValues;

	 /* local numbers of facet nodes */
	iArrayT fLocNodeNums;

	/* (global) node numbers/equations */
	iArrayT fNodes;
	iArrayT fEqnos;
};

/* inlines */

/* return the element and facet number specified for the force */
inline void Traction_CardT::Destination(int& elem_num, int& facet_num) const
{
	elem_num  = fElemNum;
	facet_num = fFacetNum;
}

inline Traction_CardT::CoordSystemT Traction_CardT::CoordSystem(void) const
{
	return fCoordSystem;
}

/* local numbers of facet nodes */
inline const iArrayT& Traction_CardT::LocalNodeNumbers(void) const
{
	return fLocNodeNums;
}

/* read/write (global) node numbers and equations */
inline iArrayT& Traction_CardT::Nodes(void) { return fNodes; }
inline iArrayT& Traction_CardT::Eqnos(void) { return fEqnos; }

/* read only */
/* read/write (global) node numbers and equations */
inline const iArrayT& Traction_CardT::Nodes(void) const { return fNodes; }
inline const iArrayT& Traction_CardT::Eqnos(void) const { return fEqnos; }

} // namespace Tahoe 
#endif /* _TRACTION_T_H_ */
