/* $Id: BEMelement.h,v 1.7 2004-07-15 08:27:56 paklein Exp $ */
/* created: AFLP (02/28/1998) */

#ifndef _BEM_ELEMENT_H_
#define _BEM_ELEMENT_H_

/* base class */
#include "ElementBaseT.h"

namespace Tahoe {

/* forward declarations */
class StringT;

class BEMelement: public ElementBaseT
{
public:

	/* constructor */
	BEMelement(const ElementSupportT& support, const FieldT& field, 
		const StringT& infile);

	/* destructor */
	virtual ~BEMelement(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* solution calls */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);
	
	/* output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);
	virtual void SendOutput(int kincode);
	 		
protected: /* for derived classes only */
		 	
	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

protected:

	/* input file name */
	const StringT& fInfile;

};

} // namespace Tahoe 
#endif /* _BEM_ELEMENT_H_ */
