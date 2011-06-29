/* $Id: SSVE_test.h,v 1.1 2009-04-23 14:38:49 tdnguye Exp $ */
/* created: TDN (5/31/2001) */
#ifndef H_SSVE_test_H
#define H_SSVE_test_H
 
#include "SSViscoelasticityT.h"

namespace Tahoe {

/** small strain linear viscoelastic constitutive law */
class SSVE_test: public SSViscoelasticityT
{
	public:

	/** constructor */
	SSVE_test(void);
			
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

};

} // namespace Tahoe 
#endif /*_SS_VISCO_H_*/
