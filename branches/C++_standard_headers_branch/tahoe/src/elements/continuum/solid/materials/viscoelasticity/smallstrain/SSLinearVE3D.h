/* $Id: SSLinearVE3D.h,v 1.4 2008-08-09 15:00:37 tdnguye Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_LINEAR_VE_3D_H_
#define _SS_LINEAR_VE_3D_H_

/* base class */
#include "SSViscoelasticityT.h"

namespace Tahoe {

/** base class for standard solid Kirchhoff St. Venant constitutive models 
 * constitutive law */
class SSLinearVE3D: public SSViscoelasticityT
{
	public:
	
	/** constructor */
	SSLinearVE3D(void);
		
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

};

}

#endif  /* _SS_LINEAR_VE_3D_H_ */
