/* $Id: DiffusionMatListT.h,v 1.7 2003-12-10 07:14:28 paklein Exp $ */
/* created: paklein (10/02/1999) */
#ifndef _DIFFUSE_MAT_LIST_T_H_
#define _DIFFUSE_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class DiffusionMatSupportT;

/** list of materials for diffusion analysis */
class DiffusionMatListT: public MaterialListT
{
public:

	/** enum defining material types */
	enum TypeT {
        kLinear = 1,
     kNonLinear = 2};

	/** constructors */
	DiffusionMatListT(int length, const DiffusionMatSupportT& support);
	DiffusionMatListT(void);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/
	
private:

	/** support for diffusion materials */
	const DiffusionMatSupportT* fDiffusionMatSupport;
};

} // namespace Tahoe 
#endif /* _DIFFUSE_MAT_LIST_T_H_ */
