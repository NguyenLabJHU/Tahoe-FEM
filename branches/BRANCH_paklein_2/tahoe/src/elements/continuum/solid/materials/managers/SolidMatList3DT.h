/* $Id: SolidMatList3DT.h,v 1.6.8.1 2002-10-28 06:49:16 paklein Exp $ */
/* created: paklein (02/14/1997) */
#ifndef _MATLIST_3D_T_H_
#define _MATLIST_3D_T_H_

/* base class */
#include "StructuralMatListT.h"
#include "MaterialT.h"

namespace Tahoe {

/* forward declaration */
class ElasticT;
class SSMatSupportT;
class FDMatSupportT;

class SolidMatList3DT: public StructuralMatListT, public MaterialT
{
public:

	/* constructors */
	SolidMatList3DT(int length, const ElasticT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:
	
	/* errror messages */
	void Error_no_small_strain(ostream& out, int matcode) const;
	void Error_no_finite_strain(ostream& out, int matcode) const;

private:

	const ElasticT&      fElementGroup;
	const SSMatSupportT* fSSMatSupport;
	const FDMatSupportT* fFDMatSupport;
};

} // namespace Tahoe 
#endif /* _MATLIST_3D_T_H_ */
