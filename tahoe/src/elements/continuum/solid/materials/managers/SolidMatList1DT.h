#ifndef _MATLIST_1D_T_H_
#define _MATLIST_1D_T_H_

/* base classes */
#include "StructuralMatListT.h"
#include "MaterialT.h"

namespace Tahoe {

/* forward declaration */
class ElasticT;
class SSMatSupportT;
class FDMatSupportT;
//class MultiScaleT;

class SolidMatList1DT: public StructuralMatListT, public MaterialT
{
public:

	/* constructor */
	SolidMatList1DT(int length, const ElasticT& element_group);

	/* read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:
	
	/* errror messages */
	void Error_no_small_strain(ostream& out, int matcode) const;
	void Error_no_finite_strain(ostream& out, int matcode) const;
	void Error_no_multi_scale(ostream& out, int matcode) const;

private:

	const ElasticT&      fElementGroup;
	const SSMatSupportT* fSSMatSupport;
	const FDMatSupportT* fFDMatSupport;
//	const MultiScaleT*   fMultiScale;
};

} // namespace Tahoe 

#endif /* _MATLIST_1D_T_H_ */
