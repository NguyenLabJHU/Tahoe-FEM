#ifndef _MATLIST_1D_T_H_
#define _MATLIST_1D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/** materials list for 1D structural analysis */
class SolidMatList1DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	SolidMatList1DT(int length, const SolidMatSupportT& support);

	/** read material data from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

private:
	
	/** \name errror messages */
	/*@{*/
	void Error_no_small_strain(ostream& out, int matcode) const;
	void Error_no_finite_strain(ostream& out, int matcode) const;
	void Error_no_multi_scale(ostream& out, int matcode) const;
	/*@}*/
};

} // namespace Tahoe 

#endif /* _MATLIST_1D_T_H_ */
