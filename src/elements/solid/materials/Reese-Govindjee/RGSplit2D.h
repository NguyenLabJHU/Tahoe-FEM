/* created: TDN (01/22/2001) */

#ifndef _RGSplit_2D_
#define _RGSplit_2D_

/* base classes */
#include "RGSplit3D.h"

namespace Tahoe {

/* forward declarations */
class PotentialT;

class RGSplit2D: public RGSplit3D
{
   public:
  
	/* constructor */
	RGSplit2D(ifstreamT& in, const FSMatSupportT& support);
	~RGSplit2D(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
	virtual const dMatrixT& C_IJKL(void);
	virtual const dSymMatrixT& S_IJ(void);

    /*compute output variables*/ 
    virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
    virtual void ComputeOutput(dArrayT& output);

   protected:
	/* return values */
	dMatrixT	fModulus2D;
	dSymMatrixT fStress2D;

   private:  
	/*work space*/
	dSymMatrixT fb2D;
};
}
#endif /* _RGSplit_2D_ */
