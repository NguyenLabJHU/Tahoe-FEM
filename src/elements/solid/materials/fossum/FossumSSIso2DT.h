/* $Id: FossumSSIso2DT.h,v 1.7 2003-01-31 09:46:54 paklein Exp $ */
#ifndef _FOSSUM_SS_ISO_2D_T_H_
#define _FOSSUM_SS_ISO_2D_T_H_
/* DEVELOPMENT */

/* base class */
#include "Material2DT.h"
#include "FossumSSIsoT.h"

#include "SSSolidMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

class FossumSSIso2DT: public FossumSSIsoT, 
                      public Material2DT
{
  public:

        /* constructor */
        FossumSSIso2DT(ifstreamT& in, const SSMatSupportT& support);

        /* initialization */
        virtual void Initialize(void);

        /* returns elastic strain (3D) */
        virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
                const ElementCardT& element, int ip);

        /* print parameters */
        virtual void Print(ostream& out) const;
        virtual void PrintName(ostream& out) const;
        
        /* modulus */
        virtual const dMatrixT& c_ijkl(void);
        virtual const dMatrixT& cdisc_ijkl(void);
        
        /* stress */
        virtual const dSymMatrixT& s_ij(void);

        /* returns the strain energy density for the specified strain */
        virtual double StrainEnergyDensity(void);

  private:
  
        /* return values */
        dSymMatrixT        fStress2D;
        dMatrixT        fModulus2D;

        /* work space */
        dSymMatrixT        fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _FOSSUM_SS_ISO_2D_T_H_ */
