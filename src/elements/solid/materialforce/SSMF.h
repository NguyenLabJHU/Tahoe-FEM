/* $Id: SSMF.h,v 1.1 2003-08-11 00:22:10 thao Exp $ */

#ifndef _SSMF_H_
#define _SSMF_H_

/* base class */
#include "SmallStrainT.h"
#include "MFSupportT.h"
#include "ofstreamT.h"
namespace Tahoe {

/* forward declarations */
class SSSolidMatT;
class SSMatSupportT;
class ifstreamT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class SSMF: public SmallStrainT, public MFSupportT
{
  public:
    /** constructor */
    SSMF(const ElementSupportT& support, const FieldT& field);

    virtual void Initialize(void);
    virtual void SetGlobalShape(void);
    
    /*accessor for displacement gradient*/
    const dMatrixT& DisplacementGradient(void) const;

    /*output*/
    /* register self for output */
    virtual void RegisterOutput(void);

    /* send output */
    virtual void WriteOutput(void);

 protected:
    /*material force evaluation*/
    void ComputeMatForce(dArray2DT& output);
    void MatForceVolMech(dArrayT& elem_val);
    void MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch);
    void MatForceSurfMech(dArrayT& global_val);

    /*utility funtions*/
    /*extrapolate element ip state variables values to nodes*/
    void Extrapolate(void);
  
 protected:	

    /*current material*/
    SSSolidMatT* fCurrSSMat;

 private:
    ArrayT<dMatrixT> fGradU_List;
    dMatrixT fEshelby;

    /*nodal and interpolated body force*/
    dArrayT fBodyForce;
    dArrayT fip_body;
    
    /*material force */
    dArrayT fMatForce;
    dArrayT fDissipForce;
    dArrayT felem_rhs;

    /*internal variables*/
    iArrayT fInternalDOF;
    double fNumInternalVal;

    dArrayT fGlobalMass;
    dArray2DT fGlobalVal;

    dArrayT felem_mass;
    dArray2DT felem_val;
    
    dArray2DT fGradInternalStrain;
};

/* inlines */
inline const dMatrixT& SSMF::DisplacementGradient(void) const
{
	return fGradU_List[CurrIP()];
}

} // namespace Tahoe 
#endif /* _SSMF_H_ */
