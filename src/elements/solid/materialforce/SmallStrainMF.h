/* $Id: SmallStrainMF.h,v 1.3 2003-04-10 22:03:07 thao Exp $ */

#ifndef _SMALL_STRAIN_MF_H_
#define _SMALL_STRAIN_MF_H_

/* base class */
#include "SmallStrainT.h"
#include "ofstreamT.h"
namespace Tahoe {

/* forward declarations */
class SSMatSupportT;
class SSSolidMatT;
class OutputSetT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class SmallStrainMF: public SmallStrainT
{
  public:
    /** constructor */
    SmallStrainMF(const ElementSupportT& support, const FieldT& field);

    /** destructor */
    ~SmallStrainMF(void);

    virtual void Initialize(void);
    virtual void SetGlobalShape(void);

    /* register self for output */
    virtual void RegisterOutput(void);

    /* send output */
    virtual void WriteOutput(void);

 protected:
    /*material force evaluation*/
    void ComputeMatForce(dArray2DT& output);
    void MatForceVolMech(dArrayT& elem_val);
    void MatForceDissip(dArrayT& elem_val, const dArrayT& statev);
    void MatForceSurfMech(dArrayT& global_val);

    /*Assemble nodal material force vectors for element group*/
    void AssembleMatForce(const dArrayT& elem_val, dArrayT& global_val,	const iArrayT& nodes);

    /* map nodal ordering of element group*/
    void MapOutput(void);

    /* returns displacement gradient*/
    const dMatrixT& DisplacementGradient(void) const;

 protected:	

    /*current material*/
    SSSolidMatT* fCurrSSMat;

    /* material force output ID */
    int fMatForceOutputID;

    /*output set and dimensions*/
    OutputSetT* fOutputSet;
    int fNumGroupNodes;
    iArrayT fMap;

    ArrayT<StringT> fNID;
    int fnumset;

    bool fopen;
    ofstreamT fout;
    StringT fsummary_file;

 private:
    ArrayT<dMatrixT> fGradU_List;
    dArrayT fgrad_viscstretch;
    dSymMatrixT fdevQ;
    dSymMatrixT fviscstress;
};

/* inlines */
inline const dMatrixT& SmallStrainMF::DisplacementGradient(void) const
{
	return fGradU_List[CurrIP()];
}

} // namespace Tahoe 
#endif /* _SMALLSTRAIN_MF_H_ */
