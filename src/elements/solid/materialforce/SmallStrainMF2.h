/* $Id: SmallStrainMF2.h,v 1.1 2003-05-15 05:15:23 thao Exp $ */

#ifndef _SMALL_STRAIN_MF2_H_
#define _SMALL_STRAIN_MF2_H_

/* base class */
#include "SmallStrainT.h"
#include "ofstreamT.h"
namespace Tahoe {

/* forward declarations */
class SSMatSupportT;
class SSSolidMatT;
class OutputSetT;
class ifstreamT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class SmallStrainMF2: public SmallStrainT
{
  public:
    /** constructor */
    SmallStrainMF2(const ElementSupportT& support, const FieldT& field);

    /** destructor */
    ~SmallStrainMF2(void);

    virtual void Initialize(void);
    virtual void SetGlobalShape(void);
    
    /*accessor for displacement gradient*/
    const dMatrixT& DisplacementGradient(void) const;

    /*output*/
    /* register self for output */
    virtual void RegisterOutput(void);

    /* send output */
    virtual void WriteOutput(void);

    /*write summary output file*/
    void WriteSummary(dArray2DT& output);


 protected:
    /*material force evaluation*/
    void ComputeMatForce(dArray2DT& output);
    void MatForceVolMech(dArrayT& elem_val);
    void MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch);
    void MatForceSurfMech(dArrayT& global_val);

    /*utility funtions*/
    /*Assemble nodal material force vectors for element group*/
    void AssembleArray(const dArrayT& elem_val, dArrayT& global_val,const iArrayT& nodes);
    void AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val,const iArrayT& nodes);

    /*set nodal values for current element in local ordering*/
    void ExtractArray2D(const dArray2DT& global_val, dArray2DT& elem_val,const iArrayT& nodes);

    /*extrapolate element ip state variables values to nodes*/
    void Extrapolate(void);

    /* map nodal ordering of element group*/
    void MapOutput(void);

    double ScalarProduct(double* pa, double* pb, const iArrayT& dims);
    /* returns displacement gradient*/
  
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


    int fhas_dissipation;
    iArrayT finternaldof;
    int fnumval;
    dArray2DT fgrad_intstrain;
    dSymMatrixT fgradstrain;
    dSymMatrixT fstress;

    dArrayT fglobal_mass;
    dArrayT felem_mass;
    dArray2DT fglobal_val;
    dArray2DT felem_val;
};

/* inlines */
inline const dMatrixT& SmallStrainMF2::DisplacementGradient(void) const
{
	return fGradU_List[CurrIP()];
}

} // namespace Tahoe 
#endif /* _SMALLSTRAIN_MF2_H_ */
