/* $Id: SSMFBaseT.h,v 1.1 2003-08-10 23:29:12 thao Exp $ */

#ifndef _SSMFBaseT_
#define _SSMFBaseT_

/* base class */
#include "ofstreamT.h"
#include "SmallStrainT.h"

namespace Tahoe {

/* forward declarations */
class SSSolidMatT;
class SSMatSupportT;
class OutputSetT;
class ifstreamT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class SSMFBaseT
{
  public:
    /** constructor */
    SSMFBaseT(const ElementSupportT& support, const FieldT& field);

    /** destructor */
    ~SSMFBaseT(void);

    void Initialize(void);
    
    /*output*/
    /* register self for output */
    virtual void RegisterOutput(void);

    /* send output */
    virtual void WriteOutput(void);

    /*write summary output file*/
    void WriteSummary(dArray2DT& output);

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

} // namespace Tahoe 
#endif /* _SSMFBaseT_ */
