/* $Id: MFSupportT.h,v 1.3 2003-11-12 19:21:19 thao Exp $ */

#ifndef _MFSupportT_
#define _MFSupportT_

/* base class */
#include "ofstreamT.h"
#include "iArrayT.h"
namespace Tahoe {

/* forward declarations */
class dArrayT;
class dArray2DT;
class ElementSupportT;
class OutputSetT;
class ifstreamT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class MFSupportT
{
  public:
    /** constructor */
    MFSupportT(const ElementSupportT& support);

    /** destructor */
    ~MFSupportT(void);

  protected:
    /* map nodal ordering of element group*/
    void MapOutput(void);

    /*write summary output file*/
    void WriteSummary(dArray2DT& output);

    /*utility funtions*/
    /*Assemble nodal material force vectors for element group*/
    void AssembleArray(const dArrayT& elem_val, dArrayT& global_val,const iArrayT& nodes);
    void AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val,const iArrayT& nodes);

    /*set nodal values for current element in local ordering*/
    void ExtractArray2D(const dArray2DT& global_val, dArray2DT& elem_val,const iArrayT& nodes);

    double ScalarProduct(double* pa, double* pb, const iArrayT& dims);
    /* returns displacement gradient*/
  
 protected:	
 
    /* material force output ID */
    int fMatForceOutputID;

    /*output set and dimensions*/
    OutputSetT* fOutputSet;
    int fNumGroupNodes;
    int fNumGroupElem;
    iArrayT fMap;

    ArrayT<StringT> fNID;
    int fnumset;

    int fhas_dissipation;

 private:
    /*fio for material support summary file*/ 
    bool fopen;
    ofstreamT fout;
    StringT fsummary_file;

    /*element support*/
    const ElementSupportT& fSupport; 
};

} // namespace Tahoe 
#endif /* _MFSupportT_ */
