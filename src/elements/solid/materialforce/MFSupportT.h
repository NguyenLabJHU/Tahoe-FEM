/* $Id: MFSupportT.h,v 1.6 2005-03-02 17:41:47 paklein Exp $ */
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

	/** set nodes over which material force is calculated and nodes on the boundary */
	void SetNodes(const ArrayT<StringT>& mf_nodes, const ArrayT<StringT>& boundary_nodes);

    /* map nodal ordering of element group*/
    void MapOutput(void);

    /*write summary output file*/
    void WriteSummary(dArray2DT& output);

    /*utility funtions*/
    /*Gather displacements for element group*/
    void GatherDisp(const dArray2DT& global_disp, dArrayT& group_disp, const iArrayT& nodes);
    /*Assemble nodal material force vectors for element group*/
    void AssembleArray(const dArrayT& elem_val, dArrayT& global_val,const iArrayT& nodes);
    void AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val,const iArrayT& nodes);

    /*set nodal values for current element in local ordering*/
    void ExtractArray2D(const dArray2DT& global_val, dArray2DT& elem_val,const iArrayT& nodes);

    double ScalarProduct(const double* pa, const double* pb, const iArrayT& dims);
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

	ArrayT<StringT> fBoundID;
    iArrayT fExclude;

    bool fhas_dissipation;

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
