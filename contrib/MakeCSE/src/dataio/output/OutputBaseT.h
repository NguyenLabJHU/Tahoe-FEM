// file: OutputBaseT.h

// created      : SAW (05/18/1999)
// last modified: PAK (03/07/2000)

/* initialization:
 *     1. construct
 *     2. SetCoordinates
 *     3. AddElementSet
 */

#include "IOBaseT.h"

#ifndef _OUTPUTBASE_T_H_
#define _OUTPUTBASE_T_H_

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "GeometryT.h"

/* forward declarations */
class dArray2DT;
class iArray2DT;
class OutputSetT;

class OutputBaseT: public IOBaseT
{
   public:
  	
  	/* constructor */
	OutputBaseT(ostream& out, const ArrayT<StringT>& outstrings);

  	/* destructor */
	~OutputBaseT(void);

	/* accessors */
	const StringT& Version(void) const;
	const StringT& CodeName(void) const;
	const StringT& Title(void) const;

	/* increment sequence, create new output file series */
	virtual void NextTimeSequence(int sequence_number); 

	/* set nodal coordinates - passing in map implies coordinates list is not
	 * complete or compact */
	virtual void SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_map);

	/* register the output for an element set. returns the output ID */
	int AddElementSet(const OutputSetT& output_set);
	const ArrayT<OutputSetT*>& ElementSets(void) const;
	int NumElements(void) const;

	void AddNodeSet(const iArrayT& nodeset, int setID);
	void AddSideSet(const iArray2DT& sideset, int setID, int group_ID); // ID from AddElementSet

	/* output functions */
	virtual void WriteGeometry(void) = 0;
	void WriteGeometryFile(const StringT& file_name, IOBaseT::IOFileType format) const;
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values, 
		const dArray2DT& e_values) = 0;

   protected:

	enum DataType {kNode = 0, 
	            kElement = 1};

	void LocalConnectivity(const iArrayT& node_map, const iArray2DT& connects, iArray2DT& local_connects) const;
   protected:

	StringT fTitle;    // title: description of problem
	StringT fCodeName; // qa_record codename and version
	StringT fVersion;  // qa_record inputfile version
	StringT fOutroot;  // root of all output files

	/* output data */
	const dArray2DT*        fCoordinates;
	const iArrayT*          fNodeMap;

	AutoArrayT<OutputSetT*> fElementSets;
	
	AutoArrayT<const iArrayT*>   fNodeSets;
	AutoArrayT<const iArray2DT*> fSideSets;
	iAutoArrayT                  fNodeSetIDs;
	iAutoArrayT                  fSideSetIDs;
	iAutoArrayT                  fSSGroupID; // fElementList group number

	int fSequence; // solution sequence number

};

inline const ArrayT<OutputSetT*>& OutputBaseT::ElementSets(void) const
{
	return fElementSets;
}

/* accessors */
inline const StringT& OutputBaseT::Version(void) const { return fVersion; }
inline const StringT& OutputBaseT::CodeName(void) const { return fCodeName; }
inline const StringT& OutputBaseT::Title(void) const { return fTitle; }

#endif /* _OUTPUTMANAGER_H_ */
