// DEVELOPMENT
#ifndef _OUTPUTLATTICE_T_H_
#define _OUTPUTLATTICE_T_H_

#include <iostream.h>
#include <fstream.h>

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "IOBaseT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class OutputBaseT;
class iArray2DT;
class dArray2DT;
class OutputSetT;

/** class to handle output of atom files in 
 *  format-independent output. Similar to IOManager*/


class OutPutLatticeT 
{

public:
  
         // constructor
        OutPutLatticeT(ostream& outfile, const StringT& program_name,
		       const StringT& version, const StringT& title, 
		       const StringT& input_file,
		       IOBaseT::FileTypeT output_format);
	// destructor
	~OutPutLatticeT(void);

	static OutputBaseT* NewOutput(const StringT& program_name,
				      const StringT& version, 
				      const StringT& title, const StringT& input_file,
				      IOBaseT::FileTypeT output_format, 
				      ostream& log);

	/* output functions */
	void SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_id);
	void SetBounds(const dArray2DT& bounds);
	void SetTypes(const iArrayT& types);
	void SetParts(const iArrayT& parts);

	void WriteGeometryFile(const StringT& file_name, 
			       IOBaseT::FileTypeT format) const;
	void WriteGeometry(void);
	void SetOutputTime(double time);

	int AddElementSet(const OutputSetT& output_set);
	const ArrayT<OutputSetT*>& ElementSets(void) const;
	void AddNodeSet(const iArrayT& nodeset, const StringT& setID);
	virtual void WriteOutput(int ID, const dArray2DT& n_values, 
				 const dArray2DT& e_values);
	const OutputSetT& OutputSet(int ID) const;

protected:

	/* output formatter */
	IOBaseT::FileTypeT fOutputFormat;
	OutputBaseT* fOutput;

	ostream& fLog;


private:

	double fOutputTime;

};// class


/* inline */
inline void OutPutLatticeT::SetOutputTime(double time) { fOutputTime = time; }


} // namespace Tahoe 

#endif

