/* $Id: ModelFileT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/15/1999)                                          */

#ifndef _MODEL_FILE_T_H_
#define _MODEL_FILE_T_H_

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class ExodusT;

class ModelFileT
{
public:

	enum ModeT {kClosed = 0,
	              kRead = 1,
	             kWrite = 2};
	
	enum StatusT {kFail = 0,
	                kOK = 1};	

	/* constructor */
	ModelFileT(void);

	/* destructor */
	~ModelFileT(void);

	/* translate */
	StatusT Translate(const ExodusT& exo_file);

	/* open file */
	StatusT OpenRead(const StringT& file_name);
	StatusT OpenWrite(const StringT& file_name, bool extern_file);
	
	/* close */
	void Close(void);

	/* title */
	StatusT PutTitle(const StringT& title);
	StatusT GetTitle(StringT& title) const;
	
	/* coordinates */
	StatusT PutCoordinates(const dArray2DT& coords);
	StatusT GetDimensions(int& num_nodes, int& dimension) const;
	StatusT GetCoordinates(dArray2DT& coords) const;

	/* element sets */
	StatusT PutElementSet(int ID, const iArray2DT& set);
	StatusT GetElementSetID(iArrayT& ID) const;
	StatusT GetElementSetDimensions(int ID, int& num_elements, int& dimension) const;
	StatusT GetElementSet(int ID, iArray2DT& set) const;

	/* node sets */
	StatusT PutNodeSet(int ID, const iArrayT& set);
	StatusT GetNodeSetID(iArrayT& ID) const;
	StatusT GetNodeSetDimensions(int ID, int& num_nodes) const;
	StatusT GetNodeSet(int ID, iArrayT& set) const;
	StatusT GetNodeSets(const iArrayT& ID, iArrayT& set) const;

	/* side sets */
	StatusT PutSideSet(int ID, int element_set_ID, const iArray2DT& set);
	StatusT GetSideSetID(iArrayT& ID) const;
	StatusT GetSideSetDimensions(int ID, int& num_sides) const;
	StatusT GetSideSet(int ID, int& element_set_ID, iArray2DT& set) const;

private:

	/* check file version */
	StatusT CheckVersion(ifstreamT& in) const;

	/* advance to line after next occurence of key */
	StatusT AdvanceStream(istream& in, const char* key) const;
	StatusT AdvanceStreamToSubsection(istream& in, const char* section,
		const char* subsection, int index) const;

	/* get set information from file */
	StatusT GetInformation(void);
	
	/* write data to file */
	void WriteFile(bool extern_file) const;

	/* return reference to external or inline stream */
	ifstreamT& OpenExternal(ifstreamT& in,  ifstreamT& in2, ostream& out,
		bool verbose, const char* fail) const;

	/* open output file */
	ostream& OpenStream(ofstream& out, const StringT& file_name) const;

	/* convert string to lower case */
	void ToLower(char* str) const;
		
private:

	ModeT   fMode;
	StringT fFileName;

	/* dimensions */
	int fNumNodes;
	int fDimension;

	/* set info lists */
	iArray2DT fElementID;
	iArray2DT fNodeSetID;
	iArray2DT fSideSetID;

	/******* only used for writing files ********/
	bool fExternFile;

	/* title */
	StringT fTitle;

	/* coordinates */
	dArray2DT fCoordinates;

	/* element sets */
	ArrayT<iArray2DT*> fElementSets;
	ArrayT<iArrayT*>   fNodeSets;
	ArrayT<iArray2DT*> fSideSets;
};

#endif /* _MODEL_FILE_T_H_ */
