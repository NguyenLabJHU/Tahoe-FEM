/*
 * File: ModelFileT.h
 */

/*
 * created      : PAK (12/15/1999)
 * last modified: PAK (12/15/1999)
 */

#ifndef _MODEL_FILE_T_H_
#define _MODEL_FILE_T_H_

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstream_x;
class ExodusT;

class ModelFileT
{
  public:

	enum Mode {kClosed = 0,
	             kRead = 1,
	            kWrite = 2};
	          
	enum Status {kFail = 0,
	               kOK = 1};	          

	/* constructor */
	ModelFileT(void);

	/* destructor */
	~ModelFileT(void);

	/* translate */
	Status Translate(const ExodusT& exo_file);

	/* open file */
	Status OpenRead(const StringT& file_name);
	Status OpenWrite(const StringT& file_name, bool extern_file);
	
	/* close */
	void Close(void);

	/* title */
	Status PutTitle(const StringT& title);
	Status GetTitle(StringT& title) const;
	
	/* coordinates */
	Status PutCoordinates(const dArray2DT& coords);
	Status GetDimensions(int& num_nodes, int& dimension) const;
	Status GetCoordinates(dArray2DT& coords) const;

	/* element sets */
	Status PutElementSet(int ID, const iArray2DT& set);
	Status GetElementSetID(iArrayT& ID) const;
	Status GetElementSetDimensions(int ID, int& num_elements, int& dimension) const;
	Status GetElementSet(int ID, iArray2DT& set) const;

	/* node sets */
	Status PutNodeSet(int ID, const iArrayT& set);
	Status GetNodeSetID(iArrayT& ID) const;
	Status GetNodeSetDimensions(int ID, int& num_nodes) const;
	Status GetNodeSet(int ID, iArrayT& set) const;
	Status GetNodeSets(const iArrayT& ID, iArrayT& set) const;

	/* side sets */
	Status PutSideSet(int ID, int element_set_ID, const iArray2DT& set);
	Status GetSideSetID(iArrayT& ID) const;
	Status GetSideSetDimensions(int ID, int& num_sides) const;
	Status GetSideSet(int ID, int& element_set_ID, iArray2DT& set) const;

  private:
  
  	/* check file version */
  	Status CheckVersion(ifstream_x& in) const;

  	/* advance to line after next occurence of key */
  	Status AdvanceStream(istream& in, const char* key) const;
  	Status AdvanceStreamToSubsection(istream& in, const char* section,
  		const char* subsection, int index) const;

	/* get set information from file */
	Status GetInformation(void);
  	
  	/* write data to file */
  	void WriteFile(bool extern_file) const;

	/* return reference to external or inline stream */
	ifstream_x& OpenExternal(ifstream_x& in,  ifstream_x& in2, ostream& out,
		bool verbose, const char* fail) const;

	/* open output file */
	ostream& OpenStream(ofstream& out, const StringT& file_name) const;

	/* convert string to lower case */
	void ToLower(char* str) const;
		
  private:

	Mode    fMode;
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
