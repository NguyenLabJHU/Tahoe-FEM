// file: AbaqusT.h

// created:       SAW (05/30/00)

/* 
   Add a variable type: 
       (1) add to VariableKey
       (2) modify WriteOutputDef
       (3) modify ReadPreliminary
       (4) modify WritePreliminary
       (5) modify WriteElementHeader
       (6) modify Record2VariableKey

   If variable type added, update key to label mapping in AbaqusInput and AbaqusOutput.

*/

#ifndef _ABAQUS_T_H_
#define _ABAQUS_T_H_

/* direct members */
#include "iosfwd.h"
#include "ArrayT.h"
#include "GeometryT.h"

/* forward declarations */
class StringT;
class iArrayT;
class dArrayT;
class dArray2DT;
class iArray2DT;

class AbaqusT
{
 public:

	enum VariableKey { aUnknown       = -1,
			   aSDV           = 5,
			   aStress        = 11,
			   aTotalStrain   = 21,
			   aUVARM         = 87,
			   aDisplacement  = 101,
			   aVelocity      = 102,
			   aAcceleration  = 103,
			   aCoordVariable = 107,
			   aPrinStress    = 401 };

	enum OutputFlag { aNoFlag      = -1,
			  aElementFlag = 0,
			  aNodalFlag   = 1,
			  aModalFlag   = 2,
			  aElemSetFlag = 3 };

	enum Analysis { aStatic = 1,
			aDynamic = 12 };
	
	AbaqusT (ostream& fMainOut, bool binary);

	void ReadVersion (istream& in, StringT& version, StringT& date, StringT& time, int& elements, int& nodes);
	void IsModal (istream& in);
	void ReadElementSets (istream& in, iArray2DT& data);
	void ReadNodeSets (istream& in, iArray2DT& data);

	void ReadElement (istream& in, int& ID, GeometryT::GeometryCode& geocode, iArrayT& nodes);
	void ReadCoordinate (istream& in, int& ID, dArrayT& nodes);

	bool NextTimeIncrement (istream& in, double& time);

	void ReadLabels (istream& in, ArrayT<AbaqusT::VariableKey>& nkeys, ArrayT<AbaqusT::VariableKey>& ekeys, bool& energydata) ;

	void ReadNodeVariables (istream& in, const ArrayT<AbaqusT::VariableKey>& nkeys, dArray2DT& nvalues);
	void ReadElementVariables (istream& in, const ArrayT<AbaqusT::VariableKey>& ekeys, dArray2DT& values, const iArrayT& elements);

	// writing/reading FORTRAN buffer header/footer
	// reset everytime you open a file
	void ResetBufferSize (int amount = 0);
	int GetBufferSize (void) const;

	void Create (ostream& out, int numelems, int numnodes, double elemsize);

	void WriteConnectivity (ostream& out, GeometryT::GeometryCode code, int startnumber, const iArray2DT& connects);
	void WriteCoordinates (ostream& out, const iArrayT& nodes_used, const dArray2DT& coords);
	void WriteElementSet (ostream& out, const StringT& name, const iArrayT& elem_map);
	void WriteNodeSet (ostream& out, const StringT& name, const iArrayT& nodes);
	void WriteActiveDOF (ostream& out, const iArrayT& activedofs);
	void WriteHeading (ostream& out, const StringT& heading);

	void WriteEndIncrement (ostream& out, bool endgeometry);
	void WriteStartIncrement (ostream& out, int step, int inc, double totaltime, double time, double timeincrement, AbaqusT::Analysis analysistype);

	void WriteNodalData (ostream& out, const ArrayT<AbaqusT::VariableKey>& key, const dArray2DT& values, const iArrayT& nodemap, GeometryT::GeometryCode code, int numnodes);
	void WriteElementData (ostream& out, const ArrayT<AbaqusT::VariableKey>& key, const dArray2DT& values);

 private:
  
	enum DoubleFormat { dprecision = 15 };

	enum KeyRecord { aElementHead  = 1,
			 aElement      = 1900,
			 aCoordinates  = 1901,
			 aDOF          = 1902,
			 aOutputDef    = 1911,
			 aVersion      = 1921,
			 aHeading      = 1922,
			 aNodeSet      = 1931,
			 aNodeSetCont  = 1932,
			 aElementSet   = 1933,
			 aElemSetCont  = 1934,
			 aModal        = 1980,
			 aEnergy       = 1999,
			 aIncrement    = 2000,
			 aEndIncrement = 2001 };

	enum ElementData { aQuadrature   = 0,
			   aCentroidal   = 1,
			   aAtNodes      = 2,
			   aRebar        = 3,
			   aNodeAveraged = 4,
			   aWholeElement = 5 };
	
	bool AdvanceStreamTo (istream& in, AbaqusT::KeyRecord key, int& length);
	bool AdvanceStreamTo (istream& in, AbaqusT::KeyRecord key, AbaqusT::KeyRecord keycont, int& length, AbaqusT::KeyRecord& keyfound, bool search);
	bool NextRecord (istream& in, AbaqusT::KeyRecord& key, int& length, int clear = 0);
	bool AdvanceToStar (istream& in) const; // for ASCII
	void ClearLine (istream& in, int n); // for binary

	void Read (istream& in, StringT& s, int n = 1);
	void Read (istream& in, int& i);
	void Read (istream& in, AbaqusT::KeyRecord& i);
	void Read (istream& in, AbaqusT::OutputFlag& i);
	void Read (istream& in, double& d);
	void CheckBufferSize (istream& in);

	bool TranslateElementName (const StringT& name, GeometryT::GeometryCode& code, int& numelemnodes) const;
	bool TranslateContinuum (const char* name, GeometryT::GeometryCode& code, int& numelemnodes) const;
	bool Translate2DName (const char* n, GeometryT::GeometryCode& code, int& numelemnodes) const;
	bool Translate3DName (const char* n, GeometryT::GeometryCode& code, int& numelemnodes) const;

	void CheckBufferSize (ostream& out);
	void StartRecord (ostream& out, int length, int i);
	void WriteASCII (ostream& out, const StringT& s);
	void Write (ostream& out, int i);
	void Write (ostream& out, double d);
	void Write (ostream& out, const StringT& s, int blocks = 1);

	void GetElementName (GeometryT::GeometryCode code, int numnodes, int& numoutnodes, StringT& name) const;

	// write before every single set of one type of variable record
	void WriteOutputDef (ostream& out, int key, const StringT& setname, GeometryT::GeometryCode code, int numnodes, int count, int& numdir, int& numshear);

	// write before every element variable record
	void WriteElementHeader (ostream& out, int key, int number, int integrationpoint, int sectionpoint, AbaqusT::ElementData flag, int numdirect, int numshear, int numdir, int numsecforc);

	// easy groupings of variable records with preliminary data
	void ReadPreliminary (istream& in, int key, int& i, int& clear);
	void WritePreliminary (ostream& out, int key, int i);
	void VariableRecordLength (const ArrayT<AbaqusT::VariableKey>& key, int index, int& length, int& count);
	AbaqusT::VariableKey Record2VariableKey (AbaqusT::KeyRecord key) const;

 private:
	ostream& fOut;    // error messaging
	bool     fBinary; // flag to write binary file
	bool     fModal;
	int      fBufferDone; // FORTRAN binary, length of buffer read/wrote
	                      // ASCII, only write 80 characters per line
	int      fBufferSize; // FORTRAN binary buffer size
};

inline void AbaqusT::ResetBufferSize (int amount) { fBufferDone = amount; }
inline int  AbaqusT::GetBufferSize (void) const { return fBufferDone; }

#endif
