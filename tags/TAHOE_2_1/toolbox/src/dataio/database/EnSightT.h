/* $Id: EnSightT.h,v 1.6 2004-07-01 16:39:40 paklein Exp $ */
/* created: sawimme (05/13/1999) */
#ifndef _ENSIGHT_T_H_
#define _ENSIGHT_T_H_

/* direct members */
#include "GeometryT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;
class ifstreamT;
class StringT;
class iArrayT;
class dArrayT;
class iArray2DT;
class dArray2DT;
template <class TYPE> class AutoArrayT;

/** interface for I/O of EnSight6 Gold format */
class EnSightT
{
public:
enum VariableTypeT {kScalarElemental,
		              kVectorElemental,
		              kScalarNodal,
		              kVectorNodal};

	EnSightT (ostream& fMainOut, bool binary, int dof);

	// geometry
	void WriteHeader (ostream& fgeo, ArrayT<StringT>& header) const;
	void WritePartInfo (ostream& fgeo, int num, StringT& desc) const;

	void WriteCoordinateHeader (ostream& fgeo, int numnodes) const;
	void WriteCoordinates (ostream& fgeo, const dArray2DT& coords) const;
	void WriteCoordinateMap (ostream& fgeo, const iArrayT& nodemap) const;

	// returns num_elem_nodes
	int WriteConnectivityHeader (ostream& fgeo, GeometryT::CodeT code, int numelems, int numelemnodes) const;
	void WriteConnectivityMap (ostream& fgeo, const iArrayT& elementmap) const;
	void WriteConnectivity (ostream& fgeo, GeometryT::CodeT code, int numelemnodes, const iArray2DT& connects) const;

	// variables
	void WriteVector (ostream& fvar, const dArray2DT& values, int i = 0) const;
	void WriteScalar (ostream& fvar, const dArray2DT& values, int i = 0) const;

	// case file
	void WriteCaseFormat (ostream& fvar) const;
	void WriteCaseGeometry (ostream& fvar, int sequence, StringT& geofile) const;
	void WriteVariableLabels (ostream& fvar, const ArrayT<StringT>& labels, const ArrayT<StringT>& filenames, const ArrayT<VariableTypeT>& t) const;
	void WriteTime (ostream& fvar, int sequence, int start, int increment, const ArrayT<double>& timesteps) const;

	void GetElementName (StringT& fElementName, int& num_output_nodes, GeometryT::CodeT geocode) const;

	// read case file
	bool CaseFile (ifstreamT& in, StringT& geofile) const;
	bool ReadVariableSection (ifstreamT& in, AutoArrayT<StringT>& nlabels, AutoArrayT<StringT>& elabels, AutoArrayT<bool>& nvector, AutoArrayT<bool>& evector, bool filename) const;
	int NumTimeSteps (ifstreamT& in) const;
	bool ReadTimeSection (ifstreamT& in, int& start, int& increment, dArrayT& timesteps) const;

	// read or skip geometry file sections
	void ReadGeometryHeader (istream& in, bool& nodemap, bool& elemmap) const;
	bool ReadPart (istream& in, int& partID) const;
	void SkipPart (istream& in, bool nodemapgiven, bool elemmapgiven, int& num_nodes, int& num_elems, int& num_elem_nodes) const;
	void SkipCoordinates (istream& in, int& num_nodes, bool nodemapgiven) const;
	  void SkipConnectivity (istream& in, int& num_elems, int& num_elem_nodes, bool elemmapgiven) const;

	void ReadCoordinates (istream& in, dArray2DT& coords, iArrayT& map, bool mapgiven) const;
	void ReadConnectivity (istream& in, iArray2DT& conn, iArrayT& map, bool mapgiven, GeometryT::CodeT& code) const;

	// returns values from variables file in array for a given part ID
	void ReadVariableHeader (istream& in, StringT& header) const;
	void ReadVariable (istream& in, dArray2DT& values) const;

private:

	/* translate element number convention between EnSight and tahoe */
	void ConvertElementNumbering(GeometryT::CodeT code, iArray2DT& conn) const;

void WritedArray2DT(ostream& out, const dArray2DT& values, int column_position) const;
void Fillto3D (ostream& out, int width, int length) const;
	void WriteString (ostream& o, const char* s) const;
	bool GeometryCode (GeometryT::CodeT& code, const StringT& name, int& num_elems) const;

	enum WidthsT { iwidth = 10, dwidth = 12, dprecision = 5};

private:
	ostream& fOut;    // error messaging
	bool     fBinary; // flag to write binary file
	int      fDOF;
};

} // namespace Tahoe 
#endif
