/* $Id: EnSightT.cpp,v 1.19 2004-07-01 16:39:40 paklein Exp $ */
/* created: sawimme (05/13/1999) */
#include "EnSightT.h"

#include <ctype.h>

#include "StringT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "ifstreamT.h"
#include "AutoArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<EnSightT::VariableTypeT>::fByteCopy = true;
} /* namespace Tahoe */

EnSightT::EnSightT (ostream& out, bool binary, int dof) :
fOut(out),
fBinary (binary),
fDOF (dof)
{
}

void EnSightT::WriteHeader (ostream& fgeo, ArrayT<StringT>& header) const
{
for (int i=0; i < header.Length(); i++)
if (fBinary) WriteString (fgeo, header[i].Pointer());
else fgeo << header[i] << '\n';
}

void EnSightT::WritePartInfo (ostream& fgeo, int num, StringT& desc) const
{
if (fBinary)
{
WriteString (fgeo, "part");
fgeo.write (reinterpret_cast<const char *> (&num), sizeof (int));
WriteString (fgeo, desc.Pointer());
}
else
{
fgeo << "part" << '\n';
fgeo << setw(iwidth) << num << '\n';
fgeo << desc << '\n';
}
}

void EnSightT::WriteCoordinateHeader (ostream& fgeo, int numnodes) const
{
if (fBinary)
{
WriteString (fgeo, "coordinates");
fgeo.write (reinterpret_cast<const char *> (&numnodes), sizeof (int));
}
else
fgeo << "coordinates\n" << setw(iwidth) << numnodes << '\n';
}

void EnSightT::WriteCoordinateMap (ostream& fgeo, const iArrayT& nodesmap) const
{
if (fBinary)
{
int itemp;
const int *pnodemap = nodesmap.Pointer();
for (int n=0; n < nodesmap.Length(); n++)
	{
	  itemp = pnodemap[n];
	  fgeo.write (reinterpret_cast<const char *> (&itemp), sizeof (int));
	}
}
else
nodesmap.WriteWithFormat(fgeo, iwidth, 0, 1, 0);
}

void EnSightT::WriteCoordinates (ostream& fgeo, const dArray2DT& coords) const
{
if (fBinary)
{
float temp; // need to convert double to float
for (int i=0; i < coords.MinorDim(); i++)
	for (int j=0; j < coords.MajorDim(); j++)
	  {
	    temp = float(coords(j,i));
	    fgeo.write (reinterpret_cast<const char *> (&temp), sizeof (float));
	  }
}
else
for (int i=0; i < fDOF; i++)
WritedArray2DT (fgeo, coords, i);

Fillto3D (fgeo, 3 - fDOF, coords.MajorDim());
}

int EnSightT::WriteConnectivityHeader (ostream& fgeo, GeometryT::CodeT code, int numelems, int numelemnodes) const
{
StringT elementname;
GetElementName (elementname, numelemnodes, code);

if (fBinary)
{
WriteString (fgeo, elementname.Pointer());
fgeo.write (reinterpret_cast<const char *> (&numelems), sizeof (int));
}
else
fgeo << elementname << '\n' << setw(iwidth) << numelems << '\n';

return numelemnodes;
}

void EnSightT::WriteConnectivityMap (ostream& fgeo, const iArrayT& elementmap) const
{
if (fBinary)
{
int itemp;
for (int nn=0; nn < elementmap.Length(); nn++)
	{
	  itemp = elementmap[nn];
	  fgeo.write (reinterpret_cast<const char *> (&itemp), sizeof (int));
	}
}
else
elementmap.WriteWithFormat (fgeo, iwidth, 0, 1, 0);
}

void EnSightT::WriteConnectivity (ostream& fgeo, GeometryT::CodeT code, int numelemnodes, const iArray2DT& connects) const
{
	/* resolve numbering convention differences between EnSight and tahoe */
	iArray2DT connects_tmp;
	connects_tmp.Alias(connects);
	ConvertElementNumbering(code, connects_tmp);

  /* do not write all columns of connectivity data,
     only write up to numelemnodes */
  if (fBinary)
    {
      for (int ic = 0; ic < connects.MajorDim(); ic++)
	for (int j=0; j < numelemnodes; j++)
	  {
	    int itemp = connects (ic,j);
	    fgeo.write (reinterpret_cast<const char *> (&itemp), sizeof (int));
	  }
    }
  else
    {
      const int *pc = connects.Pointer();
      for (int ic = 0; ic < connects.MajorDim(); ic++)
	{
	  for (int j=0; j < numelemnodes; j++)
	    fgeo << setw (iwidth) << *(pc + j);
	  fgeo << '\n';

	  pc += connects.MinorDim();
	}
    }

	/* convert back since we modified const connects */
	ConvertElementNumbering(code, connects_tmp);
}

void EnSightT::WriteVector (ostream& fvar, const dArray2DT& values, int i) const
{
/* write x-coords, then y-coords, and then z-coords */
if (fBinary)
{
float temp; // need to convert double to float
for (int m=0; m<fDOF; m++)
	for (int k = 0; k < values.MajorDim(); k++)
	  {
	    temp = float (values(k,i+m));
	    fvar.write (reinterpret_cast<const char *> (&temp), sizeof (float));
	  }
}
else
for (int k=0; k < fDOF; k++)
WritedArray2DT (fvar, values, i+k);

/* fill to 3D */
Fillto3D (fvar, 3 - fDOF, values.MajorDim());
}

void EnSightT::WriteScalar (ostream& fvar, const dArray2DT& values, int i) const
{
if (fBinary)
{
float temp; // need to convert double to float
for (int k = 0; k < values.MajorDim(); k++)
	{
	  temp = float (values(k,i));
	  fvar.write (reinterpret_cast<const char *> (&temp), sizeof (float));
	}
}
else
WritedArray2DT (fvar, values, i);
}

void EnSightT::WriteCaseFormat (ostream& fvar) const
{
fvar << "FORMAT\ntype:\t ensight gold\n\n";
}

void EnSightT::WriteCaseGeometry (ostream& fvar, int sequence, StringT& geofile) const
{
	/* remove file path */
	StringT path, file(geofile);
	path.FilePath(file);
	file.Drop(path.StringLength());
	file.ToUNIXPath();
	fvar << "GEOMETRY\nmodel:\t " << sequence << "\t " << file << "\n\n";
}

void EnSightT::WriteVariableLabels (ostream& fvar, const ArrayT<StringT>& labels, const ArrayT<StringT>& filenames, const ArrayT<EnSightT::VariableTypeT>& t) const
{
	const char *type [4] = {"scalar per element", "vector per element",
                      "scalar per node", "vector per node"};

	fvar << "VARIABLE\n";
	for (int i=0; i < labels.Length(); i++)
	{
		/* remove file path */
		StringT path, file(filenames[i]);
		path.FilePath(file);
		file.Drop(path.StringLength());
		file.ToUNIXPath();
		
		fvar << type[t[i]] << ":\t " << labels[i] << "\t " << file << "\n";
	}
	fvar << '\n';
}

void EnSightT::WriteTime (ostream& fvar, int sequence, int start, int increment, const ArrayT<double>& timesteps) const
{
int prec = fvar.precision();
fvar.setf(ios::showpoint);
fvar.setf(ios::right, ios::adjustfield);
fvar.setf(ios::scientific, ios::floatfield);
fvar.precision(dprecision);

fvar << "TIME\n";
fvar << "time set: " << sequence << '\n';
fvar << "number of steps: " << timesteps.Length() << '\n';
fvar << "filename start number: " << start << " \n";
fvar << "filename increment: " << increment << " \n";
fvar << "time values:\n";

for (int t=0; t < timesteps.Length(); t++)
{
fvar << setw(12) << timesteps[t];
if ( (t+1)%6 == 0) fvar << '\n';
}
fvar << '\n';

/* reset */
fvar.precision(prec);
}

void EnSightT::GetElementName (StringT& fElementName, int& num_output_nodes, GeometryT::CodeT geocode) const
{
int num_nodes = num_output_nodes;
switch (geocode)
{
case GeometryT::kPoint:
fElementName = "point";
num_output_nodes = 1;
//num_output_nodes not appended for this case
break;

case GeometryT::kLine:
fElementName = "bar";
if (num_nodes > 3 || num_nodes < 2)
	{
	  cout << "\nEnSightT::GetElementName cannot do bar with fNumberElementNodes = "
	       << num_nodes << "\n\n";
	  fOut << "\nEnSightT::GetElementName cannot do bar with fNumberElementNodes = "
	       << num_nodes << "\n\n";
	  throw ExceptionT::kGeneralFail;
	}
num_output_nodes = num_nodes;
fElementName.Append (num_output_nodes);
break;

case GeometryT::kQuadrilateral:
fElementName = "quad";
num_output_nodes = (num_nodes < 8) ? 4 : 8;
fElementName.Append(num_output_nodes);
break;	
	
case GeometryT::kTriangle:
fElementName = "tria";
num_output_nodes = (num_nodes < 6) ? 3 : 6;
fElementName.Append(num_output_nodes);
break;	

case GeometryT::kHexahedron:
fElementName = "hexa";
num_output_nodes = (num_nodes < 20) ? 8 : 20;
fElementName.Append(num_output_nodes);
break;	

case GeometryT::kTetrahedron:
fElementName = "tetra";
num_output_nodes = (num_nodes < 10) ? 4 : 10;
fElementName.Append(num_output_nodes);
break;	
		
case GeometryT::kPentahedron:
fElementName = "penta";
num_output_nodes = (num_nodes < 15) ? 6 : 15;
fElementName.Append(num_output_nodes);
break;

default:
cout << "\nEnSightT::GetElementName cannot find name for Geometry CodeT "
	   << geocode << "\n\n";
fOut << "\nEnSightT::GetElementName cannot find name for Geometry CodeT "
	   << geocode << "\n\n";
throw ExceptionT::kGeneralFail;
}
}

bool EnSightT::CaseFile (ifstreamT& in, StringT& geofile) const
{
// verify it is a case file
StringT word;
in >> word;
if (word != "FORMAT") return false;

// verify it is a gold format
in >> word; // geometry
in >> word; // ensight
in >> word; // gold
if (word != "gold" && word != "GOLD") return false;

// geometry file
in >> word;
if (word != "GEOMETRY") return false;
in >> word;
int timeset;
if (isdigit (in.next_char())) in >> timeset;
in >> geofile;

// do not read changing geometry at this point
if (strstr (geofile.Pointer(), "*") != NULL) return false;

return true;
}

// if filename == true, then filenames are returned instead of labels
bool EnSightT::ReadVariableSection (ifstreamT& in, AutoArrayT<StringT>& nlabels, AutoArrayT<StringT>& elabels, AutoArrayT<bool>& nvector, AutoArrayT<bool>& evector, bool filename) const
{
StringT per, name, vector, node, label;
int timeset;

in >> vector;
while (vector != "TIME" && in.good())
{
in >> per >> node;
if (isdigit (in.next_char())) in >> timeset;
in >> label >> name;

if (vector == "vector")
	{
	  for (int j=0; j < fDOF; j++)
	    {
	      StringT temp = label;
	      if (filename)
		temp = name;
	      else
		temp.Append ("_", j+1);

	      if (node == "node:")
		{
		  nlabels.Append (temp);
		  nvector.Append (true);
		}
	      else
		{
		  elabels.Append (temp);
		  evector.Append (true);
		}
	    }
	}
else if (vector == "scalar")
	{
	  StringT temp = label;
	  if (filename)
	    temp = name;

	  if (node == "node:")
	    {
	      nlabels.Append (temp);
	      nvector.Append (false);
	    }
	  else
	    {
	      elabels.Append (temp);
	      evector.Append (false);
	    }	
	}
else
	return false;

// read next line
in >> vector;
}
return true;
}

int EnSightT::NumTimeSteps (ifstreamT& in) const
{
  char c = 'c';
  int timeset, num;

  while (c != ':') in >> c;
  in >> timeset;
  c = 'c';
  
  while (c != ':') in >> c;
  in >> num;

  return num;
}

bool EnSightT::ReadTimeSection (ifstreamT& in, int& start, int& increment, dArrayT& timesteps) const
{
char c = 'c';
int timeset, num;

while (c != ':') in >> c;
in >> timeset;
c = 'c';

while (c != ':') in >> c;
in >> num;
c = 'c';

while (c != ':') in >> c;
in >> start;
c = 'c';

while (c != ':') in >> c;
in >> increment;
c = 'c';

while (c != ':') in >> c;
timesteps.Allocate (num);
in >> timesteps;
return true;
}

void EnSightT::ReadGeometryHeader (istream& in, bool& nodemap, bool& elemmap) const
{
nodemap = false;
elemmap = false;

StringT line (81);
if (fBinary)
{
in.read (line.Pointer(), sizeof (char)*80);
if (strncmp (line.Pointer(), "C Binary", 8) != 0)
	{
	  cout << "\n\nEnSight can only read C Binary\n";
	  throw ExceptionT::kGeneralFail;
	}

in.read (line.Pointer(), sizeof (char)*80); // header 1
in.read (line.Pointer(), sizeof (char)*80); // header 2
in.read (line.Pointer(), sizeof (char)*80); // node id
if (strncmp (line, "node id given", 12) == 0) nodemap = true;
in.read (line.Pointer(), sizeof (char)*80); // element id
if (strncmp (line, "element id given", 12) == 0) elemmap = true;
}
else
{
in.getline (line.Pointer(), 80, '\n'); // header 1
in.getline (line.Pointer(), 80, '\n'); // header 2
in.getline (line.Pointer(), 80, '\n'); // node id
if (strncmp (line.Pointer(), "node id given", 12) == 0) nodemap = true;
in.getline (line.Pointer(), 80, '\n'); // element id
if (strncmp (line.Pointer(), "element id given", 12) == 0) elemmap = true;
}
}

bool EnSightT::ReadPart (istream& in, int& partID) const
{
StringT line (81);
if (fBinary)
{
in.read (line.Pointer(), sizeof (char)*80); // part
if (strncmp (line.Pointer(), "part", 4) != 0) return false;
in.read (reinterpret_cast<char *> (&partID), sizeof (int));
in.read (line.Pointer(), sizeof (char)*80); // description
}
else
{
in.getline (line.Pointer(), 80, '\n'); // part
if (strncmp (line.Pointer(), "part", 4) != 0) return false;
in >> partID;
in.getline (line.Pointer(), 80, '\n'); // clear line
in.getline (line.Pointer(), 80, '\n'); // description
}

return true;
}

void EnSightT::SkipPart (istream& in, bool nodemapgiven, bool elemmapgiven, int& num_nodes, int& num_elems, int& num_elem_nodes) const
{
SkipCoordinates (in, num_nodes, nodemapgiven);
SkipConnectivity (in, num_elems, num_elem_nodes, elemmapgiven);
}

void EnSightT::SkipCoordinates (istream& in, int& num_nodes, bool nodemapgiven) const
{
StringT line (81);
int idum;
float fdum;
if (fBinary)
{
in.read (line.Pointer(), sizeof (char)*80); // coordinates
in.read (reinterpret_cast<char *> (&num_nodes), sizeof (int));
if (nodemapgiven)
	for (int i=0; i < num_nodes; i++)
	  in.read (reinterpret_cast<char *> (&idum), sizeof (int));
for (int j=0; j < num_nodes*fDOF; j++)
	in.read (reinterpret_cast<char *> (&fdum), sizeof (float));
}
else
{
in.getline (line.Pointer(), 80, '\n');
in >> num_nodes;
if (nodemapgiven)
	for (int i=0; i < num_nodes; i++)
	  in >> idum;
for (int j=0; j < num_nodes*fDOF; j++)
	in >> fdum;

in.getline (line.Pointer(), 80, '\n'); // clear endline character
}
}

void EnSightT::SkipConnectivity (istream& in, int& num_elems, int& num_elem_nodes, bool elemmapgiven) const
{
StringT line (81);
int idum;
GeometryT::CodeT code;
if (fBinary)
{
in.read (line.Pointer(), sizeof (char)*80); // element name
if (!GeometryCode (code, line, num_elem_nodes)) return;
in.read (reinterpret_cast<char *> (&num_elems), sizeof (int));
if (elemmapgiven)
	for (int k=0; k < num_elems; k++)
	  in.read (reinterpret_cast<char *> (&idum), sizeof (int));
for (int l=0; l < num_elems*num_elem_nodes; l++)
	in.read (reinterpret_cast<char *> (&idum), sizeof (int));
}
else
{
in.getline (line.Pointer(), 80, '\n'); // element name
if (!GeometryCode (code, line, num_elem_nodes)) return;
in >> num_elems;
if (elemmapgiven)
	for (int k=0; k < num_elems; k++)
	  in >> idum;
for (int l=0; l < num_elems*num_elem_nodes; l++)
	in >> idum;

in.getline (line.Pointer(), 80, '\n'); // clear endline character
}
}

void EnSightT::ReadCoordinates (istream& in, dArray2DT& coords, iArrayT& map, bool nodemapgiven) const
{
StringT line (81);
int num_nodes;
if (fBinary)
{
in.read (line.Pointer(), sizeof (char)*80); // coordinates
in.read (reinterpret_cast<char *> (&num_nodes), sizeof (int));
if (nodemapgiven)
	{
	  map.Allocate (num_nodes);
	  int *pm = map.Pointer();
	  for (int i=0; i < num_nodes; i++)
	    in.read (reinterpret_cast<char *> (pm++), sizeof (int));
	}
dArray2DT temp (fDOF, num_nodes);
double *pd = temp.Pointer();
float pf;
for (int j=0; j < num_nodes*fDOF; j++)
	{
	  in.read (reinterpret_cast<char *> (&pf), sizeof (float));
	  *pd++ = pf; // must convert to float separately from converting to double
	}
coords.Allocate (num_nodes, fDOF);
coords.Transpose (temp);
}
else
{
in.getline (line.Pointer(), 80, '\n');
in >> num_nodes;
if (nodemapgiven)
	{
	  map.Allocate (num_nodes);
	  in >> map;
	}
dArray2DT temp (fDOF, num_nodes);
in >> temp;
coords.Allocate (num_nodes, fDOF);
coords.Transpose (temp);

in.getline (line.Pointer(), 80, '\n'); // clear endline character
}
}

void EnSightT::ReadConnectivity (istream& in, iArray2DT& conn, iArrayT& map, bool elemmapgiven, GeometryT::CodeT& code) const
{
	StringT line (81);
	int num_elems, num_elem_nodes;
	if (fBinary)
	{
		in.read (line.Pointer(), sizeof (char)*80); // element name
		if (!GeometryCode (code, line, num_elem_nodes)) return;
		in.read (reinterpret_cast<char *> (&num_elems), sizeof (int));
		if (elemmapgiven)
		{
			map.Dimension (num_elems);
			int *mp = map.Pointer();
			for (int k=0; k < num_elems; k++)
			in.read (reinterpret_cast<char *> (mp++), sizeof (int));
		}

		conn.Dimension (num_elems, num_elem_nodes);
		int *cp = conn.Pointer();
		for (int l=0; l < num_elems*num_elem_nodes; l++)
			in.read (reinterpret_cast<char *> (cp++), sizeof (int));
	}
	else
	{
		in.getline (line.Pointer(), 80, '\n'); // element name
		if (!GeometryCode (code, line, num_elem_nodes)) return;
		in >> num_elems;
		if (elemmapgiven)
		{
			map.Dimension (num_elems);
			in >> map;
		}
		
		conn.Dimension (num_elems, num_elem_nodes);
		in >> conn;

		in.getline (line.Pointer(), 80, '\n'); // clear endline character
	}

	/* translate numbering conventions */
	ConvertElementNumbering(code, conn);
}

void EnSightT::ReadVariableHeader (istream& in, StringT& header) const
{
header.Allocate (81);
if (fBinary)
in.read (header.Pointer(), sizeof (char)*80); // description line
else
in.getline (header.Pointer(), 80, '\n'); // description line
}

void EnSightT::ReadVariable (istream& in, dArray2DT& values) const
{
dArray2DT temp (values.MinorDim(), values.MajorDim());
if (fBinary)
{
float pf;
double *dt = temp.Pointer();
for (int k=0; k < values.Length(); k++)
	{
	  in.read (reinterpret_cast<char *> (&pf), sizeof (float));
	  *dt++ = pf;
	}
}
else
{
in >> temp;
StringT line (81);
in.getline (line.Pointer(), 80, '\n'); // clear endline character
}

values.Transpose (temp);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* translate element number convention between EnSight and tahoe */
void EnSightT::ConvertElementNumbering(GeometryT::CodeT code, iArray2DT& conn) const
{
	/* 3-noded line elements */
	if (code == GeometryT::kLine && conn.MinorDim() == 3)
		for (int i = 0; i < conn.MajorDim(); i++)
		{
			/* swap last two nodes */
			int tmp = conn(i,1);
			conn(i,1) = conn(i,2);
			conn(i,2) = tmp;
		}
}

void EnSightT::WritedArray2DT(ostream& out, const dArray2DT& values, int column_position) const
{
/* store */
int prec = out.precision();
out.setf(ios::showpoint);
out.setf(ios::right, ios::adjustfield);
out.setf(ios::scientific, ios::floatfield);
out.precision(dprecision);
const double* pv = values.Pointer() + column_position;
for (int i=0; i < values.MajorDim(); i++)
{
out << setw(dwidth) << *pv << '\n';
pv += values.MinorDim();
}

/* reset */
out.precision(prec);
}

void EnSightT::Fillto3D (ostream& out, int width, int length) const
{
if (fBinary)
{
float temp = 0;
for (int j=0; j < width*length; j++)
	out.write (reinterpret_cast<const char *> (&temp), sizeof (float));
}
else
{
int prec = out.precision();
out.precision(dprecision);
for (int i=0; i < width*length; i++)
	out << setw(dwidth) << 0.0 << '\n';
out.precision(prec);
}
}

void EnSightT::WriteString (ostream& o, const char* s) const
{
char buffer [80];
strncpy (buffer, s, 80);
o.write (buffer, sizeof (char)*80);
}

bool EnSightT::GeometryCode(GeometryT::CodeT& code, const StringT& name, int& numelemnodes) const
{
if (strncmp (name, "point", 3) == 0)
{
code = GeometryT::kPoint;
numelemnodes = 1;
return true;
}
else if (strncmp (name, "tria", 3) == 0)
{
code = GeometryT::kTriangle;
numelemnodes = (name[4] == '3') ? 3 : 6;
return true;
}
else if (strncmp (name, "quad", 3) == 0)
{
code = GeometryT::kQuadrilateral;
numelemnodes = (name[4] == '4') ? 4 : 8;
return true;
}
else if (strncmp (name, "tetra", 3) == 0)
{
code = GeometryT::kTetrahedron;
numelemnodes = (name[5] == '4') ? 4 : 10;
return true;
}
else if (strncmp (name, "hexa", 3) == 0)
{
code = GeometryT::kHexahedron;
numelemnodes = (name[4] == '8') ? 8 : 20;
return true;
}
else if (strncmp (name, "penta", 3) == 0)
{
code = GeometryT::kPentahedron;
numelemnodes = (name[5] == '6') ? 6 : 15;
return true;
}

cout << "\n\nEnSightT::GeomtryCode element type unsupported\n";
cout << name << "." << endl;
return false;
}

