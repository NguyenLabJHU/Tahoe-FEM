/* $Id: TecPlotT.cpp,v 1.7 2003-11-21 22:41:46 paklein Exp $ */
/* created: saw (06.06.2000)                                              */
/* version 7.5                                                            */
/* rules:                                                                 */
/* 1. Maximum 32,766 zone records                                         */
/* 2. Maximum 10 custom label records (does not refer to variable labels) */
/*    (Custom labels are for axis labels, text tick marks, etc ...)       */
/*    Custom labels are not currently used.                               */
/* 3. Maximum ascii line length is 4000 characters                        */
/* 4. Don't wrap amid a "" string                                         */

#include "TecPlotT.h"
#include "ios_fwd_decl.h"
#include "StringT.h"
#include "iArray2DT.h"
#include "iArrayT.h"
#include "dArray2DT.h"


using namespace Tahoe;

TecPlotT::TecPlotT (ostream& out, bool point) :
fOut (out),
fPoint (point)
{
}

void TecPlotT::WriteHeader (ostream& out, const StringT& title, const ArrayT<StringT>& variablenames) const
{
if (title.Length() > 1)
out << "TITLE = \"" << title << "\"\n";

if (variablenames.Length() > 0)
out << "VARIABLES = ";
for (int i=0; i < variablenames.Length(); i++)
out << "\"" << variablenames[i] << "\"\n";
}

void TecPlotT::WriteIJKZone (ostream& out, const StringT& title, const iArrayT& ijk) const
{
out << "ZONE ";

if (title.Length() > 1)
out << "T=\"" << title << "\" ";

char I = 'I';
for (int i=0; i < ijk.Length() && i < 3; i++)
out << I++ << "=" << ijk[i] << " ";

if (fPoint)
out << "F=POINT\n";
else
out << "F=BLOCK\n";
}

// must use WriteConnecitivity with this
void TecPlotT::WriteFEZone (ostream& out, const StringT& title, int numnodes, int numelems, GeometryT::CodeT code, bool connectivity) const
{
out << "ZONE ";

// title
if (title.Length() > 1)
out << "T=\"" << title << "\" ";

// dimensions
out << "N=" << numnodes << " ";
if (numelems > 0)
out << "E=" << numelems << " ";

// format
if (fPoint)
out << "F=FEPOINT ";
else
out << "F=FEBLOCK ";

// element type
out << "ET=";
switch (code)
{
case GeometryT::kTriangle: out << "TRIANGLE "; break;
case GeometryT::kQuadrilateral: out << "QUADRILATERAL "; break;
case GeometryT::kHexahedron: out << "BRICK "; break;
case GeometryT::kTetrahedron: out << "TETRAHEDRON "; break;
default:
{
	cout << "\n\nTecPlotT::WriteFEZone, unknown geometry code \"" << GeometryT::ToString(code) << "\" ("
	     << code << ')' << endl;
	throw ExceptionT::kGeneralFail;
}
}

// duplicate
if (!connectivity)
out << "D=(FECONNECT) ";

out << '\n';
}

// write data can only be call once if using point format
// but may be called repeatly, in proper order, for block format
void TecPlotT::WriteData (ostream& out, const dArray2DT& data) const
{
if (fPoint)
out << data << '\n';
else
{
dArray2DT temp (data.MinorDim(), data.MajorDim());
temp.Transpose (data);

// keep row length under 4000 characters
double *pt = temp.Pointer();
for (int i=0; i < temp.Length(); i++)
	{
	  out << *pt++ << " ";
	  if ((i+1)%100 == 0 || i == temp.Length() - 1) out << '\n';
	}
}
}

// write data should not be called using the  point format
// but may be called repeatly, in proper order, for block format
void TecPlotT::WriteData (ostream& out, const ArrayT<double>& data, const int rows, const int cols) const
{
  if (fPoint)
    {
      fOut << "\n\nTecPlot::WriteData, This function should not be used to write Point format data\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    {
      for (int j=0; j < cols; j++)
	{
	  const double *p = data.Pointer(j);
	  for (int k=0; k < rows; k++)
	    {
	      out << *p << " ";
	      // keep row length under 4000 characters
	      if ((k+1)%100 == 0) 
		out << '\n';
	      p += cols;
	    }
	  out << '\n';
	}
    }
}

// only used with WriteFEZone
void TecPlotT::WriteConnectivity (ostream& out, GeometryT::CodeT code, const iArray2DT& connects) const
{
int numnodes = 0;
switch (code)
{
case GeometryT::kTriangle: numnodes = 3; break;
case GeometryT::kQuadrilateral: numnodes = 4; break;
case GeometryT::kHexahedron: numnodes = 8; break;
case GeometryT::kTetrahedron: numnodes = 4; break;
default:
{
	cout << "\n\nTecPlotT::WriteConnectivity, unknown geometry code "
	     << code << endl;
	throw ExceptionT::kGeneralFail;
}
}

for (int i=0; i < connects.MajorDim(); i++)
connects.PrintRow (i, numnodes, out);

out << '\n';
}
