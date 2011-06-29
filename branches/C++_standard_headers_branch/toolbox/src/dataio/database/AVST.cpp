/* created: saw (05.10.2001)                                              */

#include "AVST.h"
#include "ios_fwd_decl.h"
#include "ArrayT.h"
#include "StringT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"


using namespace Tahoe;

AVST::AVST (ostream& out, bool binary) :
  fOut (out),
  fBinary (binary)
{
}

void AVST::WriteHeader (ostream &out, int numpoints, int numcells, int numnodedata, int numcelldata, int nummodeldata) const
{
  if (fBinary)
    {
      fOut << "\nAVST::WriteHeader not programed for Binary\n\n";
      cout << "\nAVST::WriteHeader not programed for Binary\n\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    {
      out << setw(iwidth) << numpoints << setw(iwidth) << numcells;
      out << setw(iwidth) << numnodedata << setw(iwidth) << numcelldata;
      out << setw(iwidth) << nummodeldata << '\n';
    }
}

void AVST::WriteCoordinates (ostream &out, dArray2DT& coords, int firstID) const
{
  if (fBinary)
    {
      fOut << "\nAVST::WriteCoordinates not programed for Binary\n\n";
      cout << "\nAVST::WriteCoordinates not programed for Binary\n\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    WriteArray2DT (out, coords, firstID);
}

void AVST::WriteCells (ostream &out, GeometryT::CodeT code, iArray2DT& connects, int matid, int firstID) const
{
  int numnodes = connects.MinorDim();
  StringT name;
  GetElementName (code, name, numnodes);
  if (fBinary)
    {
      fOut << "\nAVST::WriteCellData not programed for Binary\n\n";
      cout << "\nAVST::WriteCellData not programed for Binary\n\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    {
      for (int i=0, id=firstID; i < connects.MajorDim(); i++, id++)
	{
	  out << setw(iwidth) << id << setw(iwidth) << matid;
	  out << " " << name << " ";
	  int *pc = connects(i);
	  switch (code)
	    {
	    case GeometryT::kHexahedron:
	      WriteBackward (out, pc, numnodes, 4);
	      break;
	    case GeometryT::kPentahedron:
	      WriteBackward (out, pc, numnodes, 3);
	      break;
	    default:
	      WriteForward (out, pc, numnodes);
	    }
	}
    }
}

void AVST::WriteDataHeader (ostream &out, const ArrayT<StringT>& labels) const
{
  StringT units (" ");
  if (fBinary)
    {
      fOut << "\nAVST::WriteDataHeader not programed for Binary\n\n";
      cout << "\nAVST::WriteDataHeader not programed for Binary\n\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    {
      /* write all variable as scalar for now */
      out << labels.Length(); /* number of components */
      for (int i=0; i < labels.Length(); i++)
	out << setw(iwidth) << 1; /* size of each component */
      out << '\n';
      for (int j=0; j < labels.Length(); j++)
	out << labels[j] << ", " << units << '\n';
    }
}

void AVST::WriteDataHeader (ostream &out, const ArrayT<StringT>& labels, const ArrayT<int>& dimensions) const
{
  StringT units (" ");
  if (fBinary)
    {
      fOut << "\nAVST::WriteDataHeader not programed for Binary\n\n";
      cout << "\nAVST::WriteDataHeader not programed for Binary\n\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    {
      /* write all variable as scalar for now */
      out << labels.Length(); /* number of components */
      for (int i=0; i < labels.Length(); i++)
	out << setw(iwidth) << dimensions[i]; /* size of each component */
      out << '\n';
      for (int j=0; j < labels.Length(); j++)
	out << labels[j] << ", " << units << '\n';
    }
}

void AVST::WriteData (ostream &out, const dArray2DT& data, int firstID) const
{
  if (fBinary)
    {
      fOut << "\nAVST::WriteData not programed for Binary\n\n";
      cout << "\nAVST::WriteData not programed for Binary\n\n";
      throw ExceptionT::kGeneralFail;
    }
  else
    WriteArray2DT (out, data, firstID);
}

/*************************************************************************
* Private
*************************************************************************/

void AVST::GetElementName (GeometryT::CodeT code, StringT& name, int &numnodes) const
{
  switch (code)
    {
    case GeometryT::kPoint:
      {
	name = "pt";
	numnodes = 1;
	break;
      }
    case GeometryT::kQuadrilateral:
      {
	name = "quad";
	numnodes = (numnodes < 8) ? 4 : 8;
	break;
      }
    case GeometryT::kTriangle:
      {
	name = "tri";
	numnodes = (numnodes < 6) ? 3 : 6;
	break;
      }
    case GeometryT::kHexahedron:
      {
	name = "hex";
	numnodes = (numnodes < 20) ? 8 : 20;
	break;
      }
    case GeometryT::kTetrahedron:
      {
	name = "tet";
	numnodes = (numnodes < 10) ? 4 : 10;
	break;
      }
    case GeometryT::kPentahedron:
      {
	name = "prism";
	numnodes = (numnodes < 6) ? 6 : 15;
	break;
      }
    default:
      {
	cout << "\n\n AVST::GetElementName unknown Geometry CodeT "
	     << code << "\n\n";
	fOut << "\n\n AVST::GetElementName unknown Geometry CodeT "
	      << code << "\n\n";
	throw ExceptionT::kGeneralFail;
      }
    }
}

void AVST::WriteForward (ostream &out, int *pc, int numnodes) const
{
  for (int i=0; i < numnodes; i++)
    out << setw (iwidth) << *(pc + i);
  out << '\n';
}

void AVST::WriteBackward (ostream &out, int *pc, int numnodes, int numfacenodes) const
{
  int z = numnodes / numfacenodes;
  for (int j=1, k = numfacenodes - 1; j < z; j++, k += numfacenodes)
    for (int i=0; i < numfacenodes; i++)
      out << setw (iwidth) << pc [k - i];
  out << '\n';
}

void AVST::WriteArray2DT (ostream &out, const dArray2DT& data, int firstID) const
{
  int prec = fOut.precision();
  out.precision(dprecision);
  out.setf (ios::scientific);
  const double *pc = data.Pointer();
  for (int i=0, id=firstID; i < data.MajorDim(); i++, id++)
    {
      out << setw (iwidth) << id;
      for (int j=0; j < data.MinorDim(); j++, pc++)
	out << setw(dwidth) << *pc;
      for (int k = data.MinorDim(); k < 3; k++)
	out << setw(dwidth) << 0.0;
      out << '\n';
    }
  out.precision(prec);
}
