/* $Id: ParaDynT.cpp,v 1.8 2003-11-04 01:21:10 paklein Exp $ */
#include "ParaDynT.h"

#include <ctype.h>

#include "StringT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "fstreamT.h"
#include "ios_fwd_decl.h"
#include "AutoArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe  {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParaDynT::VariableTypeT>::fByteCopy = true;
} /* namespace Tahoe */

ParaDynT::ParaDynT (ostream& out) : fOut (out)
{
}

void ParaDynT::WriteHeader (ostream& fgeo, ArrayT<StringT>& header) const
{
  for (int i=0; i < header.Length(); i++)
    fgeo << header[i] << '\n';
}

void ParaDynT::WriteCoordinateHeader (ostream& fgeo) const
{
  fgeo << " ITEM: ATOMS" << '\n';
}

/*
void ParaDynT::WriteCoordinates (ostream& fgeo, 
				 const dArray2DT& coords,
				 const iArrayT& types,
				 const iArrayT& parts)  const
{
  if(coords.MajorDim() != types.Length()) 
    throw ExceptionT::kSizeMismatch;

  if(coords.MinorDim()==2)
    { 
      for (int i=0; i < coords.MajorDim(); i++)
	fgeo << i+1  << "  "  << types[i] << "  " 
	     << float(coords(i)[0]) << "  " 
             << float(coords(i)[1]) << "  "
             << 0.0 << "\n";
    }
  else if(coords.MinorDim()==3)
    { 
      for (int i=0; i < coords.MajorDim(); i++)
	fgeo << i+1  << "  " << types[i] << "  " 
	     << float(coords(i)[0]) << "  " 
	     << float(coords(i)[1]) << "  "
	     << float(coords(i)[2]) << "\n";
    }
  else
    throw ExceptionT::kBadInputValue;
 
}
*/

/* -> For Sylvie's version of ParaDyn.....*/
void ParaDynT::WriteCoordinates (ostream& fgeo, 
				 const dArray2DT& coords,
				 const iArrayT& types,
				 const iArrayT& parts) const
{
  if(coords.MinorDim()==2)
    { 
      for (int i=0; i < coords.MajorDim(); i++)
	fgeo << i+1  << "  "  << types[i] << "  " 
	     << float(coords(i)[0]) << "  " 
             << float(coords(i)[1]) << "  " 
             << float(0.0) << "  "
	     << parts[i] << "\n";
    }
  else if(coords.MinorDim()==3)
    { 
      for (int i=0; i < coords.MajorDim(); i++)
	fgeo << i+1  << "  " << types[i] << "  " 
	     << float(coords(i)[0]) << "  " 
	     << float(coords(i)[1]) << "  "
	     << float(coords(i)[2]) << "  " 
	     << parts[i] << "  " 
	     << "\n";
    }
  else
    throw ExceptionT::kBadInputValue;
 
}


void ParaDynT::WriteTime (ostream& fvar) const
{
  fvar << " ITEM: TIME" << '\n';
  fvar << 0.0 << "\n";
}


void ParaDynT::WriteBoundHeader (ostream& fgeo) const
{
  fgeo << " ITEM: BOX BOUNDS" << '\n';
}


void ParaDynT::WriteBounds (ostream& fgeo, const dArray2DT& bounds) const
{
  if(bounds.MinorDim()==2)
    { 
      for (int i=0; i < bounds.MajorDim(); i++)
	fgeo << float(bounds(i)[0]) << "  " 
             << float(bounds(i)[1]) << "\n";
      if (bounds.MajorDim()==2)
        fgeo << float(-10000.) << "  " 
           << float(10000.) << "\n";
    }
  else if(bounds.MinorDim()==3)
    { 
      for (int i=0; i < bounds.MajorDim(); i++)
	fgeo << float(bounds(i)[0]) << "  " 
	     << float(bounds(i)[1]) << "  "
	     << float(bounds(i)[2]) << "\n";
    }
  else
    {
      cout << "Error: ParaDynT::WriteBounds\n";
      throw ExceptionT::kBadInputValue;
    }
}
