// $Id: ParameterFileManagerT.cpp,v 1.3 2002-10-25 21:02:59 paklein Exp $
#include "ParameterFileManagerT.h"
#include "ExceptionCodes.h"
#include "ifstreamT.h"

using namespace Tahoe;

ParameterFileManagerT::ParameterFileManagerT (const StringT& infile) :
  MakeCSE_IOManager (),
  fInFile (infile)
{
}

void ParameterFileManagerT::Initialize (void)
{
  ifstreamT in ('#', fInFile);
  if (!in.is_open())
    {
      cout << "ParameterFileManagerT::Initialize: Unable to open file: "
	   << fInFile << endl;
      throw eGeneralFail;
    }


  /* make sure the file has a *EOF, since we are reading an unknown number
   of StringT values after *KEYWORD, operator>>(StringT) cannot indicate 
   finding nothing, so mark the end of the file with *EOF for now */
  if (!AdvanceTo (in, "*EOF"))
    {
      cout << "\n *EOF not found in file.   *EOF is being added. \n\n";
      in.close ();

      ofstreamT o (fInFile, ios::app);
      o << "\n\n*EOF\n";
    }
}

void ParameterFileManagerT::InputFormat (IOBaseT::FileTypeT &f, StringT& s)
{
  ifstreamT in ('#', fInFile);
  if (!AdvanceTo (in, "*INPUT"))
    {
      cout << "\nParameterFileManagerT::InputFormat: No *INPUT in file.\n\n";
      throw eGeneralFail;
    }
  in >> f >> s;
}

void ParameterFileManagerT::OutputFormat (IOBaseT::FileTypeT &f, StringT& s)
{
  ifstreamT in ('#', fInFile);
  if (!AdvanceTo (in, "*OUTPUT"))
    {
      cout << "\nParameterFileManagerT::OutputFormat: No *OUTPUT in file.\n";
      cout << " Using Tahoe II with a file name of defaultoutput. \n\n";
      f = IOBaseT::kTahoeII;
      s = "defaultoutput";
    }
  else
    in >> f >> s;
}

bool ParameterFileManagerT::Verbose (void)
{
  bool b = false;
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*VERBOSE")) in >> b;
  return b;
}

void ParameterFileManagerT::Facets (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*FACET"))
    ReadIDValues (in, names);
}

void ParameterFileManagerT::Zones (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*ZONE"))
    ReadIDValues (in, names);

}

void ParameterFileManagerT::Boundaries (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*BOUNDARY"))
    ReadIDValues (in, names);
}

CSEConstants::ZoneEdgeT ParameterFileManagerT::ZoneMethod (void)
{
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*EDGETYPE"))
    {
      int i;
      in >> i;
      return int2ZoneEdgeT (i);
    }
  return CSEConstants::kSingleZE;
}

void ParameterFileManagerT::ZoneEdgeNodeSets (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*EDGENODESETS"))
    ReadIDValues (in, names);
}

void ParameterFileManagerT::Contact (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*CONTACT"))
    ReadIDValues (in, names);
}

void ParameterFileManagerT::SingleNodes (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*SINGLE"))
    ReadIDValues (in, names);
}

void ParameterFileManagerT::BlockToNode (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*BLOCKTONODE"))
    ReadIDValues (in, names);
}

void ParameterFileManagerT::NodeSetsMapped (sArrayT& names, ArrayT<CSEConstants::NodeMapMethodT>& meths)
{
  names.Free();
  iArrayT temp;
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*MAPNODE"))
    ReadID_Parameter (in, names, temp);

  meths.Dimension (temp.Length());
  for (int i=0; i < temp.Length(); i++)
    meths[i] = int2NodeMapMethodT (temp[i]);
}

void ParameterFileManagerT::SideSetsMapped (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*COPYSIDE"))
    ReadIDValues (in, names);
}

CSEConstants::RenumberMethodT ParameterFileManagerT::RenumberMethod (void)
{
  int f = CSEConstants::kNoRenumber;
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*BLOCKTONODE"))
    in >> f;
  return int2RenumberMethodT (f);
}

void ParameterFileManagerT::SplitBlocks (sArrayT& names, ArrayT<CSEConstants::SplitMethodT>& meths)
{
  names.Free();
  iArrayT temp;
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*SPLITELEMENT"))
    ReadID_Parameter (in, names, temp);

  meths.Dimension (temp.Length());
  for (int i=0; i < temp.Length(); i++)
    meths[i] = int2SplitMethodT (temp[i]);
}


/*********** PRIVATE ***************/

bool ParameterFileManagerT::AdvanceTo (ifstreamT& in, const StringT& key) const
{
  StringT test;
  while (in.good())
    {
      in >> test;
      test.ToUpper();
      if (strncmp (test.Pointer(), key.Pointer(), key.StringLength()) == 0)
	return true;
    }
  return false;
}

void ParameterFileManagerT::ReadIDValues (ifstreamT& in, sArrayT& names) const
{
  AutoArrayT<StringT> t;
  StringT temp;
  while (in.good())
    {
      in >> temp;
      if (temp[0] == '*') break;
      t.Append (temp);
    }
  names = t;
}

void ParameterFileManagerT::ReadID_Parameter (ifstreamT& in, sArrayT& name, iArrayT& params) const
{
  AutoArrayT<StringT> t;
  StringT temp;
  iAutoArrayT i;
  int itemp;
  while (in.good())
    {
      in >> temp;
      if (temp[0] == '*') break;
      in >> itemp;
      t.Append (temp);
      i.Append (itemp);
    }

  name = t;
  params.Dimension (i.Length());
  params.CopyPart (0, i, 0, i.Length());
}
