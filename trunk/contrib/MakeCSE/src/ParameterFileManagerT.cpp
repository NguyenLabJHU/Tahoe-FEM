// $Id: ParameterFileManagerT.cpp,v 1.7 2004-11-19 22:57:24 paklein Exp $
#include "ParameterFileManagerT.h"
#include "ExceptionT.h"
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
      throw ExceptionT::kGeneralFail;
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
      throw ExceptionT::kGeneralFail;
    }
  in >> f >> s;
}

void ParameterFileManagerT::OutputFormat (IOBaseT::FileTypeT &f, StringT& s)
{
  s = "defaultoutput";
  ifstreamT in ('#', fInFile);
  if (!AdvanceTo (in, "*OUTPUT"))
    {
      cout << "\nParameterFileManagerT::OutputFormat: No *OUTPUT in file.\n";
      cout << " Using Tahoe II with a file name of defaultoutput. \n\n";
      f = IOBaseT::kTahoeII;
    }
  else
    in >> f;
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

  // check for lists of ID values within the ID name
  int numcolumns = 2;
  int check = 0;
  CheckIDList (names, numcolumns, check);
}

void ParameterFileManagerT::Boundaries (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*BOUNDARY"))
    ReadIDValues (in, names);

  int numcolumns = 3;
  int check1 = 0;
  int check2 = 1;
  CheckIDList (names, numcolumns, check1);
  CheckIDList (names, numcolumns, check2);
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
  
  int numcols = 1;
  int check = 0;
  CheckIDList (names, numcols, check);
}

void ParameterFileManagerT::Contact (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*CONTACT"))
    ReadIDValues (in, names);
  
  int numcols = 1;
  int check = 0;
  CheckIDList (names, numcols, check);
}

void ParameterFileManagerT::SingleNodes (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*SINGLE"))
    ReadIDValues (in, names);
  
  int numcols = 1;
  int check = 0;
  CheckIDList (names, numcols, check);
}

void ParameterFileManagerT::BlockToNode (sArrayT& names)
{
  names.Free();
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*BLOCKTONODE"))
    ReadIDValues (in, names);
  
  int numcols = 1;
  int check = 0;
  CheckIDList (names, numcols, check);
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

void ParameterFileManagerT::SideSetsMapped(sArrayT& names)
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
  if (AdvanceTo (in, "*RENUMBER"))
    in >> f;
  return int2RenumberMethodT (f);
}

void ParameterFileManagerT::SplitBlocks (sArrayT& names, ArrayT<CSEConstants::SplitMethodT>& meths)
{
  names.Free();
  iArrayT temp;
  ifstreamT in ('#', fInFile);
  if (AdvanceTo (in, "*SPLITELEM"))
    ReadID_Parameter (in, names, temp);

  CheckIDList (names, temp);

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
		if (test == key)
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

void ParameterFileManagerT::CheckIDList (sArrayT& names, int numcols, int check) const
{
  if (check < 0 || check > numcols) throw ExceptionT::kOutOfRange;

  //* are there any?
  int numrows = names.Length()/numcols;
  bool found = false;
  int stop;
  for (int i=0, j=check; i < numrows && !found; i++, j += numcols)
    if (names[j].Tail ('-', stop))
      found = true;
  if (!found) return;

  AutoArrayT<StringT> t;
  for (int k=0; k < numrows; k++)
    {
      int offset = k*numcols;
      if (names[offset + check].Tail ('-', stop))
	{
	  int start = atoi (names[offset + check]);
	  for (int m=start; m < stop+1; m++)
	    {
	      // columns before list
	      for (int begin=0; begin < check; begin++)
		t.Append (names[offset + begin]);

	      // list column
	      StringT n;
	      n.Append (m);
	      t.Append (n);

	      // colums after
	      for (int rest=check+1; rest < numcols; rest++)
		t.Append (names[offset + rest]);
	    }
	}
      else
	{
	  for (int col=0; col < numcols; col++)
	    t.Append (names[offset + col]);
	}
    }
  names.Free ();
  names.Dimension (t.Length());
  names.CopyPart (0, t, 0, t.Length());
}

void ParameterFileManagerT::CheckIDList (sArrayT& names, iArrayT& itemp) const
{
  //* are there any?
  int numrows = names.Length();
  bool found = false;
  int stop;
  for (int i=0; i < numrows && !found; i++)
    if (names[i].Tail ('-', stop))
      found = true;
  if (!found) return;

  AutoArrayT<StringT> t;
  iAutoArrayT it2;
  for (int k=0; k < numrows; k++)
    {
      if (names[k].Tail ('-', stop))
	{
	  int start = atoi (names[k]);
	  for (int m=start; m < stop+1; m++)
	    {
	      // append name
	      StringT n;
	      n.Append (m);
	      t.Append (n);

	      // append parameter
	      it2.Append (itemp[k]);
	    }
	}
      else
	{
	  t.Append (names[k]);
	  it2.Append (itemp[k]);
	}
    }
  names.Free ();
  names.Dimension (t.Length());
  names.CopyPart (0, t, 0, t.Length());

  itemp.Free ();
  itemp.Dimension (it2.Length());
  itemp.CopyPart (0, it2, 0, it2.Length());
}
