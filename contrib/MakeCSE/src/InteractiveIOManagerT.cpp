// $Id: InteractiveIOManagerT.cpp,v 1.6 2004-11-19 22:57:24 paklein Exp $
#include "InteractiveIOManagerT.h"
#include "ExceptionT.h"
#include "ModelManagerT.h"

using namespace Tahoe;

InteractiveIOManagerT::InteractiveIOManagerT(const ModelManagerT& model):
	fModel(model)
{

}

void InteractiveIOManagerT::Initialize (void)
{
  /* echoing an input file */
  StringT filename (81);
  StringT answer (81);
  cout << "\n Do you want to write an input file? (1 or y)? ";
  cin.getline (answer.Pointer(), 80, '\n');
  if (answer[0] == 'y' || answer[0] == 'Y' || answer[0] == '1') 
    {
      cout << " Enter file name for input file: ";
      cin.getline (filename.Pointer(), 80, '\n');
      fEchoInput.open (filename);
      fEchoInput << "%\n";
      fEcho = true;
    }
  else {
    fEcho = false;
    cout << endl;
   }
  
  Method ();
}

void InteractiveIOManagerT::InputFormat (IOBaseT::FileTypeT &f, StringT& s)
{
  if (fEcho) 
    fEchoInput << "*INPUT " << f << "\n" << s << "\n";
}

void InteractiveIOManagerT::OutputFormat (IOBaseT::FileTypeT &f, StringT& s)
{
  IOBaseT temp (cout);
  cout << "\n\n";
  temp.OutputFormats (cout);
  cout << "\n Enter the Output Format: ";
  cin >> f;
  cout << "\n Enter the root of the output files: ";
  cin >> s;
  
  if (fEcho) fEchoInput << "*OUTPUT " << f << '\n';
}

bool InteractiveIOManagerT::Verbose (void)
{
  if (fEcho) fEchoInput << "*VERBOSE 0" << endl;
  return false;
}

void InteractiveIOManagerT::Facets (sArrayT& names)
{
  if (fMethod == CSEConstants::kFacet)
    names = fMethodSets;
  else
    names.Free();
}

void InteractiveIOManagerT::Zones (sArrayT& names)
{
  if (fMethod == CSEConstants::kZone)
    names = fMethodSets;
  else
    names.Free();
}

void InteractiveIOManagerT::Boundaries (sArrayT& names)
{
  if (fMethod == CSEConstants::kBoundary)
    names = fMethodSets;
  else
    names.Free();
}

CSEConstants::ZoneEdgeT InteractiveIOManagerT::ZoneMethod (void)
{
  int z;
  StringT answer (81);
  cout << "\n" << CSEConstants::kSingleZE << ". All single noded.\n";
  cout << CSEConstants::kDoubleZE << ". All double noded.\n";
  cout << CSEConstants::kMixSingZE << ". Mixed, I'll specify the node sets for single noding.\n";
  cout << CSEConstants::kMixDoubZE << ". Mixed, I'll specify the node sets for double noding.\n";
  cout << "\n How do you want your edge zone nodes done: ";
  cin >> z;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (fEcho) fEchoInput << "*EDGETYPE " << z << "\n";
  
  return int2ZoneEdgeT (z);
}

void InteractiveIOManagerT::ZoneEdgeNodeSets (sArrayT& names)
{
  int number;
  StringT answer (81);
  cout << "\n Enter the number of node sets: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*EDGENODESET\n";
      names.Dimension (number);
      sArrayT n (1);
      n[0] = "NodeSet";
      ReadIDValues (n, names);
    }
}

void InteractiveIOManagerT::Contact (sArrayT& names)
{
  int number;
  StringT answer (81);
  sArrayT q (1);

  cout << "\n Enter the number of Output CSE Blocks to save side/node set data about: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0)  
    {
      if (fEcho) fEchoInput << "*CONTACT\n";
      q[0] = "CSE Output Block";
      names.Dimension (number);
      ReadIDValues (q, names);
    }
}

void InteractiveIOManagerT::SingleNodes (sArrayT& names)
{
  int number;
  StringT answer (81);
  sArrayT q (1);

  cout << "\n Enter the number of Single Node node sets: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*SINGLE\n";
      q[0] = "Node Set";
      names.Dimension (number);
      ReadIDValues (q, names);
    }
}

void InteractiveIOManagerT::BlockToNode (sArrayT& names)
{
  int number;
  StringT answer (81);
  sArrayT q (1);

  cout << "\n Enter the number of Element Blocks to save as node sets: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*BLOCKTONODE\n";
      q[0] = "Element Block";
      names.Dimension (number);
      ReadIDValues (q, names);
    }
}

void InteractiveIOManagerT::NodeSetsMapped (sArrayT& names, ArrayT<CSEConstants::NodeMapMethodT>& meths)
{
  int number;
  StringT answer (81);
  sArrayT q (1);

	/* list node set ids */
	cout << "\n Node sets:\n";
	const ArrayT<StringT>& ns_ID = fModel.NodeSetIDs();
	for (int i = 0; i < ns_ID.Length(); i++)
		cout << setw(5) << i+1 << ": " << ns_ID[i] << '\n';

  cout << "\n Enter the number of Node Sets to transfer: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0)
    {
      cout << "\n" << CSEConstants::kSurface1 << ". Surface 1 Transfer Method\n";
      cout << CSEConstants::kSurface2 << ". Surface 2 Transfer Method\n";
      cout << CSEConstants::kMap << ". Mapping (Surface 2) Transfer Method\n";
      cout << CSEConstants::kSplit << ". Splitting Transfer Method\n\n";
      if (fEcho) fEchoInput << "*MAPNODE\n";
      q.Dimension (2);
      q[0] = "Node Set";
      q[1] = "Transfer Method";
      names.Dimension (number);
      meths.Dimension (number);
      iArrayT temp (number);
      ReadID_Parameter (q, names, temp);

      for (int i=0; i < number; i++)
	meths[i] = int2NodeMapMethodT (temp[i]);
    }
}

void InteractiveIOManagerT::SideSetsMapped (sArrayT& names)
{
  int number;
  StringT answer (81);
  sArrayT q (1);

	/* list side set ids */
	cout << "\n Side sets:\n";
	const ArrayT<StringT>& ss_ID = fModel.SideSetIDs();
	for (int i = 0; i < ss_ID.Length(); i++) {
		const StringT& block = fModel.SideSetGroupID(ss_ID[i]);
		cout << setw(5) << i+1 << ": " << ss_ID[i] << " in element block " << block << '\n';
	}

  cout << "\n Enter the number of Side Sets to transfer: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*COPYSIDE\n";
      q.Dimension (2);
      q[0] = "Side Set";
      q[1] = "Element Block";
      names.Dimension (number);
      ReadIDValues (q, names);
    }
}

CSEConstants::RenumberMethodT InteractiveIOManagerT::RenumberMethod (void)
{
  int number;
  StringT answer (81);
  cout << "\n" << CSEConstants::kNoRenumber << ". No Renumbering\n";
  cout << CSEConstants::kRenumberAdded << ". Renumber Added Nodes\n";
  cout << CSEConstants::kRenumberAll << ". Renumber All Nodes\n";
  cout << "\n Enter renumbering option: ";
  cin >> number;
  if (fEcho) fEchoInput << "*RENUMBER " << number << "\n";
  cin.getline (answer.Pointer(), 80, '\n'); // clear line  
  return int2RenumberMethodT (number);
}

void InteractiveIOManagerT::SplitBlocks (sArrayT& names, ArrayT<CSEConstants::SplitMethodT>& meths)
{
  StringT answer (81);
  cout << "\n Did you want to split elements (y or 1, n or 0)? ";
  cin.getline (answer.Pointer(), 80, '\n');
  if (answer[0] == 'y' || answer[0] == 'Y' || answer[0] == '1')
    {
      int num;
      cout << "\n Enter the number of element blocks to convert: ";
      cin >> num;
      cout << "\n" << CSEConstants::kXMethod << ". X-Method for Quad 4 to Tria 3\n";
      cout << CSEConstants::kSlashMethod << ". Slash Method for Quad 4 to Tria 3\n";
      cout << CSEConstants::kBackSlashMethod << ". Back Slash Method for Quad 4 to Tria 3\n";
      cout << CSEConstants::kStarMethod << ". Star Method for Quad8 to Tria 3\n\n";
      if (fEcho) fEchoInput << "*SPLITELEMENT\n";
      sArrayT q (2);
      q[0] = "Element Block";
      q[1] = "Split Method";
      names.Dimension (num);
      meths.Dimension (num);
      iArrayT temp (num);
      ReadID_Parameter (q, names, temp);

      for (int i=0; i < num; i++)
	meths[i] = int2SplitMethodT (temp[i]);
    }  
}

//***************** private **************

void InteractiveIOManagerT::Method (void)
{
  int m;
  StringT answer (81);
  cout << "\n " << CSEConstants::kFacet << ". Facet Method\n";
  cout << " " << CSEConstants::kZone << ". Zone Method\n";
  cout << " " << CSEConstants::kBoundary << ". Boundary Method\n";
  cout << "\n How would you like to specify the facets: ";
  cin >> m;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line

  int num;
  sArrayT qnames;
  fMethod = int2CSEMethodT (m);
  switch (m)
    {
    case CSEConstants::kFacet:
      {
	cout << "\n Enter the number of side sets (" << fModel.NumSideSets() << "): ";
	cin >> num;
	cin.getline (answer.Pointer(), 80, '\n'); // clear line
	if (fEcho) fEchoInput << "*FACET\n";
	fMethodSets.Dimension (num*3);
	qnames.Dimension (3);
	qnames[0] = "Side Set";
	qnames[1] = "Element Block";
	qnames[2] = "Output Element Block";
	break;
      }
    case CSEConstants::kZone:
      {
	cout << "\n Enter the number of element blocks (" << fModel.NumElementGroups() << "): ";
	cin >> num;
	cin.getline (answer.Pointer(), 80, '\n'); // clear line
	if (fEcho) fEchoInput << "*ZONE\n";
	fMethodSets.Dimension (num*2);
	qnames.Dimension (2);
	qnames[0] = "Element Block";
	qnames[1] = "Output Element Block";
	break;
      }
    case CSEConstants::kBoundary:
      {
	cout << "\n Enter the number of cohesive element blocks to create: ";
	cin >> num;
	cin.getline (answer.Pointer(), 80, '\n'); // clear line
	if (fEcho) fEchoInput << "*BOUNDARY\n";
	fMethodSets.Dimension (num*3);
	qnames.Dimension (3);
	qnames[0] = "Element Block 1";
	qnames[1] = "Element Block 2";
	qnames[2] = "Output Element Block";
	break;
      }
    }
  ReadIDValues (qnames, fMethodSets);
}

void  InteractiveIOManagerT::ReadIDValues (const sArrayT& q, sArrayT& names)
{
  int num = names.Length() / q.Length();
  StringT answer (81);
  for (int i=0, k=0; i < num; i++)
    {
      for (int j=0; j < q.Length(); j++, k++)
	{
	  cout << "(" << i+1 <<") Enter " << q[j] << ": ";
	  cin >> names[k];
	  cin.getline (answer.Pointer(), 80, '\n'); // clear line  
	  
	  if (fEcho)
	    fEchoInput << names[k] << "   ";
	}

      if (fEcho)
	fEchoInput << '\n';
    }
}

void InteractiveIOManagerT::ReadID_Parameter (const sArrayT& q, sArrayT& names, iArrayT& vals)
{
  if (vals.Length() != names.Length()) ExceptionT::SizeMismatch("InteractiveIOManagerT::ReadID_Parameter");

  int num = vals.Length();
  StringT answer (81);
  for (int i=0; i < num; i++)
    {
      cout << " Enter " << q[0] << " " << i+1 << ": ";
      cin >> names[i];
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
      
      cout << " Enter " << q[1] << " " << i+1 << ": ";
      cin >> vals[i];
      cin.getline (answer.Pointer(), 80, '\n'); // clear line
      
      if (fEcho) fEchoInput << names[i] << " " << vals[i] << '\n';
    }
}
