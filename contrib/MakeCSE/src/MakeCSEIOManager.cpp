// file: MakeCSEIOManager.cpp

// created: SAW 12/21/99

// this class redirects the Reader variable to MakeCSEReader

#include "MakeCSEIOManager.h"
#include "MakeCSE.h"
#include "ifstreamT.h"
#include "Quad2Tri.h"
#include "FEManager.h"

#include "ExodusInputT.h"
#include "TahoeInputT.h"
#include "EnSightInputT.h"
#include "AbaqusInputT.h"
#include "ExodusOutputT.h"
#include "FE_ASCIIT.h"
#include "EnSightOutputT.h"
#include "AbaqusOutputT.h"
#include "TecPlotOutputT.h"

using namespace Tahoe;

MakeCSEIOManager::MakeCSEIOManager (ostream& out) :
  ModelManagerT (out),
  fVerbose (0),
  fRenumber (0),
  fZoneEdge (MakeCSE::kSingleZE),
  fData (kNumDataTypes)
{
}

void MakeCSEIOManager::Interactive (void)
{
  StringT filename (81);
  StringT answer (81);
  if (fEchoInput) fEchoInput.close();
  cout << "\n Do you want to write an input file? (1 or y)? ";
  cin.getline (answer.Pointer(), 80, '\n');
  if (answer[0] == 'y' || answer[0] == 'Y' || answer[0] == '1') 
    {
      cout << " Enter file name for input file: ";
      cin.getline (filename.Pointer(), 80, '\n');
      fEchoInput.open (filename);
      fEcho = true;
    }
  else
    fEcho = false;

  InteractiveIO ();

  fEchoInput.flush();

  if (fEcho) fEchoInput << "*VERBOSE " << fVerbose << endl;

  InteractiveCSE ();

  InteractiveSplitElement ();

  fEchoInput.flush ();
}

void MakeCSEIOManager::ReadParameters (ifstreamT& in, bool interactive, const StringT& program, const StringT& version)
{
  fMessage << "\n Welcome to: " << program << " " << version;
  fMessage << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";
  StringT database (81);
  if (interactive)
    {
      fMessage << " No input file, interactive session.\n\n";
      Interactive ();
      cout << "\n Enter root for output files: ";
      cin.getline (database.Pointer(), 80, '\n');
    }

  else 
    {
      fMessage << " Reading from Input file . . . . . . . . . . . . = " 
	  << in.filename() << '\n';
      ReadInputFile (in);

      // change input file name to make slightly different root
      // do not wnat to write over input database
      database = in.filename();
      database.Root();
      database.Append ("_db.in");
    }

  // extension will be stripped off, so need one to strip
  if (!strchr (database.Pointer(), '.')) database.Append (".in");
  SetOutput(program, version, fTitle, database, fOutputFormat);
}

void MakeCSEIOManager::InputData (int& data, int key) const
{
  switch (key)
    {
    case kExecution: data = fVerbose; break;
    case kRenumber: data = fRenumber; break;
    case kZoneEdge: data = fZoneEdge; break;
    default:
      {
	cout << "\nMakeCSEIOManager::InputData cannot find integer value.\n";
	cout << "Key = " << key << endl << endl;
	throw eGeneralFail;
      }
    }
}

void MakeCSEIOManager::InputData (iArrayT& data, int key) const
{
  if (key > -1 && key < kNumDataTypes)
    data = fData[key];
  else
      {
	cout << "\nMakeCSEIOManager::InputData cannon find iArray value.\n";
	cout << "Key = " << key << endl << endl;
	throw eGeneralFail;
      }
}

/* output functions */
void MakeCSEIOManager::WriteGeometry(void)
{
	fOutput->WriteGeometry();
}


//***************** private **************

void MakeCSEIOManager::InteractiveIO (void)
{
  StringT filename, answer (81);
  IOBaseT temp (cout);
  temp.InputFormats (cout);

  bool opened = false;
  while (!opened)
    {
      // read input format
      cout << "\n Enter database type to start with: ";
      cin >> fFormat;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line

      // read database file name
      cout << "\n Enter database file name: ";
      cin >> fInputName;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line

      try { SetInput(); }
      catch (int errorcode) 
	{ if (errorcode != eBadInputValue) throw errorcode; }
      opened = true;
    }
  if (fEcho) fEchoInput << "*INPUT " << fFormat << "\n"
			<< fInputName << "\n";

  // read output format
  cout << "\n Enter output format: ";
  cin >> fOutputFormat;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (fEcho) fEchoInput << "*OUTPUT " << fOutputFormat << endl;

  // read if Tahoe II files are external/internal
  if (fOutputFormat == IOBaseT::kTahoeII)
    {
      cout << "\n Enter 1 for external files, or any other key for inline: ";
      cin.getline (answer.Pointer(), 80, '\n');
      if (answer[0] == '1') fExternTahoeII = true;
      if (fEcho) fEchoInput << " " << fExternTahoeII << endl;
    }

  // set output from read parameters
}

// scan parameter file
void MakeCSEIOManager::ReadInputFile(ifstreamT& in)
{
  StringT word1;
  try
    {
      while (in.good())
	{
	  if (!ReadWord1 (in, word1)) return;
	  if (strncmp (word1, "EOF", 3) == 0) return;
	  Parse (in, word1);
	}
    }
 catch (int j)
   {
     cout << '\n';
     if (j > 0) cout << "   Invalid Parameter = " << j << endl;
     cout << "   Last keyword = " << word1 << endl;
     throw eBadInputValue;
   } 

 catch (StringT& j)
   {
     cout << "\n   Invalid Parameter = " << j << endl;
     cout << "   Last keyword = " << word1 << endl;
     throw eBadInputValue;
   } 
}

void MakeCSEIOManager::Parse (ifstreamT& in, StringT& word1)
{
  if (strncmp (word1, "OUTPUT", 6) == 0)
    ReadOutputFormat (in);

  else if (strncmp (word1, "INPUT", 5) == 0)
    ReadInputFormat (in);

  else if (strncmp (word1, "TITLE", 5) == 0)
    {
      fTitle.GetLineFromStream (in);
      fMessage << "\n Title: " << fTitle << "\n\n";
    }

  else if (strncmp (word1, "VERBOSE", 7) == 0)
    in >> fVerbose;

  else if (strncmp (word1, "FACET", 5) == 0)
    ReadIDColumnal (in, kFacet, 2);

  else if (strncmp (word1, "ZONE", 4) == 0)
    ReadIDColumnal (in, kZone);

  else if (strncmp (word1, "EDGETYPE", 8) == 0)
    in >> fZoneEdge;

  else if (strncmp (word1, "EDGENODESET", 8) == 0)
    ReadMultiID (in, kZoneEdgeNSets);

  else if (strncmp (word1, "BOUNDARY", 7) == 0)
    Read2IDColumnal (in, kBoundary);

  else if (strncmp (word1, "SINGLE", 6) == 0)
    ReadMultiID (in, kSingleNodes);

  else if (strncmp (word1, "MAPNODE", 7) == 0)
    ReadIDColumnal (in, kMapNodes);

  else if (strncmp (word1, "COPYSIDE", 8) == 0)
    ReadIDColumnal (in, kCopySide);

  else if (strncmp (word1, "SPLIT", 5) == 0)
    ReadIDColumnal (in, kElementSplit);

  else if (strncmp (word1, "CONTACT", 6) == 0)
    ReadMultiID (in, kContactData);

  else if (strncmp (word1, "RENUMBER", 7) == 0)
    in >> fRenumber;

  else if (strncmp (word1, "BLOCKTONODE", 10) == 0)
    ReadMultiID (in, kBlockToNode);

  else
    MakeCSEIOManager::Parse (in, word1);
}

bool MakeCSEIOManager::ReadWord1 (ifstreamT& in, StringT& word1) const
{
  // read string instead of character, to assist in user debugging
  StringT word;
  in >> word;
  word.ToUpper();
  if (word.Length() <= 1) return false;
  if (word[0] != '*') 
    {
      cout << "\n\nError Reading Input File: Encountered " << word 
	   << " when expecting *KEYWORD.\n";
      throw -1;
    }

  word1 = word.Pointer(1);
  return true;
}

void MakeCSEIOManager::ReadIDColumnal (ifstreamT& in, int key, int numoptions)
{
  int start, stop;
  AutoArrayT<int> data;
  while (ReadID (in, start, stop))
    {
      iArrayT options (numoptions);
      for (int j=0; j < numoptions; j++)
	in >> options[j];

      for (int i=start; i < stop + 1; i++)
	{
	  data.Append (i);
	  data.Append (options);
	}
    }
  fData[key].Allocate (data.Length());
  fData[key].CopyPart (0, data, 0, data.Length());
}

void MakeCSEIOManager::Read2IDColumnal (ifstreamT& in, int key, int numoptions)
{
  int start1, stop1;
  int start2, stop2;
  iAutoArrayT data;
  while (ReadID (in, start1, stop1)) // read first ID values
    {
      ReadID (in, start2, stop2); // read second ID values

      iArrayT options (numoptions); // read options
      for (int j=0; j < numoptions; j++)
	in >> options[j];

      // expand start, stop
      for (int i=start1; i < stop1 + 1; i++)
	  {
	    data.Append (i);
	    data.Append (start2);
	    data.Append (stop2);
	    data.Append (options);
	  }
    }
  fData[key].Allocate (data.Length());
  fData[key].CopyPart (0, data, 0, data.Length());
}

bool MakeCSEIOManager::ReadID (ifstreamT& in, int& start, int& stop) const
{
  start = -1;
  char next;
  next = in.next_char();
  if (!isdigit(next)) return false;
  if (!in.good()) return false;

  in >> start;
  stop = start;
  next = in.next_char();
  if (next == '-')
    {
      in >> next;
      in >> stop;
    }
  return true;
}

void MakeCSEIOManager::ReadMultiID (ifstreamT& in, int key)
{
  int start, stop;
  iAutoArrayT ids;
  while (ReadID (in, start, stop))
    for (int i=start; i < stop + 1; i++)
      ids.Append (i);
  fData[key].Allocate (ids.Length());
  fData[key].CopyPart (0, ids, 0, ids.Length());
}

void MakeCSEIOManager::InteractiveCSE (void)
{
  int method;
  StringT answer (81);
  cout << "\n 1. Facet Method\n";
  cout << " 2. Zone Method\n";
  cout << " 3. Boundary Method\n";
  cout << "\n How would you like to specify the facets: ";
  cin >> method;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line

  switch (method)
    {
    case 1:
      {
	int num;
	cout << "\n Enter the number of side sets: ";
	cin >> num;
	cin.getline (answer.Pointer(), 80, '\n'); // clear line
	if (fEcho) fEchoInput << "*FACET\n";
	Read3D ("Side Set ID", "Element Block ID", "Output Element Block ID", kFacet, num);
	break;
      }
    case 2:
      {
	int num;
	cout << "\n Enter the number of element blocks: ";
	cin >> num;
	cin.getline (answer.Pointer(), 80, '\n'); // clear line
	if (fEcho) fEchoInput << "*ZONE\n";
	Read2D ("Element Block ID", "Output Element Block ID", kZone, num);
	break;
      }
    case 3:
      {
	int num;
	cout << "\n Enter the number of element blocks: ";
	cin >> num;
	cin.getline (answer.Pointer(), 80, '\n'); // clear line
	if (fEcho) fEchoInput << "*BOUNDARY\n";
	Read3D ("Element Block ID 1", "Element Block ID 2", "Output Element Block ID", kBoundary, num);
	break;
      }
    default:
      {
	cout << "\n Invalid Method \n";
	throw eGeneralFail;
      }
    }

  // see if user does not want edge single noded
  int number;
  if (method == 2)
    {
      cout << "\n" << MakeCSE::kSingleZE << ". All single noded.\n";
      cout << MakeCSE::kDoubleZE << ". All double noded.\n";
      cout << MakeCSE::kMixSingZE << ". Mixed, I'll specify the node sets for single noding.\n";
      cout << MakeCSE::kMixDoubZE << ". Mixed, I'll specify the node sets for double noding.\n";
      cout << "\n How do you want your edge zone nodes done: ";
      cin >> fZoneEdge;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line
      if (fEcho) fEchoInput << "*EDGETYPE " << fZoneEdge << "\n";

      if (fZoneEdge ==  MakeCSE::kMixSingZE || fZoneEdge ==  MakeCSE::kMixDoubZE)
	{
	  cout << "\n Enter the number of node sets: ";
	  cin >> number;
	  cin.getline (answer.Pointer(), 80, '\n'); // clear line
	  if (number > 0) 
	    {
	      if (fEcho) fEchoInput << "*EDGENODESET\n";
	      Read ("Node Set ID", kZoneEdge, number);
	    }
	}
    }

  // contact data
  cout << "\n Enter the number of Output CSE Blocks to save side/node set data about: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0)  
    {
      if (fEcho) fEchoInput << "*CONTACT\n";
      Read ("CSE Output Block ID", kContactData, number);
    }

  cout << "\n Enter the number of Single Node node sets: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*SINGLE\n";
      Read ("Node Set ID", kSingleNodes, number);
    }

  cout << "\n Enter the number of Element Blocks to save as node sets: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*BLOCKTONODE\n";
      Read ("ElementBlock ID", kBlockToNode, number);
    }

  cout << "\n Enter the number of Node Sets to transfer: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0)
    {
      cout << "\n" << NodeManagerPrimitive::kSurface1 << ". Surface 1 Transfer Method\n";
      cout << NodeManagerPrimitive::kSurface2 << ". Surface 2 Transfer Method\n";
      cout << NodeManagerPrimitive::kMap << ". Mapping (Surface 2) Transfer Method\n";
      cout << NodeManagerPrimitive::kSplit << ". Splitting Transfer Method\n\n";
      if (fEcho) fEchoInput << "*MAPNODE\n";
      Read2D ("Node Set ID", "Transfer Method", kMapNodes, number);
    }

  cout << "\n Enter the number of Side Sets to transfer: ";
  cin >> number;
  cin.getline (answer.Pointer(), 80, '\n'); // clear line
  if (number > 0) 
    {
      if (fEcho) fEchoInput << "*COPYSIDE\n";
      Read2D ("Side Set ID", "Element Block ID", kCopySide, number);
    }

  cout << "\n" << FEManager::kNoRenumber << ". No Renumbering\n";
  cout << FEManager::kRenumberAdded << ". Renumber Added Nodes\n";
  cout << FEManager::kRenumberAll << ". Renumber All Nodes\n";
  cout << "\n Enter renumbering option: ";
  cin >> fRenumber;
  if (fEcho) fEchoInput << "*RENUMBER " << fRenumber << "\n";
  cin.getline (answer.Pointer(), 80, '\n'); // clear line  
}

void MakeCSEIOManager::InteractiveSplitElement (void)
{
  StringT answer (81);
  cout << "\n Did you want to split elements (y or 1, n or 0)? ";
  cin.getline (answer.Pointer(), 80, '\n');
  if (answer[0] == 'y' || answer[0] == 'Y' || answer[0] == '1')
    {
      int num;
      cout << "\n Enter the number of element blocks to convert: ";
      cin >> num;
      cout << "\n" << Quad2Tri::kXMethod << ". X-Method for Quad 4 to Tria 3\n";
      cout << Quad2Tri::kSlashMethod << ". Slash Method for Quad 4 to Tria 3\n";
      cout << Quad2Tri::kBackSlashMethod << ". Back Slash Method for Quad 4 to Tria 3\n";
      cout << Quad2Tri::kStarMethod << ". Star Method for Quad8 to Tria 3\n\n";
      if (fEcho) fEchoInput << "*SPLITELEMENT\n";
      Read2D ("Element Block ID", "Split Method", kElementSplit, num);
    }  
}

void MakeCSEIOManager::Read (const char* first, int key, int num)
{
  StringT answer (81);
  fData[key].Allocate (num);
  int *p = fData[key].Pointer();
  for (int i=0, j=0; i < num; i++)
    {
      cout << " Enter " << first << " " << i+1 << ": ";
      cin >> *p++;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
    }
  if (fEcho) fEchoInput << fData[key] << '\n';
}

void MakeCSEIOManager::Read2D (const char* first, const char* second, int key, int num)
{
  StringT answer (81);
  fData[key].Allocate (num*2);
  int *p = fData[key].Pointer();
  for (int i=0, j=0; i < num; i++)
    {
      cout << " Enter " << first << " " << i+1 << ": ";
      cin >> *p++;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
      cout << " Enter " << second << " " << i+1 << ": ";
      cin >> *p++;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
    }
  if (fEcho) 
      fData[key].WriteWrapped (fEchoInput, 2);
}

void MakeCSEIOManager::Read3D (const char* first, const char* second, const char* third, int key, int num)
{
  StringT answer (81);
  fData[key].Allocate (num*3);
  int *p = fData[key].Pointer();
  for (int i=0, j=0; i < num; i++)
    {
      cout << " Enter " << first << " " << i+1 << ": ";
      cin >> *p++;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
      cout << " Enter " << second << " " << i+1 << ": ";
      cin >> *p++;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
      cout << " Enter " << third << " " << i+1 << ": ";
      cin >> *p++;
      cin.getline (answer.Pointer(), 80, '\n'); // clear line  
    }
  if (fEcho) 
    fData[key].WriteWrapped (fEchoInput, 3);
}

void MakeCSEIOManager::ReadOutputFormat (ifstreamT& in)
{
  in >> fOutputFormat;
  fMessage << " Output format . . . . . . . . . . . . . . . . . = " 
       << fOutputFormat << '\n';
  IOBaseT temp (cout);
  temp.OutputFormats (fMessage);

  // if TahoeII, choose external or inline data
  if (fOutputFormat == IOBaseT::kTahoeII)
    {
      int tahoetype;
      in >> tahoetype;
      if (tahoetype == 1) fExternTahoeII = true;
      fMessage << " External Files. . . . . . . . . . . . . . . . . = " 
	   << fExternTahoeII << '\n';
    }
}

void MakeCSEIOManager::ReadInputFormat (ifstreamT& in)
{
  int InputFormat;
  in >> fFormat;

  fMessage << " Input format. . . . . . . . . . . . . . . . . . = " 
       << fFormat  << '\n';
  IOBaseT temp (cout);
  temp.InputFormats (fMessage);

  // read database name
  if (fFormat > IOBaseT::kTahoe)
    {
      in >> fInputName;
      if (fInputName[0] == '*') throw fInputName;
    }
  SetInput ();
}

void MakeCSEIOManager::SetInput(void)
{
  cout << " Reading data from: " << fInputName << endl;
  fMessage << " Reading data from: " << fInputName << endl;
	switch (fFormat)
	{
		case IOBaseT::kExodusII:
		{
		  fInput = new ExodusInputT (fMessage);
		  fInput->Open (fInputName);
			break;
		}
	        case IOBaseT::kTahoeII:
		{
			fInput = new TahoeInputT (fMessage);
			fInput->Open (fInputName);
			break;		 
		}
	        case IOBaseT::kEnSight:
		{
		        fInput = new EnSightInputT (fMessage, false);
			fInput->Open (fInputName);
			break;
		}
	        case IOBaseT::kEnSightBinary:
		{
		        fInput = new EnSightInputT (fMessage, true);
			fInput->Open (fInputName);
			break;
		}
	        case IOBaseT::kAbaqus:
	        case IOBaseT::kAbaqusBinary:
		{
		        fInput = new AbaqusInputT (fMessage);
			fInput->Open (fInputName);
			break;
		}
	        case IOBaseT::kTahoe:
		{
		  cout << "\nIOManager::SetInput Tahoe I formatted input no longer accepted. Please create a Tahoe II geometry file.\n";
		}  
		default:
		{
			cout << "\nIOManager::SetInput Wrong Input File Type.";
			fMessage << "\nIOManager::SetInput Wrong Input File Type.";
			throw eBadInputValue;
      	}
	}
	
	if (!fInput) throw eOutOfMemory;
}

void MakeCSEIOManager::SetOutput(const StringT& program_name, 
	const StringT& version, const StringT& title, const StringT& input_file, 
	IOBaseT::FileTypeT output_format)
{	
	ArrayT<StringT> outstrings (4);
	outstrings[0] = input_file;
	outstrings[1] = title;
	outstrings[2] = program_name;
	outstrings[3] = version;

	const int kdigits = 4;
	switch (output_format)
	  {
	  case IOBaseT::kExodusII:
	    fOutput = new ExodusOutputT(fMessage, outstrings);
	    break;
	  case IOBaseT::kTahoe:
	  case IOBaseT::kTahoeII:
	    fOutput = new FE_ASCIIT(fMessage, fExternTahoeII, outstrings);
	    break;
	  case IOBaseT::kEnSight:
	    fOutput = new EnSightOutputT (fMessage, outstrings, kdigits, false);
	    break;
	  case IOBaseT::kEnSightBinary:
	    fOutput = new EnSightOutputT (fMessage, outstrings, kdigits, true);
	    break;
	  case IOBaseT::kAbaqus:
	    fOutput = new AbaqusOutputT (fMessage, outstrings, false);
	    break;
	  case IOBaseT::kAbaqusBinary:
	    fOutput = new AbaqusOutputT (fMessage, outstrings, true);
	    break;

	  case IOBaseT::kTecPlot:
	    fOutput = new TecPlotOutputT (fMessage, outstrings, kdigits);
	    break;
	    
	  default:
	    {			
	      cout << "\n IOManager::SetOutput unknown output format:" 
		   << output_format << endl;
	      fMessage  << "\n IOManager::SetOutput unknown output format:" 
		    << output_format << endl;
	      throw eBadInputValue;
	    } 
	  }	
	if (!fOutput) throw eOutOfMemory;
}

