/* $Id: IOManager.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: sawimme (10/12/1999)                                          */
/* this class creates InputBaseT and OutputBaseT pointers                 */

#include "IOManager.h"

#include "fstreamT.h"
#include "OutputSetT.h"
#include "dArrayT.h"

// input
#include "InputBaseT.h"
#include "ExodusInputT.h"
#include "TextInput.h"
#include "EnSightInputT.h"
#include "AbaqusInputT.h"

// output
#include "FE_ASCIIT.h"
#include "ExodusOutputT.h"
#include "EnSightOutputT.h"
#include "AbaqusOutputT.h"
#include "TecPlotOutputT.h"

IOManager::IOManager(ostream& outfile, const StringT& program_name,
	const StringT& version, const StringT& title, const StringT& input_file,
	IOBaseT::FileTypeT output_format):
	fLog(outfile),
	fOutputFormat(output_format),
	fOutput(NULL),
	fInputFormat(IOBaseT::kTahoe),
	fInput(NULL),
	fEcho (false),
	fOutput_tmp(NULL)
{
	/* construct output formatter */
	fOutput = SetOutput(program_name, version, title, input_file, fOutputFormat);
}

IOManager::IOManager(ifstreamT& in, const IOManager& io_man):
	fLog(io_man.fLog),
	fOutputFormat(io_man.fOutputFormat),
	fOutput(NULL),
	fInputFormat(io_man.fInputFormat),
	fInput(NULL),
	fEcho (false),
	fOutput_tmp(NULL)
{
	/* construct output formatter */
	fOutput = SetOutput((io_man.fOutput)->CodeName(),
				(io_man.fOutput)->Version(),
				(io_man.fOutput)->Title(), in.filename(), fOutputFormat);
}

// constructor to use in conjunction with ReadParameters
IOManager::IOManager (ostream& out) :
	fLog (out),
	fOutput(NULL),
	fInput(NULL),
	fOutputFormat (IOBaseT::kExodusII),
	fInputFormat (IOBaseT::kTahoe),
	fEcho (false)
{

}

IOManager::~IOManager(void)
{
	delete fInput;
	fInput = NULL;

	/* in case output is diverted */
	RestoreOutput();
	
	delete fOutput;
	fOutput = NULL;
}

/*********** OUTPUT **************/


void IOManager::NextTimeSequence(int sequence_number)
{
	fOutput->NextTimeSequence(sequence_number);
}

/* set model coordinates */
void IOManager::SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_map)
{
	fOutput->SetCoordinates(coordinates, node_map);
}
	
/* register the output for an element set. returns the output ID */
int IOManager::AddElementSet(const OutputSetT& output_set)
{
	return fOutput->AddElementSet(output_set);
}

const ArrayT<OutputSetT*>& IOManager::ElementSets(void) const
{
	return fOutput->ElementSets();
}

void IOManager::AddNodeSet(const iArrayT& nodeset, int setID)
{
	fOutput->AddNodeSet(nodeset, setID);
}

void IOManager::AddSideSet(const iArray2DT& sideset, int setID, int group_ID)
{
	fOutput->AddSideSet(sideset, setID, group_ID);
}

/* output functions */
void IOManager::WriteGeometry(void)
{
	fOutput->WriteGeometry();
}

void IOManager::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT format) const
{
	if (!fOutput)
	{
		cout << "\n IOManager::WriteGeometryFile: output must be configured" << endl;
		throw eGeneralFail;		
	}

	fOutput->WriteGeometryFile(file_name, format);
}

void IOManager::WriteOutput(double time, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	fOutput->WriteOutput(time, ID, n_values, e_values);
}

/* (temporarily) re-route output */
void IOManager::DivertOutput(const StringT& outfile)
{
	/* can only divert once */
	if (fOutput_tmp != NULL)
		cout << "\n IOManager::DivertOutput: cannot divert output to \""
		     << outfile << "\".\n"
		     <<   "     Output is already diverted to \"" <<
		     fOutput->OutputRoot() << "\"" << endl;
	else
	{
		/* store main out */
		fOutput_tmp = fOutput;
	
		/* construct temporary output formatter */
		StringT tmp(outfile);
		tmp.Append(".ext"); //OutputBaseT takes root of name passed in
		fOutput = SetOutput(fOutput_tmp->CodeName(), fOutput_tmp->Version(),
			fOutput_tmp->Title(), tmp, fOutputFormat);
		
		/* add all output sets */
		const ArrayT<OutputSetT*>& element_sets = fOutput_tmp->ElementSets();
		for (int i = 0; i < element_sets.Length(); i++)
			AddElementSet(*(element_sets[i]));
			
		/* set coordinate data */
		fOutput->SetCoordinates(fOutput_tmp->Coordinates(), fOutput_tmp->NodeMap());		
	}
}

void IOManager::RestoreOutput(void)
{
	if (fOutput_tmp != NULL)
	{
		/* delete temp output formatter */
		delete fOutput;
	
		/* restore main out */
		fOutput = fOutput_tmp;
		fOutput_tmp = NULL;
	}
}

const OutputSetT& IOManager::OutputSet(int ID) const
{
	return fOutput->OutputSet(ID);
}

/*********** INPUT **************/

int IOManager::NumElementGroups(void) const
{
	return fInput->NumElementGroups();
}

int IOManager::NumSideSets(void) const
{
	return fInput->NumSideSets();
}

int IOManager::NumNodeSets(void) const
{
	return fInput->NumNodeSets();
}

void IOManager::GroupNumbers(iArrayT& groupnums) const
{
	fInput->GroupNumbers(groupnums);
}

void IOManager::SideSetNumbers(iArrayT& sidenums) const
{
	fInput->SideSetNumbers(sidenums);
}

void IOManager::NodeSetNumbers(iArrayT& nodenums) const
{
	fInput->NodeSetNumbers(nodenums);
}

void IOManager::ReadConnectivity(const int group, GeometryT::CodeT & geocode, iArray2DT& connects, iArrayT& elementmap)
{
	fInput->ReadConnectivity(group, geocode, connects, elementmap);
}

void IOManager::ReadCoordinates(dArray2DT& coords, iArrayT& nodemap)
{
	fInput->ReadCoordinates(coords, nodemap);
}

void IOManager::ReadSideSet(const int set_num, iArray2DT& sides, const int global)
{
	if (global)
		fInput->ReadSideSetGlobal(set_num, sides);
	else
		fInput->ReadSideSet(set_num, sides);
}

void IOManager::ReadNodeSet(const int set_num, iArrayT& nodes)
{
	fInput->ReadNodeSet(set_num, nodes);
}

void IOManager::ReadTimeSteps (dArrayT& steps)
{
fInput->ReadTimeSteps (steps);
}

void IOManager::ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id)
{
fInput->ReadLabels (nlabels, elabels, group_id);
}

void IOManager::ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues)
{
fInput->ReadVariables (step, group_id, nvalues, evalues);
}

void IOManager::ReadParameters (ifstreamT& in, bool interactive, const StringT& program, const StringT& version)
{
fLog << "\n Welcome to: " << program << " " << version;
fLog << "\n\n Build Date: " << __DATE__ " " << __TIME__ << "\n\n";
StringT database (81);
if (interactive)
{
fLog << " No input file, interactive session.\n\n";
Interactive ();
cout << "\n Enter root for output files: ";
cin.getline (database.Pointer(), 80, '\n');
}

else
{
fLog << " Reading from Input file . . . . . . . . . . . . = "
	  << in.filename() << '\n';
ReadInputFile (in);

// change input file name to make slightly different root
// do not wnat to write over input database
database.DefaultName (in.filename(), ".in", "_db.in", -1);
}

// extension will be stripped off, so need one to strip
if (!strchr (database.Pointer(), '.')) database.Append (".in");
fOutput = SetOutput(program, version, fTitle, database, fOutputFormat);
}

// gather parameter data from user interactively
void IOManager::Interactive (void)
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
}

// scan parameter file
void IOManager::ReadInputFile(ifstreamT& in)
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

void IOManager::InputData (int& data, int key) const
{
#pragma unused (data)
#pragma unused (key)
}

void IOManager::InputData (iArrayT& data, int key) const
{
#pragma unused (data)
#pragma unused (key)
}

void IOManager::SetInput(void)
{
cout << " Reading data from: " << fInDatabase << endl;
fLog << " Reading data from: " << fInDatabase << endl;
	switch (fInputFormat)
	{
		case IOBaseT::kExodusII:
		{
			fInput = new ExodusInputT (fLog, fInDatabase);
			break;
		}
	        case IOBaseT::kTahoeII:
		{
			fInput = new TextInput (fLog, fInDatabase);
			break;		
		}
	        case IOBaseT::kEnSight:
		{
		        fInput = new EnSightInputT (fLog, false, fInDatabase);
			break;
		}
	        case IOBaseT::kEnSightBinary:
		{
		        fInput = new EnSightInputT (fLog, true, fInDatabase);
			break;
		}
	        case IOBaseT::kAbaqus:
		{
		        fInput = new AbaqusInputT (fLog, false, fInDatabase);
			break;
		}
	        case IOBaseT::kAbaqusBinary:
		{
		        fInput = new AbaqusInputT (fLog, true, fInDatabase);
			break;
		}
	        case IOBaseT::kTahoe:
		{
		  cout << "\nIOManager::SetInput Tahoe I formatted input no longer accepted. Please create a Tahoe II geometry file.\n";
		}
		default:
		{
			cout << "\nIOManager::SetInput Wrong Input File Type.";
			fLog << "\nIOManager::SetInput Wrong Input File Type.";
			throw eBadInputValue;
	}
	}
	
	if (!fInput) throw eOutOfMemory;
}

// call this function to translate between fInput and fOutput
void IOManager::Translate (void)
{
// coordinates
iArrayT nodemap;
dArray2DT coords;
ReadCoordinates (coords, nodemap);
SetCoordinates (coords, &nodemap);
fLog << "\n\nNumber of Nodes: " << coords.MajorDim() << endl;

// node sets
int num_sets = NumNodeSets ();
iArrayT ids (num_sets);
NodeSetNumbers (ids);
ArrayT<iArrayT> nodeset (num_sets);
for (int n=0; n < num_sets; n++)
{
ReadNodeSet (ids[n], nodeset[n]);
AddNodeSet (nodeset[n], ids[n]);
fLog << "NodeSet " << ids[n] << "\n   Length: " << nodeset[n].Length() << endl;
}

// element groups
int num_groups = NumElementGroups ();
iArrayT eids (num_groups), outputID (num_groups);
ArrayT<iArray2DT> connects (num_groups);
GroupNumbers (eids);
for (int e=0; e < num_groups; e++)
{
ArrayT<StringT> n_labels (0), e_labels (0);
iArrayT map;
GeometryT::CodeT geocode;
ReadConnectivity (eids[e], geocode, connects[e], map);
ReadLabels (n_labels, e_labels, eids[e]);
OutputSetT output_set (eids[e], geocode, connects[e], n_labels, e_labels, false);
outputID[e] = AddElementSet (output_set);
fLog << "ElementSet " << eids[e]
	   << "\n   Length: " << connects[e].MajorDim()
	   << "\n   Nodes: " << connects[e].MinorDim()
	   << "\n   Geometry: " << geocode << endl;
for (int z=0; z < n_labels.Length(); z++)
	fLog << "   Node Variable " << z+1 << ": " << n_labels[z] << endl;
for (int y=0; y < e_labels.Length(); y++)
	fLog << "   Element Variable " << y+1 << ": " << e_labels[y] << endl;
}

// side sets
num_sets = NumSideSets ();
SideSetNumbers (ids);
ArrayT<iArray2DT> sideset (num_sets);
for (int s=0; s < num_sets; s++)
{
int group = -1;
ReadSideSet (ids[s], sideset[s], 1);
int *elem = sideset[s].Pointer();
for (int i=0, total = 0; i < num_groups && group < 0; i++)
	{
	  total += connects[i].MajorDim();
	  if (*elem < total) group = i;
	}
ReadSideSet (ids[s], sideset[s], 0);
AddSideSet (sideset[s], ids[s], group);
fLog << "SideSet " << ids[s]
	   << "\n   Length: " << sideset[s].MajorDim()
	   << "\n   Element Group: " << eids[group] << endl;
}

// variable data
dArrayT timesteps;
ReadTimeSteps (timesteps);
if (timesteps.Length() > 0)
for (int t=0; t < timesteps.Length(); t++)
{
	cout << "Time Step " << t+1 << ": " << timesteps[t] << endl;
	for (int g=0; g < num_groups; g++)
	  {
	    fLog << "   Element Group: " << eids[g] << endl;
	    dArray2DT nvalues (0,0), evalues (0,0);
	    ReadVariables (t, eids[g], nvalues, evalues);
	    WriteOutput (timesteps[t], outputID[g], nvalues, evalues);
	  }
}
else
WriteGeometry ();
}

/************************************************************************
* Protected
************************************************************************/

bool IOManager::ReadWord1 (ifstreamT& in, StringT& word1) const
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

void IOManager::Parse (ifstreamT& in, StringT& word1)
{
if (strncmp (word1, "OUTPUT", 6) == 0)
ReadOutputFormat (in);

else if (strncmp (word1, "INPUT", 5) == 0)
ReadInputFormat (in);

else if (strncmp (word1, "TITLE", 5) == 0)
{
fTitle.GetLineFromStream (in);
fLog << "\n Title: " << fTitle << "\n\n";
}

else
{
cout << "\n\nUnknown keyword encountered in input file.\n";
throw -1;
}
}

/************************************************************************
* Private
************************************************************************/

void IOManager::ReadOutputFormat (ifstreamT& in)
{
in >> fOutputFormat;
fLog << " Output format . . . . . . . . . . . . . . . . . = "
<< fOutputFormat << '\n';
PrintFormat (fLog);

// if TahoeII, choose external or inline data
if (fOutputFormat == IOBaseT::kTahoeII)
{
int tahoetype;
in >> tahoetype;
if (tahoetype == 1) fExternTahoeII = true;
fLog << " External Files. . . . . . . . . . . . . . . . . = "
	   << fExternTahoeII << '\n';
}
}

void IOManager::ReadInputFormat (ifstreamT& in)
{
in >> fInputFormat;

fLog << " Input format. . . . . . . . . . . . . . . . . . = "
<< fInputFormat  << '\n';
PrintFormat (fLog);

// read database name
if (fInputFormat > IOBaseT::kTahoe)
{
in >> fInDatabase;
if (fInDatabase[0] == '*') throw fInDatabase;
}
SetInput ();
}

void IOManager::PrintFormat (ostream &log) const
{
log << "    eq. " << IOBaseT::kTahoe         << ", Tahoe\n";
log << "    eq. " << IOBaseT::kTahoeII       << ", Tahoe II\n";
log << "    eq. " << IOBaseT::kTecPlot       << ", TecPlot 7.5\n";
log << "    eq. " << IOBaseT::kEnSight       << ", Ensight 6 Gold ASCII\n";
log << "    eq. " << IOBaseT::kEnSightBinary << ", Ensight 6 Gold Binary\n";
log << "    eq. " << IOBaseT::kExodusII      << ", Exodus II\n";
log << "    eq. " << IOBaseT::kAbaqus        << ", Abaqus ASCII (.fin)\n";
log << "    eq. " << IOBaseT::kAbaqusBinary  << ", Abaqus Binary (.fil)\n";
}

void IOManager::InteractiveIO (void)
{
StringT filename, answer (81);
PrintFormat (cout);

bool opened = false;
while (!opened)
{
// read input format
cout << "\n Enter database type to start with: ";
cin >> fInputFormat;
cin.getline (answer.Pointer(), 80, '\n'); // clear line

// read database file name
cout << "\n Enter database file name: ";
cin >> fInDatabase;
cin.getline (answer.Pointer(), 80, '\n'); // clear line

try { SetInput(); }
catch (int errorcode)
	{ if (errorcode != eBadInputValue) throw errorcode; }
opened = true;
}
if (fEcho) fEchoInput << "*INPUT " << fInputFormat << "\n"
			<< fInDatabase << "\n";

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

/* construct and return new output formatter */
OutputBaseT* IOManager::SetOutput(const StringT& program_name,
	const StringT& version, const StringT& title, const StringT& input_file,
	IOBaseT::FileTypeT output_format)
{
	ArrayT<StringT> outstrings (4);
	outstrings[0] = input_file;
	outstrings[1] = title;
	outstrings[2] = program_name;
	outstrings[3] = version;

	const int kdigits = 4;
	OutputBaseT* output = NULL;
	switch (output_format)
	  {
	  case IOBaseT::kExodusII:
	    output = new ExodusOutputT(fLog, outstrings);
	    break;
	  case IOBaseT::kTahoe:
	  case IOBaseT::kTahoeII:
	    output = new FE_ASCIIT(fLog, fExternTahoeII, outstrings);
	    break;
	  case IOBaseT::kEnSight:
	    output = new EnSightOutputT (fLog, outstrings, kdigits, false);
	    break;
	  case IOBaseT::kEnSightBinary:
	    output = new EnSightOutputT (fLog, outstrings, kdigits, true);
	    break;
	  case IOBaseT::kAbaqus:
	    output = new AbaqusOutputT (fLog, outstrings, false);
	    break;
	  case IOBaseT::kAbaqusBinary:
	    output = new AbaqusOutputT (fLog, outstrings, true);
	    break;
	  case IOBaseT::kTecPlot:
	    output = new TecPlotOutputT (fLog, outstrings, kdigits);
	    break;
	  default:
	    {			
	      cout << "\n IOManager::SetOutput unknown output format:"
		   << output_format << endl;
	      fLog  << "\n IOManager::SetOutput unknown output format:"
		    << output_format << endl;
	      throw eBadInputValue;
	    }
	  }	
	if (!output) throw eOutOfMemory;
	return output;
}

