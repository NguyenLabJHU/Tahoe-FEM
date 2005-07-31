/* $Id: AbaqusOutputT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/31/2000)                                          */

#include "AbaqusOutputT.h"
#include "AbaqusT.h"
#include "OutputSetT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "ArrayT.h"

AbaqusOutputT::AbaqusOutputT(ostream& out, const ArrayT<StringT>& out_strings, bool binary):
OutputBaseT(out, out_strings),
fBinary (binary),
fBufferWritten (0),
fOldTime (0)
{

}

/* print geometry from multiple element groups to one file */
void AbaqusOutputT::WriteGeometry(void)
{
AbaqusT aba (fout, fBinary);
ofstream out;
CreateResultsFile (0, aba, out);

// write node sets
for (int i=0; i < fNodeSets.Length(); i++)
{
StringT name ("NodeSet ");
name.Append (fNodeSetIDs[i]);
iArrayT& set = *((iArrayT*) fNodeSets[i]);
set++;
aba.WriteNodeSet (out, name, set);
set--;
}
}

void AbaqusOutputT::WriteOutput(double time, int ID, const dArray2DT& n_values,
const dArray2DT& e_values)
{
	/* inherited */
	OutputBaseT::WriteOutput(time, ID, n_values, e_values);
	if (fElementSets[ID]->NumNodes() == 0) return;

	// open file
	AbaqusT aba (fout, fBinary);
	ofstream out;
	if (fElementSets[ID]->PrintStep() == 0)
	  CreateResultsFile (ID, aba, out);
	else
	  {
	    StringT filename;
	    FileName (ID, filename);
	    out.open (filename, ios::app);
	    aba.ResetBufferSize (fBufferWritten);
	  }

	// start time increment
	int inc = fElementSets[ID]->PrintStep() + 1;
	int step = fSequence + 1;
	double timeincrement = time - fOldTime;
	AbaqusT::Analysis analysistype = AbaqusT::aDynamic;
	aba.WriteStartIncrement (out, step, inc, time, time, timeincrement, analysistype);
	fOldTime = time;

	// gather variable data
	const ArrayT<StringT>& n_labels = fElementSets[ID]->NodeOutputLabels();	
	const ArrayT<StringT>& e_labels = fElementSets[ID]->ElementOutputLabels();
	const iArray2DT& connects = fElementSets[ID]->Connectivities();
	int numelemnodes = connects.MinorDim();

	// write variable data
	if (n_labels.Length() > 0)
	  {
	    ArrayT<AbaqusT::VariableKeyT> nkeys (n_labels.Length());
	    SetRecordKey (n_labels, nkeys);
	    iArrayT nodes_used;
	    nodes_used.Alias(fElementSets[ID]->NodesUsed());
	    nodes_used++;
	    aba.WriteNodalData (out, nkeys, n_values, nodes_used, fElementSets[ID]->Geometry(), numelemnodes);
	    nodes_used--;
	  }

	if (e_labels.Length() > 0)
	  {
	    ArrayT<AbaqusT::VariableKeyT> ekeys (e_labels.Length());
	    SetRecordKey (e_labels, ekeys);
	    //aba.WriteElementData (out, ekeys, e_values);
	  }

	// write end increment
	aba.WriteEndIncrement (out, false);

	// store amount of buffer written
	fBufferWritten = aba.GetBufferSize ();
}

/*************************************************************************
* Protected
*************************************************************************/

/*************************************************************************
* Private
*************************************************************************/

/* generate database file name for the given ID */
void AbaqusOutputT::FileName(int ID, StringT& filename) const
{
	/* root */
	filename = fOutroot;

	/* tack on sequence number */
	if (fSequence > 0) filename.Append(".sq", fSequence + 1);

	/* I/O ID */
	filename.Append(".io", ID);

	/* changing geometry */
	if (fElementSets[ID]->Changing())
		filename.Append(".ps", fElementSets[ID]->PrintStep() + 1);

	/* extension */
	if (fBinary)
	  filename.Append(".fil");
	else
	  filename.Append(".fin");
}

bool AbaqusOutputT::OpenFile (ofstream& out, const StringT& filename, AbaqusT& aba)
{
if (out) out.close();
aba.ResetBufferSize (fBufferWritten);
out.open (filename);
if (out) return true;
else return false;
}

void AbaqusOutputT::CreateResultsFile (int ID, AbaqusT& aba, ofstream& out)
{
// file name
StringT filename;
FileName (ID, filename);

fBufferWritten = 0;
if (!OpenFile (out, filename, aba)) return;

// create new file
int numelems = fElementSets[ID]->NumElements();
int numnodes = fElementSets[ID]->NumNodes();
double elemsize = 1; // default for now
aba.Create (out, numelems, numnodes, elemsize);

// write element records
int startelemnum = 1;
const iArray2DT& connects = fElementSets[ID]->Connectivities();
iArray2DT unoffset = connects;
unoffset++;
aba.WriteConnectivity (out, fElementSets[ID]->Geometry(), startelemnum, unoffset);

// write coordinate records
iArrayT nodes_used;
nodes_used.Alias(fElementSets[ID]->NodesUsed());
nodes_used++;
aba.WriteCoordinates (out, nodes_used, *fCoordinates);
nodes_used--;

// write element set
iArrayT elem_map (numelems);
elem_map.SetValueToPosition ();
elem_map++;
StringT name = "Grp ";
name.Append (fElementSets[ID]->ID());
aba.WriteElementSet (out, name, elem_map);

// write active dof
iArrayT activedof (fCoordinates->MinorDim());
activedof.SetValueToPosition();
activedof++;
aba.WriteActiveDOF (out, activedof);

// write heading
aba.WriteHeading (out, fTitle);

// write end increment
aba.WriteEndIncrement (out, true);
}

void AbaqusOutputT::SetRecordKey (const ArrayT<StringT>& labels, ArrayT<AbaqusT::VariableKeyT>& keys) const
{
keys = AbaqusT::aUnknown;
for (int i=0; i < labels.Length(); i++)
{
const char* l = labels[i];
if (strncmp (l, "D_X", 3) == 0 ||
	  strncmp (l, "D_Y", 3) == 0 ||
	  strncmp (l, "D_Z", 3) == 0 )
	keys[i] = AbaqusT::aDisplacement;
else if (strncmp (l, "s11", 3) == 0 ||
	       strncmp (l, "s22", 3) == 0 ||
	       strncmp (l, "s33", 3) == 0 ||
	       strncmp (l, "s12", 3) == 0 ||
	       strncmp (l, "s13", 3) == 0 ||
	       strncmp (l, "s23", 3) == 0 )
	keys[i] = AbaqusT::aStress;
else if (strncmp (l, "e11", 3) == 0 ||
	       strncmp (l, "e22", 3) == 0 ||
	       strncmp (l, "e33", 3) == 0 ||
	       strncmp (l, "e12", 3) == 0 ||
	       strncmp (l, "e13", 3) == 0 ||
	       strncmp (l, "e23", 3) == 0 )
	keys[i] = AbaqusT::aTotalStrain;
else if (strncmp (l, "V_X", 3) == 0 ||
	       strncmp (l, "V_Y", 3) == 0 ||
	       strncmp (l, "V_Z", 3) == 0 )
	keys[i] = AbaqusT::aVelocity;
else if (strncmp (l, "A_X", 3) == 0 ||
	       strncmp (l, "A_Y", 3) == 0 ||
	       strncmp (l, "A_Z", 3) == 0 )
	keys[i] = AbaqusT::aAcceleration;
else if (strncmp (l, "x_X", 3) == 0 ||
	       strncmp (l, "x_Y", 3) == 0 ||
	       strncmp (l, "x_Z", 3) == 0 )
	keys[i] = AbaqusT::aCoordVariable;
else if (strncmp (l, "p1", 2) == 0 ||
	       strncmp (l, "p2", 2) == 0 ||
	       strncmp (l, "p3", 2) == 0)
	keys[i] = AbaqusT::aPrinStress;
else
	{
	  keys[i] = AbaqusT::aUVARM;
	  if (fElementSets[0]->PrintStep() == 0)
	    fout << "AbaqusOutputT::SetLabelName, unknown label "
	       << labels[i] << ", writing as UVARM\n";
	}
}
}


