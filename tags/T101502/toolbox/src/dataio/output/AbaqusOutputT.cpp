/* $Id: AbaqusOutputT.cpp,v 1.7 2002-07-26 18:43:57 sawimme Exp $ */
/* created: sawimme (05/31/2000)                                          */

#include "AbaqusOutputT.h"

#include <fstream.h>

#include "AbaqusResultsT.h"
#include "OutputSetT.h"
#include "ios_fwd_decl.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "ofstreamT.h"


using namespace Tahoe;

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
  AbaqusResultsT aba (fout);
  CreateResultsFile (0, aba);
  
  // write node sets
  for (int i=0; i < fNodeSets.Length(); i++)
    {
       iArrayT& set = *((iArrayT*) fNodeSets[i]);
      set++;
      aba.WriteNodeSet (fNodeSetNames[i], set);
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
  AbaqusResultsT aba (fout);
  if (fElementSets[ID]->PrintStep() == 0)
    CreateResultsFile (ID, aba);
  else
    {
      StringT filename;
      FileName (ID, filename);
      aba.OpenWrite (filename, fBinary, fBufferWritten);
    }
  
  // start time increment
  int inc = fElementSets[ID]->PrintStep() + 1;
  int step = fSequence + 1;
  double timeincrement = time - fOldTime;
  AbaqusResultsT::AnalysisTypeT analysistype = AbaqusResultsT::kDynamic;
  aba.WriteStartIncrement (step, inc, time, time, timeincrement, analysistype);
  fOldTime = time;
  
  // gather variable data
  const ArrayT<StringT>& n_labels = fElementSets[ID]->NodeOutputLabels();	
  const ArrayT<StringT>& e_labels = fElementSets[ID]->ElementOutputLabels();
  
  // if stress or strain variable, define num direct and num shear components 
  int numdir = 0;
  int numshear = 0;
  if (fCoordinates->MinorDim() == 2)
    {
      numdir = 2;
      numshear = 1;
    }
  else
    {
      numdir = 3;
      numshear = 3;
    }

  // break up variables by block
  const ArrayT<StringT>& blockids = fElementSets[ID]->BlockID();
  int elstart = 1;
  for (int ec=0; ec < fElementSets[ID]->NumBlocks(); ec++)
    {
      const iArray2DT* connects = fElementSets[ID]->Connectivities(blockids[ec]);
      int numelemnodes = connects->MinorDim();
      StringT setname = "Grp";
      setname.Append (blockids[ec]);
      GeometryT::CodeT code = fElementSets[ID]->Geometry();
      
      // nodal variables
      if (n_labels.Length() > 0)
	{
	  // set labels
	  iArrayT nkeys (n_labels.Length()); 
	  SetRecordKey (aba, n_labels, nkeys);

	  // collect node numbers and break up group into block
	  iArrayT nodes_used;
	  dArray2DT blockvals;
	  NodalBlockValues (ID, ec, n_values, blockvals, nodes_used);

	  nodes_used++;
	  for (int n=0; n < nkeys.Length(); n++)
	    {
	      aba.WriteOutputDefinition (nkeys[n], setname, code, numelemnodes);
	      
	      // if stress or strain variable, define num direct and num shear components 
	      if (nkeys[n] == aba.VariableKey ("S") || nkeys[n] == aba.VariableKey ("E"))
		aba.WriteNodeVariables (n, nkeys, blockvals, nodes_used, numdir, numshear);
	      else
		{
		  int index = aba.VariableKeyIndex (nkeys[n]);
		  int dimension = aba.VariableDimension (index);
		  aba.WriteNodeVariables (n, nkeys, blockvals, nodes_used, dimension, 0);
		}
	    }
	  nodes_used--;
	}
      
      if (e_labels.Length() > 0)
	{
	  // set labels
	  iArrayT ekeys (e_labels.Length());
	  SetRecordKey (aba, e_labels, ekeys);
	  
	  // collect by block
	  dArray2DT blockvals (fElementSets[ID]->NumBlockElements(blockids[ec]), e_values.MinorDim());
	  ElementBlockValues (ID, ec, e_values, blockvals);

	  // collect element IDs used
	  iArrayT els_used (fElementSets[ID]->NumBlockElements(blockids[ec]));
	  els_used.SetValueToPosition();
	  els_used += elstart;

	  for (int e=0; e < ekeys.Length(); e++)
	    {
	      aba.WriteOutputDefinition (ekeys[e], setname, code, numelemnodes);

	      // if stress or strain variable, define num direct and num shear components 
	      if (ekeys[e] == aba.VariableKey ("S") || ekeys[e] == aba.VariableKey ("E"))
		aba.WriteElementVariables (e, ekeys, blockvals, els_used, numdir, numshear);
	      else
		{
		  int index = aba.VariableKeyIndex (ekeys[e]);
		  int dimension = aba.VariableDimension (index);
		  aba.WriteElementVariables (e, ekeys, blockvals, els_used, dimension, 0);
		}
	    }
	}
      elstart += fElementSets[ID]->NumBlockElements(blockids[ec]);
    }
  
  // write end increment
  aba.WriteEndIncrement ();
  
  // store amount of buffer written
  fBufferWritten = aba.Close ();
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

void AbaqusOutputT::CreateResultsFile (int ID, AbaqusResultsT& aba)
{
  // file name
  StringT filename;
  FileName (ID, filename);
  fBufferWritten = 0;

  // create new file
  int numelems = fElementSets[ID]->NumElements();
  int numnodes = fElementSets[ID]->NumNodes();
  double elemsize = 1; // default for now
  aba.Create (filename, fBinary, numelems, numnodes, elemsize);

  const ArrayT<StringT>& blockids = fElementSets[ID]->BlockID();
  
  // write element records
  int startnum = 1;
  for (int ic=0; ic < fElementSets[ID]->NumBlocks(); ic++)
    {
      const iArray2DT* connects = fElementSets[ID]->Connectivities(blockids[ic]);
      iArray2DT unoffset = *connects;
      unoffset++;
      aba.WriteConnectivity (fElementSets[ID]->Geometry(), startnum, unoffset);
      startnum += connects->MajorDim();
      unoffset--;
    }

  // write coordinate records
  iArrayT nodes_used;
  nodes_used.Alias(fElementSets[ID]->NodesUsed());
  nodes_used++;
  aba.WriteCoordinates (nodes_used, *fCoordinates);
  nodes_used--;
  
  // write element sets
  startnum = 1;
  for (int ec=0; ec < fElementSets[ID]->NumBlocks(); ec++)
    {
      iArrayT elem_map (fElementSets[ID]->NumBlockElements(blockids[ec]));
      elem_map.SetValueToPosition ();
      elem_map += startnum;
      aba.WriteElementSet (fElementSets[ID]->ID(), elem_map);
      startnum += elem_map.Length();
    }

  // write node_used as node sets
  for (int nc=0; nc < fElementSets[ID]->NumBlocks(); nc++)
    {
      iArrayT nodemap;
      nodemap.Alias(fElementSets[ID]->BlockNodesUsed(blockids[nc]));
      nodemap++;
      aba.WriteNodeSet (fElementSets[ID]->ID(), nodemap);
      nodemap--;
    }
  
  // write active dof
  iArrayT activedof (6);
  activedof = 0;
  switch (fCoordinates->MinorDim())
    {
    case 1:
      activedof[0] = 1;
      break;
    case 2:
      activedof[0] = 1;
      activedof[1] = 2;
      activedof[5] = 3;
      break;
    case 3:
      activedof.SetValueToPosition ();
      activedof++;
      break;
    }
  aba.WriteActiveDOF (activedof);
  
  // write heading
  aba.WriteHeading (fTitle);
  
  // write end increment
  aba.WriteEndIncrement ();
}

void AbaqusOutputT::SetRecordKey (AbaqusResultsT& aba, const ArrayT<StringT>& labels, iArrayT& keys) const
{
  keys = -1;
  for (int i=0; i < labels.Length(); i++)
    {
      const char* l = labels[i];
      if (strncmp (l, "D_X", 3) == 0 ||
	  strncmp (l, "D_Y", 3) == 0 ||
	  strncmp (l, "D_Z", 3) == 0 )
	keys[i] = aba.VariableKey ("U");
      else if (strncmp (l, "U_1", 3) == 0 ||
	       strncmp (l, "U_2", 3) == 0 ||
	       strncmp (l, "U_3", 3) == 0 ||
	       strncmp (l, "U_4", 3) == 0 ||
	       strncmp (l, "U_5", 3) == 0 ||
	       strncmp (l, "U_6", 3) == 0 )
	keys[i] = aba.VariableKey ("U");
      else if (strncmp (l, "s11", 3) == 0 ||
	       strncmp (l, "s22", 3) == 0 ||
	       strncmp (l, "s33", 3) == 0 ||
	       strncmp (l, "s12", 3) == 0 ||
	       strncmp (l, "s13", 3) == 0 ||
	       strncmp (l, "s23", 3) == 0 )
	keys[i] = aba.VariableKey ("S");
      else if (strncmp (l, "e11", 3) == 0 ||
	       strncmp (l, "e22", 3) == 0 ||
	       strncmp (l, "e33", 3) == 0 ||
	       strncmp (l, "e12", 3) == 0 ||
	       strncmp (l, "e13", 3) == 0 ||
	       strncmp (l, "e23", 3) == 0 )
	keys[i] = aba.VariableKey ("E");
      else if (strncmp (l, "V_X", 3) == 0 ||
	       strncmp (l, "V_Y", 3) == 0 ||
	       strncmp (l, "V_Z", 3) == 0 )
	keys[i] = aba.VariableKey ("V");
      else if (strncmp (l, "A_X", 3) == 0 ||
	       strncmp (l, "A_Y", 3) == 0 ||
	       strncmp (l, "A_Z", 3) == 0 )
	keys[i] = aba.VariableKey ("A");
      else if (strncmp (l, "x_X", 3) == 0 ||
	       strncmp (l, "x_Y", 3) == 0 ||
	       strncmp (l, "x_Z", 3) == 0 )
	keys[i] = aba.VariableKey ("COORD");
      else if (strncmp (l, "p1", 2) == 0 ||
	       strncmp (l, "p2", 2) == 0 ||
	       strncmp (l, "p3", 2) == 0)
	keys[i] = aba.VariableKey ("SP");
      else
	{
	  keys[i] = aba.VariableKey ("UVARM");
	  if (fElementSets[0]->PrintStep() == 0)
	    fout << "AbaqusOutputT::SetLabelName, unknown label "
		 << labels[i] << ", writing as UVARM\n";
	}
    }
}


