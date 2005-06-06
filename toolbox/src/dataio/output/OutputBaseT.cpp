/* $Id: OutputBaseT.cpp,v 1.22 2005-06-06 06:36:42 paklein Exp $ */
/* created: sawimme (05/18/1999) */
#include "OutputBaseT.h"
#include "OutputSetT.h"

/* database types */
#include "ExodusT.h"
#include "ModelFileT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "AutoArrayT.h"
#include "InverseMapT.h"

using namespace Tahoe;

OutputBaseT::OutputBaseT(ostream& out, const ArrayT<StringT>& out_strings):
	IOBaseT(out),
	fCoordinates(NULL),
	fNodeID(NULL),
	fBounds(NULL),
	fTypes(NULL),
	fParts(NULL),
	fSequence(0)
{
	if (out_strings.Length() > 3)
	{
		fOutroot.Root(out_strings[0]);
		fTitle    = out_strings[1];
		fCodeName = out_strings[2];
		fVersion  = out_strings[3];
	}
	else // general default
	{
		fOutroot  = "default";
		fTitle    = "not given";
		fCodeName = "tahoe";
		fVersion  = "unknown";
	}	
}

/* destructor */
OutputBaseT::~OutputBaseT(void)
{
	for (int i = 0; i < fElementSets.Length(); i++)
		delete fElementSets[i];
	fElementSets.Dimension(0);
}

const OutputSetT& OutputBaseT::OutputSet(int ID) const
{
	return *(fElementSets[ID]);
}

/* return the array of nodes used */
const iArrayT& OutputBaseT::NodesUsed(int ID) const
{
	return fElementSets[ID]->NodesUsed();
}

void OutputBaseT::NextTimeSequence(int sequence_number)
{
	fSequence  = sequence_number;
	for (int i = 0; i < fElementSets.Length(); i++)
		fElementSets[i]->ResetPrintStep();
}

/* set nodal coordinates */
void OutputBaseT::SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_id)
{
	fCoordinates = &coordinates;
	fNodeID = node_id;
	
	/* id list check */
	if (fNodeID && fNodeID->Length() != fCoordinates->MajorDim()) {
		cout << "\n OutputBaseT::SetCoordinates: id list length " << fNodeID->Length() << " doesn't\n"
		     <<   "     match the number of nodes " << fCoordinates->MajorDim() << endl;
		throw ExceptionT::kSizeMismatch;
	}
}

/* SA: add definitions of special arrays for ParaDyn Format Output */
void OutputBaseT::SetBounds(const dArray2DT& bounds)
{
  fBounds = &bounds;
}

void OutputBaseT::SetTypes(const iArrayT& types)
{
  fTypes = &types;
}

void OutputBaseT::SetParts(const iArrayT& parts)
{
  fParts = &parts;
}
/* end added definitions */

int OutputBaseT::AddElementSet(const OutputSetT& output_set)
{
	OutputSetT* copy = new OutputSetT(output_set);
	if (!copy) throw ExceptionT::kOutOfMemory;

	/* ID is just position in array */
	StringT ID;
	ID.Append(fElementSets.Length());
	copy->SetID(ID);
	fElementSets.Append(copy);
	return fElementSets.Length() - 1;
}

int OutputBaseT::NumElements(void) const
{
	int count = 0;
	for (int i = 0; i < fElementSets.Length(); i++)
		count += fElementSets[i]->NumElements();
	return count;
}

void OutputBaseT::AddNodeSet(const iArrayT& nodeset, const StringT& setID)
{
	fNodeSets.Append(&nodeset);
	fNodeSetNames.Append (setID);
}

void OutputBaseT::AddSideSet(const iArray2DT& sideset, const StringT& setID, const StringT& group_ID)
{
	fSideSets.Append(&sideset);
	fSideSetNames.Append (setID);
	fSSGroupNames.Append(group_ID);
}


/* output functions */
void OutputBaseT::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT format) const
{
	if (!fCoordinates)
	{
		cout << "\n OutputBaseT::WriteGeometryFile: pointer to coordinates not set" << endl;
		throw ExceptionT::kGeneralFail;
	}

	if (format == IOBaseT::kTahoeII)
	{
		/* database file */
		ModelFileT tahoeII;
		tahoeII.OpenWrite(file_name, true);

		/* coordinate data */
		tahoeII.PutCoordinates(*fCoordinates);

		/* element set data */
		for (int i = 0, id=1; i < fElementSets.Length(); i++)
		{
		  const ArrayT<StringT>& blockIDs = fElementSets[i]->BlockID();
		  for (int b=0; b < fElementSets[i]->NumBlocks(); b++, id++)
		    {
		      const iArray2DT* c = fElementSets[i]->Connectivities(blockIDs[b]);
		      iArray2DT conn = *c;
		      
		      iArrayT tmp(conn.Length(), conn.Pointer());
		      tmp++;
		      tahoeII.PutElementSet (id, conn);
		      tmp--;
		    }
		}
	}
	else if (format == IOBaseT::kExodusII)
	{
		/* total number of element blocks */
		int num_elem_blocks = 0;
		for (int i = 0; i < fElementSets.Length(); i++)
			num_elem_blocks += fElementSets[i]->NumBlocks();
	
		/* database file */
		ExodusT exo(cout);
		ArrayT<StringT> nothing;
		exo.Create(file_name, fTitle, nothing, nothing, fCoordinates->MinorDim(),
			fCoordinates->MajorDim(), NumElements(), num_elem_blocks, 0, 0);

		/* coordinate data */
		exo.WriteCoordinates(*fCoordinates);

		/* element set data */
		for (int i = 0, id=1; i < fElementSets.Length(); i++)
		{
			/* write connectivities */
		  const ArrayT<StringT>& blockIDs = fElementSets[i]->BlockID();
		  for (int b=0; b < fElementSets[i]->NumBlocks(); b++, id++)
		    {
				/* cast away const-ness to we can shift numbers */
				iArray2DT& connects = *((iArray2DT*) fElementSets[i]->Connectivities(blockIDs[b]));
	
				/* write to file */
				connects++;
				exo.WriteConnectivities(id, fElementSets[i]->Geometry(), connects);
				connects--;
		    }

		}
	}
	else
	{
		cout << "\n OutputBaseT::WriteGeometryFile: output format must be "
		     << IOBaseT::kTahoeII << " or " << IOBaseT::kExodusII << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void OutputBaseT::WriteOutput(double time, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
#pragma unused(time)
#pragma unused(ID)
#pragma unused(n_values)
#pragma unused(e_values)

	const char caller[] = "OutputBaseT::WriteOutput";

	/* checks */
	if (ID < 0 || ID >= fElementSets.Length())
		ExceptionT::GeneralFail(caller, "%d out of range: 0 < ID < %d", ID, fElementSets.Length());

	/* updtae output set */
	const OutputSetT& set = OutputSet(ID);
	if (set.Changing())
	{
		/* update nodes used */
		set.NodesUsed();
		
		/* update nodes used in blocks */
		const ArrayT<StringT>& block_ID = set.BlockID();
		for (int i = 0; i < block_ID.Length(); i++)
			set.BlockNodesUsed(block_ID[i]);
	}
	
	/* dimension checks */
	if (n_values.MinorDim() > 0 && set.NumNodes() != n_values.MajorDim()) 
		ExceptionT::SizeMismatch(caller, "expecting %d not %d nodes in set ID %d", set.NumNodes(), n_values.MajorDim(), ID);
	if (set.NumNodeValues() != n_values.MinorDim()) 
		ExceptionT::SizeMismatch(caller, "expecting %d not %d nodal values in set ID %d", set.NumNodes(), n_values.MinorDim(), ID);
	if (e_values.MinorDim() > 0 && set.NumElements() != e_values.MajorDim()) 
		ExceptionT::SizeMismatch(caller, "expecting %d not %d elements in set ID %d", set.NumElements(), e_values.MajorDim(), ID);
	if (set.NumElementValues() != e_values.MinorDim()) 
		ExceptionT::SizeMismatch(caller, "expecting %d not %d element values in set ID %d", set.NumElementValues(), e_values.MinorDim(), ID);

	/* set block ID to string names if possible 
	   else to the global block index position */
	CreateElementBlockIDs ();

	if (!fCoordinates)
		ExceptionT::GeneralFail(caller, "pointer to coordinates not set");

	/* increment the print step */
	fElementSets[ID]->IncrementPrintStep();
}

/* send data for output */
void OutputBaseT::WriteOutput(double time, int ID, const ArrayT<int>& nodes, const dArray2DT& n_values, 
		const ArrayT<int>& elements, const dArray2DT& e_values)
{
	/* get the output set */
	const OutputSetT& set = OutputSet(ID);
	
	/* assembly map */
	InverseMapT node_to_index;
	node_to_index.SetMap(set.NodesUsed());

	/* collect nodal values */
	dArray2DT n_out(set.NumNodes(), set.NumNodeValues());
	n_out = 0.0;
	for (int i = 0; i < nodes.Length(); i++)
		n_out.SetRow(node_to_index.Map(nodes[i]), n_values(i));

	/* collect element values */
	dArray2DT e_out(set.NumElements(), set.NumElementValues());
	e_out = 0.0;
	e_out.RowCollect(elements, e_values);
	
	/* write output */
	WriteOutput(time, ID, n_values, e_values);
}

/*************************************************************************
* Protected
*************************************************************************/

void OutputBaseT::LocalConnectivity(const iArrayT& node_map,
	const iArray2DT& connects, iArray2DT& local_connects) const
{
	/* sizes must match */
	if (connects.MajorDim() != local_connects.MajorDim() ||
	    connects.MinorDim() != local_connects.MinorDim()) throw ExceptionT::kSizeMismatch;

	/* quick exit - nothing to do */
	if (connects.MajorDim() == 0) return;

	/* generate inverse map */
	InverseMapT inv_node_map;
	inv_node_map.SetMap(node_map);

	/* generate local connects */
	int length = local_connects.Length();
	int* p_loc = local_connects.Pointer();
	const int* p_glb = connects.Pointer();
	for (int j = 0; j < length; j++)
		*p_loc++ = inv_node_map.Map(*p_glb++);
}

void OutputBaseT::ElementBlockValues (int ID, int block, const dArray2DT& allvalues, dArray2DT& blockvalues) const
{
  int length = fElementSets[ID]->NumBlockElements(fElementSets[ID]->BlockID(block));
  if (blockvalues.MajorDim() != length ||
      blockvalues.MinorDim() != allvalues.MinorDim()) throw ExceptionT::kSizeMismatch;

  /* find start point */
  int start = 0;
  for (int s=0; s < block; s++)
    start += fElementSets[ID]->NumBlockElements(fElementSets[ID]->BlockID(s));

  /* set row tags */
  iArrayT rows (length);
  rows.SetValueToPosition ();
  rows += start;

  /* copy certain rows */
  blockvalues.RowCollect (rows, allvalues);
}

void OutputBaseT::NodalBlockValues(int ID, int block, const dArray2DT& allvalues, dArray2DT& blockvalues, iArrayT& block_nodes) const
{
	const iArrayT& group_nodes = fElementSets[ID]->NodesUsed();
	block_nodes.Alias(fElementSets[ID]->BlockNodesUsed(fElementSets[ID]->BlockID(block)));
	if (block_nodes.Length() == group_nodes.Length()) /* nodes set nodes used by block */
	{
		blockvalues.Alias(allvalues);
	}
	else
	{
		/* block index to set index map */
		const iArrayT& index_map = fElementSets[ID]->BlockIndexToSetIndexMap(fElementSets[ID]->BlockID(block));
		
		/* collect block values */
		blockvalues.Dimension(index_map.Length(), allvalues.MinorDim());
		blockvalues.RowCollect(index_map, allvalues);	
    }
}

void OutputBaseT::ElementGroupBlockIndex (const StringT& n, int& g, int& b) const
{
  for (int i=0; i < fElementSets.Length(); i++)
    {
      const ArrayT<StringT>& ids = fElementSets[i]->BlockID ();
      for (int j=0; j < ids.Length(); j++)
	if (strncmp (ids[j].Pointer(), n.Pointer(), n.StringLength()) == 0)
	  {
	    g = i;
	    b = j;
	    return;
	  }
    }
}

void OutputBaseT::CreateElementBlockIDs (void)
{
  bool unique = true;
  int numblocks = 0;
  fElementBlockIDs.Dimension (fElementSets.Length());
  for (int i=0; i < fElementSets.Length(); i++)
    {
      const ArrayT<StringT>& bids = fElementSets[i]->BlockID();
      iArrayT& bints = fElementBlockIDs[i];
      bints.Dimension (bids.Length());
      if (unique)
	{
	  for (int j=0; j < bids.Length() && unique; j++)
	    {
	      bints[j] = atoi (bids[j]);

	      /* ZEROES ARE NOT ALLOWED BY SOME PROGRAMS */
	      if (bints[j] == 0) 
		bints[j]++;
	      
	      /* ensure uniqueness within the block */
	      for (int k=0; k < j && unique; k++)
		{
		  if (bints[k] == bints[j])
		    unique = false;
		}

	      /* ensure uniqueness to previous group blocks */
	      for (int g=0; g < i && unique; g++)
		{
		  iArrayT& gints = fElementBlockIDs[g];
		  for (int h=0; h < gints.Length(); h++)
		    {
		      if (gints[h] == bints[j]) 
			unique = false;
		    }
		}
	    }
	}
      
      if (!unique)
	{
	  bints.SetValueToPosition();
	  bints += 1 + numblocks;
	}

      numblocks += bids.Length();
    }
}

void OutputBaseT::String2IntIDs (const ArrayT<StringT>& s, iArrayT& i) const
{
  i.Dimension (s.Length());
  for (int j=0; j < s.Length(); j++)
    {
      i[j] = atoi (s[j]);

      /* ensure uniqueness */
      for (int k=0; k < j; k++)
	if (i[k] == i[j])
	  {
	    i.SetValueToPosition ();
	    i++;
	    return;
	  }
    }
}
