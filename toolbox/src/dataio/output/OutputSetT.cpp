/* $Id: OutputSetT.cpp,v 1.3.2.2 2001-10-28 23:40:53 paklein Exp $ */
/* created: paklein (03/07/2000)                                          */

#include "OutputSetT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

/* array behavior */
const bool ArrayT<OutputSetT*>::fByteCopy = true;

/* constructor */
OutputSetT::OutputSetT(int ID, GeometryT::CodeT geometry_code,
	const iArrayT& block_ID, const ArrayT<const iArray2DT*> connectivities, 
	const ArrayT<StringT>& n_labels, const ArrayT<StringT>& e_labels, bool changing):
	fPrintStep(-1),
	fID(ID),
	fChanging(changing),
	fGeometry(geometry_code),
	fBlockID(block_ID),
	fConnectivities (connectivities.Length())
{
	for (int i=0; i < fConnectivities.Length(); i++)
	        fConnectivities[i] = connectivities[i];

	fNodeOutputLabels.Allocate(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = n_labels[i];

	fElementOutputLabels.Allocate(e_labels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
		fElementOutputLabels[j] = e_labels[j];

	if (fConnectivities.Length() != fBlockID.Length()) 
	  throw eSizeMismatch;

	/* set the nodes used array */
	if (!fChanging) SetNodesUsed();
}

OutputSetT::OutputSetT(const OutputSetT& source):
	fPrintStep(-1),
	fID(source.fID),
	fChanging(source.fChanging),
	fGeometry(source.fGeometry),
	fBlockID(source.fBlockID),
	fNodesUsed(source.fNodesUsed),
	fConnectivities (source.NumBlocks())
{
	for (int i=0; i < fConnectivities.Length(); i++)
	        fConnectivities[i] = source.fConnectivities[i];

	fNodeOutputLabels.Allocate(source.fNodeOutputLabels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = source.fNodeOutputLabels[i];

	fElementOutputLabels.Allocate(source.fElementOutputLabels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
		fElementOutputLabels[j] = source.fElementOutputLabels[j];

	if (fConnectivities.Length() != fBlockID.Length()) 
	  throw eSizeMismatch;
}

/* dimensions */
int OutputSetT::NumBlockElements (int index) const
{
  if (index < 0 || index >= fConnectivities.Length())
    throw eOutOfRange;
	return fConnectivities[index]->MajorDim();
}

int OutputSetT::NumElements(void) const 
{
  int num = 0;
  for (int i=0; i < fConnectivities.Length(); i++)
    num += fConnectivities[i]->MajorDim();
  return num; 
}

const iArray2DT* OutputSetT::Connectivities(int index) const
{
  if (index < 0 || index >= fConnectivities.Length())
    throw eOutOfRange;
	return fConnectivities[index];
}

//TEMP - used to write all set connectivities at once
#if 0
void  OutputSetT::AllConnectivities (iArray2DT& connects) const
{
  if (fConnectivities.Length() == 1)
    connects = *(fConnectivities[0]);
  else
    {
      int row = 0;
      for (int i=0; i < fConnectivities.Length(); i++)
	{
	  connects.BlockRowCopyAt (*(fConnectivities[i]), row);
	  row += fConnectivities[i]->MajorDim();
	}
    }
}
#endif

/* set nodes used */
void OutputSetT::SetNodesUsed(void)
{
	/* quick exit */
  //if (fConnectivities.Length() == 0) return;

       int num_blocks = fConnectivities.Length();

       iArrayT mins (num_blocks);
       iArrayT maxes (num_blocks);
       for (int i=0; i < fConnectivities.Length(); i++)
	 {
	   mins[i] = fConnectivities[i]->Min();
	   maxes[i] = fConnectivities[i]->Max();
	 }

	/* compressed number range */
	int min   = mins.Min();
	int max   = maxes.Max();
	int range = max - min + 1;

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int b=0; b < num_blocks; b++)
	  {
	    const iArray2DT* conn = fConnectivities[b];
	    int *pc = conn->Pointer();
	    for (int i = 0; i < conn->Length(); i++)
	      node_map[*pc++ - min] = 1;
	  }

	/* collect list */
	fNodesUsed.Allocate(node_map.Count(1));
	int dex = 0;
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) fNodesUsed[dex++] = j + min;
}
