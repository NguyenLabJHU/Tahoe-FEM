/* $Id: OutputSetT.cpp,v 1.4 2001-12-16 23:57:06 paklein Exp $ */
/* created: paklein (03/07/2000) */

#include "OutputSetT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

/* array behavior */
const bool ArrayT<OutputSetT*>::fByteCopy = true;

/* constructor */
OutputSetT::OutputSetT(int ID, GeometryT::CodeT geometry_code,
		       const iArrayT& block_ID, 
		       const ArrayT<const iArray2DT*> connectivities, 
		       const ArrayT<StringT>& n_labels, 
		       const ArrayT<StringT>& e_labels, bool changing):
	fPrintStep(-1),
	fID(ID),
	fChanging(changing),
	fGeometry(geometry_code),
	fBlockID(block_ID),
	fConnectivities(connectivities.Length()),
	fBlockNodesUsed(fConnectivities.Length()),
	fBlockIndexToSetIndexMap(fConnectivities.Length())
{
	if (fConnectivities.Length() != fBlockID.Length()) 
	  {
	    cout << "\n\nOutputSetT::OutputSetT size mismatch: \n";
	    cout << " fConnectivities.Length = " << fConnectivities.Length();
	    cout << "\n    fBlockID.Length = " << fBlockID.Length() << "\n\n";
	    throw eSizeMismatch;
	  }

	for (int i=0; i < fConnectivities.Length(); i++)
	        fConnectivities[i] = connectivities[i];

	fNodeOutputLabels.Allocate(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = n_labels[i];

	fElementOutputLabels.Allocate(e_labels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
		fElementOutputLabels[j] = e_labels[j];

	/* set the nodes used array */
	fChanging = true; // force calculation of nodes used
	NodesUsed();
	for (int i = 0; i < fConnectivities.Length(); i++)
		BlockNodesUsed(i);
	fChanging = changing; // reset
}

OutputSetT::OutputSetT(const OutputSetT& source):
	fPrintStep(-1),
	fID(source.fID),
	fChanging(source.fChanging),
	fGeometry(source.fGeometry),
	fBlockID(source.fBlockID),
	fNodesUsed(source.fNodesUsed),
	fConnectivities(source.NumBlocks()),
	fBlockNodesUsed(fConnectivities.Length()),
	fBlockIndexToSetIndexMap(source.fBlockIndexToSetIndexMap)
{
	for (int i=0; i < fConnectivities.Length(); i++)
	        fConnectivities[i] = source.fConnectivities[i];

	fNodeOutputLabels.Allocate(source.fNodeOutputLabels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
		fNodeOutputLabels[i] = source.fNodeOutputLabels[i];

	fElementOutputLabels.Allocate(source.fElementOutputLabels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
		fElementOutputLabels[j] = source.fElementOutputLabels[j];

	if (fConnectivities.Length() != fBlockID.Length()) throw eSizeMismatch;
	
	/* set nodes used by blocks */
	if (fConnectivities.Length() == 1)
		fBlockNodesUsed[0].Alias(fNodesUsed);	
	else
		fBlockNodesUsed = source.fBlockNodesUsed;
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

const iArrayT& OutputSetT::BlockNodesUsed(int index)
{
	if (fChanging) /* need to reset data */
	{
		/* just one set */
		if (fBlockNodesUsed.Length() == 1)
		{
			fBlockNodesUsed[index].Alias(fNodesUsed);
			fBlockIndexToSetIndexMap[0].Allocate(fNodesUsed.Length());
			fBlockIndexToSetIndexMap[0].SetValueToPosition();		
		}
		else /* more than one block */
		{
			/* determine nodes used by block */
			SetNodesUsed(*fConnectivities[index], fBlockNodesUsed[index]);

			/* block to set index map */
			iArrayT& map = fBlockIndexToSetIndexMap[index];
			iArrayT& used = fBlockNodesUsed[index];
			map.Allocate(used.Length());
		
			/* range of nodes numbers */
			int min, max;
			fNodesUsed.MinMax(min, max);
			int range = max - min + 1;
			
			/* nodes used sequence */
			iArrayT sequence(range);
			sequence = -1;
			
			/* mark sequence */
			int dex = 0;
			for (int i = 0; i < fNodesUsed.Length(); i++)
				sequence[fNodesUsed[i] - min] = dex++;
					
			/* collect index list */
			for (int i = 0; i < map.Length(); i++)
			{
				int dex = sequence[used[i] - min];
				if (dex < 0) {
					cout << "\n OutputSetT::BlockNodesUsed: ERROR: block node used " << used[i]+1 
					     << " is not marked as used by the set" << endl;
					throw eGeneralFail;
				}
				else
					map[i] = dex;
			}
		}
	}
	
	return fBlockNodesUsed[index];
}

/* determine the nodes used */
void OutputSetT::SetNodesUsed(const iArray2DT& connects, iArrayT& nodes_used)
{
	nodes_used.Union(connects);
}

/* determine the nodes used */
void OutputSetT::SetNodesUsed(const ArrayT<const iArray2DT*>& connects_list, 
	iArrayT& nodes_used)
{
	/* compiler won't cast array type */
	ArrayT<const nArrayT<int>*> tmp(connects_list.Length(), 
		(const nArrayT<int>**) connects_list.Pointer());

	nodes_used.Union(tmp);
}
