/* $Id: OutputSetT.cpp,v 1.25 2005-06-06 06:38:24 paklein Exp $ */
/* created: paklein (03/07/2000) */
#include "OutputSetT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

using namespace Tahoe;

namespace Tahoe {
/* array behavior */
DEFINE_TEMPLATE_STATIC const bool ArrayT<OutputSetT*>::fByteCopy = true;
}

/* constructor */
OutputSetT::OutputSetT(GeometryT::CodeT geometry_code,
	const ArrayT<StringT>& block_ID, 
	const ArrayT<const iArray2DT*>& connectivities, 
	const ArrayT<StringT>& n_labels, 
	const ArrayT<StringT>& e_labels, bool changing):
	fMode(kElementBlock),
	fPrintStep(-1),
	fID("1"), /* dummy ID */
	fChanging(changing),
	fGeometry(geometry_code),
	fBlockID(block_ID),
	fConnectivities(connectivities),
	fBlockNodesUsed(fConnectivities.Length()),
	fBlockIndexToSetIndexMap(fConnectivities.Length()),
	fPoints(NULL)
{
	if (fConnectivities.Length() != fBlockID.Length()) 
		ExceptionT::SizeMismatch("OutputSetT::OutputSetT",
			"fConnectivities.Length = %d, fBlockID.Length = %d",
			fConnectivities.Length(), fBlockID.Length());

	fNodeOutputLabels.Dimension(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
	  {
		fNodeOutputLabels[i] = n_labels[i];
		fNodeOutputLabels[i].Replace (' ', '_');
	  }

	fElementOutputLabels.Dimension(e_labels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
	  {
		fElementOutputLabels[j] = e_labels[j];
		fElementOutputLabels[j].Replace (' ', '_');
	  }

	/* initialize memory managers */
	fNodesUsed_man.SetWard(0, fNodesUsed);
	if (fConnectivities.Length() > 1) {
		fBlockNodesUsed_man.Dimension(fBlockNodesUsed.Length());
		for (int i = 0; i < fBlockNodesUsed_man.Length(); i++)
			fBlockNodesUsed_man[i].SetWard(0, fBlockNodesUsed[i]);
	}
	fBlockIndexToSetIndexMap_man.Dimension(fBlockIndexToSetIndexMap.Length());
	for (int i = 0; i < fBlockIndexToSetIndexMap_man.Length(); i++)
		fBlockIndexToSetIndexMap_man[i].SetWard(0, fBlockIndexToSetIndexMap[i]);

	/* set the nodes used arrays */
	fChanging = true; // force calculation of nodes used
	NodesUsed();
	if (fConnectivities.Length() > 1) {
		for (int i = 0; i < fConnectivities.Length(); i++)
			BlockNodesUsed(fBlockID[i]);
	} 
	else 
		fBlockNodesUsed[0].Alias(fNodesUsed);
	fChanging = changing; // reset
}

OutputSetT::OutputSetT(GeometryT::CodeT geometry_code,
	const iArray2DT& connectivities, const ArrayT<StringT>& n_labels, bool changing):
	fMode(kFreeSet),
	fPrintStep(-1),
	fID("1"), /* dummy ID */
	fChanging(false),
	fGeometry(geometry_code),
	fBlockID(1),
	fConnectivities(1),
	fBlockNodesUsed(1),
	fBlockIndexToSetIndexMap(1),
	fPoints(NULL)
{
	/* keep reference to connectivities */
	fConnectivities[0] = &connectivities;
	fBlockID[0] = fID; /* must give connectivities a reasonable ID for compatibility
	                    * with the output classes */

	/* copy node labels */
	fNodeOutputLabels.Dimension(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
	  {
		fNodeOutputLabels[i] = n_labels[i];
		fNodeOutputLabels[i].Replace (' ', '_');
	  }

	/* initialize memory managers */
	fNodesUsed_man.SetWard(0, fNodesUsed);
	fBlockIndexToSetIndexMap_man.Dimension(1);
	fBlockIndexToSetIndexMap_man[0].SetWard(0, fBlockIndexToSetIndexMap[0]);

	/* set the nodes used array */
	fChanging = true; // force calculation of nodes used
	NodesUsed();
	fBlockNodesUsed[0].Alias(fNodesUsed);
	fChanging = changing;
}

/* output data record for a set of points */
OutputSetT::OutputSetT(const iArrayT& points, const ArrayT<StringT>& n_labels, bool changing):
	fMode(kFreeSet),
	fPrintStep(-1),
	fID("1"), /* dummy ID */
	fChanging(false),
	fGeometry(GeometryT::kPoint),
	fBlockID(1),
	fConnectivities(1),
	fBlockNodesUsed(1),
	fBlockIndexToSetIndexMap(1),
	fPoints(&points),
	fConnects2D(fPoints->Length(), 1, fPoints->Pointer())
{
	/* keep reference to connectivities */
	fConnectivities[0] = &fConnects2D;
	fBlockID[0] = fID; /* must give connectivities a reasonable ID for compatibility
	                    * with the output classes */

	/* copy node labels */
	fNodeOutputLabels.Dimension(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
	  {
		fNodeOutputLabels[i] = n_labels[i];
		fNodeOutputLabels[i].Replace (' ', '_');
	  }

	/* initialize memory managers */
	fNodesUsed_man.SetWard(0, fNodesUsed);
	fBlockIndexToSetIndexMap_man.Dimension(1);
	fBlockIndexToSetIndexMap_man[0].SetWard(0, fBlockIndexToSetIndexMap[0]);

	/* set the nodes used array */
	fChanging = true; // force calculation of nodes used
	NodesUsed();
	fBlockNodesUsed[0].Alias(fNodesUsed);
	fChanging = changing;
}

OutputSetT::OutputSetT(GeometryT::CodeT geometry_code,
		       const ArrayT<StringT>& block_ID, const ArrayT<StringT>& sideset_ID,
	const ArrayT<const iArray2DT*>& connectivities, 
	const ArrayT<StringT>& n_labels, 
	const ArrayT<StringT>& e_labels, bool changing):
	fMode(kElementFromSideSet),
	fPrintStep(-1),
	fID("1"), /* dummy ID */
	fChanging(changing),
	fGeometry(geometry_code),
	fBlockID(block_ID),
	fSSID(sideset_ID),
	fConnectivities(connectivities),
	fBlockNodesUsed(fConnectivities.Length()),
	fBlockIndexToSetIndexMap(fConnectivities.Length()),
	fPoints(NULL)
{
	if (fConnectivities.Length() != fBlockID.Length() &&
		fBlockID.Length() != fSSID.Length())
		ExceptionT::SizeMismatch("OutputSetT::OutputSetT",
			"fConnectivities.Length = %d, fBlockID.Length = %d, fSSID.Length = %d",
			fConnectivities.Length(), fBlockID.Length(), fSSID.Length());

	fNodeOutputLabels.Dimension(n_labels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
	  {
		fNodeOutputLabels[i] = n_labels[i];
		fNodeOutputLabels[i].Replace (' ', '_');
	  }

	fElementOutputLabels.Dimension(e_labels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
	  {
		fElementOutputLabels[j] = e_labels[j];
		fElementOutputLabels[j].Replace (' ', '_');
	  }

	/* initialize memory managers */
	fNodesUsed_man.SetWard(0, fNodesUsed);
	if (fConnectivities.Length() > 1) {
		fBlockNodesUsed_man.Dimension(fBlockNodesUsed.Length());
		for (int i = 0; i < fBlockNodesUsed_man.Length(); i++)
			fBlockNodesUsed_man[i].SetWard(0, fBlockNodesUsed[i]);
	}
	fBlockIndexToSetIndexMap_man.Dimension(fBlockIndexToSetIndexMap.Length());
	for (int i = 0; i < fBlockIndexToSetIndexMap_man.Length(); i++)
		fBlockIndexToSetIndexMap_man[i].SetWard(0, fBlockIndexToSetIndexMap[i]);

	/* set the nodes used array */
	fChanging = true; // force calculation of nodes used
	NodesUsed();
	if (fConnectivities.Length() > 1) {
		for (int i = 0; i < fConnectivities.Length(); i++)
			BlockNodesUsed(fBlockID[i]);
	} 
	else 
		fBlockNodesUsed[0].Alias(fNodesUsed);
	fChanging = changing; // reset
}

OutputSetT::OutputSetT(const OutputSetT& source):
	fMode(source.fMode),
	fPrintStep(-1),
	fID(source.fID),
	fChanging(source.fChanging),
	fGeometry(source.fGeometry),
	fBlockID(source.fBlockID),
	fSSID(source.fSSID),
	fConnectivities(source.NumBlocks()),
	fBlockNodesUsed(fConnectivities.Length()),
	fBlockIndexToSetIndexMap(fConnectivities.Length()),
	fPoints(source.fPoints)
{
	if (!fPoints) {
		for (int i=0; i < fConnectivities.Length(); i++)
			fConnectivities[i] = source.fConnectivities[i];
	} else { /* data over list of points */
		if (fConnectivities.Length() > 1) 
			ExceptionT::GeneralFail("OutputSetT::OutputSetT",
				"expecting 1 block not %d", fConnectivities.Length());
		fConnects2D.Alias(fPoints->Length(), 1, fPoints->Pointer());
		fConnectivities[0] = &fConnects2D;
	}

	fNodeOutputLabels.Dimension(source.fNodeOutputLabels.Length());
	for (int i = 0; i < fNodeOutputLabels.Length(); i++)
	  {
		fNodeOutputLabels[i] = source.fNodeOutputLabels[i];
		fNodeOutputLabels[i].Replace (' ', '_');
	  }

	fElementOutputLabels.Dimension(source.fElementOutputLabels.Length());
	for (int j = 0; j < fElementOutputLabels.Length(); j++)
	  {
		fElementOutputLabels[j] = source.fElementOutputLabels[j];
		fElementOutputLabels[j].Replace (' ', '_');
	  }

	if (fMode == kElementBlock &&
	    fConnectivities.Length() != fBlockID.Length()) ExceptionT::SizeMismatch("OutputSetT::OutputSetT");

	/* initialize memory managers */
	fNodesUsed_man.SetWard(0, fNodesUsed);
	fNodesUsed_man.SetLength(source.fNodesUsed.Length(), false);
	fNodesUsed = source.fNodesUsed;

	if (fConnectivities.Length() > 1) {
		fBlockNodesUsed_man.Dimension(fBlockNodesUsed.Length());
		for (int i = 0; i < fBlockNodesUsed_man.Length(); i++) {
			fBlockNodesUsed_man[i].SetWard(0, fBlockNodesUsed[i]);
			fBlockNodesUsed_man[i].SetLength(source.fBlockNodesUsed[i].Length(), false);
			fBlockNodesUsed[i] = source.fBlockNodesUsed[i];
		}
	}
	else
		fBlockNodesUsed[0].Alias(fNodesUsed);

	fBlockIndexToSetIndexMap_man.Dimension(fBlockIndexToSetIndexMap.Length());
	for (int i = 0; i < fBlockIndexToSetIndexMap_man.Length(); i++) {
		fBlockIndexToSetIndexMap_man[i].SetWard(0, fBlockIndexToSetIndexMap[i]);
		fBlockIndexToSetIndexMap_man[i].SetLength(source.fBlockIndexToSetIndexMap[i].Length(), false);
		fBlockIndexToSetIndexMap[i] = source.fBlockIndexToSetIndexMap[i];
	}
}

/* dimensions */
int OutputSetT::NumBlockElements(const StringT& ID) const
{
	return fConnectivities[BlockIndex(ID)]->MajorDim();
}

int OutputSetT::NumElements(void) const 
{
  int num = 0;
  for (int i=0; i < fConnectivities.Length(); i++)
    num += fConnectivities[i]->MajorDim();
  return num; 
}

void OutputSetT::SetID(const StringT& id)
{
	fID = id;
	/* must give connectivities a reasonable ID for compatibility with the 
	 * output classes */
	if (Mode() == kFreeSet) {
		fBlockID[0] = fID;
	}
}

const iArray2DT* OutputSetT::Connectivities(const StringT& ID) const
{
	return fConnectivities[BlockIndex(ID)];
}

const iArrayT& OutputSetT::BlockNodesUsed(const StringT& ID) const
{
	int index = BlockIndex(ID);
	if (fChanging) /* need to reset data */
	{
		/* not so const */
		OutputSetT* non_const_this = (OutputSetT*) this;

		/* reset alias */
		if (fPoints) (non_const_this->fConnects2D).Alias(fPoints->Length(), 1, fPoints->Pointer());
	
		/* just one set fBlockNodesUsed[0] is alias to fNodesUsed */
		if (fBlockNodesUsed.Length() == 1)
		{
			/* reset fNodesUsed */
			non_const_this->NodesUsed();
			non_const_this->fBlockNodesUsed[0].Alias(fNodesUsed);
			(non_const_this->fBlockIndexToSetIndexMap_man)[0].SetLength(fNodesUsed.Length(), false);
			(non_const_this->fBlockIndexToSetIndexMap)[0].SetValueToPosition();		
		}
		else /* more than one block */
		{
			/* determine nodes used by block */
			SetNodesUsed(*fConnectivities[index], (non_const_this->fBlockNodesUsed)[index],
				(non_const_this->fBlockNodesUsed_man)[index]);
			
			/* block could be empty */
			if (fBlockNodesUsed[index].Length() > 0)
			{
				/* block to set index map */
				iArrayT& map = (non_const_this->fBlockIndexToSetIndexMap)[index];
				const iArrayT& used = fBlockNodesUsed[index];
				(non_const_this->fBlockIndexToSetIndexMap_man)[index].SetLength(used.Length(), false);
		
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
					if (dex < 0)
						ExceptionT::GeneralFail("OutputSetT::BlockNodesUsed",
							"block node used %d is not marked as used by the set", used[i]+1);
					else
						map[i] = dex;
				}
			}
		}
	}
	
	return fBlockNodesUsed[index];
}

/* determine the nodes used */
void OutputSetT::SetNodesUsed(const iArray2DT& connects, iArrayT& nodes_used, 
	VariArrayT<int>& nodes_used_man) const
{
	/* collect union */
	iArrayT tmp;
	tmp.Union(connects);
	
	/* copy */
	nodes_used_man.SetLength(tmp.Length(), false);
	nodes_used = tmp;
}

/* determine the nodes used */
void OutputSetT::SetNodesUsed(const ArrayT<const iArray2DT*>& connects_list, 
	iArrayT& nodes_used, VariArrayT<int>& nodes_used_man)
{
	/* compiler won't cast array type */
	ArrayT<const nArrayT<int>*> tmp(connects_list.Length(), 
		(const nArrayT<int>**) connects_list.Pointer());

	/* collect union */
	iArrayT i_tmp;
	i_tmp.Union(tmp);
	
	/* copy */
	nodes_used_man.SetLength(i_tmp.Length(), false);
	nodes_used = i_tmp;
}

/*************************************************************************
* Private
*************************************************************************/

/* returns the index for the element block for the given ID */
int OutputSetT::BlockIndex(const StringT& ID) const
{
	if (fMode == kFreeSet)
		return 0;
	else
	{
		int index = -1;
		for (int i = 0; index == -1 && i < fBlockID.Length(); i++)
			if (fBlockID[i] == ID)
				index = i;

		if (index == -1)
			ExceptionT::GeneralFail("OutputSetT::BlockIndex",
				"block ID %s not found", ID.Pointer());

		return index;
	}
}

/* returns the index for the side set for the given ID */
/*int OutputSetT::SideSetIndex(const StringT& ID) const
{
	if (fMode == kFreeSet || fMode == kElementBlock)
		return 0;
	else
	{
		int index = -1;
		for (int i = 0; index == -1 && i < fSSID.Length(); i++)
			if (fSSID[i] == ID)
				index = i;

		if (index == -1) {
			cout << "\n OutputSetT::SideSetIndex: side set ID not found: " << ID << endl;
			throw ExceptionT::kGeneralFail;
		}
		return index;
	}
}*/
