/* $Id: OutputBaseT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1999)                                          */

#include "OutputBaseT.h"

#include "OutputSetT.h"

/* database types */
#include "ExodusT.h"
#include "ModelFileT.h"

#include "dArray2DT.h"
#include "iArray2DT.h"
#include "AutoArrayT.h"

OutputBaseT::OutputBaseT(ostream& out, const ArrayT<StringT>& out_strings):
	IOBaseT(out),
	fCoordinates(NULL),
	fNodeMap(NULL),
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
	fElementSets.Allocate(0);
}

const OutputSetT& OutputBaseT::OutputSet(int ID) const
{
	return *(fElementSets[ID]);
}

void OutputBaseT::NextTimeSequence(int sequence_number)
{
	fSequence  = sequence_number;
	for (int i = 0; i < fElementSets.Length(); i++)
		fElementSets[i]->ResetPrintStep();
}

/* set nodal coordinates */
void OutputBaseT::SetCoordinates(const dArray2DT& coordinates,
	const iArrayT* node_map)
{
	fCoordinates = &coordinates;
	fNodeMap = node_map;
}

int OutputBaseT::AddElementSet(const OutputSetT& output_set)
{
	OutputSetT* copy = new OutputSetT(output_set);
	if (!copy) throw eOutOfMemory;
	fElementSets.Append(copy);

	/* ID is just position in array */
	return fElementSets.Length() - 1;
}

int OutputBaseT::NumElements(void) const
{
	int count = 0;
	for (int i = 0; i < fElementSets.Length(); i++)
		count += fElementSets[i]->NumElements();
	return count;
}

void OutputBaseT::AddNodeSet(const iArrayT& nodeset, int setID)
{
	fNodeSets.Append(&nodeset);
	fNodeSetIDs.Append (setID);
}

void OutputBaseT::AddSideSet(const iArray2DT& sideset, int setID, int group_ID)
{
	/* check */
	if (group_ID < 0 || group_ID >= fElementSets.Length())
	{
		cout << "\n OutputBaseT::AddSideSet: group ID out of range: " << group_ID << endl;
		throw eOutOfRange;
	}

	fSideSets.Append(&sideset);
	fSideSetIDs.Append (setID);
	fSSGroupID.Append(group_ID);
}


/* output functions */
void OutputBaseT::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT format) const
{
	if (!fCoordinates)
	{
		cout << "\n OutputBaseT::WriteGeometryFile: pointer to coordinates not set" << endl;
		throw eGeneralFail;
	}

	if (format == IOBaseT::kTahoeII)
	{
		/* database file */
		ModelFileT tahoeII;
		tahoeII.OpenWrite(file_name, true);

		/* coordinate data */
		tahoeII.PutCoordinates(*fCoordinates);

		/* element set data */
		for (int i = 0; i < fElementSets.Length(); i++)
		{
			const iArray2DT& connects = fElementSets[i]->Connectivities();
			iArrayT tmp(connects.Length(), connects.Pointer());
			tmp++;
			tahoeII.PutElementSet(i, connects);
			tmp--;
		}
	}
	else if (format == IOBaseT::kExodusII)
	{
		/* database file */
		ExodusT exo(cout);
		ArrayT<StringT> nothing;
		exo.Create(file_name, fTitle, nothing, nothing, fCoordinates->MinorDim(),
			fCoordinates->MajorDim(), NumElements(), fElementSets.Length(), 0, 0);

		/* coordinate data */
		exo.WriteCoordinates(*fCoordinates);

		/* element set data */
		for (int i = 0; i < fElementSets.Length(); i++)
		{
			/* write connectivities */
			const iArray2DT& connects = fElementSets[i]->Connectivities();
			iArrayT tmp(connects.Length(), connects.Pointer());
			tmp++;
			exo.WriteConnectivities(i+1, fElementSets[i]->Geometry(), connects); // ID cannot be 0
			tmp--;
		}
	}
	else
	{
		cout << "\n OutputBaseT::WriteGeometryFile: output format must be "
		     << IOBaseT::kTahoeII << " or " << IOBaseT::kExodusII << endl;
		throw eGeneralFail;
	}
}

void OutputBaseT::WriteOutput(double time, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
#pragma unused(time)
#pragma unused(ID)
#pragma unused(n_values)
#pragma unused(e_values)

	/* checks */
	if (ID < 0 || ID >= fElementSets.Length())
	{
		cout << "\n OutputBaseT::WriteOutput: set ID " << ID
		     << " is out of range {" << 0 << "," << fElementSets.Length()
		     << "}" << endl;
		throw eOutOfRange;
	}

	if (!fCoordinates)
	{
		cout << "\n OutputBaseT::WriteOutput: pointer to coordinates not set" << endl;
		throw eGeneralFail;
	}

	/* increment the print step */
	fElementSets[ID]->IncrementPrintStep();
}

/*************************************************************************
* Protected
*************************************************************************/

void OutputBaseT::LocalConnectivity(const iArrayT& node_map,
	const iArray2DT& connects, iArray2DT& local_connects) const
{
	/* sizes must match */
	if (connects.MajorDim() != local_connects.MajorDim() ||
	    connects.MinorDim() != local_connects.MinorDim()) throw eSizeMismatch;

	/* compressed number range */
	int max, shift;
	node_map.MinMax(shift, max);
	int range = max - shift + 1;

	/* generate inverse map */
	iArrayT inv_node_map(range);
	inv_node_map = -1;
	
	int nmap = node_map.Length();
	int* pmap = node_map.Pointer();
	for (int i = 0; i < nmap; i++)
		inv_node_map[*pmap++ - shift] = i;

	/* generate local connects */
	int length = local_connects.Length();
	int* p_loc = local_connects.Pointer();
	int* p_glb = connects.Pointer();
	for (int j = 0; j < length; j++)
		*p_loc++ = inv_node_map[*p_glb++ - shift];
}

