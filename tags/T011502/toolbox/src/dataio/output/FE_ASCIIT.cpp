/* $Id: FE_ASCIIT.cpp,v 1.4 2002-01-09 12:15:28 paklein Exp $ */
/* created: sawimme (05/20/1999) */

#include "FE_ASCIIT.h"

#include "GeometryT.h"
#include "OutputSetT.h"
#include "ModelFileT.h"
#include "ofstreamT.h"

#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

/* constructor */
FE_ASCIIT::FE_ASCIIT(ostream& out, bool external, const ArrayT<StringT>& out_strings):
	OutputBaseT(out, out_strings),
	fExternTahoeII(external),
	fInitGeom(false),
	fInitRun(false)
{

}

/* increment sequence, create new output file series */
void FE_ASCIIT::NextTimeSequence(int sequence_number)
{
	/* inherited */
	OutputBaseT::NextTimeSequence(sequence_number);
	
	/* reset flags */
	fInitGeom = false;
	fInitRun  = false;
}

//NOTE: to write a geometry definition file
void FE_ASCIIT::WriteGeometry(void)
{
	ModelFileT mf;
	StringT filename = fOutroot;

	/* changing geometry */
	bool change = false;
	for (int j=0; j < fElementSets.Length() && !change; j++)
	  if (fElementSets[j]->Changing()) change = true;
	if (change)
	  filename.Append(".ps", fElementSets[0]->PrintStep() + 1);
	filename.Append (".geom");
	mf.OpenWrite (filename, fExternTahoeII);
	
	mf.PutTitle (fTitle);
	mf.PutCoordinates (*fCoordinates);
	
	for (int e=0; e < fElementSets.Length(); e++)
	  {
	    const iArrayT& blockIDs = fElementSets[e]->BlockID();
	    for (int b=0; b < fElementSets[e]->NumBlocks(); b++)
	      {
		const iArray2DT* c = fElementSets[e]->Connectivities(b);
		iArray2DT conn = *c;
		
		iArrayT tmp(conn.Length(), conn.Pointer());
		tmp++;
		mf.PutElementSet (blockIDs[b], conn);
		tmp--;
	      }
	  }
	
	for (int n=0; n < fNodeSets.Length(); n++)
	  {
	    iArrayT& set = *((iArrayT*) fNodeSets[n]);
	    set++;
	    mf.PutNodeSet (fNodeSetIDs[n], set);
	    set--;
	  }

	for (int s=0; s < fSideSets.Length(); s++)
	  {
	    iArray2DT& set = *((iArray2DT*) fSideSets[s]);
	    set++;
	    int block_ID = fElementSets[fSSGroupID[s]]->ID();
	    mf.PutSideSet (fSideSetIDs[s], block_ID, set);
	    set--;
	  }

	mf.Close(); // datafile actually written here
}

void FE_ASCIIT::WriteOutput(double time, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	/* inherited */
	OutputBaseT::WriteOutput(time, ID, n_values, e_values);

	/* geometry data */
	if (fElementSets[ID]->PrintStep() == 0 ||
		fElementSets[ID]->Changing())
	{	
		/* file name */
		StringT geom_file(fOutroot);
		if (fSequence > 0) geom_file.Append(".seq", fSequence + 1);
		geom_file.Append(".geo");

		/* open stream */
		ofstreamT out;
		SetStreamPrefs(out);
		if (!fInitGeom)
		{
			/* initialize geometry file */
			out.open(geom_file);
			out << "\n G E O M E T R Y   D A T A:\n\n";
	
			/* set flag */
			fInitGeom = true;
		}
		else /* re-open file */
			out.open_append(geom_file);

		/* check */
		if (!out.is_open())
		{
			cout << "\n FE_ASCIIT::WriteOutput: error opening file: " << geom_file << endl;
			throw eGeneralFail;
		}

		/* ID information */
		out << " Group number. . . . . . . . . . . . . . . . . . = "
		     << fElementSets[ID]->ID() << '\n';
		out << " Output ID . . . . . . . . . . . . . . . . . . . = "
		     << ID << '\n';
		out << " Number of Blocks. . . . . . . . . . . . . . . . = "
		     << fElementSets[ID]->NumBlocks() << '\n';
		if (fElementSets[ID]->Changing())
		{
			out << " Print Step. . . . . . . . . . . . . . . . . . . = "
			    << fElementSets[ID]->PrintStep() << '\n';
			out << " Time. . . . . . . . . . . . . . . . . . . . . . = "
			    << time << '\n';
		}
	
		/* write geometry */
		WriteGeometryData(out, ID);
	}

	/* file name */
	StringT dat_file(fOutroot);
	if (fSequence > 0) dat_file.Append(".seq", fSequence + 1);
	dat_file.Append(".run");

	/* open stream */
	ofstreamT out;
	SetStreamPrefs(out);
	if (!fInitRun)
	{
		/* initialize geometry file */
		out.open(dat_file);
		out << "\n O U T P U T   D A T A :\n\n";
	
		/* set flag */
		fInitRun = true;
	}
	else /* re-open file */
		out.open_append(dat_file);

	/* check */
	if (!out.is_open())
	{
		cout << "\n FE_ASCIIT::WriteOutput: error opening file: " << dat_file << endl;
		throw eGeneralFail;
	}

	/* print header */
	out << "\n Group number. . . . . . . . . . . . . . . . . . = "
	    << fElementSets[ID]->ID() << '\n';
	out << " Output ID . . . . . . . . . . . . . . . . . . . = "
	    << ID << '\n';
	out << " Print Step. . . . . . . . . . . . . . . . . . . = "
	    << fElementSets[ID]->PrintStep() << '\n';
	out << " Time. . . . . . . . . . . . . . . . . . . . . . = "
	    << time << '\n';
	out << " Number of Blocks. . . . . . . . . . . . . . . . = "
	    << fElementSets[ID]->NumBlocks() << '\n';

	/* write data */
	WriteOutputData(out, ID, n_values, e_values);
}

/*************************************************************************
* Private
*************************************************************************/

void FE_ASCIIT::WriteGeometryData(ostream& out, int ID)
{
	/* dimensions */
	int nsd = fCoordinates->MinorDim();

	/* generate coordinate labels */
	ArrayT<StringT> coord_labels(nsd);
	for (int i = 0; i < nsd; i++)
	{
		coord_labels[i] = "x";
		coord_labels[i].Append(i+1);
	}

	const iArrayT& blockIDs = fElementSets[ID]->BlockID();
	for (int b=0; b < fElementSets[ID]->NumBlocks(); b++)
	  {
	    /* write coordinates */
	    const iArrayT& nodes_used = fElementSets[ID]->BlockNodesUsed(b);
	    dArray2DT local_coordinates(nodes_used.Length(), nsd);
	    local_coordinates.RowCollect(nodes_used, *fCoordinates);
	
	    out << "\n Nodal coordinates:\n";	
	    out << " Block number . . . .  . . . . . . . . . . . . . = "
		<< blockIDs[b] << '\n';
	    WriteNodeHeader(out, local_coordinates.MajorDim(), coord_labels);
	    WriteNodeValues(out, nodes_used, local_coordinates);
	
	    /* write connectivities */
	    const iArray2DT* c = fElementSets[ID]->Connectivities(b);
	    out << "\n Connectivities:\n";
	    out << " Block number . . . .  . . . . . . . . . . . . . = "
		<< blockIDs[b] << '\n';
	    out << " Number of elements .  . . . . . . . . . . . . . = "
		<< c->MajorDim() << '\n';
	    out << " Number of element nodes . . . . . . . . . . . . = "
		<< c->MinorDim() << '\n';
	    out << " Geometry code . . . . . . . . . . . . . . . . . = "
		<< fElementSets[ID]->Geometry() << "\n\n";
	    out << setw(kIntWidth) << "element";
	    for (int j = 0; j < c->MinorDim(); j++)
	      out << setw(kIntWidth - 1) << "n" << j+1;
	    out << '\n';
			
	    /* correct offset for output */
	    iArray2DT connects_temp;
	    connects_temp.Alias(*c);	
	    connects_temp++;
	    connects_temp.WriteNumbered(out);
	    connects_temp--;
	  }
	out.flush();
}

void FE_ASCIIT::WriteOutputData(ostream& out, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	const iArrayT& blockIDs = fElementSets[ID]->BlockID();
	for (int b=0; b < fElementSets[ID]->NumBlocks(); b++)
	  {
	    /* write node header */
	    out << "\n Nodal data:\n";	
	    out << " Block number . . . .  . . . . . . . . . . . . . = "
		<< blockIDs[b] << '\n';
	    const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();

	    /* write node vars */
	    if (n_values.MajorDim () > 0)
		{
			iArrayT nodes_used;
			dArray2DT blockvals;
			NodalBlockValues (ID, b, n_values, blockvals, nodes_used);

			WriteNodeHeader(out, nodes_used.Length(), node_labels);
			WriteNodeValues(out, nodes_used, blockvals);
		}
		else
			WriteNodeHeader(out, 0, node_labels);


	    /* write element header */
	    out << "\n Element data:\n";
	    out << " Block number . . . .  . . . . . . . . . . . . . = "
		<< blockIDs[b] << '\n';
	    const ArrayT<StringT>& elem_labels = fElementSets[ID]->ElementOutputLabels();


	    /* write element values */
	    if (e_values.MajorDim() > 0)
		{
			/* collect block values */
			dArray2DT local_vals(fElementSets[ID]->NumBlockElements (b), e_values.MinorDim());
			ElementBlockValues(ID, b, e_values, local_vals);
			
			/* write */
			WriteElementHeader(out, local_vals.MajorDim(), elem_labels);
			WriteElementValues(out, local_vals);
		}
		else
			WriteElementHeader(out, 0, elem_labels);
	}
	out.flush();
}

void FE_ASCIIT::WriteNodeHeader(ostream& out, int num_output_nodes,
	const ArrayT<StringT>& labels) const
{
	double* junk = NULL;
	int d_width = OutputWidth(out, junk);

	out << " Number of nodal points. . . . . . . . . . . . . = "
	    << num_output_nodes << '\n';
	out << " Number of values. . . . . . . . . . . . . . . . = "
		<< labels.Length() << "\n\n";

	if (labels.Length())
	{
		out << setw(kIntWidth) << "node";
		for (int i = 0; i < labels.Length(); i++)
			out << setw(d_width) << labels[i];
		out << '\n';
	}
}

void FE_ASCIIT::WriteElementHeader(ostream& out, int num_output_elems,
	const ArrayT<StringT>& labels) const
{
	double* junk = NULL;
	int d_width = OutputWidth(out, junk);

	out << " Number of elements. . . . . . . . . . . . . . . = "
	    << num_output_elems << '\n';
	out << " Number of values. . . . . . . . . . . . . . . . = "
	    << labels.Length() << "\n\n";

	if (labels.Length() > 0)
	{
		out << setw(kIntWidth) << "element";
		for (int i = 0; i < labels.Length(); i++)
			out << setw(d_width) << labels[i];
		out << '\n';
	}
}

void FE_ASCIIT::WriteNodeValues(ostream& out, const iArrayT& node_numbers,
	const dArray2DT& values) const
{
	/* no values */
	if (values.Length() == 0) return;

	/* check */
	if (node_numbers.Length() != values.MajorDim()) throw eSizeMismatch;

	/* write */
	for (int i = 0; i < node_numbers.Length(); i++)
	{
		out << setw(kIntWidth) << node_numbers[i] + 1;
		values.PrintRow(i, out);	
	}
}

void FE_ASCIIT::WriteElementValues(ostream& out, const dArray2DT& values) const
{
	/* no values */
	if (values.Length() == 0) return;

	/* write */
	for (int i = 0; i < values.MajorDim(); i++)
	{
		out << setw(kIntWidth) << i + 1;
		values.PrintRow(i, out);	
	}
}

