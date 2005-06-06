/* $Id: TextOutputT.cpp,v 1.5 2005-06-06 06:35:33 paklein Exp $ */
/* created: sawimme (05/20/1999) */
#include "TextOutputT.h"

#include "GeometryT.h"
#include "OutputSetT.h"
#include "ModelFileT.h"
#include "ofstreamT.h"

#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* constructor */
TextOutputT::TextOutputT(ostream& out, bool external, const ArrayT<StringT>& out_strings):
	OutputBaseT(out, out_strings),
	fExternTahoeII(external)
{

}

/* register the output for an element set. returns the output ID */
int TextOutputT::AddElementSet(const OutputSetT& output_set)
{
	/* set flags */
	fInitGeom.Append(false);
	fInitRun.Append(false);

	/* inherited */
	return OutputBaseT::AddElementSet(output_set);
}

/* increment sequence, create new output file series */
void TextOutputT::NextTimeSequence(int sequence_number)
{
	/* inherited */
	OutputBaseT::NextTimeSequence(sequence_number);
	
	/* reset flags */
	fInitGeom = false;
	fInitRun  = false;
}

//NOTE: to write a geometry definition file
void TextOutputT::WriteGeometry(void)
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
	
	CreateElementBlockIDs ();
	for (int e=0; e < fElementSets.Length(); e++)
	  {
	    const ArrayT<StringT>& blockIDs = fElementSets[e]->BlockID();
	    for (int b=0; b < fElementSets[e]->NumBlocks(); b++)
	      {
		const iArray2DT* c = fElementSets[e]->Connectivities(blockIDs[b]);
		iArray2DT conn = *c;
		
		iArrayT tmp(conn.Length(), conn.Pointer());
		tmp++;
		mf.PutElementSet (fElementBlockIDs[e][b], conn);
		tmp--;
	      }
	  }
	
	/* create integer ID values from string values */
	String2IntIDs (fNodeSetNames, fNodeSetIntIDs);
	String2IntIDs (fSideSetNames, fSideSetIntIDs);

	for (int n=0; n < fNodeSets.Length(); n++)
	  {
	    iArrayT& set = *((iArrayT*) fNodeSets[n]);
	    set++;
	    mf.PutNodeSet (fNodeSetIntIDs[n], set);
	    set--;
	  }

	for (int s=0; s < fSideSets.Length(); s++)
	  {
	    /* search for group name */
	    StringT& gname = fSSGroupNames [s];
	    int gindex, bindex;
	    ElementGroupBlockIndex (gname, gindex, bindex);
	    int block_ID = fElementBlockIDs [gindex][bindex];

	    iArray2DT& set = *((iArray2DT*) fSideSets[s]);
	    set++;
	    mf.PutSideSet (fSideSetIntIDs[s], block_ID, set);
	    set--;
	  }

	mf.Close(); // datafile actually written here
}

void TextOutputT::WriteOutput(double time, int ID, const dArray2DT& n_values,
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
		geom_file.Append(".io", ID);
		if (fElementSets[ID]->Changing()) /* changing - assuming no more than 1000 output steps */
			geom_file.Append(".ps", fElementSets[ID]->PrintStep() + 1, 4);		
		geom_file.Append(".geo");

		/* open stream */
		ofstreamT out;
		SetStreamPrefs(out);
		if (!fInitGeom[ID])
		{
			/* initialize geometry file */
			out.open(geom_file);
			out << "\n G E O M E T R Y   D A T A:\n\n";
	
			/* set flag */
			fInitGeom[ID] = true;
		}
		else /* re-open file */
			out.open_append(geom_file);

		/* check */
		if (!out.is_open())
		{
			cout << "\n TextOutputT::WriteOutput: error opening file: " << geom_file << endl;
			throw ExceptionT::kGeneralFail;
		}

		/* ID information */
		out << " Group number. . . . . . . . . . . . . . . . . . = "
		     << fElementSets[ID]->ID() << '\n';
		out << " Output ID . . . . . . . . . . . . . . . . . . . = "
		     << ID << '\n';
		out << " Number of blocks. . . . . . . . . . . . . . . . = "
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

	/* toc file name */
	StringT toc_file(fOutroot);
	if (fSequence > 0) toc_file.Append(".seq", fSequence + 1);
	toc_file.Append(".io", ID);
	if (fElementSets[ID]->Changing()) /* changing - assuming no more than 1000 output steps */
		toc_file.Append(".ps", fElementSets[ID]->PrintStep() + 1, 4);		
	toc_file.Append(".run");

	/* open stream */
	ofstreamT toc;
	SetStreamPrefs(toc);
	if (!fInitRun[ID] || fElementSets[ID]->Changing())
	{
		/* initialize toc file */
		toc.open(toc_file);
		InitResultsFile(toc, ID);
	
		/* set flag */
		fInitRun[ID] = true;
	}
	else /* re-open file */
		toc.open_append(toc_file);

	/* data file name */
#ifndef __USE_CONSTANT_SUFFIX__
	StringT dat_file(toc_file);
	dat_file.Append(".ps", fElementSets[ID]->PrintStep(), 4);
#else // for Windows, it's nice to have a suffix to associate with output files
	StringT dat_file(fOutroot);
	if (fSequence > 0) dat_file.Append(".seq",fSequence + 1);
	dat_file.Append(".io", ID);
	dat_file.Append(".ps", fElementSets[ID]->PrintStep(), 4);
	dat_file.Append(".run");
#endif
	/* write toc entry - drop the file path */
	StringT file_path;
	file_path.FilePath(dat_file);
	StringT dat_file_relative(dat_file);
	dat_file_relative.Drop(file_path.StringLength());
	toc << dat_file_relative << '\n';
	toc.close();
		
	/* open data file */
	ofstreamT out(dat_file);
	if (!out.is_open())
		ExceptionT::GeneralFail("TextOutputT::WriteOutput", "error opening file: %s", dat_file.Pointer());
	SetStreamPrefs(out);

	/* print header */
	out << "\n Group number. . . . . . . . . . . . . . . . . . = "
        << fElementSets[ID]->ID() << '\n';
	out << " Output ID . . . . . . . . . . . . . . . . . . . = "
		<< ID << '\n';	
	out << " Print Step. . . . . . . . . . . . . . . . . . . = "
	    << fElementSets[ID]->PrintStep() << '\n';
	out << " Time. . . . . . . . . . . . . . . . . . . . . . = "
	    << time << '\n';
	out << " Number of blocks. . . . . . . . . . . . . . . . = "
	    << fElementSets[ID]->NumBlocks() << '\n';

	/* write data */
	WriteOutputData(out, ID, n_values, e_values);
}

/*************************************************************************
* Private
*************************************************************************/

/* initialize the results file */
void TextOutputT::InitResultsFile(ostream& out, int ID)
{
	/* output set */
	OutputSetT& set = *fElementSets[ID];

	/* dimensions section */
	out << "\n S U M M A R Y :\n\n";
	out << " Group number. . . . . . . . . . . . . . . . . . = "
        << fElementSets[ID]->ID() << '\n';
	out << " Output ID . . . . . . . . . . . . . . . . . . . = "
		<< ID << '\n';	
	out << " Number of nodal values. . . . . . . . . . . . . = "
	    << set.NumNodeValues() << '\n';
	if (set.NumNodeValues() > 0)
	{
		out << " Labels:\n";
		//int count = 0;
		const ArrayT<StringT>& n_labels = set.NodeOutputLabels();
		for (int i = 0; i < n_labels.Length(); i++)
			out << '\t' << n_labels[i] << '\n';
	}
	out << " Number of element values. . . . . . . . . . . . = "
	    << set.NumElementValues() << '\n';
	if (set.NumElementValues() > 0)
	{
		out << " Labels:\n";
		//int count = 0;
		const ArrayT<StringT>& e_labels = set.ElementOutputLabels();
		for (int i = 0; i < e_labels.Length(); i++)
			out << '\t' << e_labels[i] << '\n';
	}
	out << " Number of blocks. . . . . . . . . . . . . . . . = "
	    << set.NumBlocks() << '\n';
	out << " Changing geometry . . . . . . . . . . . . . . . = "
	    << set.Changing() << '\n';

	/* write block information */
	out << '\n';
	for (int i = 0; i < fElementSets[ID]->NumBlocks(); i++)
	{
		/* block ID */
		const StringT& block_ID = set.BlockID(i);		
		out << " Block ID. . . . . . . . . . . . . . . . . . . . = " << block_ID << '\n';
		out << " Number of nodes . . . . . . . . . . . . . . . . = " << set.BlockNodesUsed(block_ID).Length() << '\n';
		out << " Number of elements. . . . . . . . . . . . . . . = " << set.NumBlockElements(block_ID) << '\n';
	}

	/* result section header */
	out << "\n O U T P U T   D A T A :  T O C\n";
}

void TextOutputT::WriteGeometryData(ostream& out, int ID)
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

	/* collect set coordinates */
	const iArrayT& nodes_used = fElementSets[ID]->NodesUsed();
	dArray2DT set_coords(nodes_used.Length(), fCoordinates->MinorDim());
	set_coords.RowCollect(nodes_used, *fCoordinates);

	/* collect set ids */
	iArrayT node_id;
	if (!fNodeID) {
		node_id = nodes_used;
		node_id++;
	}
	else {
		node_id.Dimension(nodes_used.Length());
		node_id.Collect(nodes_used, *fNodeID);
	}

	/* write coords */
	out << "\n Nodal coordinates:\n";	
	WriteNodeHeader(out, set_coords.MajorDim(), coord_labels);
	WriteNodeValues(out, node_id, set_coords);
	
	const ArrayT<StringT>& blockIDs = fElementSets[ID]->BlockID();
	for (int b=0; b < fElementSets[ID]->NumBlocks(); b++)
	  {
	    /* write connectivities */
	    const iArray2DT* c = fElementSets[ID]->Connectivities(blockIDs[b]);
	    out << "\n Connectivities:\n";
	    out << " Block ID . . . .  . . . . . . . . . . . . . . . = "
		<< blockIDs[b] << '\n';
	    out << " Number of elements .  . . . . . . . . . . . . . = "
		<< c->MajorDim() << '\n';
	    out << " Number of element nodes . . . . . . . . . . . . = "
		<< c->MinorDim() << '\n';
	    out << " Geometry code . . . . . . . . . . . . . . . . . = "
		<< fElementSets[ID]->Geometry() << "\n\n";
	    out << setw(kIntWidth) << "index" << setw(kIntWidth) << "element";
	    for (int j = 0; j < c->MinorDim(); j++)
	      out << setw(kIntWidth - 1) << "n" << j+1;
	    out << '\n';
			
	    /* map to local numbering */
	    iArray2DT connects_temp;
	    connects_temp.Dimension(*c);
	    LocalConnectivity(nodes_used, *c, connects_temp);
	
		/* write */
	    connects_temp++;
		for (int i = 0; i < connects_temp.MajorDim(); i++)
		{
			out << setw(kIntWidth) << i + 1  // index
				<< setw(kIntWidth) << i + 1; // element number
			connects_temp.PrintRow(i, out);	
		}
	    connects_temp--;

	  }
	out.flush();
}

void TextOutputT::WriteOutputData(ostream& out, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values)
{
	/* write node header */
	out << "\n Nodal data:\n";	
	const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();

	/* write node vars */
	if (n_values.MajorDim () > 0)
	{
		/* collect set ids */
		const iArrayT& nodes_used = fElementSets[ID]->NodesUsed();
		iArrayT node_id;
		if (!fNodeID) {
			node_id = nodes_used;
			node_id++;
		}
		else {
			node_id.Dimension(nodes_used.Length());
			node_id.Collect(nodes_used, *fNodeID);
		}

		/* write values */
		WriteNodeHeader(out, node_id.Length(), node_labels);
		WriteNodeValues(out, node_id, n_values);
	}
	else
		WriteNodeHeader(out, 0, node_labels);

	const ArrayT<StringT>& blockIDs = fElementSets[ID]->BlockID();
	for (int b=0; b < fElementSets[ID]->NumBlocks(); b++)
	  {
	    /* write element header */
	    out << "\n Element data:\n";
	    out << " Block ID . . . . . .  . . . . . . . . . . . . . = "
		<< blockIDs[b] << '\n';
	    const ArrayT<StringT>& elem_labels = fElementSets[ID]->ElementOutputLabels();


	    /* write element values */
	    if (e_values.MajorDim() > 0)
		{
			/* collect block values */
			dArray2DT local_vals(fElementSets[ID]->NumBlockElements(blockIDs[b]), e_values.MinorDim());
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

void TextOutputT::WriteNodeHeader(ostream& out, int num_output_nodes,
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
		out << setw(kIntWidth) << "index" << setw(kIntWidth) << "node";
		for (int i = 0; i < labels.Length(); i++)
			out << setw(d_width) << labels[i];
		out << '\n';
	}
}

void TextOutputT::WriteElementHeader(ostream& out, int num_output_elems,
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
		out << setw(kIntWidth) << "index" << setw(kIntWidth) << "element";
		for (int i = 0; i < labels.Length(); i++)
			out << setw(d_width) << labels[i];
		out << '\n';
	}
}

void TextOutputT::WriteNodeValues(ostream& out, const ArrayT<int>& node_numbers,
	const dArray2DT& values)
{
	/* no values */
	if (values.Length() == 0) return;

	/* check */
	if (node_numbers.Length() != values.MajorDim()) 
		ExceptionT::SizeMismatch("TextOutputT::WriteNodeValues");

	/* write */
	for (int i = 0; i < node_numbers.Length(); i++)
	{
		out << setw(kIntWidth) << i+1
		    << setw(kIntWidth) << node_numbers[i];
		values.PrintRow(i, out);	
	}
}

void TextOutputT::WriteElementValues(ostream& out, const dArray2DT& values)
{
	/* no values */
	if (values.Length() == 0) return;

	/* write */
	for (int i = 0; i < values.MajorDim(); i++)
	{
		out << setw(kIntWidth) << i + 1  // index
		    << setw(kIntWidth) << i + 1; // element number
		values.PrintRow(i, out);	
	}
}
