// file: FEASCIIT.cpp

// created      : SAW(05/20/1999)
// last modified: PAK(11/08/1999)

#include "FE_ASCIIT.h"

#include "GeometryT.h"
#include "OutputSetT.h"
#include "ModelFileT.h"

#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

/* constructor */
FE_ASCIIT::FE_ASCIIT(ostream& out, bool external, const ArrayT<StringT>& out_strings):
	OutputBaseT(out, out_strings),
	fExternTahoeII (external)
{
}

/* destructor */
FE_ASCIIT::~FE_ASCIIT(void)
{
	/* close streams */
	if (IsOpen(fmovie)) fmovie.close();
	if (IsOpen(fgeo)) fgeo.close();
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
      const iArray2DT& conn = fElementSets[e]->Connectivities();
      iArrayT tmp(conn.Length(), conn.Pointer());
      tmp++;
      mf.PutElementSet (fElementSets[e]->ID(), conn);
      tmp--;
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

	/* open streams */
	if (!IsOpen(fmovie) || !IsOpen(fgeo)) OpenStreams();

	/* write geometry */
	fgeo << " Time. . . . . . . . . . . . . . . . . . . . . . = ";
	fgeo << time << '\n';
	WriteGeometryData(ID);
	fgeo << '\n';

	/* Print header */
	fmovie << " Time. . . . . . . . . . . . . . . . . . . . . . = ";
	fmovie << time << '\n';
	WriteOutputData(ID, n_values, e_values);
	fmovie << '\n';
}

/*************************************************************************
 * Protected
 *************************************************************************/

/*************************************************************************
 * Private
 *************************************************************************/

void FE_ASCIIT::OpenStreams(void)
{
	if (IsOpen(fmovie)) fmovie.close();
	if (IsOpen(fgeo)) fgeo.close();

	StringT var(fOutroot);
	StringT geom(fOutroot);

	/* tack on sequence number */
	if (fSequence > 0)
    {
		var.Append("_seq", fSequence + 1);
		geom.Append("_seq", fSequence + 1);
	}

	/* tack on group number extension */
	var.Append(".run");
	geom.Append(".geo");

	fmovie.open(var);
	SetStreamPrefs(fmovie);
	fmovie << "\n O U T P U T   D A T A :\n\n";

	fgeo.open(geom);
	SetStreamPrefs(fgeo);
	fgeo << "\n G E O M E T R Y   D A T A:\n\n";
}

void FE_ASCIIT::WriteGeometryData(int ID)
{
	/* do not need to rewrite geometry */
	if (fElementSets[ID]->PrintStep() > 0 && 
	   !fElementSets[ID]->Changing()) return;

	/* dimensions */
	int nsd = fCoordinates->MinorDim();

	/* generate coordinate labels */
	ArrayT<StringT> coord_labels(nsd);
	for (int i = 0; i < nsd; i++)
	{
		coord_labels[i] = "x";
		coord_labels[i].Append(i+1);
	}

	/* ID information */
	fgeo << " Group number. . . . . . . . . . . . . . . . . . = "
	     << fElementSets[ID]->ID() << '\n';
	fgeo << " Output ID . . . . . . . . . . . . . . . . . . . = "
	     << ID << '\n';
	fgeo << " Print Step. . . . . . . . . . . . . . . . . . . = ";
	fgeo << fElementSets[ID]->PrintStep() << '\n'; 

	/* write coordinates */
	const iArrayT& nodes_used = fElementSets[ID]->NodesUsed();
	dArray2DT local_coordinates(nodes_used.Length(), nsd);
	local_coordinates.RowCollect(nodes_used, *fCoordinates);
	
	fgeo << "\n Nodal coordinates:\n";	
	WriteNodeHeader(fgeo, local_coordinates.MajorDim(), coord_labels);
	WriteNodeValues(fgeo, nodes_used, local_coordinates);
	
	/* write connectivities */
	const iArray2DT& connectivities = fElementSets[ID]->Connectivities();
	fgeo << "\n Connectivities:\n";
	fgeo << " Number of elements .  . . . . . . . . . . . . . = " 
         << connectivities.MajorDim() << '\n';
	fgeo << " Number of element nodes . . . . . . . . . . . . = " 
         << connectivities.MinorDim() << "\n\n";
	fgeo << setw(kIntWidth) << "element";
	for (int j = 0; j < connectivities.MinorDim(); j++)
		fgeo << setw(kIntWidth - 1) << "n" << j+1;
	fgeo << '\n';
			
	/* correct offset for output */
	iArray2DT connects_temp;
	connects_temp.ShallowCopy(connectivities);	
	connects_temp++;
	connects_temp.WriteNumbered(fgeo);
	connects_temp--;
	fgeo << endl;
}

void FE_ASCIIT::WriteOutputData(int ID, const dArray2DT& n_values, 
	const dArray2DT& e_values)
{
	fmovie << " Group number. . . . . . . . . . . . . . . . . . = "
	       << fElementSets[ID]->ID() << '\n';
	fmovie << " Output ID . . . . . . . . . . . . . . . . . . . = "
	       << ID << '\n';
	fmovie << " Print Step. . . . . . . . . . . . . . . . . . . = ";
	fmovie << fElementSets[ID]->PrintStep() << '\n'; 

	fmovie << "\n Nodal data:\n";	
	const ArrayT<StringT>& node_labels = fElementSets[ID]->NodeOutputLabels();
	WriteNodeHeader(fmovie, n_values.MajorDim(), node_labels);
	WriteNodeValues(fmovie, fElementSets[ID]->NodesUsed(), n_values);
	     
	fmovie << "\n Element data:\n";
	const ArrayT<StringT>& elem_labels = fElementSets[ID]->ElementOutputLabels();
	WriteElementHeader(fmovie, e_values.MajorDim(), elem_labels);
	WriteElementValues(fmovie, e_values);
}

void FE_ASCIIT::WriteNodeHeader(ostream& out, int num_output_nodes,
	const ArrayT<StringT>& labels) const
{
	int d_width = out.precision() + kDoubleExtra;

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
	int d_width = out.precision() + kDoubleExtra;

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

/********** delete ***********/
#if 0
void FE_ASCIIT::PrintSingleTextFile(StringT& file, const iArray2DT& a) const
{
  ofstream mout(file);
  mout << a.MajorDim() << "  " << a.MinorDim() << endl;
  a.WriteNumbered(mout);
}

void FE_ASCIIT::PrintSingleTextFile(StringT& file, const iArrayT& a) const
{
  ofstream mout(file);
  mout << a.Length() << endl << a;
}
#endif
/********** delete ***********/

