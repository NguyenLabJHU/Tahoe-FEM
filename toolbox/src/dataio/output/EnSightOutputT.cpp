/* $Id: EnSightOutputT.cpp,v 1.13 2004-07-01 16:39:43 paklein Exp $ */
/* created: sawimme (05/18/1999) */

#include "EnSightOutputT.h"

#include "EnSightT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "OutputSetT.h"
#include "AutoArrayT.h"
#include <fstream.h>


using namespace Tahoe;

EnSightOutputT::EnSightOutputT (ostream& out, const ArrayT<StringT>& out_strings, int numdigs, bool binary) :
  OutputBaseT (out, out_strings),
  fBinary(binary),
  fNumDigits (numdigs),
  fTimeValues (20)
{
}

void EnSightOutputT::WriteGeometry (void)
{
  int dof = fCoordinates->MinorDim();
  EnSightT ens (fout, fBinary, dof);
  
  // write geometry file
  ofstream geo;
  StringT geocase = OpenGeometryFile (ens, geo, -1);
  if (geo)
    {
      CreateElementBlockIDs ();
      int partID = 1;
      for (int i=0; i < fElementSets.Length(); i++)
	{
	  WritePart (geo, ens, i);
	  if (partID < fElementBlockIDs[i][0])
	    partID = fElementBlockIDs[i][0] + 1;
	}
      // node sets
      for (int n=0; n < fNodeSets.Length(); n++)
	{
	  StringT description (fOutroot);
	  description.Append (" NodeSet ", fNodeSetNames[n]);
	  ens.WritePartInfo (geo, partID++, description);
	  WriteCoordinates (geo, ens, *fNodeSets[n]);
	}
    }
  
  // write case file
  StringT label = "case";
  StringT casefile = CreateFileName (label, kNoIncFile, -1);
  ofstream out (casefile);
  ens.WriteCaseFormat (out);
  ens.WriteCaseGeometry (out, fSequence + 1, geocase);
}

void EnSightOutputT::WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
  //cout << "***** writing " << ID << endl;
  OutputBaseT::WriteOutput (time, ID, n_values, e_values);
  
  // save time value
  if (fElementSets[ID]->PrintStep() == fTimeValues.Length())
    fTimeValues.Append (time);
  
  int dof = fCoordinates->MinorDim();
  EnSightT ens (fout, fBinary, dof);
  
  // write geometry file
  ofstream geo;
  StringT geocase = OpenGeometryFile (ens, geo, ID);
  WritePart (geo, ens, ID);
  
  // variable files
  AutoArrayT<StringT> names (20), files (20);
  AutoArrayT<EnSightT::VariableTypeT> vtypes;
  const ArrayT<StringT>& n_labels = fElementSets[ID]->NodeOutputLabels();
  /* Original:  const ArrayT<StringT>& e_labels = fElementSets[ID]->ElementOutputLabels();
     changed by TDN: 5/30/2002*/
  const ArrayT<StringT>& e_labels = fElementSets[ID]->ElementOutputLabels();
  if (n_values.MajorDim() > 0)
    WriteVariable (ens, true, ID, n_values, n_labels, names, files, vtypes);
  if (e_values.MajorDim() > 0)
    WriteVariable (ens, false, ID, e_values, e_labels, names, files, vtypes);
  
  // write case file
  StringT label = "case";
  StringT casefile = CreateFileName (label, kNoIncFile, fElementBlockIDs[ID][0]);
  ofstream out (casefile);
  ens.WriteCaseFormat (out);
  ens.WriteCaseGeometry (out, fSequence + 1, geocase);
  if (names.Length() > 0)
    ens.WriteVariableLabels (out, names, files, vtypes);
  if (fTimeValues.Length() > 0)
    {
      int startinc = 0;
      int increment = 1;
      ens.WriteTime (out, fSequence + 1, startinc, increment, fTimeValues);
    }
}

// *************** PRIVATE ********************

StringT EnSightOutputT::OpenGeometryFile (EnSightT& ens, ofstream& geo, int ID) const
{
  StringT label = "geo";
  StringT geofile, geocase;
  bool change = false;
  for (int j=0; j < fElementSets.Length() && !change; j++)
    if (fElementSets[j]->Changing()) 
      change = true;

  int group = (ID > 0) ? fElementBlockIDs[ID][0] : -1;

  if (change)
    {
      geofile = CreateFileName (label, fElementSets[0]->PrintStep(), group);
      geocase = CreateFileName (label, kWildFile, group);
    }
  else
    {
      geofile = CreateFileName (label, kNoIncFile, group);
      geocase = geofile;
    }
  
  // is it necessary to write a geometry file
  if (!change && fElementSets[0]->PrintStep() > 0) return geocase;
  geo.open (geofile);
  
  // header
  int h = 0;
  ArrayT<StringT> header;
  if (fBinary)
    {
      header.Allocate (5);
      header[h++] = "C Binary";
    }
  else
    header.Allocate (4);
  header[h] = fOutroot;
  header[h++].Append (" ", fTitle.Pointer());
  header[h] = fCodeName;
  header[h++].Append (" ", fVersion.Pointer());
  header[h++] = "node id assign";
  header[h++] = "element id assign";
  ens.WriteHeader (geo, header);
  
return geocase;
}

StringT EnSightOutputT::CreateFileName (const StringT& Label, int increment, int groupnumber) const
{
  StringT var (fOutroot);
  
  /* tack on sequence number */
  if (fSequence > 0) var.Append(".sq", fSequence + 1);
  
  /* tack on group number */
  if (groupnumber > 0)
    var.Append (".gp", groupnumber);
  
  /* tack on print increment number or wildcards */
  if (increment == kWildFile)
    {
      var.Append (".ps");
      for (int i = 0; i < fNumDigits; i++)
	var.Append ("*");
    }
  else if (increment > kNoIncFile)
    var.Append (".ps", increment, fNumDigits);
  
  /* tack on extension */
  var.Append (".");
  var.Append (Label);
  return var;
}

void EnSightOutputT::WritePart (ostream& geo, EnSightT& ens, int index) const
{
  const ArrayT<StringT>& blockIDs = fElementSets[index]->BlockID();
  for (int b=0; b < fElementSets[index]->NumBlocks(); b++)
    {
      StringT description = fOutroot;
      description.Append (" Grp ", fElementSets[index]->ID());
      description.Append (".", blockIDs[b]);
      ens.WritePartInfo (geo, fElementBlockIDs[index][b], description);
  
      const iArrayT& nodes_used = fElementSets[index]->BlockNodesUsed(blockIDs[b]);
      WriteCoordinates (geo, ens, nodes_used);
      WriteConnectivity (geo, ens, nodes_used, index, b);
    }
}

void EnSightOutputT::WriteCoordinates (ostream& geo, EnSightT& ens, const iArrayT& nodes) const
{
  dArray2DT local (nodes.Length(), fCoordinates->MinorDim());
  for (int i=0; i < nodes.Length(); i++)
    local.SetRow (i, (*fCoordinates)(nodes[i]));
  ens.WriteCoordinateHeader (geo, nodes.Length());
  ens.WriteCoordinates (geo, local);
}

void EnSightOutputT::WriteConnectivity (ostream& geo, EnSightT& ens, const iArrayT& nodes_used, int i, int block) const
{
  const iArray2DT* connects = fElementSets[i]->Connectivities(fElementSets[i]->BlockID(block));
  int outputnodes = ens.WriteConnectivityHeader (geo, fElementSets[i]->Geometry(), 
						 connects->MajorDim(), connects->MinorDim());
  
  iArray2DT localconn (connects->MajorDim(), connects->MinorDim());
  LocalConnectivity (nodes_used, *connects, localconn);
  localconn++;
  ens.WriteConnectivity (geo, fElementSets[i]->Geometry(), outputnodes, localconn);
}

void EnSightOutputT::WriteVariable (EnSightT& ens, bool nodal, int ID,
				    const dArray2DT& values, const ArrayT<StringT>& labels,
				    AutoArrayT<StringT>& names, AutoArrayT<StringT>& files,
				    AutoArrayT<EnSightT::VariableTypeT>& vtypes) const
{
  // extract block values from output set
  const ArrayT<StringT>& block_ID = fElementSets[ID]->BlockID();
  ArrayT<dArray2DT> blockvalues(fElementSets[ID]->NumBlocks());
  for (int block=0; block < fElementSets[ID]->NumBlocks(); block++)
    {
      if (nodal)
	{
	  iArrayT nodes_used;
	  NodalBlockValues (ID, block, values, blockvalues[block], nodes_used);
	}
      else
	{
	  blockvalues[block].Allocate (fElementSets[ID]->NumBlockElements(block_ID[block]), values.MinorDim());
	  ElementBlockValues (ID, block, values, blockvalues[block]);
	}
    }

  // print each variable to a separate file
  int dof = fCoordinates->MinorDim();
  for (int n=0; n < values.MinorDim(); n++)
    {

      // determine variable name
      StringT extension;
      bool vector = IsVector (labels, n, extension, dof);
      names.Append (extension);
      
      // create variable file name
      StringT varfile = CreateFileName (extension, fElementSets[ID]->PrintStep(), fElementBlockIDs[ID][0]);
      StringT varcase = CreateFileName (extension, kWildFile, fElementBlockIDs[ID][0]);
      
      // account for only one time step being written. 
      if (fTimeValues.Length() < 2)
	files.Append (varfile);
      else
	files.Append (varcase);
      
      // open file
      ofstream var (varfile);
      ArrayT<StringT> header (1);
      header[0] = extension;
      header[0].Append (" ", fElementSets[ID]->PrintStep());
      ens.WriteHeader (var, header);
      
      // write this variable's values for each block to the same file
      for (int block=0; block < fElementSets[ID]->NumBlocks(); block++)
	{
	  // write part information
	  StringT name = "coordinates";
	  if (!nodal)
	    {
	      const iArray2DT* connects = fElementSets[ID]->Connectivities(block_ID[block]);
	      int numelemnodes = connects->MinorDim();
	      ens.GetElementName (name, numelemnodes, fElementSets[ID]->Geometry());
	    }
	  ens.WritePartInfo (var, fElementBlockIDs[ID][block], name);

	  // write values
	  if (vector)
	    ens.WriteVector (var, blockvalues[block], n);
	  else
	    ens.WriteScalar (var, blockvalues[block], n);
	}

      /* save variable types for case file and increment n if vector */
      if (vector)
	{
	  n += dof - 1;
	  if (nodal)
	    vtypes.Append (EnSightT::kVectorNodal);
	  else
	    vtypes.Append (EnSightT::kVectorElemental);
	}
      else
	{
	  if (nodal)
	    vtypes.Append (EnSightT::kScalarNodal);
	  else
	    vtypes.Append (EnSightT::kScalarElemental);
	}
    }
}

bool EnSightOutputT::IsVector (const ArrayT<StringT>& inlabels, int index, StringT& extension, int dof) const
{
	extension = inlabels[index];

	// weed out scalar variables
	if ((strstr ((const char*) inlabels[index], "_x")) == NULL &&
	    (strstr ((const char*) inlabels[index], "_X")) == NULL)
	  return false;

	// ensure that 2D and 3D variables have a Y component
	if (dof >= 2)
	  {
	    if (inlabels.Length() < index+1)
	      return false;
	    if ( (strstr ((const char*) inlabels[index+1], "_y")) == NULL &&
		 (strstr ((const char*) inlabels[index+1], "_Y")) == NULL)
	      return false;
	  }

	// ensure that 3D variables have a Z component
	if (dof == 3)
	  {
	    if (inlabels.Length() < index+2)
	      return false;
	    if ( (strstr ((const char*) inlabels[index+2], "_z")) == NULL &&
		 (strstr ((const char*) inlabels[index+2], "_Z")) == NULL)
	      return false;
	  }

	// we have a vector, now create the proper extension
	extension.Drop(-2);
	return true;
}

