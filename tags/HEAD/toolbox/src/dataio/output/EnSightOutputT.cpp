/* $Id: EnSightOutputT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1999)                                          */

#include "EnSightOutputT.h"

#include "EnSightT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "OutputSetT.h"
#include "AutoArrayT.h"
#include <fstream.h>

EnSightOutputT::EnSightOutputT (ostream& out, const ArrayT<StringT>& out_strings, int numdigs, bool binary) :
OutputBaseT (out, out_strings),
fNumDigits (numdigs),
fBinary (binary),
fTimeValues (20, true)
{
}

void EnSightOutputT::WriteGeometry (void)
{
int dof = fCoordinates->MinorDim();
EnSightT ens (fout, fBinary, dof);

// write geometry file
ofstream geo;
StringT geocase = OpenGeometryFile (ens, geo);
if (geo)
{
int partID = 1;
for (int i=0; i < fElementSets.Length(); i++)
	{
	  WritePart (geo, ens, i);
	  if (partID < fElementSets[i]->ID())
	    partID = fElementSets[i]->ID() + 1;
	}
// node sets
for (int n=0; n < fNodeSets.Length(); n++)
	{
	  StringT description (fOutroot);
	  description.Append (" NodeSet ", fNodeSetIDs[n]);
	  ens.WritePartInfo (geo, partID++, description);
	  WriteCoordinates (geo, ens, *fNodeSets[n]);
	}
}

// write case file
StringT label = "case";
StringT casefile = CreateFileName (label, kNoIncFile);
ofstream out (casefile);
ens.WriteCaseFormat (out);
ens.WriteCaseGeometry (out, fSequence + 1, geocase);
}

void EnSightOutputT::WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
OutputBaseT::WriteOutput (time, ID, n_values, e_values);

// save time value
fTimeValues.Append (time);

int dof = fCoordinates->MinorDim();
EnSightT ens (fout, fBinary, dof);

// write geometry file
ofstream geo;
StringT geocase = OpenGeometryFile (ens, geo);
WritePart (geo, ens, ID);

// variable files
AutoArrayT<StringT> names (20, false), files (20, false);
AutoArrayT<EnSightT::VariableTypeT> vtypes;
const ArrayT<StringT>& n_labels = fElementSets[ID]->NodeOutputLabels();
const ArrayT<StringT>& e_labels = fElementSets[ID]->NodeOutputLabels();
WriteVariable (ens, true, ID, n_values, n_labels, names, files, vtypes);
WriteVariable (ens, false, ID, e_values, e_labels, names, files, vtypes);

// write case file
StringT label = "case";
StringT casefile = CreateFileName (label, kNoIncFile);
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

StringT EnSightOutputT::OpenGeometryFile (EnSightT& ens, ofstream& geo) const
{
StringT label = "geo";
StringT geofile, geocase;
bool change = false;
for (int j=0; j < fElementSets.Length() && !change; j++)
if (fElementSets[j]->Changing()) change = true;
if (change)
{
geofile = CreateFileName (label, fElementSets[0]->PrintStep());
geocase = CreateFileName (label, kWildFile);
}
else
{
geofile = CreateFileName (label, kNoIncFile);
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

StringT EnSightOutputT::CreateFileName (const StringT& Label, int increment) const
{
StringT var (fOutroot);

/* tack on sequence number */
if (fSequence > 0) var.Append(".sq", fSequence + 1);

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
StringT description = fOutroot;
description.Append (" Grp ", fElementSets[index]->ID());
ens.WritePartInfo (geo, fElementSets[index]->ID(), description);

iArrayT nodes_used;
nodes_used.Alias (fElementSets[index]->NodesUsed());
WriteCoordinates (geo, ens, nodes_used);
WriteConnectivity (geo, ens, nodes_used, index);
}

void EnSightOutputT::WriteCoordinates (ostream& geo, EnSightT& ens, const iArrayT& nodes) const
{
dArray2DT local (nodes.Length(), fCoordinates->MinorDim());
for (int i=0; i < nodes.Length(); i++)
local.SetRow (i, (*fCoordinates)(nodes[i]));
ens.WriteCoordinateHeader (geo, nodes.Length());
ens.WriteCoordinates (geo, local);
}

void EnSightOutputT::WriteConnectivity (ostream& geo, EnSightT& ens, const iArrayT& nodes_used, int i) const
{
const iArray2DT& connects = fElementSets[i]->Connectivities();
int numelems = connects.MajorDim();
int numelemnodes = connects.MinorDim();
int outputnodes = ens.WriteConnectivityHeader (geo, fElementSets[i]->Geometry(), numelems, numelemnodes);

iArray2DT localconn (numelems, numelemnodes);
LocalConnectivity (nodes_used, connects, localconn);
localconn++;
ens.WriteConnectivity (geo, outputnodes, localconn);
}

void EnSightOutputT::WriteVariable (EnSightT& ens, bool nodal, int ID,
	const dArray2DT& values, const ArrayT<StringT>& labels,
	AutoArrayT<StringT>& names, AutoArrayT<StringT>& files,
	AutoArrayT<EnSightT::VariableTypeT>& vtypes) const
{
int dof = fCoordinates->MinorDim();
for (int n=0; n < values.MinorDim(); n++)
{
// determine variable name
StringT extension;
bool vector = IsVector (labels, n, extension, dof);
names.Append (extension);

// open variable file
StringT varfile = CreateFileName (extension, fElementSets[ID]->PrintStep());
StringT varcase = CreateFileName (extension, kWildFile);
files.Append (varcase);

// open file
ofstream var (varfile);
ArrayT<StringT> header (1);
header[0] = extension;
header[0].Append (" ", fElementSets[ID]->PrintStep());
ens.WriteHeader (var, header);

// write part information
StringT name = "coordinates";
if (!nodal)
	{
	  const iArray2DT& connects = fElementSets[ID]->Connectivities();
	  int numelemnodes = connects.MinorDim();
	  ens.GetElementName (name, numelemnodes, fElementSets[ID]->Geometry());
	}
ens.WritePartInfo (var, fElementSets[ID]->ID(), name);

// write values
if (vector)
	{
	  ens.WriteVector (var, values, n);
	  n += dof - 1;
	  if (nodal)
	    vtypes.Append (EnSightT::kVectorNodal);
	  else
	    vtypes.Append (EnSightT::kVectorElemental);
	}
else
	{
	  ens.WriteScalar (var, values, n);
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

if ( (strstr ((const char*) inlabels[index], "_x")) == NULL &&
	 (strstr ((const char*) inlabels[index], "_X")) == NULL)
return false;

if (dof >= 2)
{
	if (inlabels.Length() < index+1)
	  return false;
	if ( (strstr ((const char*) inlabels[index+1], "_y")) == NULL &&
	     (strstr ((const char*) inlabels[index+1], "_Y")) == NULL)
	  return false;
}
if (dof == 3)
{
	if (inlabels.Length() < index+2)
	  return false;
	if ( (strstr ((const char*) inlabels[index+2], "_z")) == NULL &&
	     (strstr ((const char*) inlabels[index+2], "_Z")) == NULL)
	  return false;
}

// create vector extension
if ( (strstr ((const char*) inlabels[index], "_x")) != NULL)
extension.DefaultName (inlabels[index], "_x", "", -1);
else
extension.DefaultName (inlabels[index], "_X", "", -1);

return true;
}

