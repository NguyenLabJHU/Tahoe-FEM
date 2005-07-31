/* $Id: AbaqusInputT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1998)                                          */

#include "AbaqusInputT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iAutoArrayT.h"

AbaqusInputT::AbaqusInputT (ostream& out, bool binary, const char* filename) :
InputBaseT (out),
	fData (out, binary),
	fResultFile (filename),
	fNumElements (0),
	fNumNodes (0),
	fDimensions (0),
	fNumElementSets (0),
	fNumNodeSets (0),
	fEnergyData (false)
{
Initialize ();
if (fNumElements > 0) ReadElementSets ();
if (fNumNodes > 0) ReadNodeSets ();
}

int AbaqusInputT::NumElementGroups (void) const
{
iArrayT temp;
GroupNumbers (temp);
return temp.Length();
}

int AbaqusInputT::NumNodeSets (void) const
{
iArrayT ids;
NodeSetNumbers (ids);
return ids.Length();
}

void AbaqusInputT::GroupNumbers (iArrayT& groupnums) const
{
iAutoArrayT temp;
// must account for some element sets having zero elements ??
// only return those element sets that appear in fElementData
for (int i=0; i < fNumElements; i++)
temp.AppendUnique (fElementData (i, 2));

groupnums.Allocate (temp.Length());
groupnums.CopyPart (0, temp, 0, temp.Length());
}

void AbaqusInputT::NodeSetNumbers (iArrayT& nodenums) const
{
iAutoArrayT temp;
// must account for some element sets having zero elements ??
// only return those element sets that appear in fElementData
for (int i=0; i < fNumNodes; i++)
if (fCoordinateData (i, 1) > -1)
temp.AppendUnique (fCoordinateData (i, 1));

nodenums.Allocate (temp.Length());
nodenums.CopyPart (0, temp, 0, temp.Length());
}

void AbaqusInputT::ReadNodeSet (int set_num, iArrayT& nodes) const
{
int count = 0;
for (int i=0; i < fNumNodes; i++)
if (fCoordinateData (i,1) == set_num) count++;

nodes.Allocate (count);
int *pn = nodes.Pointer();
for (int j=0; j < fNumNodes; j++)
if (fCoordinateData (j,1) == set_num)
*pn++ = j; //change to global consecutive numbering, offset to zero
}

void AbaqusInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
coords.Allocate (fNumNodes, fDimensions);
nodemap.Allocate (fNumNodes);
dArrayT x;
int ID;
ifstream in;
if (!OpenFile (in)) throw eGeneralFail;
for (int i=0; i < fNumNodes; i++)
{
fData.ReadCoordinate (in, ID, x);
coords.SetRow (i, x);
nodemap [i] = ID;
}
}

void AbaqusInputT::ReadConnectivity (int group, GeometryT::CodeT& geocode, iArray2DT& connects, iArrayT& elementmap)
{
int count = 0;
int numelemnodes;
for (int i=0, j=0; i < fNumElements; i++)
if (fElementData (i, 2) == group)
{
	count++;
	if (j++ == 0) numelemnodes = fElementData (i, 3);
}

connects.Allocate (count, numelemnodes);
elementmap.Allocate (count);
int ID;
GeometryT::CodeT code;
iArrayT nodes, internal;
ifstream in;
if (!OpenFile (in)) throw eGeneralFail;
for (int k=0, m=0; k < fNumElements; k++)
if (fElementData (k,2) == group)
{
	fData.ReadElement (in, ID, code, nodes);

	// change to consecutive numbering, offset from zero
	internal.Allocate (nodes.Length());
	int dex;
	for (int n=0; n < nodes.Length(); n++)
	  {
	    if (!fCoordinateData.ColumnHasValue (0, nodes[n], dex))
	      throw eGeneralFail;
	    internal[n] = dex;
	  }
	if (m == 0) geocode = code;
	else if (geocode != code)
	  {
	    cout << "\n\nAbaqusInputT::ReadConnectivity, mismatched geo. code\n";
	    throw eGeneralFail;
	  }

	connects.SetRow (m, internal);
	elementmap[m] = ID;
	m++;
}
}

void AbaqusInputT::Close (void)
{
}

void AbaqusInputT::QARecords (ArrayT<StringT>& records) const
{
records.Allocate (4);
records[0] = "Abaqus";
records[1] = fVersion;
records[2] = fDate;
records[3] = fTime;
}

void AbaqusInputT::ReadTimeSteps (dArrayT& steps)
{
ifstream in;
if (!OpenFile (in)) throw eGeneralFail;

double time;
AutoArrayT<double> temp;
while (fData.NextTimeIncrement (in, time))
temp.Append (time);

steps.Allocate (temp.Length());
steps.CopyPart (0, temp, 0, temp.Length());

// check time steps to see if duplicate time value is used
for (int i=1; i < steps.Length(); i++)
if (steps[i] == steps[i-1])
{
	double shift;
	if (i+1 != steps.Length())
	  shift = (steps[i+1]-steps[i])/1.e-5;
	else
	  shift = steps[i]/1.e-5;
	steps[i] += shift;
}
}

void AbaqusInputT::ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id)
{
#pragma unused(group_id)
ifstream in;
if (!OpenFile (in)) throw eGeneralFail;
double time;
if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
ArrayT<AbaqusT::VariableKeyT> nkeys, ekeys;
fData.ReadLabels (in, nkeys, ekeys, fEnergyData);
nlabels.Allocate (nkeys.Length());
elabels.Allocate (ekeys.Length());

StringT newlabel;
char incrementor = '1';
int unknown = 0;
for (int n=0; n < nkeys.Length(); n++)
{
if (n == 0 || nkeys[n] != nkeys[n-1]) incrementor = '1';
SetLabelName (nkeys[n], unknown, nlabels[n], incrementor);
incrementor++;
}
for (int e=0; e < ekeys.Length(); e++)
{
if (e == 0 || ekeys[e] != ekeys[e-1]) incrementor = '1';
SetLabelName (ekeys[e], unknown, elabels[e], incrementor);
incrementor++;
}
}

void AbaqusInputT::ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues)
{
// determine number of variables
ifstream in;
if (!OpenFile (in)) throw eGeneralFail;
ArrayT<AbaqusT::VariableKeyT> nkeys, ekeys;
double time;
for (int i=0; i < step + 1; i++)
if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
fData.ReadLabels (in, nkeys, ekeys, fEnergyData);

// read nodal variables
if (nkeys.Length() > 0)
{
if (!OpenFile (in)) throw eGeneralFail; // rewind file
for (int i1=0; i1 < step + 1; i1++)
	if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;
nvalues.Allocate (fNumNodes, nkeys.Length());
fData.ReadNodeVariables (in, nkeys, nvalues);
}

// read elemental variables
if (ekeys.Length() > 0)
{
if (!OpenFile (in)) throw eGeneralFail; // rewind file
for (int i2=0; i2 < step + 1; i2++)
	if (!fData.NextTimeIncrement (in, time)) throw eGeneralFail;

iAutoArrayT elems;
for (int i=0, j=0; i < fNumElements; i++)
	if (fElementData (i, 2) == group_id)
	  elems.Append (fElementData (i, 0));

evalues.Allocate (elems.Length(), ekeys.Length());
iArrayT els (elems.Length());
els.Set (elems.Length(), elems.Pointer());
fData.ReadElementVariables (in, ekeys, evalues, els);
}
}

/*************************************************************************
*
* Private
*
*************************************************************************/

bool AbaqusInputT::OpenFile (ifstream& in)
{
if (in) in.close();
fData.ResetBufferSize ();
in.open (fResultFile);
if (in) return true;
else
{
cout << "\n\nAbaqusInputT::OpenFile unable to open file: "
	   << fResultFile << endl;
return false;
}
}

void AbaqusInputT::Initialize (void)
{
ifstream in;
if (!OpenFile (in))
{
cout << "\n\nAbaqusInputT, unable to open file: " << fResultFile << endl;
throw eGeneralFail;
}

// read num elements, num dimensions, etc.
fData.ReadVersion (in, fVersion, fDate, fTime, fNumElements, fNumNodes);

// set up element dat
if (!OpenFile (in)) throw eGeneralFail; // rewind file
int ID;
GeometryT::CodeT geocode;
iArrayT nodes;
fElementData.Allocate (fNumElements, 4);
fElementData = -1;
int *p = fElementData.Pointer();
for (int e=0; e < fNumElements; e++)
{
fData.ReadElement (in, ID, geocode, nodes);
*p++ = ID;
*p++ = geocode;
*p++;
*p++ = nodes.Length();
}

// set up coord data
if (!OpenFile (in)) throw eGeneralFail; // rewind file
fCoordinateData.Allocate (fNumNodes, 2);
fCoordinateData = -1;
dArrayT x;
p = fCoordinateData.Pointer();
for (int n=0; n < fNumNodes; n++)
{
fData.ReadCoordinate (in, ID, x);
*p++ = ID;
*p++;
}
fDimensions = x.Length();

// determine if modal analysis
if (!OpenFile (in)) throw eGeneralFail; // rewind file
fData.IsModal (in);
}

void AbaqusInputT::ReadElementSets (void)
{
ifstream (in);
if (!OpenFile (in)) throw eGeneralFail;
fData.ReadElementSets (in, fElementData);
}

void AbaqusInputT::ReadNodeSets (void)
{
ifstream (in);
if (!OpenFile (in)) throw eGeneralFail;
fData.ReadNodeSets (in, fCoordinateData);
}

void AbaqusInputT::SetLabelName (AbaqusT::VariableKeyT key, int& unknown, StringT& name, char incrementor) const
{
switch (key)
{
case AbaqusT::aSDV:           name = "SDV_"; break;
case AbaqusT::aStress:        name = "S_";   break;
case AbaqusT::aTotalStrain:   name = "E_";   break;
case AbaqusT::aDisplacement:  name = "D_";   break;
case AbaqusT::aVelocity:      name = "V_";   break;
case AbaqusT::aAcceleration:  name = "A_";   break;
case AbaqusT::aCoordVariable: name = "X_";   break;
case AbaqusT::aPrinStress:    name = "PS_";  break;
default:
name = "unk";
name.Append (unknown++);
name.Append ("_");
}
name.Append (incrementor);
}

