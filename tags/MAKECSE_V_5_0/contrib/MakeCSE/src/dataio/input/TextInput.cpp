// file: TextInput.cpp

// created: SAW (08/11/99)

/* base class */
#include "TextInput.h"
#include "dArrayT.h"

TextInput::TextInput (ostream& out, const StringT& filename) :
  InputBaseT (out)
{
  fData.OpenRead (filename);
}

int  TextInput::NumElementGroups (void) const
{
  iArrayT ids;
  GroupNumbers (ids);
  return ids.Length();
}

int  TextInput::NumSideSets (void) const
{
  iArrayT ids;
  SideSetNumbers (ids);
  return ids.Length();
}

int  TextInput::NumNodeSets (void) const
{
  iArrayT ids;
  NodeSetNumbers (ids);
  return ids.Length();
}

void TextInput::GroupNumbers (iArrayT& groupnums) const
{
  if (!fData.GetElementSetID (groupnums)) throw eBadInputValue;
}

void TextInput::SideSetNumbers (iArrayT& sidenums) const
{
  if (!fData.GetSideSetID (sidenums)) throw eBadInputValue;
}

void TextInput::NodeSetNumbers (iArrayT& nodenums) const
{
  if (!fData.GetNodeSetID (nodenums)) throw eBadInputValue;
}


void TextInput::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
#pragma unused (nodemap)
  fData.GetCoordinates (coords);
}

void TextInput::ReadConnectivity (int group, GeometryT::GeometryCode& geocode, iArray2DT& connects, iArrayT& elementmap)
{
#pragma unused (elementmap)
  geocode = GeometryT::kNone;
  if (!fData.GetElementSet (group, connects)) throw eBadInputValue;

  connects += -1;

  int numnodes, numdims;
  if (!fData.GetDimensions (numnodes, numdims)) throw eBadInputValue;
  if (numdims == 2)
    switch (connects.MinorDim())
      {
      case 6:
      case 3: geocode = GeometryT::kTriangle; break;
      case 8: 
      case 4: geocode = GeometryT::kQuadrilateral; break;
      }
  else if (numdims == 3)
    switch (connects.MinorDim())
      {
      case 4:
      case 10: geocode = GeometryT::kTetrahedron; break;
      case 8:
      case 20: geocode = GeometryT::kHexahedron; break;
      case 6:
      case 15: geocode = GeometryT::kPentahedron; break;
      }
}

void TextInput::ReadNodeSet (int set_num, iArrayT& nodes) const
{
  if (!fData.GetNodeSet (set_num, nodes)) throw eBadInputValue;

  nodes += -1;
}

void TextInput::ReadSideSet (int set_num, iArray2DT& sides) const
{
  int block;
  if (!fData.GetSideSet (set_num, block, sides)) throw eBadInputValue;

  sides += -1;
}
                
void TextInput::ReadSideSetGlobal (int set_num, iArray2DT& sides) const
{
  int block;
  if (!fData.GetSideSet (set_num, block, sides)) throw eBadInputValue;

  iArrayT ids;
  GroupNumbers (ids);
  int num_elems, dim, offset = 0;
  for (int i=0; i < ids.Length(); i++)
    {
      if (ids[i] == block) break;
      fData.GetElementSetDimensions (ids[i], num_elems, dim);
      offset += num_elems;
    }

  int *pelem = sides.Pointer();
  for (int j=0; j < sides.MajorDim(); j++, pelem += 2)
    *pelem += offset;

  sides += -1;
}
                
void TextInput::Close (void)
{
  fData.Close();
}
                
void TextInput::QARecords (ArrayT<StringT>& records) const
{
  records.Allocate (4);
  records[0] = "Tahoe II";
  records[1] = "v1.0";
  records[2] = "";
  records[3] = "";
}

void TextInput::ReadTimeSteps (dArrayT& steps)
{
  steps.Allocate (0);
}

void TextInput::ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id)
{
#pragma unused (group_id)
  nlabels.Allocate (0);
  elabels.Allocate (0);
}

void TextInput::ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues)
{
#pragma unused (group_id)
#pragma unused (step)
  nvalues.Allocate (0,0);
  evalues.Allocate (0,0);
}

