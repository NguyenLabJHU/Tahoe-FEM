
#include "ExtractNode.h"

ExtractNode::ExtractNode (ostream& out, istream& in, bool write) :
  ExtractIOManager (out, in, write)
{
}

/**************** PROTECTED **********************/

void ExtractNode::Initialize (void)
{
  InitializeNodeVariables ();
  if (fNumNV < 1)
    {
      fMessage << "\n No nodal variables found.";
      return;
    }

  iArrayT nodes;
  InitializeNodePoints(nodes, fItemIndex);
  fNumItems = nodes.Length();
  if (fNumItems < 1)
    {
      fMessage << "\n No node points found.";
      return;
    }

  fItemNames.Allocate (fNumItems);
  for (int i=0; i < fNumItems; i++)
    fItemNames[i].Append (nodes[i]);
}

void ExtractNode::TranslateVariables (void)
{
  PrepFiles (fNVUsed, fNodeLabels);

  int numnodes, numdims;
  fModel.CoordinateDimensions (numnodes, numdims);
  fCoordinates = fModel.Coordinates ();

  fVarData.Allocate (numnodes, fNumNV);
  for (int t=0; t < fNumTS; t++)
    {
      fModel.AllNodeVariables (fTimeIncs[t], fVarData);
      WriteVarData (fNVUsed, t);
    }
}

