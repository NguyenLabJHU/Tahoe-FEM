
#include "ExtractQuad.h"

using namespace Tahoe;

ExtractQuad::ExtractQuad (ostream& out, istream& in, bool write) :
  ExtractIOManager (out, in, write)
{
}

/**************** PROTECTED **********************/

void ExtractQuad::Initialize (void)
{
  InitializeQuadVariables ();
  if (fNumQV < 1)
    {
      fMessage << "\n No quadrature variables found.";
      return;
    }

  int numelems, numelemnodes;
  InitializeElements(fElementGroup, fElementName);
  fModel.ElementGroupDimensions (fElementName, numelems, numelemnodes);
  if (numelems < 1)
    {
      fMessage << "\n No elements found.";
      return;
    }

  // change from elements to quad points
  int numquadpts = fModel.NumElementQuadPoints (fElementName);
  iArrayT elementmap (numelems);
  fModel.ElementIDs (fElementName, elementmap);

  fNumItems = numelems * numquadpts;
  fItemNames.Allocate (fNumItems);
  fItemIndex.Allocate (fNumItems);

  fItemIndex.SetValueToPosition ();
  for (int i=0, j=0; i < numelems; i++)
    for (int k=0; k < numquadpts; k++, j++)
      {
	fItemNames [j] = "";
	fItemNames [j].Append (elementmap[i]);
	fItemNames [j].Append ("_qp", k+1);
      }
}

void ExtractQuad::TranslateVariables (void)
{
  PrepFiles (fQVUsed, fQuadratureLabels);

  // need to allocate correct number of quad points
  int numelems, numelemnodes;
  fModel.ElementGroupDimensions (fElementName, numelems, numelemnodes);
  int numquadpts = fModel.NumElementQuadPoints (fElementName);

  fVarData.Allocate (numelems * numquadpts, fNumQV);
  for (int t=0; t < fNumTS; t++)
    {
      fModel.QuadratureVariables (fTimeIncs[t], fElementName, fVarData);
      WriteVarData (fQVUsed, t);
    }
}

