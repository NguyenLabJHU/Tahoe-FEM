
#include "ExtractQuad.h"

ExtractQuad::ExtractQuad (ostream& out) :
  ExtractIOManager (out)
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
  fModel.ElementGroupDimensions (fElementGroup, numelems, numelemnodes);
  if (fNumItems < 1)
    {
      fMessage << "\n No elements found.";
      return;
    }

  // change fItems and fItemIndex from elements to quad points
  int numquadpts = fModel.NumElementQuadPoints (fElementName);
  iArrayT elementmap (numelems);
  fModel.ElementMap (fElementName, elementmap);

  fNumItems = numelems * numquadpts;
  fItemNames.Allocate (fNumItems);
  fItemIndex.Allocate (fNumItems);

  fItemIndex.SetValueToPosition ();
  for (int i=0, j=0; i < numelems; i++)
    for (int k=0; k < numquadpts; k++)
      {
	fItemNames [j].Append (elementmap[i]);
	fItemNames [j].Append ("_", k+1);
      }
}

void ExtractQuad::TranslateVariables (void)
{
  PrepFiles (fQVUsed, fQuadratureLabels);

  // need to allocate correct number of quad points
  int numelems, numelemnodes;
  fModel.ElementGroupDimensions (fElementGroup, numelems, numelemnodes);
  int numquadpts = fModel.NumElementQuadPoints (fElementName);

  fVarData.Allocate (numelems * numquadpts, fNumQV);
  for (int t=0; t < fNumTS; t++)
    {
      fModel.QuadratureVariables (fTimeIncs[t], fElementName, fVarData);
      WriteVarData (fQVUsed, t);
    }
}

