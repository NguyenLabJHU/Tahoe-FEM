// $Id: Quad2Tri.cpp,v 1.7 2004-04-08 23:45:35 paklein Exp $
// created: SAW 12/21/99
#include "Quad2Tri.h"

#include "ExceptionT.h"
#include "dArrayT.h"
#include "MakeCSE_FEManager.h"

const char* fMethodName [] = { "X-Method",
			      "Slash Method",
			      "Back Slash Method",
			      "Star Method" };

using namespace Tahoe;

Quad2Tri::Quad2Tri (ostream& fMainOut, NodeManagerPrimitive& NMP, CSEConstants::SplitMethodT method, const StringT& ID) :
  MakeCSE_ElementBaseT (fMainOut, ID),
  theNodes (&NMP),
  fMethod (method)
{
}

// *********** PROTECTED *************

void Quad2Tri::EchoConnectivity (ModelManagerT& theInput)
{
  // read quadrilateral data
  GeometryT::CodeT geocode;
  ReadConnectivity (theInput, geocode, fConn);

  // translate fConn into fNodeNums, by splitting quadrilaterals
  out << "\n  Translating Element Group ID " 
      << MakeCSE_ElementBaseT::GroupNumber()
      << " from Quad to Tri by " << fMethodName[fMethod] << "\n";
  Translate ();

  // initialize triangle data
  fGeometryCode = GeometryT::kTriangle;
  InitializeConnectivity ();
}

void Quad2Tri::EchoSideSets (ModelManagerT& model, MakeCSE_IOManager& theInput)
{
  ArrayT<iArray2DT> sidesets;
  ReadSideSetData (model, theInput, sidesets);

  switch (fMethod)
    {
    case CSEConstants::kXMethod:
      XMethodSideSets (sidesets);
      break;
    case CSEConstants::kSlashMethod:
      SlashSideSets (sidesets);
      break;
    case CSEConstants::kBackSlashMethod:
      BackSlashSideSets (sidesets);
      break;
    case CSEConstants::kStarMethod:
      StarSideSets (sidesets);
      break;
    default:
      cout << "Quad2Tri::EchoSideSets, unknown method: " 
	   << fMethod << endl;
      throw ExceptionT::kGeneralFail;
    }

  CheckAllSideSets ();    
}

// *********** PRIVATE *************

void Quad2Tri::Translate (void)
{
  // allocate space
  int numQuadNodes = fConn.MinorDim();
  Allocate (numQuadNodes);

  // translate
  const dArray2DT& coords = theNodes->InitialCoordinates();
  int *quad = fConn.Pointer();
  int count = 0;
  for (int i=0; i < fConn.MajorDim(); i++, quad += numQuadNodes)
    {
      switch (fMethod)
	{
	case CSEConstants::kXMethod:
	  {
	    int newnode = ElementCentroid (quad, numQuadNodes, coords);
	    XMethodNumbering (count, newnode, quad);
	    break;
	  }
	case CSEConstants::kSlashMethod:
	  SlashNumbering (count, quad);
	  break;
	case CSEConstants::kBackSlashMethod:
	  BackSlashNumbering (count, quad);
	  break;
	case CSEConstants::kStarMethod:
	  {
	    int newnode = ElementCentroid (quad, numQuadNodes, coords);
	    StarNumbering (count, newnode, quad);
	    break;
	  }
	}
    }

  // check
  if (count != fNodeNums.MajorDim())
    {
      cout << "\n\nQuad2Tri::Translate: error numbering tris\n"
	   << count << " created\n"
	   << fNodeNums.MajorDim() << " initialized\n\n";
      throw ExceptionT::kGeneralFail;
    }

  // free temp space
  fConn.Free();
}

void Quad2Tri::Allocate (int numQuadNodes)
{
  int quadnodesrequired;
  int numTriNodes;
  int numCreated;
  // determine necessary space
  switch (fMethod)
    {
    case CSEConstants::kXMethod:
      numTriNodes = 3;
      numCreated = 4;
      quadnodesrequired = 4;
      break;
    case CSEConstants::kSlashMethod:
    case CSEConstants::kBackSlashMethod:
      numTriNodes = 3;
      numCreated = 2;
      quadnodesrequired = 4;
      break;
    case CSEConstants::kStarMethod:
      numTriNodes = 3;
      numCreated = 8;
      quadnodesrequired = 8;
      break;
    default:
      cout << "Quad2Tri::Allocate, unknown method: " << fMethod << endl;
      throw ExceptionT::kGeneralFail;
    }

  // check
  if (numQuadNodes != quadnodesrequired)
    {
      cout << "\n\nQuad2Tri::Allocate, Only quad" 
	   << quadnodesrequired << " are supported for "
	   << fMethodName[fMethod] << " " << numQuadNodes<< "\n";
      throw ExceptionT::kBadInputValue;
    }

  // allocate space
  fNodeNums.Allocate (fConn.MajorDim() *numCreated, numTriNodes);
  fNodeNums = CSEConstants::kNotSet;
}

/* find centroid of element */
int Quad2Tri::ElementCentroid (int* quad, int numQuadNodes, const dArray2DT& coords)
{
      int dof = coords.MinorDim();
      dArrayT xbar (dof);
      xbar = 0.0;
      for (int k=0; k < numQuadNodes; k++)
	for (int y=0; y < dof; y++)
	  xbar[y] += coords (quad[k], y);
      for (int h=0; h < dof; h++)
	xbar[h] = xbar[h]/numQuadNodes;

      return theNodes->AddCoord (xbar);
}

void Quad2Tri::XMethodNumbering (int& count, int newnode, int* quad)
{
      fNodeNums (count, 0) = quad[0];
      fNodeNums (count, 1) = quad[1];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[1];
      fNodeNums (count, 1) = quad[2];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[2];
      fNodeNums (count, 1) = quad[3];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[3];
      fNodeNums (count, 1) = quad[0];
      fNodeNums (count, 2) = newnode;
      count++;
}

void Quad2Tri::SlashNumbering (int& count, int* quad)
{
      fNodeNums (count, 0) = quad[0];
      fNodeNums (count, 1) = quad[1];
      fNodeNums (count, 2) = quad[2];
      count++;
      fNodeNums (count, 0) = quad[2];
      fNodeNums (count, 1) = quad[3];
      fNodeNums (count, 2) = quad[0];
      count++;
}

void Quad2Tri::BackSlashNumbering (int& count, int* quad)
{
      fNodeNums (count, 0) = quad[3];
      fNodeNums (count, 1) = quad[0];
      fNodeNums (count, 2) = quad[1];
      count++;
      fNodeNums (count, 0) = quad[1];
      fNodeNums (count, 1) = quad[2];
      fNodeNums (count, 2) = quad[3];
      count++;
}

void Quad2Tri::StarNumbering (int& count, int newnode, int* quad)
{
      fNodeNums (count, 0) = quad[0];
      fNodeNums (count, 1) = quad[4];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[4];
      fNodeNums (count, 1) = quad[1];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[1];
      fNodeNums (count, 1) = quad[5];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[5];
      fNodeNums (count, 1) = quad[2];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[2];
      fNodeNums (count, 1) = quad[6];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[6];
      fNodeNums (count, 1) = quad[3];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[3];
      fNodeNums (count, 1) = quad[7];
      fNodeNums (count, 2) = newnode;
      count++;
      fNodeNums (count, 0) = quad[7];
      fNodeNums (count, 1) = quad[0];
      fNodeNums (count, 2) = newnode;
      count++;
}

// with this method original quad facets are always a triangle 0 facet
void Quad2Tri::XMethodSideSets (ArrayT<iArray2DT>& sidesets)
{
  fSideSetData.Allocate (sidesets.Length());
  for (int s=0; s < sidesets.Length(); s++)
    {
      fSideSetData[s].Allocate (sidesets[s].MajorDim(), 2);
      fSideSetData[s] = 0;
      int *e = sidesets[s].Pointer();
      int *f = sidesets[s].Pointer(1);
      int *en = fSideSetData[s].Pointer();
      for (int i=0; i < sidesets[s].MajorDim(); i++)
	{
	  *en = *e *4 + *f;
	  e += 2;
	  f += 2;
	  en += 2;
	}
    }
}

void Quad2Tri::SlashSideSets (ArrayT<iArray2DT>& sidesets)
{
  fSideSetData.Allocate (sidesets.Length());
  for (int s=0; s < sidesets.Length(); s++)
    {
      fSideSetData[s].Allocate (sidesets[s].MajorDim(), 2);
      fSideSetData[s] = 0;
      int *e = sidesets[s].Pointer();
      int *f = sidesets[s].Pointer(1);
      int *enew = fSideSetData[s].Pointer();
      int *fnew = fSideSetData[s].Pointer(1);
      for (int i=0; i < sidesets[s].MajorDim(); i++)
	{
	  *enew = *e * 2 + ((*f < 2) ? 0 : 1);
	  *fnew = *f - ((*f < 2) ? 0 : 2);
	  e += 2;
	  f += 2;
	  enew += 2;
	  fnew += 2;
	}
    }
}

void Quad2Tri::BackSlashSideSets (ArrayT<iArray2DT>& sidesets)
{
  fSideSetData.Allocate (sidesets.Length());
  for (int s=0; s < sidesets.Length(); s++)
    {
      fSideSetData[s].Allocate (sidesets[s].MajorDim(), 2);
      fSideSetData[s] = 0;
      int *e = sidesets[s].Pointer();
      int *f = sidesets[s].Pointer(1);
      int *enew = fSideSetData[s].Pointer();
      int *fnew = fSideSetData[s].Pointer(1);
      for (int i=0; i < sidesets[s].MajorDim(); i++)
	{
	  *enew = *e * 2;
	  if (*f == 1 || *f == 2) 
	    {
	      *enew += 1;
	      *fnew = *f - 1;
	    }
	  else if (*f == 0) *fnew = 1;
	  else *fnew = 0;
	  e += 2;
	  f += 2;
	  enew += 2;
	  fnew += 2;
	}
    }
}

// with this method original quad facets are always a triangle 0 facet
void Quad2Tri::StarSideSets (ArrayT<iArray2DT>& sidesets)
{
  fSideSetData.Allocate (sidesets.Length());
  for (int s=0; s < sidesets.Length(); s++)
    {
      fSideSetData[s].Allocate (2*sidesets[s].MajorDim(), 2);
      fSideSetData[s] = 0;
      int *e = sidesets[s].Pointer();
      int *f = sidesets[s].Pointer(1);
      int *enew = fSideSetData[s].Pointer();
      for (int i=0; i < sidesets[s].MajorDim(); i++)
	{
	  *enew = *e *8 + *f *2;
	  enew += 2;
	  *enew = *e *8 + *f *2 + 1;
	  enew += 2;
	  e += 2;
	  f += 2;
	}
    }
}
