
#include "BandT.h"
#include "SSEnhLocCraigT.h"

using namespace Tahoe;

/* constructor */
BandT::BandT(const dArrayT normal, const dArrayT slipDir, const dArrayT perpSlipDir, dArrayT &coord, SSEnhLocCraigT *element):
fNormal(normal),
fSlipDir(slipDir),
fPerpSlipDir(perpSlipDir),
fJump(0.0),
fJumpIncrement(0.0),
currentElement(element)
{
  kNSD = normal.Length();

  if (kNSD > 2)
    {
      cout << "BandT::BandT, not ready for 3 dimensions.\n";
      throw ExceptionT::kGeneralFail;
    }

  ActivateNodes(coord);
  SetEndPoints(coord);
}

const iAutoArrayT& BandT::ActiveNodes() const
{ 
  return fActiveNodes;
}

const dArrayT& BandT::Normal() const
{
  return fNormal;
}

const dArrayT& BandT::SlipDir() const
{
  return fSlipDir;
}

const dArrayT& BandT::PerpSlipDir() const
{
  return fPerpSlipDir;
}

double BandT::Jump() const
{
  return fJump;
}

void BandT::IncrementJump (const double increment)
{
  fJump += increment;
}


/*---------------------------------------------------------------------
private
-----------------------------------------------------------------------*/

void BandT::ActivateNodes(dArrayT& coord)
{
  //ElementCardT element = CurrentElement();
  //iArrayT nodes = element.NodesX();
  LocalArrayT nodalCoords = currentElement->InitialCoordinates();
  int nen = currentElement->NumElementNodes();

  dArrayT nodalCoord;

  for (int i = 0; i < nen; i++)
    {
      for (int j = 0; j < nodalCoords.MinorDim(); j++)
	nodalCoord [j] = nodalCoords(i,j);

      nodalCoord -= coord;
      // if dot product is greater than one, then nodes
      if (nodalCoord.Dot(nodalCoord, fNormal) > 0.0)
	fActiveNodes.Append(i);

    }

}

void BandT::SetEndPoints(dArrayT& coord)
{
  //implement
}
