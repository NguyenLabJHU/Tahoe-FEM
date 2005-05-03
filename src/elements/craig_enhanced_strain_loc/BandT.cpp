#include "BandT.h"
#include "SSEnhLocCraigT.h"

#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
BandT::BandT(const dArrayT normal, const dArrayT slipDir, const dArrayT
perpSlipDir, dArrayT &coord, double h_delta, double residCohesion, ArrayT<dSymMatrixT> stressList, SSEnhLocCraigT *element):
fNormal(normal),
fSlipDir(slipDir),
fPerpSlipDir(perpSlipDir),
fCoords(coord), 
fH_delta(h_delta),
fEffectiveSoftening(h_delta),
fResidualCohesion(residCohesion),
fJump(0.0),
fJumpIncrement(0.0),
fIsBandActive(true),
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

  /* get stresses */

  /* should create deep copy */
  fStress_List = stressList;

  /*
  int numIP = currentElement->NumIP();
  fStress_List.Dimension(numIP);
  //ShapeFunctionT* shapes = currentElement->fShapes;
  //shapes -> TopIP();
  //currentElement -> fShapes -> TopIP();
  for (int i = 0; i < numIP; i++)
    {
      //shapes -> NextIP();
      //currentElement -> fShapes -> TopIP();
      fStress_List[i].Dimension(currentElement->NumSD());
      fStress_List[i] = stressList[i];
    }
  */

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

double BandT::H_delta()const
{
  return fH_delta;
}

double BandT::ResidualCohesion() const
{
  return fResidualCohesion;
}


double BandT::Jump() const
{
  return fJump;
}

double BandT::JumpIncrement() const
{
  return fJumpIncrement;
}

/*
void BandT::IncrementJump ()
{
  fJump += fJumpIncrement;
}
*/

void BandT::StoreJumpIncrement(double increment)
{
  fJumpIncrement = increment;
}

void BandT::SetEffectiveSoftening(double effectiveSoftening)
{
  fEffectiveSoftening = effectiveSoftening;
}

double BandT::EffectiveSoftening()
{
  //cout << "EffectiveSoftening = " << fEffectiveSoftening;
  return fEffectiveSoftening;
}

void BandT::SetActive(bool active)
{
  fIsBandActive = active;
}

bool BandT::IsActive()
{
  return fIsBandActive;
}

dArrayT& BandT::Coords()
{
 return fCoords;
}

/*
void BandT::CloseStep()
{
  	  IncrementJump();
	  cout << "JumpIncrement = " << fBand -> JumpIncrement() << endl;
	  cout << "Jump = " << fBand -> Jump() << endl; 

	  ShapefunctionT* shapes = currentElement->fShapes;

	  shapes->TopIP();
	  while (shapes->NextIP())
	    {
	      dSymMatrixT strainIncr = currentElement->fStrain_List [currentElement->CurrIP()];
	      strainIncr -= currentElement->fStrain_last_List [currentElement->CurrIP()]; 
	      dSymMatrixT stressIncr(currentElement->NumSD());
	      stressIncr.A_ijkl_B_kl(currentElement-> fCurrMaterial -> ce_ijkl(), strainIncr);
	      fStress_List[currentElement->CurrIP()] += stressIncr; 
	    }
}
*/


dSymMatrixT BandT::Stress_List(int ip)
{
  //put in check to see if ip is valid?

  return fStress_List[ip];
}

void BandT::IncrementStress(dSymMatrixT stressIncr, int ip)
{
  fStress_List [ip] += stressIncr;
}

/*
void BandT::UpdateCohesion()
{
  fResidualCohesion += fH_delta * fabs(fJumpIncrement);
  cout << "fResidualCohesion = " << fResidualCohesion << endl;

  if(fResidualCohesion < 0.0)
    {
      fResidualCohesion = 0.0;
      fH_delta = 0.0; //no more softening possible
    }
}
*/

void BandT::FlipSlipDir()
{
  /* keep component parallel to normal same, multiply component 
     paraellel to perpSlipDir by -1.0 */
  fSlipDir.AddScaled(-2.0*fSlipDir.Dot(fSlipDir, fPerpSlipDir), fPerpSlipDir);
  fPerpSlipDir *= -1.0;

  /* NOTE: if there is reverse loading after the first step, fJump becomes
  a cumulative slip vector, but not a value in the direction. It does not
  play a critical role solving the BVP (since the everything is solving in
  terms of the incremental jump, but the strains that are output
  are not the regular strains, as they would be for a monotonic loading 
  problem.
   */
}

void BandT::CloseStep()
{
  //update jump
  fJump += fJumpIncrement;

  //update cohesion
  fResidualCohesion += fH_delta * fabs(fJumpIncrement);
  //cout << "fResidualCohesion = " << fResidualCohesion << endl;

  if(fResidualCohesion <= 0.0)
    {
      fResidualCohesion = 0.0;
      fH_delta = 0.0; //no more softening possible
    }

  //update stress list


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
  cout << "fNormal = " << fNormal << endl;
  cout << "fSlipDir = " << fSlipDir << endl << endl; 

  dArrayT nodalCoord(kNSD);

  for (int i = 0; i < nen; i++)
    {
      for (int j = 0; j < nodalCoords.MinorDim(); j++)
	nodalCoord [j] = nodalCoords(i,j);

      nodalCoord -= coord;
      // if dot product is greater than one, then nodes
      if (nodalCoord.Dot(nodalCoord, fNormal) > 0.0)
	{
	  fActiveNodes.Append(i);
	  cout << "node " << i << " is active.Coords =\n" << nodalCoord << endl; 
	}
    }

}

void BandT::SetEndPoints(dArrayT& coord)
{
  //implement
}
