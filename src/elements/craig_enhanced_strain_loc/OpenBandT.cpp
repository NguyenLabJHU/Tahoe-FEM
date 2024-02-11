#include "OpenBandT.h"
#include "SSEnhLocCraigT.h"

#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
OpenBandT::OpenBandT(const dArrayT normal, dArrayT shearDir, dArrayT &coord, double residCohesion,
 ArrayT<dSymMatrixT> stressList, SSEnhLocOpenT *element):
fNormal(normal),
//fSlipDir(slipDir),
//fPerpSlipDir(perpSlipDir),
fShearDir(shearDir),
fCoords(coord), 
//fH_delta(h_delta),
//fEffectiveSoftening(h_delta),
fInitialCohesion(residCohesion),
fLastCohesion(residCohesion),
fCohesion(residCohesion),
//fJump(0.0),
//fLastJump(0.0),
fBandState(2), //kDamage
currentElement(element)
{

  kNSD = normal.Length();

  if (kNSD > 2)
    {
      cout << "OpenBandT::OpenBandT, not ready for 3 dimensions.\n";
      throw ExceptionT::kGeneralFail;
    }
	
 fJump.Dimension(kNSD);
 fJump = 0.0;
 fLastJump.Dimension(kNSD);
 fLastJump = 0.0;
		
 //fShearDir.Dimension(2);
 //fShearDir(0) = fNormal(1);
 //fShearDir(1) = -1.0 * fNormal(0);

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

const iAutoArrayT& OpenBandT::ActiveNodes() const
{ 
  return fActiveNodes;
}

const dArrayT& OpenBandT::Normal() const
{
  return fNormal;
}

const dArrayT& OpenBandT::ShearDir() const
{
  return fShearDir;
}

/*
const dArrayT& OpenBandT::SlipDir() const
{
  return fSlipDir;
}

const dArrayT& OpenBandT::PerpSlipDir() const
{
  return fPerpSlipDir;
}
*/

/*
double OpenBandT::H_delta()const
{
  return fH_delta;
}
*/

double OpenBandT::InitialCohesion() const
{
 return fInitialCohesion;
}

/*
double LastCohesion() const;
{
	return fLastCohesion;
}
*/

double OpenBandT::Cohesion() const
{
	return fCohesion;
}

void OpenBandT::UpdateCohesion(double cohesion)
{
	fCohesion = cohesion;
}

/*
double OpenBandT::ResidualCohesion() const
{
  return fResidualCohesion;
}
*/

dArrayT OpenBandT::Jump() const
{
  return fJump;
}

dArrayT OpenBandT::LastJump() const
{
  return fLastJump;
}

/*
double OpenBandT::JumpIncrement() const
{
  return fJumpIncrement;
}
*/

void OpenBandT::UpdateJump()
{
  fLastJump = fJump;
}


void OpenBandT::StoreJump(dArrayT jump)
{
  fJump = jump;
}


/*
void OpenBandT::SetEffectiveSoftening(double effectiveSoftening)
{
  fEffectiveSoftening = effectiveSoftening;
}

double OpenBandT::EffectiveSoftening()
{
  //cout << "EffectiveSoftening = " << fEffectiveSoftening;
  return fEffectiveSoftening;
}
*/

void OpenBandT::SetBandState(int state)
{
  cout << "band state = " << state << endl;
  fBandState = state;
}

int OpenBandT::BandState()
{
  return fBandState;
}

dArrayT& OpenBandT::Coords()
{
 return fCoords;
}


void OpenBandT::CloseStep()
{
	UpdateJump();
	
	/*update cohesion */
	fLastCohesion = fCohesion;
	
}



dSymMatrixT OpenBandT::Stress_List(int ip)
{
  //put in check to see if ip is valid?
  //cout << "ip = " << ip << endl;
  //cout << "fStress_List[ip] = " << fStress_List[ip] << endl;

  return fStress_List[ip];
}

void OpenBandT::IncrementStress(dSymMatrixT stressIncr, int ip)
{
  fStress_List [ip] += stressIncr;
}

/*
void OpenBandT::UpdateCohesion()
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

void OpenBandT::FlipSlipDir()
{
  /* keep component parallel to normal same, multiply component 
     paraellel to perpSlipDir by -1.0 */
  //fSlipDir.AddScaled(-2.0*fSlipDir.Dot(fSlipDir, fPerpSlipDir), fPerpSlipDir);
  fShearDir *= -1.0;

  /* NOTE: if there is reverse loading after the first step, fJump becomes
  a cumulative slip vector, but not a value in the direction. It does not
  play a critical role solving the BVP (since the everything is solving in
  terms of the incremental jump, but the strains that are output
  are not the regular strains, as they would be for a monotonic loading 
  problem.
   */
}



/*---------------------------------------------------------------------
private
-----------------------------------------------------------------------*/

void OpenBandT::ActivateNodes(dArrayT& coord)
{
  //ElementCardT element = CurrentElement();
  //iArrayT nodes = element.NodesX();
  LocalArrayT nodalCoords = currentElement -> InitialCoordinates();
  int nen = currentElement->NumElementNodes();

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
	  cout << "node " << i << " is active. Coords =\n" << nodalCoord << endl; 
	}
    }

}

void OpenBandT::SetEndPoints(dArrayT& coord)
{
  //implement
}
