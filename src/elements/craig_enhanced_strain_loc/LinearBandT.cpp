#include "LinearBandT.h"
#include "SSEnhLocLinearT.h"

#include "ShapeFunctionT.h"


using namespace Tahoe;

/* constructor */
LinearBandT::LinearBandT(const dArrayT normal, const dArrayT slipDir, const dArrayT
perpSlipDir, dArrayT &coord, double h_delta, ArrayT<dSymMatrixT> stressList,
SSEnhLocLinearT *element):
fNormal(normal),
//fSlipDir(slipDir),
//fPerpSlipDir(perpSlipDir),
fCoords(coord), 
//fH_delta(h_delta),
//fEffectiveSoftening(h_delta), //check - this is dArrayT now, can it be initialized w/ double?
//fResidualCohesion(residCohesion),
//fJump(0.0),
//fJumpIncrement(0.0),
//fDistanceBetweenSurfaceIPs(distanceBetweenSurfaceIPs),
fIsBandActive(true), //check - as above
fCurrentElement(element)
{
  kNSD = normal.Length();

  if (kNSD > 2)
  {
      cout << "LinearBandT::LinearBandT, not ready for 3 dimensions.\n";
      throw ExceptionT::kGeneralFail;
  }
  
  fH_delta.Dimension(kNSD);
  for (int i = 0; i < kNSD; i++)
	fH_delta [i] = h_delta;
	
  fEffectiveSoftening = fH_delta;

  /* get stresses */ 
  /* should create deep copy */
  fStress_List = stressList;
  //fIPBandCoords = IPBandCoords;
  //fSurfaceIPTractions = surfaceIPTractions;

  fJump.Dimension(kNumSurfaceIPs);
  fJumpIncrement.Dimension(kNumSurfaceIPs);
  fJump = 0.0;
  fJumpIncrement = 0.0;  
  fSlipDir.Dimension(kNumSurfaceIPs);
  fPerpSlipDir.Dimension(kNumSurfaceIPs);
  
  for (int i=0; i < kNumSurfaceIPs; i++)
  {
	/*
	(fJump[i]).Dimension(kNSD);
	fJump[i] = 0.0;
	(fJumpIncrement[i]).Dimension(kNSD);
	fJumpIncrement[i] = 0.0; */
	
	fSlipDir[i] = slipDir;
	fPerpSlipDir[i] = perpSlipDir;
  }
  
  ActivateNodes(coord);
  ActivateBulkIPs(coord);
  SetEndPoints(coord);
}

const iAutoArrayT& LinearBandT::ActiveNodes() const
{ 
  return fActiveNodes;
}

const dArrayT& LinearBandT::Normal() const
{
  return fNormal;
}

const dArrayT& LinearBandT::SlipDir(int bandIP) const
{
  return fSlipDir[bandIP];
}

const dArrayT& LinearBandT::PerpSlipDir(int bandIP) const
{
  return fPerpSlipDir[bandIP];
}

double LinearBandT::H_delta(int bandIP)const
{
  return fH_delta[bandIP];
}

double LinearBandT::ResidualCohesion(int ip) const
{
  return fResidualCohesion[ip];
}


double LinearBandT::Jump(int ip) const
{
  return fJump[ip];
}

double LinearBandT::JumpIncrement(int ip) const
{
  //cout << "fJumpIncrement = \n" << fJumpIncrement << endl; 
  return fJumpIncrement[ip];
}


void LinearBandT::IncrementJump (int ip)
{
  fJump[ip] += fJumpIncrement[ip];
}

/* check - should this be done at once or by ip? */
void LinearBandT::StoreJumpIncrement(double increment, int bandIP)
{
  fJumpIncrement[bandIP] = increment;
  //cout << "fJumpIncrement = " << fJumpIncrement << endl;
}

void LinearBandT::SetEffectiveSoftening(double effectiveSoftening, int bandIP)
{
  fEffectiveSoftening [bandIP] = effectiveSoftening;
}

double LinearBandT::EffectiveSoftening(int bandIP)
{
  //cout << "EffectiveSoftening = " << fEffectiveSoftening;
  return fEffectiveSoftening[bandIP];
}

void LinearBandT::SetActive(int ip, bool active)
{
  fIsBandActive[ip] = active;
}

bool LinearBandT::IsActive(int ip)
{
  return fIsBandActive[ip];
}

dArrayT& LinearBandT::Coords()
{
 return fCoords;
}

/*
void LinearBandT::CloseStep()
{
  	  IncrementJump();
	  cout << "JumpIncrement = " << fBand -> JumpIncrement() << endl;
	  cout << "Jump = " << fBand -> Jump() << endl; 

	  ShapefunctionT* shapes = fCurrentElement->fShapes;

	  shapes->TopIP();
	  while (shapes->NextIP())
	    {
	      dSymMatrixT strainIncr = fCurrentElement->fStrain_List [fCurrentElement->CurrIP()];
	      strainIncr -= fCurrentElement->fStrain_last_List [fCurrentElement->CurrIP()]; 
	      dSymMatrixT stressIncr(fCurrentElement->NumSD());
	      stressIncr.A_ijkl_B_kl(fCurrentElement-> fCurrMaterial -> ce_ijkl(), strainIncr);
	      fStress_List[fCurrentElement->CurrIP()] += stressIncr; 
	    }
}
*/


dSymMatrixT LinearBandT::Stress_List(int ip)
{
  //put in check to see if ip is valid?

  return fStress_List[ip];
}

double LinearBandT::IPBandCoord(int ip)
{
	return fIPBandCoords[ip];
}

void LinearBandT::IncrementStress(dSymMatrixT stressIncr, int ip)
{
  fStress_List [ip] += stressIncr;
}

void LinearBandT::IncrementTractionAtBandIP(dArrayT increment, int bandIP)
{
	fSurfaceIPTractions [bandIP] += increment;
}


/*
void LinearBandT::UpdateCohesion()
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


void LinearBandT::FlipSlipDir(int ip)
{
  /* keep component parallel to normal same, multiply component 
     paraellel to perpSlipDir by -1.0 */
  (fSlipDir[ip]).AddScaled(-2.0*fSlipDir[ip].Dot(fSlipDir[ip], fPerpSlipDir[ip]), fPerpSlipDir[ip]);
  fPerpSlipDir[ip] *= -1.0;

  /* NOTE: if there is reverse loading after the first step, fJump becomes
  a cumulative slip vector, but not a value in the direction. It does not
  play a critical role solving the BVP (since the everything is solving in
  terms of the incremental jump, but the strains that are output
  are not the regular strains, as they would be for a monotonic loading 
  problem.
   */
}


void LinearBandT::CloseStep()
{
  //update jump
  fJump += fJumpIncrement;

  //update cohesion
  for (int i = 0; i < kNumSurfaceIPs; i++)
  {
	fResidualCohesion[i] += fH_delta[i] * fabs(fJumpIncrement[i]);

	if(fResidualCohesion[i] <= 0.0)
    {
      fResidualCohesion[i] = 0.0;
      fH_delta[i] = 0.0; //no more softening possible
    }
  }
  //update stress list - currently done separately
}

/*---------------------------------------------------------------------
private
-----------------------------------------------------------------------*/

void LinearBandT::ActivateBulkIPs(dArrayT& coord)
{
	fIsBulkIPActive.Dimension(fCurrentElement -> NumIP()); 
    
	for (int i = 0; i < fCurrentElement -> NumIP(); i++)
	{
		dArrayT ipCoord(kNSD);
		fCurrentElement -> ShapeFunction().IPCoords(ipCoord,i); 
		dArrayT vector = ipCoord;
		vector -= coord;
	
		if(vector.Dot(vector, fNormal) > 0.0)
			fIsBulkIPActive = true;
		else 
			fIsBulkIPActive = false;
	}
}

void LinearBandT::ActivateNodes(dArrayT& coord)
{
  //ElementCardT element = CurrentElement();
  //iArrayT nodes = element.NodesX();
  LocalArrayT nodalCoords = fCurrentElement->InitialCoordinates();
  int nen = fCurrentElement->NumElementNodes();
  //cout << "fNormal = " << fNormal << endl;
  //cout << "fSlipDir = " << fSlipDir << endl << endl; 

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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int numSides;
switch (fCurrentElement -> GeometryCode())
{
 case GeometryT::kQuadrilateral:
   {
     numSides = 4;
     break;
   }
 case GeometryT::kTriangle: 
   {
     numSides = 3;
     break;
   }
 default:
   {
     cout << "SSEnhLocLinearT::AddNewEdgeElements, geometry not implemented. \n" << flush;
     throw ExceptionT::kGeneralFail;
   }
}

  ArrayT<dArrayT> endPoints(2);//2 okay for 2D, fix for 3D
  for (int i=0; i < 2; i++)
	(endPoints[i]).Dimension(kNSD);
  
  dArrayT nodalCoord1(kNSD), nodalCoord2(kNSD);
  dArrayT nodalCords = fCurrentElement -> InitialCoordinates();

  int sideNumber = 0;
  /* loop over nodes to find side between active and passive nodes */  
   for(int i = 0; i < numSides; i++)
	if ((fActiveNodes.HasValue((i+1) % numSides) && !(fActiveNodes.HasValue(i)))
	|| (!(fActiveNodes.HasValue((i+1) % numSides)) && fActiveNodes.HasValue(i)))
	{
	  /* set nodal coordinates of side */
		for (int j = 0; j < nodalCoords.MinorDim(); j++)
		{		
			nodalCoord1 [j] = nodalCoords(i,j);
			nodalCoord2 [j] = nodalCoords((i+1) % numSides, j); //replace numSides w/ num Ele nodes!
		}
	  
	  /*find endpoints on that side*/
	  //assumes straight sides
		dArrayT sideVector = nodalCoord2;
		sideVector -= nodalCoord1;
		
		double alpha = sideVector[0] * (fCoords[1] - nodalCoord1[1]) -
		sideVector[1] * (fCoords[0] - nodalCoord1[0]);
		alpha /= sideVector[1] * fPerpSlipDir[0][0] - sideVector[0] *
		fPerpSlipDir[0][1];
		
		endPoints[sideNumber] = fCoords;
		endPoints[sideNumber].AddScaled(alpha, fPerpSlipDir[0]);
	  
	  sideNumber++;
	}
	
	/* reset fCoords to center of band */
	fCoords = endPoints[0];
	fCoords += endPoints[1];
	fCoords /= 2.0;
	
	/* find band length */
	dArrayT bandVector = endPoints[1];
	bandVector -= endPoints[0];
	
	double bandLength = bandVector.Magnitude();
	fDistanceBetweenSurfaceIPs = bandLength/sqrt(3.0);
	
	/* set coords of surface IPs */
	fSurfaceIPCoords.Dimension(kNumSurfaceIPs);
	for (int i=0; i < kNumSurfaceIPs; i++)
		(fSurfaceIPCoords[i]).Dimension(kNSD);
	
	fSurfaceIPCoords[0] = fCoords;
	fSurfaceIPCoords[0].AddScaled(fDistanceBetweenSurfaceIPs/-2.0, fPerpSlipDir[0]);
	
	fSurfaceIPCoords[1] = fCoords;
	fSurfaceIPCoords[1].AddScaled(fDistanceBetweenSurfaceIPs/2.0, fPerpSlipDir[0]);
	
	/* determine surface coordinate of bulk IPs */
	fIPBandCoords.Dimension(fCurrentElement -> NumIP());
	fCurrentElement -> fShapes -> TopIP();
	//fShapes -> TopIP();
	for (int i = 0; i < fCurrentElement -> NumIP(); i++)
	{
		fCurrentElement -> fShapes -> NextIP();
		
		/* get ip coordinate */
		dArrayT ipCoord(kNSD);
		fCurrentElement -> fShapes -> IPCoords(ipCoord);	
		
		/* subtract band center */
		ipCoord -= fCoords;
		
		fIPBandCoords[i] = dArrayT::Dot(ipCoord, fPerpSlipDir[0]);
	}

	/*extrapolate bulk ip stresses to nodal values */
	ArrayT<dSymMatrixT> nodalStressList = NodalStressList();
	fSurfaceIPTractions.Dimension(kNumSurfaceIPs);	
	
	/* interpolate nodal stresses to surface ip's */
	for (int i = 0; i < kNumSurfaceIPs; i ++)
	{

		/* transform coordinates to Parent Domain (Local Coordinates)*/
		dArrayT localCoords(kNSD);
		if (! fCurrentElement -> fShapes -> ParentDomain().MapToParentDomain(fCurrentElement->InitialCoordinates(),
				fSurfaceIPCoords[i], localCoords))
		{
			cout << "LinearBandT::StrainAtCoord, failed to map to local coordinates\n " << flush;
			throw ExceptionT::kGeneralFail;
		}
	
		/*evaluate shape functions */
		dMatrixT grad_U_local(kNSD); //grad_U_global(kNSD); 
		dArray2DT DNa; //DNa - shape function gradients wrt local coords
		dArrayT Na;	// Na - shape functions
	
		fCurrentElement -> fShapes -> 
			GradU(fCurrentElement -> Displacements(), grad_U_local, fSurfaceIPCoords[i], Na, DNa); 
		
		/* evaluate stress */
	    dSymMatrixT bandIPStress(kNSD);
		bandIPStress = 0.0;
		int numNodes = fCurrentElement -> ElementSupport().NumNodes();
		
		for (int j = 0; j < numNodes; j++)
			for (int k = 0; k < dSymMatrixT::NumValues(kNSD); k++)
				bandIPStress [k] += Na [j] * nodalStressList[j] [k];	
	  
		  
		/* calculate and record tractions */
		fSurfaceIPTractions[i].Dimension(kNSD);
		bandIPStress.Multx(fNormal, fSurfaceIPTractions[i]);
	}
	
	/* determine residual cohesion at surface IPs */
	fResidualCohesion.Dimension(kNumSurfaceIPs);
	
	for (int i = 0; i < kNumSurfaceIPs; i++)
	{
		double normalTraction = dArrayT::Dot(fSurfaceIPTractions[i], fNormal);
		double shearTraction = dArrayT::Dot(fSurfaceIPTractions[i], fPerpSlipDir[i]);
		
		if (shearTraction < 0.0)
		{
			FlipSlipDir(i);
			shearTraction *= -1.0;
		}
		
		if (normalTraction < 0.0)
			fResidualCohesion[i] = shearTraction +
			 (fCurrentElement -> fLocalizedFrictionCoeff) * normalTraction;
		else
			fResidualCohesion[i] = shearTraction;
			
	}
	
}

ArrayT<dSymMatrixT> LinearBandT::NodalStressList()
{
	int numNodes = fCurrentElement -> ElementSupport().NumNodes();
	ArrayT<dSymMatrixT> nodalStressList(numNodes);
	
	for(int j = 0; j < numNodes; j++)
		(nodalStressList [j]).Dimension(kNSD);
	
	int numIP = fStress_List.Length();
	int stressLength = dSymMatrixT::NumValues(kNSD);
	dArrayT stressComponentAtIP(numIP), stressComponentAtNode(numNodes);
	
	for (int i = 0; i < stressLength; i++)
	{
	    stressComponentAtIP = 0.0;
		for (int j = 0; j < numIP; j ++)
			stressComponentAtIP[j] = fStress_List[i] [j];
		
		/*extrapolate*/
		fCurrentElement->IP_ExtrapolateAll(stressComponentAtIP, stressComponentAtNode);
		
		/* put into nodal array */
		for (int j = 0; j < numNodes; j ++)
			nodalStressList [i] [j] = stressComponentAtNode [j];
	}

	return nodalStressList;
}

void LinearBandT::SetEndPoints(dArrayT& coord)
{

  
	
}
