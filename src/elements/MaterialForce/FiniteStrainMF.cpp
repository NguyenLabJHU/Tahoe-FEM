/* $Id: FiniteStrainMF.cpp,v 1.1 2003-02-12 18:37:41 thao Exp $ */
#include "FiniteStrainMF.h"

#include "ScheduleT.h"
#include "Traction_CardT.h"
#include "ShapeFunctionT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"

/* materials lists */
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"

using namespace Tahoe;

/* constructor */
FiniteStrainMF::FiniteStrainMF(const ElementSupportT& support, const FieldT& field):
	FiniteStrainT(support, field),
	fNumGroupNodes(0),
	fMatForceOutputID(-1){}

/*destructor*/
FiniteStrainMF::~FiniteStrainMF(void) { delete fMap };

/* register self for output */
void FiniteStrainMF::RegisterOutput(void)
{
	/* inherited */
	SolidElementT::RegisterOutput();

	ArrayT<StringT> n_labels(2*NumSD());
	ArrayT<StringT> e_labels;
	
	StringT mf_label = "mF";
	StringT mfd_label = "mF_dissip";
	const char suffix[] = {"_X", "_Y", "_Z"};
	dex = 0;
	for (int i = 0; i < NumSD(); i++)
		n_labels[dex++].Append(mf_label, suffix[i]);
	for (int i = 0; i < NumSD(); i++)
		n_labels[dex++].Append(mfd_label, suffix[i]);

	/* collect ID's of the element blocks in the group */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* set output specifier */
	OutputSetT fOutputSet(fGeometryCode, block_ID, fConnectivities, n_labels, e_labels, false);
		
	/* register and get output ID */
	fMatForceOutputID = ElementSupport().RegisterOutput(fOutputSet);
}

/* send output */
void FiniteStrainMF::WriteOutput(void)
{
	/* inherited */
	SolidElementT::WriteOutput();

	/* calculate output values */
	dArray2DT n_values; 
	dArray2DT e_values;

	MapOuput();
	ComputeMatForce(n_values);
	
	/* send to output */
	ElementSupport().WriteOutput(fMatForceOutputID, n_values, e_values);
}

void FiniteStrainMF::ComputeMatForce(dArray2DT& output)
{
  const char caller[] = "FiniteStrainT::ComputeMatForce";

  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nen = NumElementNodes();
  int nsd = NumSD();
  int nmf = nnd*nsd;

  output.Dimesion(2,nmf);
  output = 0;
  double* pmat_force = output(0);
  double* pmat_fdissip = output(1);
  dArrayT mat_force(nmf,pmat_force);
  dArrayT mat_fdissip(nmf,pmat_fdissip);
  
  dArrayT elem_data(nsd*nen);
   
  Top();
  while (NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrFSMat = dynamic_cast<FSSolidMatT*>(pmat);
    if (!fCurrFSMat) ExceptionT::GeneralFail(caller);
    
    /*Set Global Shape Functions for current element*/
    SetGlobalShape();

    MatForceVolMech(elem_data);
    AssembleMatForce(elem_data, mat_force, CurrentElement().NodesX());
    if (fCurrFSMat->HasDissipVar()) 
    {
      MatForceDissip(elem_data, CurrentElement().DoubleData());
      AssembleMatForce(elem_data, mat_fdissip, CurrentElement().NodesX());
    }
  }
  MatForceSurfMech(mat_force);
  mat_force += mat_fdissip;
}

void FiniteStrainMF::AssembleMatForce(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& nodes)
{
  int nnd = nodes.Length();
  int nsd = NumSD();
  if (nnd*nsd != elem_val.Length()) 
    ExceptionT::GeneralFail('FiniteStrainMF::AssembleMatForce');

  int pos;
  for(int i = 0; i < nnd; i++)
  {
    pos = fMap[nodes[i]];
    double* p = global_val.Pointer()+(pos*nsd);
    for (int j = 0; j<nsd; j++)
      p[j] += elem_val[i*nsd+j];
  }
}

void FiniteStrainMF::MapOutput(void);
{
  const iArrayT& nodes_used = fOutputSet.NodesUsed();
  fNumGroupNodes = nodes_used.Length();

  /*find maximum node number)]*/
  int max=0;
  const int* node = nodes_used.Pointer(); 
  for (int i = 0; i<fNumGroupNodes; i++)
  {
     if (*node > max) max = *nodes;
     nodes++;
  }

  /*map ordering*/
  fMap.Dimension(max);
  for (int i = 0; i<fNumGroupNodes; i++)
  {
    fMap[nodes_used[i]]=i;
  }
}

void FiniteStrainMF::MatForceVolMech(dArrayT& elem_val)
{
  /*obtain dimensions*/
  int nsd = NumSD();
  int nen = NumElementNodes();

  /*initialize workspaces*/
  //element vectors: F=[F1x|F1y||F2x|F2y||...]
  dArrayT bodyforce(nsd*nen);
  dMatrixT Eshelby(nsd);
  dSymMatrixT C(nsd);
  dArrayT ip_body(nsd);
 
  elem_val = 0;
  bodyforce = 0.0;
  
  /*get density*/
  double density = fCurrFSMat->Density();
  
  if (fLocAcc.IsRegistered())
    SetLocalU(fLocAcc);
  else
  	fLocAcc = 0.0;

  /*copy acceleration and body force data into body force vector*/
  double* pbody = bodyforce.Pointer();    
  if (fBodySchedule)
  {
    double loadfactor = fBodySchedule->Value();
    for (int i = 0; i<nen; i++)
    {
      double* pLocAcc = fLocAcc.Pointer();
      pLocAcc += i;
      for (int j = 0; j<nsd; j++)
      {
	*pbody++ = (*pLocAcc - fBody[j]*loadfactor)*density;
	pLocAcc += nen;
      }
    }
  }
  else 
  {
    for (int i = 0; i<nen; i++)
    {
      double* pLocAcc = fLocAcc.Pointer();
      pLocAcc += i;
      for (int j = 0; j<nsd; j++)
      {
	*pbody++ = (*pLocAcc)*density;
	pLocAcc += nen;
      }
    }
  }
    
  /*intialize shape function data*/
  const double* jac = fShapes->IPDets();
  const double* weight = fShapes->IPWeights();
  fShapes->TopIP();
  while(fShapes->NextIP())
  {
    /*get shape function and derivatives at integration point*/
    const double* pQa = fShapes->IPShapeX();
    const double* pQaU = fShapes->IPShapeU();
    const dArray2DT& DQa = fShapes->Derivatives_X();
    
    /*gather material ip data*/
    double energy = fCurrFSMat->StrainEnergyDensity();
    const dMatrixT& F = fCurrFSMat->F_mechanical();
    Eshelby = 0;
    const dSymMatrixT& S = fCurrFSMat->S_IJ();
    double* pbody = bodyforce.Pointer();
    double* pforce = elem_val.Pointer(); 
    if (nsd == 2)
    {
      /*interpolate to ip*/
      ip_body = 0;
      for (int i= 0; i<nen; i++)
      {
	ip_body[0] += (*pQaU) * (*pbody++);
	ip_body[1] += (*pQaU++) * (*pbody++);
      }	 
      /*form right stretch tensor*/
      C[0] = F[0]*F[0] + F[1]*F[1];
      C[1] = F[2]*F[2] + F[3]*F[3];
      C[2] = F[0]*F[2] + F[1]*F[3];
      /*form negative of Eshelby stress SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
      Eshelby[0] = C[0]*S[0]+C[2]*S[2]-energy; //(0,0) 
      Eshelby[1] = C[2]*S[0]+C[1]*S[2]; //(1,0)
      Eshelby[2] = C[0]*S[2]+C[2]*S[1]; //(0,1)
      Eshelby[3] = C[2]*S[2]+C[1]*S[1]-energy; //(1,1)
      
      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      
      for (int j = 0; j<nen; j++)
      {
	/*add Eshelby volume integral contribution*/
	*(pforce++) += ( Eshelby[0]*(*pDQaX)+Eshelby[2]*(*pDQaY)
		     + (F[0]*ip_body[0]+F[1]*ip_body[1])*(*pQa) )
	             *(*jac)*(*weight);
	*(pforce++) += ( Eshelby[1]*(*pDQaX++)+Eshelby[3]*(*pDQaY++)
		     + (F[2]*ip_body[0]+F[3]*ip_body[1])*(*pQa++) )
	             *(*jac)*(*weight);
      }
    }
    else if (nsd ==3)
    {
      /*interpolate ip values*/
      ip_body = 0;
      for (int i= 0; i<nen; i++)
      {
	ip_body[0] += (*pQaU) * (*pbody++);
	ip_body[1] += (*pQaU) * (*pbody++);
	ip_body[2] += (*pQaU++) * (*pbody++);
      }
	  
      /*form Eshelby stress SIG_IJ = Psi I_IJ - C_IK S_KJ*/
      C[0] = F[0]*F[0]+F[1]*F[1]+F[2]*F[2];
      C[1] = F[3]*F[3]+F[4]*F[4]+F[5]*F[5];
      C[2] = F[6]*F[6]+F[7]*F[7]+F[8]*F[8];
      C[3] = F[3]*F[6]+F[4]*F[7]+F[5]*F[8];
      C[4] = F[0]*F[6]+F[1]*F[7]+F[2]*F[8];
      C[5] = F[0]*F[3]+F[1]*F[4]+F[2]*F[5];
      
      Eshelby[0] = C[0]*S[0]+C[5]*S[5]+C[4]*S[4]-energy; //(0,0)
      Eshelby[4] = C[5]*S[5]+C[1]*S[1]+C[3]*S[3]-energy; //(1,1)
      Eshelby[8] = C[4]*S[4]+C[3]*S[3]*C[2]*S[2]-energy; //(2,2)
      Eshelby[7] = C[5]*S[4]+C[1]*S[3]+C[3]*S[2]; //(1,2)
      Eshelby[5] = C[4]*S[5]+C[3]*S[1]+C[2]*S[3]; //(2,1)
      Eshelby[6] = C[0]*S[4]+C[5]*S[3]+C[4]*S[2]; //(0,2)
      Eshelby[2] = C[4]*S[0]+C[3]*S[5]+C[2]*S[4]; //(2,0)
      Eshelby[3] = C[0]*S[5]+C[5]*S[1]*C[4]*S[3]; //(0,1)
      Eshelby[1] = C[5]*S[0]+C[1]*S[5]+C[3]*S[4]; //(1,0)
      
      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      double* pDQaZ = DQa(2);
      
      for (int j = 0; j<nen; j++)
      {
	/*add Eshelby volume integral contribution*/
	*(pforce++) += (Eshelby[0]*(*pDQaX)+Eshelby[3]*(*pDQaY)
	             +  Eshelby[6]*(*pDQaZ)
                     + (F[0]*ip_body[0]+F[1]*ip_body[1]
                     +  F[2]*ip_body[2])*(*pQa) )*(*jac)*(*weight);
	*(pforce++) += (Eshelby[1]*(*pDQaX)+Eshelby[4]*(*pDQaY)
                     +  Eshelby[7]*(*pDQaZ)
                     + (F[3]*ip_body[0]+F[4]*ip_body[1]
		     +	F[5]*ip_body[2])*(*pQa) )*(*jac)*(*weight);
	*(pforce++) += (Eshelby[2]*(*pDQaX++)+Eshelby[5]*(*pDQaY++)
		     +  Eshelby[8]*(*pDQaZ++)
		     + (F[6]*ip_body[0]+F[7]*ip_body[1]
                     +  F[8]*ip_body[2])*(*pQa++))*(*jac)*(*weight);
      }
    }
    weight++;
    jac++;
  }
}

void FiniteStrainMF::MatForceDissip(dArrayT& elem_val, const dArrayT& statev)
{
  /*obtain dimensions*/
  int nsd = NumSD();
  int numstress = nsd*(nsd+1)*0.5;
  int nen = NumElementNodes();
  
  /*initialize workspaces*/
  dSymMatrixT ExtrapMatrix(nen);
  nArrayT<dSymMatrixT> Nodal_iInStretch(nen);
  nArrayT<dSymMatrixT> IP_iInStretch(nen);
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& ip_val = IP_iInStretch[i];
    dSymMatrixT& nodal_val = Nodal_iInStretch[i];
    ip_val.Allocate(nsd);
    nodal_val.Allocate(nsd);
  }
  dArrayT Grad_iInStretch(nsd*numstress);

  elem_val = 0;
    
  /*get shape function data*/
  int nip = fShapes->NumIP();
  const double* jac = fShapes->IPDets();;
  const double* weight = fShapes->IPWeights();

  /*set inelastic stretch to current values*/
  int nstatev = statev.Length()/nip;
  double* pstatev = statev.Pointer();
  double* pInStretch = pstatev;
  pstatev += numstress;
  /*skip over previous value of inelastic stretch
    and set inelastic stress measure*/
  pstatev += numstress;
  double* pInStress = pstatev;

  /*extrapolate iInStretch to nodes*/
  ExtrapMatrix = 0;
  /*get shape function and derivatives at integration point*/
  fShapes->TopIP();
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& ip_val = IP_iInStretch[i];
    ip_val = 0;
  }
  while(fShapes->NextIP())
  {
    dSymMatrixT InStretch(numstress, pInStretch);
    dSymMatrixT iInStretch = InStretch;
    iInStretch.Inverse();

    const double* pQbU = fShapes->IPShapeU();
    for (int i=0; i<nen; i++)
    {
      const double* pQaU = fShapes->IPShapeU();
      for (int j = i; j<nen; i++)
      {
	//Note: The extrapolation matrix is symmetric
	ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
      }
      
      dSymMatrixT& ip_val = IP_iInStretch[i];
      ip_val.AddCombination((*pQbU++)*(*jac)*(*weight), iInStretch);
    }
    pInStretch += nstatev;
    weight++;
    jac++;
  }
  ExtrapMatrix.Inverse();
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& nodal_val = Nodal_iInStretch[i];
    nodal_val = 0;
    for (int j = 0; j<nen; j++)
      nodal_val.AddCombination(ExtrapMatrix[i+j*nen],IP_iInStretch[j]);
  }
  
  double* pelem_val = elem_val.Pointer();
  jac = fShapes->IPDets();
  weight = fShapes->IPWeights();
  fShapes->TopIP();
  while(fShapes->NextIP())
  {
    const double* pQa = fShapes->IPShapeX();
    const dArray2DT& DQa = fShapes->Derivatives_X();
    if (nsd ==2)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);
      Grad_iInStretch = 0;

      /*Interpolate grad of iverse inelastic stretch to ip*/
      double* pGradX = Grad_iInStretch.Pointer();
      double* pGradY = pGradX+numstress;
      for (int i = 0; i<nen; i++)
      {
	dSymMatrixT& nodal_val = Nodal_iInStretch[i];

        pGradX[0] += (*pDQaX)*nodal_val[0];
	pGradX[1] += (*pDQaX)*nodal_val[1];
	pGradX[2] += (*pDQaX++)*nodal_val[2];
	
	pGradY[0] += (*pDQaY)*nodal_val[0];
	pGradY[1] += (*pDQaY)*nodal_val[1];
	pGradY[2] += (*pDQaY++)*nodal_val[2];
      }
      /*get inelastic stress*/
      dSymMatrixT InStress(numstress, pInStress);
      
      /*integrate material force*/
      for (int i = 0; i<nen; i++)
      {
	(*pelem_val++) +=0.5*(InStress[0]*pGradX[0]+InStress[1]*pGradX[1]
			+2.0* InStress[2]*pGradX[2])*(*pQa)*(*jac)*(*weight);
	(*pelem_val++) +=0.5*(InStress[0]*pGradY[0]+InStress[1]*pGradY[1]
		        +2.0* InStress[2]*pGradY[2])*(*pQa++)*(*jac)*(*weight);
      }
    }
    else if (nsd ==3)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);
      const double* pDQaZ = DQa(2);
      Grad_iInStretch = 0;
      /*Interpolate grad of iverse inelastic stretch to ip*/
      double* pGradX = Grad_iInStretch.Pointer();
      double* pGradY = pGradX+numstress;
      double* pGradZ = pGradY+numstress;
      for (int i = 0; i<nen; i++)
      {	
	dSymMatrixT& nodal_val = Nodal_iInStretch[i];

	pGradX[0] += (*pDQaX)*nodal_val[0];
	pGradX[1] += (*pDQaX)*nodal_val[1];
	pGradX[2] += (*pDQaX)*nodal_val[2];
	pGradX[3] += (*pDQaX)*nodal_val[3];
	pGradX[4] += (*pDQaX)*nodal_val[4];
	pGradX[5] += (*pDQaX++)*nodal_val[5];
	
	pGradY[0] += (*pDQaY)*nodal_val[0];
	pGradY[1] += (*pDQaY)*nodal_val[1];
	pGradY[2] += (*pDQaY)*nodal_val[2];
	pGradY[3] += (*pDQaY)*nodal_val[3];
	pGradY[4] += (*pDQaY)*nodal_val[4];
	pGradY[5] += (*pDQaY++)*nodal_val[5];
      }
      /*get inelastic stress*/
      dSymMatrixT InStress(numstress, pInStress);
      
      /*integrate material force*/
      for (int i = 0; i<nen; i++)
      {
	(*pelem_val++) += 0.5*(InStress[0]*pGradX[0]+InStress[1]*pGradX[1]
 		           +   InStress[2]*pGradX[2]+2.0*InStress[3]*pGradX[3]
                          +2.0*InStress[4]*pGradX[4]+2.0*InStress[5]*pGradX[5])
                           *(*pQa)*(*jac)*(*weight);
	(*pelem_val++) += 0.5*(InStress[0]*pGradY[0]+InStress[1]*pGradY[1]
	                   +   InStress[2]*pGradY[2]+2.0*InStress[3]*pGradY[3]
                          +2.0*InStress[4]*pGradY[4]+2.0*InStress[5]*pGradY[5])
                           *(*pQa)*(*jac)*(*weight);
	(*pelem_val++) += 0.5*(InStress[0]*pGradZ[0]+InStress[1]*pGradZ[1]
	                   +   InStress[2]*pGradZ[2]+2.0*InStress[3]*pGradZ[3]
                          +2.0*InStress[4]*pGradZ[4]+2.0*InStress[5]*pGradZ[5])
                           *(*pQa++)*(*jac)*(*weight);
      }
    }
    jac++;
    weight++;
  }
}

void FiniteStrainMF::MatForceSurfMech(dArrayT& global_val)
{
  if(fTractionList.Length() > 0)
  {
    /*obtain dimensions*/
    int nen = NumElementNodes();
    int nsd = NumSD();

    /*initialize workspace*/
    dSymMatrixT ExtrapMatrix(nen);
    nArrayT<dMatrixT> IP_F(nen);
    for (int i = 0; i<nen; i++)
    {
      dMatrixT& ip_val = IP_F[i];
      ip_val.Allocate(nsd);
    }

    dArrayT IP_energy(nen);      

    dMatrixT jacobian(nsd,nsd-1);
    dMatrixT ip_F(nsd);
    dArrayT ip_tract(nsd);
    dArrayT ip_mattract(nsd);
    dArrayT global_tract(nsd);
    double ip_energy;
    dMatrixT Q(nsd);
 
    for (int k = 0; k<fTractionList.Length(); k++)
    {
      /*retrieve ith traction card*/
      const Traction_CardT& BC_card = fTractionList[k];

      /*BC destination*/
      int elem, facet;  //element # and facet # of traction bc
      BC_card.Destination(elem,facet);

      /*dimension of ith traction card*/
      /*global node numbers of facet nodes*/
      const iArrayT& surf_nodes = BC_card.Nodes();
      int nsn = surf_nodes.Length();

      /*get local node numbers of facet nodes in local ordering of the facet*/
      iArrayT loc_surf_nodes(nsn);
      fShapes->NodesOnFacet(facet,loc_surf_nodes);

      /*retrieve coordinates of nodes in local ordering*/
      LocalArrayT surf_coords(LocalArrayT::kInitCoords, nsn, nsd); //nsnxnsd
      ElementSupport().RegisterCoordinates(surf_coords);
      surf_coords.SetLocal(surf_nodes);    

      /*allocate work space*/
      nArrayT<dMatrixT> Nodal_F(nsn);
      for (int i = 0; i<nsn; i++)
      {
		dMatrixT& nodal_val = Nodal_F[i];
		nodal_val.Allocate(nsd);
      }
      dArrayT nodal_energy(nsn);  

      /*initialize mat force vector*/
      dArrayT elem_val(nsn*nsd);
      /*get local nodal tractions scaled by LTf*/
      LocalArrayT tract(LocalArrayT::kUnspecified, nsn, nsd);
      BC_card.CurrentValue(tract);

      /*Extrapolate ip energies and def grad to nodes*/
      const ElementCardT& elem_card = fElementCards[elem];
      ContinuumMaterialT* pmat = (*fMaterialList)[elem_card.MaterialNumber()];
      FSSolidMatT* fCurrFSMat = dynamic_cast<FSSolidMatT*>(pmat);
      if (!fCurrFSMat) ExceptionT::GeneralFail();

      /*get element coords*/
      const iArrayT& elem_nodes=elem_card.NodesX();
      fLocInitCoords.SetLocal(elem_nodes);
      LocalArrayT& elem_coords = fLocInitCoords;

      /*Set element shapefunctions*/
      fShapes->SetDerivatives(); 
      const double* jac = fShapes->IPDets();
      const double* weight = fShapes->IPWeights();

      /*Form Extrapolation Matrix*/
      ExtrapMatrix = 0;
      IP_energy = 0;
      for (int i = 0; i<nen; i++)
      {
	dMatrixT& ip_val = IP_F[i];
	ip_val = 0;
      }
      fShapes->TopIP();
      while (fShapes->NextIP())
      {
	const dMatrixT& F = fCurrFSMat->F_mechanical();
	double energy = fCurrFSMat->StrainEnergyDensity();

	const double* pQbU = fShapes->IPShapeU();
	for (int i=0; i<nen; i++)
	{
	  const double* pQaU = fShapes->IPShapeU();
	  for (int j = i; j<nen; i++)
	  {
	    //Note: The extrapolation matrix is symmetric
	    ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
	  }

	  IP_energy[i] += (*pQbU)*energy*(*jac)*(*weight);
	  dMatrixT& ip_val = IP_F[i];
	  ip_val.AddCombination((*pQbU++)*(*jac)*(*weight), F);
	}
	weight++;
	jac++;
      }
      ExtrapMatrix.Inverse();
      //Q:Assumes that nodes are in local ordering
      for (int i = 0; i<nsn; i++)
      {
	int node = loc_surf_nodes[i];
	dMatrixT& nodal_val = Nodal_F[i];
	nodal_val = 0;
	for (int j = 0; j<nen; j++)
	{
	  nodal_val.AddCombination(ExtrapMatrix[node+j*nen],IP_F[j]);
	  nodal_energy[i] += ExtrapMatrix[node+j*nen]*IP_energy[j];
	}
      }

      /*get surface shape function*/
      const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
      int nip = surf_shape.NumIP();       
      double* penerg = nodal_energy.Pointer();
       
      if (nsd ==2)
      {
	elem_val = 0.0;
	int thickness = 1.0; //Q:Hardwired thickness;
	const double* ip_w = surf_shape.Weight();
	for (int j = 0; j < nip; j++)
	{
	  surf_shape.DomainJacobian(surf_coords,j,jacobian);
	  double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	  
	  /*interpolate to surface ip*/
	  double* ptract_X = tract(0);
	  double* ptract_Y = tract(1);
	  ip_F = 0;
	  ip_tract = 0;
	  ip_energy = 0;
	  
	  const double* pQaU = surf_shape.Shape(j);
	  for (int i= 0; i<nsn; i++)
	  {
	    dMatrixT& nodal_val = Nodal_F[i];
	    ip_F.AddCombination(*pQaU,nodal_val);
	    ip_energy += (*pQaU)*(*penerg++);
	    ip_tract[0] += (*pQaU)*(*ptract_X++);
	    ip_tract[1] += (*pQaU++)*(*ptract_Y++);
	  }
	  /*surface normal*/
	  double* n = Q(1);
	  if (BC_card.CoordSystem() == Traction_CardT::kLocal)
	  {
	    /*rotate traction from local to global coords*/
	    global_tract[0] = Q[0]*ip_tract[0]+Q[2]*ip_tract[1];
	    global_tract[1] = Q[1]*ip_tract[0]+Q[3]*ip_tract[1];
	  }
	  else if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    global_tract[0] = ip_tract[0];
	    global_tract[1] = ip_tract[1];
	  }
	  ip_mattract[0] = (ip_F[0]*global_tract[0]
			    +ip_F[1]*global_tract[1])*thickness-ip_energy*n[0];
	  ip_mattract[1] = (ip_F[2]*global_tract[0]
			    +ip_F[3]*global_tract[1])*thickness-ip_energy*n[1];
	  
	  /*integrate material force*/
	  const double* pQa = surf_shape.Shape(j);
	  double * pelem_val = elem_val.Pointer();
	  for (int i = 0; i<nsn; i++)
	  {
	    (*pelem_val++) +=(*pQa)*ip_mattract[0]*detj*(*ip_w);
	    (*pelem_val++) +=(*pQa++)*ip_mattract[1]*detj*(*ip_w);
	  }
	  ip_w++;
	}
      }
      else if (nsd == 3)
      {
	const double* ip_w = surf_shape.Weight();
	for (int j = 0; j < nip; j++)
	{
	  surf_shape.DomainJacobian(surf_coords,j,jacobian);
	  double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	  
	  /*interpolate to surface ip*/
	  double* ptract_X = tract(0);
	  double* ptract_Y = tract(1);
	  double* ptract_Z = tract(2);
	  ip_F = 0;
	  ip_tract = 0;
	  ip_energy = 0;
	  const double* pQaU = surf_shape.Shape(j);
	  for (int i= 0; i<nsn; i++)
	  {
	    dMatrixT& nodal_val = Nodal_F[i];
	    ip_F.AddCombination(*pQaU,nodal_val);
	    ip_energy += (*pQaU)*(*penerg++);
	    
	    ip_tract[0] += (*pQaU)*(*ptract_X++);
	    ip_tract[1] += (*pQaU)*(*ptract_Y++);
	    ip_tract[2] += (*pQaU++)*(*ptract_Z++);
	  }
	  /*surface normal*/
	  double* n = Q(2);
	  if (BC_card.CoordSystem() == Traction_CardT::kLocal)
	  {
	    /*rotate traction from local to global coords*/
	    global_tract[0]=Q[0]*ip_tract[0]+Q[3]*ip_tract[1]+Q[6]*ip_tract[2];
	    global_tract[1]=Q[1]*ip_tract[0]+Q[4]*ip_tract[1]+Q[7]*ip_tract[2];
	    global_tract[2]=Q[2]*ip_tract[0]+Q[5]*ip_tract[1]+Q[8]*ip_tract[2];
	  }
	  else if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    global_tract[0] = ip_tract[0];
	    global_tract[1] = ip_tract[1];
	    global_tract[2] = ip_tract[2];
	  }
	  
	  ip_mattract[0] = (ip_F[0]*ip_tract[0]+ip_F[1]*ip_tract[1]
			  + ip_F[2]*ip_tract[2]) - ip_energy*n[0];
	  ip_mattract[1] = (ip_F[3]*ip_tract[0]+ip_F[4]*ip_tract[1]
			  + ip_F[5]*ip_tract[2]) - ip_energy*n[1];
	  ip_mattract[1] = (ip_F[6]*ip_tract[0]+ip_F[7]*ip_tract[1]
			  + ip_F[8]*ip_tract[2]) - ip_energy*n[2];

	  /*integrate material force*/
	  const double* pQa = surf_shape.Shape(j);
	  double * pelem_val = elem_val.Pointer();
	  for (int i = 0; i<nsn; i++)
	  {
	    (*pelem_val++) +=(*pQa)*ip_mattract[0]*detj*(*ip_w);
	    (*pelem_val++) +=(*pQa)*ip_mattract[1]*detj*(*ip_w);
	    (*pelem_val++) +=(*pQa++)*ip_mattract[2]*detj*(*ip_w);
	  }
	  ip_w++;
	}
      }
      AssembleMatForce(elem_val, global_val, surf_nodes);
    }
  }
}
