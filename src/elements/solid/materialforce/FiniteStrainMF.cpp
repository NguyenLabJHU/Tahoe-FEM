/* $Id: FiniteStrainMF.cpp,v 1.3 2003-03-21 06:30:50 thao Exp $ */
#include "FiniteStrainMF.h"

#include "OutputSetT.h"
#include "ScheduleT.h"
#include "Traction_CardT.h"
#include "GeometryT.h"
#include "ModelManagerT.h"
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

FiniteStrainMF::~FiniteStrainMF(void) 
{
    delete fOutputSet;
}

/* register self for output */
void FiniteStrainMF::RegisterOutput(void)
{
	/* inherited */
	SolidElementT::RegisterOutput();

	ArrayT<StringT> n_labels(3*NumSD());
	ArrayT<StringT> e_labels;
	
	StringT mf_label = "mF";
	StringT mfd_label = "mF_dissip";
	StringT d_label = "D";
	
	const char* suffix[3] = {"_X", "_Y", "_Z"};
	int dex = 0;
	for (int i = 0; i < NumSD(); i++)
		n_labels[dex++].Append(mf_label, suffix[i]);
	for (int i = 0; i < NumSD(); i++)
		n_labels[dex++].Append(mfd_label, suffix[i]);
	for (int i = 0; i < NumSD(); i++)
		n_labels[dex++].Append(d_label, suffix[i]);

	/* collect ID's of the element blocks in the group */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* set output specifier */
	fOutputSet = new OutputSetT(GeometryCode(), block_ID, fConnectivities, n_labels, e_labels, false);
		
	/* register and get output ID */
	fMatForceOutputID = ElementSupport().RegisterOutput(*fOutputSet);
}

/* send output */
void FiniteStrainMF::WriteOutput(void)
{
	/* inherited */
	SolidElementT::WriteOutput();

	/* calculate output values */
	dArray2DT n_values; 
	dArray2DT e_values;

	MapOutput();
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

  output.Dimension(nnd,3*nsd);
  dArrayT mat_force(nmf);
  dArrayT mat_fdissip(nmf);
  mat_force = 0;
  mat_fdissip = 0;
  
  const dArray2DT& disp = Field()[0];
  if (disp.MajorDim() != output.MajorDim()) throw ExceptionT::kGeneralFail;

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

  /*assemble material forces and displacements into output array*/
  double* pout_force = output.Pointer();
  double* pout_dissip = output.Pointer(nsd);
  double* pout_disp = output.Pointer(2*nsd);
  double* pmat_force = mat_force.Pointer();
  double* pmat_fdissip = mat_fdissip.Pointer();
  const iArray2DT& eqno = Field().Equations();
  for (int i = 0; i<nnd; i++)
  {
    for (int j = 0; j<nsd; j++)
    {
      if (eqno[i*nsd+j]<1)
      {
	*pout_force++ = 0.0;
	*pout_dissip++ = 0.0;
	pmat_force++;
	pmat_fdissip++;
      }
      else
      {
	*pout_force++ = (*pmat_force++) + (*pmat_fdissip);
	*pout_dissip++ = (*pmat_fdissip++);
      }
      *pout_disp++ = disp[i*nsd+j];
    }
    pout_force += 2*nsd;
    pout_dissip += 2*nsd;
    pout_disp += 2*nsd;
  }
}

void FiniteStrainMF::AssembleMatForce(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& elem_nodes)
{
  int nen = elem_nodes.Length();
  int nsd = NumSD();
  if (elem_val.Length() != nen*nsd) throw ExceptionT::kGeneralFail;
  int index;
  double* p = global_val.Pointer();
  for(int i = 0; i < nen; i++)
  {
    index=fMap[elem_nodes[i]]*nsd;
    
    for (int j = 0; j<nsd; j++)
    {
    	index += j;
      	p[index] += elem_val[i*nsd+j];
    }
  }
}

void FiniteStrainMF::MapOutput(void)
{
  const iArrayT& nodes_used = fOutputSet->NodesUsed();
  fNumGroupNodes = nodes_used.Length();

  /*find maximum node number)]*/
  int max=0;
  const int* pnode = nodes_used.Pointer(); 
  for (int i = 0; i<fNumGroupNodes; i++)
  {
     if (*pnode > max) max = *pnode;
     pnode++;
  }

  /*map ordering*/
  fMap.Dimension(max+1);
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
  dMatrixT nEshelby(nsd);
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

  /*********************print check***************
   *  bool print = false;
   *  if (fElementCards.Position() == 635) print = true;
   *******************************************************/

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
    nEshelby = 0;
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
      /*form negative of Eshelby stress -SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
      nEshelby[0] = C[0]*S[0]+C[2]*S[2]-energy; //(0,0) 
      nEshelby[1] = C[2]*S[0]+C[1]*S[2]; //(1,0)
      nEshelby[2] = C[0]*S[2]+C[2]*S[1]; //(0,1)
      nEshelby[3] = C[2]*S[2]+C[1]*S[1]-energy; //(1,1)

      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);

      for (int j = 0; j<nen; j++)
      {
	    /*add nEshelby volume integral contribution*/
	    *(pforce++) += ( nEshelby[0]*(*pDQaX)+nEshelby[2]*(*pDQaY)
		        + (F[0]*ip_body[0]+F[1]*ip_body[1])*(*pQa) )*(*jac)*(*weight);
	    *(pforce++) += ( nEshelby[1]*(*pDQaX++)+nEshelby[3]*(*pDQaY++)
		        + (F[2]*ip_body[0]+F[3]*ip_body[1])*(*pQa++) )*(*jac)*(*weight);
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
      
      nEshelby[0] = C[0]*S[0]+C[5]*S[5]+C[4]*S[4]-energy; //(0,0)
      nEshelby[4] = C[5]*S[5]+C[1]*S[1]+C[3]*S[3]-energy; //(1,1)
      nEshelby[8] = C[4]*S[4]+C[3]*S[3]+C[2]*S[2]-energy; //(2,2)
      nEshelby[7] = C[5]*S[4]+C[1]*S[3]+C[3]*S[2]; //(1,2)
      nEshelby[5] = C[4]*S[5]+C[3]*S[1]+C[2]*S[3]; //(2,1)
      nEshelby[6] = C[0]*S[4]+C[5]*S[3]+C[4]*S[2]; //(0,2)
      nEshelby[2] = C[4]*S[0]+C[3]*S[5]+C[2]*S[4]; //(2,0)
      nEshelby[3] = C[0]*S[5]+C[5]*S[1]*C[4]*S[3]; //(0,1)
      nEshelby[1] = C[5]*S[0]+C[1]*S[5]+C[3]*S[4]; //(1,0)
      
      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      double* pDQaZ = DQa(2);
      
      for (int j = 0; j<nen; j++)
      {
	    /*add Eshelby volume integral contribution*/
	    *(pforce++) += (nEshelby[0]*(*pDQaX)+nEshelby[3]*(*pDQaY)+nEshelby[6]*(*pDQaZ)
	            +(F[0]*ip_body[0]+F[1]*ip_body[1]+F[2]*ip_body[2])*(*pQa) )*(*jac)*(*weight);
	    *(pforce++) += (nEshelby[1]*(*pDQaX)+nEshelby[4]*(*pDQaY)+nEshelby[7]*(*pDQaZ)
	            +(F[3]*ip_body[0]+F[4]*ip_body[1]+F[5]*ip_body[2])*(*pQa) )*(*jac)*(*weight);
	    *(pforce++) += (nEshelby[2]*(*pDQaX++)+nEshelby[5]*(*pDQaY++)+nEshelby[8]*(*pDQaZ++)
	            +(F[6]*ip_body[0]+F[7]*ip_body[1]+F[8]*ip_body[2])*(*pQa++))*(*jac)*(*weight);
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
  int numstress = nsd*(nsd+1)/2;
  int nen = NumElementNodes();
  
  /*initialize workspaces*/
  dMatrixT ExtrapMatrix(nen);
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
  /*initialize workspace*/
  fShapes->TopIP();
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& ip_val = IP_iInStretch[i];
    ip_val = 0;
  }
  while(fShapes->NextIP())
  {
    dSymMatrixT InStretch(nsd, pInStretch);
    dSymMatrixT iInStretch = InStretch;
    iInStretch.Inverse();

    const double* pQbU = fShapes->IPShapeU();
    for (int i=0; i<nen; i++)
    {
      const double* pQaU = fShapes->IPShapeU();
      for (int j = 0; j<nen; j++)
      {
	    //Note: The extrapolation matrix is symmetric
	    ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
      }
      
      dSymMatrixT& ip_val = IP_iInStretch[i];
      for (int cnt = 0; cnt < numstress; cnt++)
    	  ip_val[cnt] += (*pQbU)*(*jac)*(*weight)*iInStretch[cnt];
      pQbU++;  
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
    {
        dSymMatrixT& ip_val = IP_iInStretch[j];
        double M = ExtrapMatrix(i,j);
        for (int cnt = 0; cnt < numstress; cnt++)
            nodal_val[cnt] += M*ip_val[cnt];
	}
  }
  
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
      dSymMatrixT InStress(nsd, pInStress);
      
      /*integrate material force*/
      double* pelem_val = elem_val.Pointer();
      for (int i = 0; i<nen; i++)
      {
	    (*pelem_val++) -=0.5*(InStress[0]*pGradX[0]+InStress[1]*pGradX[1]
	        +2.0* InStress[2]*pGradX[2])*(*pQa)*(*jac)*(*weight);
	    (*pelem_val++) -=0.5*(InStress[0]*pGradY[0]+InStress[1]*pGradY[1]
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
      dSymMatrixT InStress(nsd, pInStress);
      
      /*integrate material force*/
      double* pelem_val = elem_val.Pointer();
      for (int i = 0; i<nen; i++)
      {
	    (*pelem_val++) -= 0.5*(InStress[0]*pGradX[0]+InStress[1]*pGradX[1]
	            +InStress[2]*pGradX[2]+2.0*InStress[3]*pGradX[3]
	            +2.0*InStress[4]*pGradX[4]+2.0*InStress[5]*pGradX[5])*(*pQa)*(*jac)*(*weight);
	    (*pelem_val++) -= 0.5*(InStress[0]*pGradY[0]+InStress[1]*pGradY[1]
	            +InStress[2]*pGradY[2]+2.0*InStress[3]*pGradY[3]
	            +2.0*InStress[4]*pGradY[4]+2.0*InStress[5]*pGradY[5])*(*pQa)*(*jac)*(*weight);
	    (*pelem_val++) -= 0.5*(InStress[0]*pGradZ[0]+InStress[1]*pGradZ[1]
	            +InStress[2]*pGradZ[2]+2.0*InStress[3]*pGradZ[3]
	            +2.0*InStress[4]*pGradZ[4]+2.0*InStress[5]*pGradZ[5])*(*pQa++)*(*jac)*(*weight);
      }
    }
    jac++;
    weight++;
  }
}

void FiniteStrainMF::MatForceSurfMech(dArrayT& global_val)
{
  if (fTractionList.Length() > 0)
  {
    
    int nen = NumElementNodes();
    int nsd = NumSD();
    int numval = nsd*nsd;
    
    /*initialize work space*/
    dMatrixT ExtrapMatrix(nen);
    dMatrixT jacobian(nsd,nsd-1);
    dMatrixT Q(nsd);
  
    /*loop through traction card*/
    for (int k = 0; k<fTractionList.Length(); k++)
    {
      const Traction_CardT& BC_card = fTractionList[k];
      const iArrayT& surf_nodes = BC_card.Nodes();
      int nfn = surf_nodes.Length();
      int elem, facet;
      BC_card.Destination(elem, facet);

      /*get traction*/
      LocalArrayT tract(LocalArrayT::kUnspecified, nfn, nsd);
      BC_card.CurrentValue(tract);      

      /*material force vector for facet*/
      dArrayT elem_val(nfn*nsd);

      /*retrieve nodal coordinates of facets in local ordering*/
      LocalArrayT surf_coords(LocalArrayT::kInitCoords, nfn, nsd); 
      ElementSupport().RegisterCoordinates(surf_coords);
      surf_coords.SetLocal(surf_nodes);    
    
      /*retrieve element information*/
      const ElementCardT& elem_card = fElementCards[elem];
      fElementCards.Current(elem);
      SetGlobalShape();
      ContinuumMaterialT* pmat = (*fMaterialList)[elem_card.MaterialNumber()];
      fCurrFSMat = dynamic_cast<FSSolidMatT*>(pmat);
      if (!fCurrFSMat) ExceptionT::GeneralFail();
 
      iArrayT loc_surf_nodes(nfn);
      fShapes->NodesOnFacet(facet,loc_surf_nodes);

      /*Set element shapefunctions*/
      fShapes->SetDerivatives(); 
      const double* jac = fShapes->IPDets();
      const double* weight = fShapes->IPWeights();

      /*project elem ip vals to elem nodal vals*/
      nArrayT<dMatrixT> IP_F(nen);
      for (int i = 0; i<nen; i++)
      {
        dMatrixT& pval = IP_F[i];
        pval.Allocate(nsd);
        pval = 0;
      }
      dArrayT IP_Energy(nen);
      IP_Energy = 0;
      ExtrapMatrix = 0;      
      fShapes->TopIP();
      while(fShapes->NextIP())
      {
        const dMatrixT& F = fCurrFSMat->F_mechanical();
        double energy = fCurrFSMat->StrainEnergyDensity();
        const double* pQbU = fShapes->IPShapeU();
        for (int i = 0; i<nen; i++)
        {
          const double* pQaU = fShapes->IPShapeU();
          for (int j = 0; j<nen; j++)
	        ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
          dMatrixT& pval = IP_F[i];
          for (int cnt = 0; cnt <numval; cnt++)
            pval[cnt] += (*pQbU)*(*jac)*(*weight)*F[cnt];
          IP_Energy[i] += (*pQbU++)*(*jac)*(*weight)*energy;
        } 
        jac++;
        weight++;
      }
      ExtrapMatrix.Inverse();
      nArrayT<dMatrixT> Nodal_F(nfn);
      dArrayT Nodal_Energy(nfn);
      for (int i = 0; i<nfn; i++)
      {
        int loc_node=loc_surf_nodes[i];
        dMatrixT& nval = Nodal_F[i];
        nval.Allocate(nsd);
        nval = 0;
        Nodal_Energy[i] = 0;
        for (int j = 0; j<nen; j++)
        {
          dMatrixT& pval = IP_F[j];
          double M = ExtrapMatrix(loc_node, j);
          for (int cnt = 0; cnt<numval; cnt++)
            nval[cnt] += M*pval[cnt];
          Nodal_Energy[i] += M*IP_Energy[j];
        }
      }
      
      /*get surface shape function*/
      const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
      int nip = surf_shape.NumIP();
      dArrayT ip_tract(nsd);
      dArrayT global_tract(nsd);
      dArrayT ip_eshelby(nsd);     
      dMatrixT ip_F(nsd);
      if (nsd == 2)
      {
        elem_val = 0.0;
        double thickness = 1.0;
        const double* ip_w = surf_shape.Weight();
        for (int j = 0; j<nip; j++)
        {
		  surf_shape.DomainJacobian(surf_coords,j,jacobian);
	  	  double detj = surf_shape.SurfaceJacobian(jacobian, Q);

          /*interpolate to surface ip*/
          double* ptract_X = tract(0);
          double* ptract_Y = tract(1);
          ip_F = 0;
          ip_tract = 0;
          double ip_energy = 0;
          const double* pQaU = surf_shape.Shape(j);
	      for (int i= 0; i<nfn; i++)
	      {
	        dMatrixT& nval = Nodal_F[i];
	        ip_F[0] += *pQaU*nval[0];
	        ip_F[1] += *pQaU*nval[1];
	        ip_F[2] += *pQaU*nval[2];
	        ip_F[3] += *pQaU*nval[3];
	        ip_energy +=  *pQaU*Nodal_Energy[i];
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
	      ip_eshelby[0] = (ip_energy*n[0]-ip_F[0]*global_tract[0]
	            -ip_F[1]*global_tract[1])*thickness;
	      ip_eshelby[1] = (ip_energy*n[1]-ip_F[2]*global_tract[0]
	            -ip_F[3]*global_tract[1])*thickness;
	      /*integrate material force*/
	      const double* pQa = surf_shape.Shape(j);
	      double * pelem_val = elem_val.Pointer();
	      for (int i = 0; i<nfn; i++)
	      {
	        (*pelem_val++) +=(*pQa)*ip_eshelby[0]*detj*(*ip_w);
	        (*pelem_val++) +=(*pQa++)*ip_eshelby[1]*detj*(*ip_w);
	      }
	      ip_w++;
        }
      }
      else if (nsd == 3)
      { 
        elem_val = 0.0;
        const double* ip_w = surf_shape.Weight();
        for (int j = 0; j<nip; j++)
        {
		  surf_shape.DomainJacobian(surf_coords,j,jacobian);
	  	  double detj = surf_shape.SurfaceJacobian(jacobian, Q);

          /*interpolate to surface ip*/
          double* ptract_X = tract(0);
          double* ptract_Y = tract(1);
          double* ptract_Z = tract(2);
          ip_F = 0;
          ip_tract = 0;
          double ip_energy = 0;
          const double* pQaU = surf_shape.Shape(j);
	      for (int i= 0; i<nfn; i++)
	      {
	        dMatrixT& nval = Nodal_F[i];
	        ip_F[0] += *pQaU*nval[0];
	        ip_F[1] += *pQaU*nval[1];
	        ip_F[2] += *pQaU*nval[2];
	        ip_F[3] += *pQaU*nval[3];
	        ip_F[4] += *pQaU*nval[4];
	        ip_F[5] += *pQaU*nval[5];
	        ip_F[6] += *pQaU*nval[6];
	        ip_F[7] += *pQaU*nval[7];
	        ip_F[8] += *pQaU*nval[8];	        
	        ip_energy +=  *pQaU*Nodal_Energy[i];
	        ip_tract[0] += (*pQaU)*(*ptract_X++);
	        ip_tract[1] += (*pQaU)*(*ptract_Y++);
	        ip_tract[2] += (*pQaU++)*(*ptract_Z++);
	      }
	      /*surface normal*/
	      double* n = Q(2);
	      if (BC_card.CoordSystem() == Traction_CardT::kLocal)
	      {
	        /*rotate traction from local to global coords*/
	        global_tract[0] = Q[0]*ip_tract[0]+Q[3]*ip_tract[1]+Q[6]*ip_tract[2];
	        global_tract[1] = Q[1]*ip_tract[0]+Q[4]*ip_tract[1]+Q[7]*ip_tract[2];
	        global_tract[2] = Q[2]*ip_tract[0]+Q[5]*ip_tract[1]+Q[8]*ip_tract[2];
     	  }
	      else if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	      {
	        global_tract[0] = ip_tract[0];
	        global_tract[1] = ip_tract[1];
	        global_tract[2] = ip_tract[2];
	      }
	      ip_eshelby[0] = (ip_energy*n[0]-ip_F[0]*global_tract[0]
	            -ip_F[1]*global_tract[1]-ip_F[2]*global_tract[2]);
	      ip_eshelby[1] = (ip_energy*n[1]-ip_F[3]*global_tract[0]
	            +ip_F[4]*global_tract[1]-ip_F[5]*global_tract[2]);
	      ip_eshelby[2] = (ip_energy*n[2]-ip_F[6]*global_tract[0]
	            +ip_F[7]*global_tract[1]-ip_F[8]*global_tract[2]);
	      /*integrate material force*/
	      const double* pQa = surf_shape.Shape(j);
	      double * pelem_val = elem_val.Pointer();
	      for (int i = 0; i<nfn; i++)
	      {
	        (*pelem_val++) +=(*pQa)*ip_eshelby[0]*detj*(*ip_w);
	        (*pelem_val++) +=(*pQa)*ip_eshelby[1]*detj*(*ip_w);
	        (*pelem_val++) +=(*pQa++)*ip_eshelby[2]*detj*(*ip_w);
	      }
	      ip_w++;
        }
      }
      AssembleMatForce(elem_val, global_val, surf_nodes);
    }
  }
}
