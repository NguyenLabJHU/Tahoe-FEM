/* $Id: SSQ2P1MF.cpp,v 1.7 2004-07-15 08:28:31 paklein Exp $ */
#include "SSQ2P1MF.h"

#include "OutputSetT.h"
#include "Traction_CardT.h"
#include "ScheduleT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"
#include "GeometryT.h"
#include "SSMatSupportT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

/* materials lists */
#include "MaterialListT.h"
#if 0
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"
#endif

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
using namespace Tahoe;

/* constructor */
SSQ2P1MF::SSQ2P1MF(const ElementSupportT& support, const FieldT& field):
  SmallStrainQ2P1(support, field),
  MFSupportT(support) {}
	
void SSQ2P1MF::Initialize(void)
{
ExceptionT::GeneralFail("SSQ2P1MF::Initialize", "out of date");
#if 0	
  SmallStrainQ2P1::Initialize();
  
  int nsd = (NumSD() == 2) ? 4 : 3;
   
  /*dimension workspace*/
  fEshelby.Dimension(NumSD());
  fCauchy.Dimension(nsd);

  fBodyForce.Dimension(NumSD()*NumElementNodes());
  fip_body.Dimension(NumSD());
    
  fGradU_List.Dimension(NumIP());
  for (int i = 0; i< NumIP(); i++)
    fGradU_List[i].Dimension(NumSD());
#endif
}

void SSQ2P1MF::SetGlobalShape(void)
{
  SmallStrainQ2P1::SetGlobalShape();
  for (int i = 0; i < NumIP(); i++)
    fShapes->GradU(fLocDisp, fGradU_List[i],i);
}

/***************************outputs managers***********************************/
/* register self for output */
void SSQ2P1MF::RegisterOutput(void)
{
  /* inherited */
  SolidElementT::RegisterOutput();

  ArrayT<StringT> n_labels(3*NumSD());
  ArrayT<StringT> e_labels;
  
  StringT mf_label = "mF";
  StringT mfd_label = "mF_dissip";
  StringT disp_label = "D";
  const char* suffix[3] = {"_X", "_Y", "_Z"};
  int dex = 0;
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(mf_label, suffix[i]);
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(mfd_label, suffix[i]);
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(disp_label, suffix[i]);
  
  /* collect ID's of the element blocks in the group */
  ArrayT<StringT> block_ID(fBlockData.Length());
  for (int i = 0; i < block_ID.Length(); i++)
    block_ID[i] = fBlockData[i].ID();
  
  /* set output specifier */
  fOutputSet = new OutputSetT(GeometryCode(), block_ID, fConnectivities, 
		   n_labels, e_labels, false);
  
  /* register and get output ID */
  fMatForceOutputID = ElementSupport().RegisterOutput(*fOutputSet);
}

/* send output */
void SSQ2P1MF::WriteOutput(void)
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

/**********************Material Force Driver**************/
/*driver*/
void SSQ2P1MF::ComputeMatForce(dArray2DT& output)
{
  const char caller[] = "SSQ2P1MF::ComputeMatForce";

  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nen = NumElementNodes();
  int nmf = nnd*NumSD();

  /*dimension output array and workspace*/
  output.Dimension(nnd,3*NumSD());
  
  const dArray2DT& disp = Field()[0];
  if (disp.MajorDim() != output.MajorDim()) throw ExceptionT::kGeneralFail;

  fMatForce.Dimension(nmf);
  fDissipForce.Dimension(nmf);
  felem_rhs.Dimension(NumSD()*nen);  

  fMatForce = 0.0;
  fDissipForce = 0.0;

  /*if internal dissipation vars exists, extrapolate from element ip to nodes*/
  if (fhas_dissipation) { 
    /*assume all materials within elemblock have the same internal dissipation variables*/
    ContinuumMaterialT* pmat0 = (*fMaterialList)[ElementCard(0).MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat0);
    if (!fCurrSSMat) throw ExceptionT::kGeneralFail;

    fInternalDOF = fCurrSSMat->InternalDOF();
    int varsets = fInternalDOF.Length();
    fNumInternalVal = 0;
    for (int i = 0; i<varsets; i++)
        fNumInternalVal += fInternalDOF[i];

    fGradInternalStrain.Dimension(NumSD(),fNumInternalVal);

    fGlobalMass.Dimension(nnd);
    fGlobalVal.Dimension(nnd, fNumInternalVal);
    felem_mass.Dimension(nen);
    felem_val.Dimension(nen, fNumInternalVal);

    fGlobalVal = 0.0;
    fGlobalMass = 0.0;
  
    Extrapolate();
  }

  /*evaluate volume contributions to material and dissipation force*/
  Top();
  while (NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
    if (!fCurrSSMat) ExceptionT::GeneralFail(caller);
    
    /*Set Global Shape Functions for current element*/
    SetGlobalShape();
    MatForceVolMech(felem_rhs);
    AssembleArray(felem_rhs, fMatForce, CurrentElement().NodesX());
    if (fhas_dissipation) 
    {
      felem_val.Free();
      ExtractArray2D(fGlobalVal, felem_val, CurrentElement().NodesX());
      MatForceDissip(felem_rhs, felem_val);
      AssembleArray(felem_rhs, fDissipForce, CurrentElement().NodesX());
    }
  }

  /*add surface contribution*/
  MatForceSurfMech(fMatForce);

  /*assemble material forces and displacements into output array*/
  double* pout_force = output.Pointer();
  double* pout_dissip = output.Pointer(NumSD());
  double* pout_disp = output.Pointer(2*NumSD());
  double* pmat_force = fMatForce.Pointer();
  double* pmat_fdissip = fDissipForce.Pointer();
  const iArray2DT& eqno = Field().Equations();
  for (int i = 0; i<nnd; i++)
  {
    for (int j = 0; j<NumSD(); j++)
    {
      /*material force set to zero for kinematically constrained nodes*/
      if(eqno[i*NumSD()+j] < 1)
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
      *pout_disp++ = disp[i*NumSD()+j];
    }
    pout_force += 2*NumSD();
    pout_dissip += 2*NumSD();
    pout_disp += 2*NumSD();
  }
  WriteSummary(output);
}

void SSQ2P1MF::MatForceVolMech(dArrayT& elem_val)
{   
  const char caller[] = "SSQ2P1MF::MatForceVolMech";
  int nen = NumElementNodes();
  int elem = CurrElementNumber();
  elem_val = 0.0;
  
  /*get density*/
  double density = fCurrSSMat->Density();
  
  if (fLocAcc.IsRegistered())
  {
    SetLocalU(fLocAcc);
    fLocAcc.ReturnTranspose(fBodyForce);
  }  
  else
    fBodyForce = 0.0;
  fBodyForce *= density;
    
  /*copy acceleration and body force data into body force vector*/
  double* pbody = fBodyForce.Pointer();    
  if (fBodySchedule)
  {
    double loadfactor = fBodySchedule->Value();
    for (int i = 0; i<nen; i++)
      for (int j = 0; j<NumSD(); j++)
        *pbody++ -= fBody[j]*loadfactor*density;
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

    /*gather material data*/
    double energy = fCurrSSMat->StrainEnergyDensity();
    fCauchy = fCurrSSMat->s_ij();
    const dMatrixT& gradU = DisplacementGradient();

    double* pbody = fBodyForce.Pointer();
    double* pforce = elem_val.Pointer(); 
    if (NumSD() == 2)
    {
      /*interpolate to ip*/
      fip_body = 0.0;
      for (int i= 0; i<nen; i++)
      {
	    fip_body[0] += (*pQaU) * (*pbody++);
	    fip_body[1] += (*pQaU++) * (*pbody++);
      }	 

      /*form negative of Eshelby stress -SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
     
      fEshelby(0,0) = gradU(0,0)*fCauchy(0,0) + gradU(1,0)*fCauchy(1,0) - energy;
      fEshelby(0,1) = gradU(0,0)*fCauchy(0,1) + gradU(1,0)*fCauchy(1,1);
      fEshelby(1,0) = gradU(0,1)*fCauchy(0,0) + gradU(1,1)*fCauchy(1,0);
      fEshelby(1,1) = gradU(0,1)*fCauchy(0,1) + gradU(1,1)*fCauchy(1,1) - energy;

      const double* pDQaX = DQa(0); 
      const double* pDQaY = DQa(1);
      
      for (int j = 0; j<nen; j++)
      {
	/*add nEshelby volume integral contribution*/
       	*(pforce++) += (fEshelby[0]*(*pDQaX) + fEshelby[2]*(*pDQaY)
	       + (gradU[0]*fip_body[0]+gradU[1]*fip_body[1])*(*pQa))*(*jac)*(*weight);
	    *(pforce++) += (fEshelby[1]*(*pDQaX++) + fEshelby[3]*(*pDQaY++)
	       + (gradU[2]*fip_body[0]+gradU[3]*fip_body[1])*(*pQa++))*(*jac)*(*weight); 
      }
    }
    else if (NumSD() ==3)
    {
      /*interpolate ip values*/
      fip_body = 0.0;
      for (int i= 0; i<nen; i++)
      {
		fip_body[0] += (*pQaU) * (*pbody++);
		fip_body[1] += (*pQaU) * (*pbody++);
		fip_body[2] += (*pQaU++) * (*pbody++);
      }
	        
      /*form negative of Eshelby stress -SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
      
      fEshelby(0,0) = gradU(0,0)*fCauchy(0,0) + gradU(1,0)*fCauchy(1,0)
        + gradU(2,0)*fCauchy(2,0) - energy;
      fEshelby(0,1) = gradU(0,0)*fCauchy(0,1) + gradU(1,0)*fCauchy(1,1)
        + gradU(2,0)*fCauchy(2,1);
      fEshelby(0,2) = gradU(0,0)*fCauchy(0,2) + gradU(1,0)*fCauchy(1,2)
        + gradU(2,0)*fCauchy(2,2);      
      fEshelby(1,0) = gradU(0,1)*fCauchy(0,0) + gradU(1,1)*fCauchy(1,0)
        + gradU(2,1)*fCauchy(2,0);      
      fEshelby(1,1) = gradU(0,1)*fCauchy(0,1) + gradU(1,1)*fCauchy(1,1)
        + gradU(2,1)*fCauchy(2,1) - energy;      
      fEshelby(1,2) = gradU(0,1)*fCauchy(0,2) + gradU(1,1)*fCauchy(1,2)
        + gradU(2,1)*fCauchy(2,2);      
      fEshelby(2,0) = gradU(0,2)*fCauchy(0,0) + gradU(1,2)*fCauchy(1,0)
        + gradU(2,2)*fCauchy(2,0);     
      fEshelby(2,1) = gradU(0,2)*fCauchy(0,1) + gradU(1,2)*fCauchy(1,1)
        + gradU(2,2)*fCauchy(2,1);      
      fEshelby(2,2) = gradU(0,2)*fCauchy(0,2) + gradU(1,2)*fCauchy(1,2)
        + gradU(2,2)*fCauchy(2,2)-energy;

      const double* pDQaX = DQa(0); 
      const double* pDQaY = DQa(1);
      const double* pDQaZ = DQa(2);
      
      for (int j = 0; j<nen; j++)
      {
	    /*add Eshelby volume integral contribution*/
	    *(pforce++) += (fEshelby[0]*(*pDQaX) + fEshelby[3]*(*pDQaY) + fEshelby[6]*(*pDQaZ) 
	        +(gradU[0]*fip_body[0] + gradU[1]*fip_body[1] + gradU[2]*fip_body[2])
	        *(*pQa) )*(*jac)*(*weight);

    	*(pforce++) += (fEshelby[1]*(*pDQaX) + fEshelby[4]*(*pDQaY) + fEshelby[7]*(*pDQaZ)
	        +(gradU[3]*fip_body[0] + gradU[4]*fip_body[1] + gradU[5]*fip_body[2])
	        *(*pQa) )*(*jac)*(*weight);

	    *(pforce++) += (fEshelby[2]*(*pDQaX++) + fEshelby[5]*(*pDQaY++) + 
	        fEshelby[8]*(*pDQaZ++) + (gradU[6]*fip_body[0]+gradU[7]*fip_body[1]
	        +gradU[8]*fip_body[2])*(*pQa++))*(*jac)*(*weight);
      }
    }
    weight++;
    jac++;
  }
}

void SSQ2P1MF::MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch)
{
  
  /*obtain dimensions*/
  int nen = NumElementNodes();
  int nip = fShapes->NumIP();
  int varsets = fInternalDOF.Length();

  elem_val = 0.0;

  const double* jac = fShapes->IPDets();
  const double* weight = fShapes->IPWeights();
  fShapes->TopIP();

  while(fShapes->NextIP())
  { 
    fGradInternalStrain = 0;
    const double* pQa = fShapes->IPShapeX();
    const dArray2DT& DQa = fShapes->Derivatives_X();
    if (NumSD() ==2)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);

      /*Interpolate grad of inverse inelastic stretch to ip*/
      double* pGradX = fGradInternalStrain(0);
      double* pGradY = fGradInternalStrain(1);
      for (int i = 0; i<nen; i++)
      {
        for (int cnt = 0; cnt < fNumInternalVal; cnt++)
        {
          pGradX[cnt] += (*pDQaX)*internalstretch(i,cnt);
          pGradY[cnt] += (*pDQaY)*internalstretch(i,cnt);
	    }
    	pDQaX++;
    	pDQaY++;
      }
      //      cout << "\nfGradInternalStrain: "<<fGradInternalStrain;
      /*integrate material force*/
      const dArrayT& internalstress = fCurrSSMat->InternalStressVars();
      const double* pstress = internalstress.Pointer();
      //      cout << "\ninternalstress: "<<internalstress;
      double xval = ScalarProduct(pstress, pGradX, fInternalDOF);
      double yval = ScalarProduct(pstress, pGradY, fInternalDOF);
      //      cout << "\nxval: "<<xval;
      //      cout << "\nyval: "<<yval;

      double* pelem_val = elem_val.Pointer();
      for (int i = 0; i<nen; i++)
      {
    	(*pelem_val++) += xval*(*pQa)*(*jac)*(*weight);
	    (*pelem_val++) += yval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    else if (NumSD() ==3)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);
      const double* pDQaZ = DQa(2);

      /*Interpolate grad of iverse inelastic stretch to ip*/
      double* pGradX = fGradInternalStrain(0);
      double* pGradY = fGradInternalStrain(1);
      double* pGradZ = fGradInternalStrain(2);
      for (int i = 0; i<nen; i++)
      {	
        for (int cnt = 0; cnt < fNumInternalVal; cnt++)
        {
          pGradX[cnt] += (*pDQaX)*internalstretch(i,cnt);
	      pGradY[cnt] += (*pDQaY)*internalstretch(i,cnt);
	      pGradZ[cnt] += (*pDQaZ)*internalstretch(i,cnt);
    	}
	    pDQaX++;
	    pDQaY++;
	    pDQaZ++;
      }
            
      /*integrate material force*/
      const dArrayT& internalstress = fCurrSSMat->InternalStressVars();
      const double* pstress = internalstress.Pointer();
      double xval = ScalarProduct(pstress, pGradX, fInternalDOF);
      double yval = ScalarProduct(pstress, pGradY, fInternalDOF);
      double zval = ScalarProduct(pstress, pGradZ, fInternalDOF);

      double* pelem_val = elem_val.Pointer();
      for (int i = 0; i<nen; i++)
      {
	    (*pelem_val++) += xval*(*pQa)*(*jac)*(*weight);
	    (*pelem_val++) += yval*(*pQa)*(*jac)*(*weight);      
	    (*pelem_val++) += zval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    jac++;
    weight++;
  }
}

void SSQ2P1MF::MatForceSurfMech(dArrayT& global_val)
{

  if (fTractionList.Length() > 0)
  {
    
    int nen = NumElementNodes();
    fGradInternalStrain = 0.0;
    int numval = NumSD()*NumSD();

    /*initialize work space*/
    dMatrixT ExtrapMatrix(nen);
    dMatrixT jacobian(NumSD(),NumSD()-1);
    dMatrixT Q(NumSD());
 
    /*loop through traction card*/
    for (int k = 0; k<fTractionList.Length(); k++)
    {
      const Traction_CardT& BC_card = fTractionList[k];
      const iArrayT& surf_nodes = BC_card.Nodes();
      int nfn = surf_nodes.Length();
      int elem, facet;
      BC_card.Destination(elem, facet);

      /*get traction*/
      LocalArrayT tract(LocalArrayT::kUnspecified, nfn, NumSD());
      BC_card.CurrentValue(tract);      

      /*material force vector for facet*/
      dArrayT elem_val(nfn*NumSD());

      /*retrieve nodal coordinates of facets in local ordering*/
      LocalArrayT surf_coords(LocalArrayT::kInitCoords, nfn, NumSD()); 
      ElementSupport().RegisterCoordinates(surf_coords);
      surf_coords.SetLocal(surf_nodes);    
    
      /*retrieve element information*/
      const ElementCardT& elem_card = fElementCards[elem];
      fElementCards.Current(elem);
      SetGlobalShape();
      ContinuumMaterialT* pmat = (*fMaterialList)[elem_card.MaterialNumber()];
      fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
      if (!fCurrSSMat) ExceptionT::GeneralFail();
 
      iArrayT loc_surf_nodes(nfn);
      fShapes->NodesOnFacet(facet,loc_surf_nodes);

      /*Set element shapefunctions*/
      fShapes->SetDerivatives(); 
      const double* jac = fShapes->IPDets();
      const double* weight = fShapes->IPWeights();

      /*project elem ip vals to elem nodal vals*/
      nArrayT<dMatrixT> ExtrapGradU(nen);
      for (int i = 0; i<nen; i++)
      {
        dMatrixT& pval = ExtrapGradU[i];
        pval.Allocate(NumSD());
        pval = 0;
      }
      dArrayT ExtrapEnergy(nen);
      ExtrapEnergy = 0.0;
      ExtrapMatrix = 0;      
      fShapes->TopIP();
      while(fShapes->NextIP())
      {
        double energy = fCurrSSMat->StrainEnergyDensity();
        const dMatrixT& GradU = DisplacementGradient();
        const double* pQbU = fShapes->IPShapeU();
        for (int i = 0; i<nen; i++)
        {
          const double* pQaU = fShapes->IPShapeU();
          for (int j = 0; j<nen; j++)
	    ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
          dMatrixT& pval = ExtrapGradU[i];
          for (int cnt = 0; cnt <numval; cnt++)
            pval[cnt] += (*pQbU)*(*jac)*(*weight)*GradU[cnt];
          ExtrapEnergy[i] += (*pQbU++)*(*jac)*(*weight)*energy;
        } 
        jac++;
        weight++;
      }
      ExtrapMatrix.Inverse();
      nArrayT<dMatrixT> NodalGradU(nfn);
      dArrayT NodalEnergy(nfn);
      for (int i = 0; i<nfn; i++)
      {
        int loc_node=loc_surf_nodes[i];
        dMatrixT& nval = NodalGradU[i];
        nval.Allocate(NumSD());
        nval = 0;
        NodalEnergy[i] = 0;
        for (int j = 0; j<nen; j++)
        {
          dMatrixT& pval = ExtrapGradU[j];
          double M = ExtrapMatrix(loc_node, j);
          for (int cnt = 0; cnt<numval; cnt++)
            nval[cnt] += M*pval[cnt];
          NodalEnergy[i] += M*ExtrapEnergy[j];
        }
      }
      
      /*get surface shape function*/
      const ParentDomainT& surf_shape = 
      ShapeFunction().FacetShapeFunction(facet);
      int nip = surf_shape.NumIP();
      dArrayT ip_tract(NumSD());
      dArrayT global_tract(NumSD());
      dArrayT ip_eshelby(NumSD());     
      dMatrixT ip_gradU(NumSD());
      if (NumSD() == 2)
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
          ip_gradU = 0;
          ip_tract = 0.0;
          double ip_energy = 0;
          const double* pQaU = surf_shape.Shape(j);
	  for (int i= 0; i<nfn; i++)
	  {
	    dMatrixT& nval = NodalGradU[i];
	    ip_gradU[0] += *pQaU*nval[0];
	    ip_gradU[1] += *pQaU*nval[1];
	    ip_gradU[2] += *pQaU*nval[2];
	    ip_gradU[3] += *pQaU*nval[3];
	    ip_energy += *pQaU*NodalEnergy[i];
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
	  ip_eshelby[0] = (ip_energy*n[0]-ip_gradU[0]*global_tract[0]
		          -ip_gradU[1]*global_tract[1])*thickness;
	  ip_eshelby[1] = (ip_energy*n[1]-ip_gradU[2]*global_tract[0]
	                  -ip_gradU[3]*global_tract[1])*thickness;
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
      else if (NumSD() == 3)
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
          ip_gradU = 0;
          ip_tract = 0.0;
          double ip_energy = 0;
          const double* pQaU = surf_shape.Shape(j);
	  for (int i= 0; i<nfn; i++)
	  {
	    dMatrixT& nval = NodalGradU[i];
	    ip_gradU[0] += *pQaU*nval[0];
	    ip_gradU[1] += *pQaU*nval[1];
	    ip_gradU[2] += *pQaU*nval[2];
	    ip_gradU[3] += *pQaU*nval[3];
	    ip_gradU[4] += *pQaU*nval[4];
	    ip_gradU[5] += *pQaU*nval[5];
	    ip_gradU[6] += *pQaU*nval[6];
	    ip_gradU[7] += *pQaU*nval[7];
	    ip_gradU[8] += *pQaU*nval[8];	        
	    ip_energy +=  *pQaU*NodalEnergy[i];
	    ip_tract[0] += (*pQaU)*(*ptract_X++);
	    ip_tract[1] += (*pQaU)*(*ptract_Y++);
	    ip_tract[2] += (*pQaU++)*(*ptract_Z++);
	  }
	  /*surface normal*/
	  double* n = Q(2);
	  if (BC_card.CoordSystem() == Traction_CardT::kLocal)
	  {
	    /*rotate traction from local to global coords*/
	    global_tract[0] = 
Q[0]*ip_tract[0]+Q[3]*ip_tract[1]+Q[6]*ip_tract[2];
	    global_tract[1] = 
Q[1]*ip_tract[0]+Q[4]*ip_tract[1]+Q[7]*ip_tract[2];
	    global_tract[2] = 
Q[2]*ip_tract[0]+Q[5]*ip_tract[1]+Q[8]*ip_tract[2];
     	  }
	  else if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    global_tract[0] = ip_tract[0];
	    global_tract[1] = ip_tract[1];
	    global_tract[2] = ip_tract[2];
	  }
	  ip_eshelby[0] = (ip_energy*n[0]-ip_gradU[0]*global_tract[0]
			   -ip_gradU[1]*global_tract[1]-ip_gradU[2]*global_tract[2]);
	  ip_eshelby[1] = (ip_energy*n[1]-ip_gradU[3]*global_tract[0]
			   +ip_gradU[4]*global_tract[1]-ip_gradU[5]*global_tract[2]);
	  ip_eshelby[2] = (ip_energy*n[2]-ip_gradU[6]*global_tract[0]
			   +ip_gradU[7]*global_tract[1]-ip_gradU[8]*global_tract[2]);
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
      AssembleArray(elem_val, global_val, surf_nodes);
    }
  }
}

/****************utitlity function******************************/
void SSQ2P1MF::Extrapolate(void)
{
  const char caller[] = "SSQ2P1MF::Extrapolate";   
  
  Top();
  while (NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
    if (!fCurrSSMat) throw ExceptionT::kGeneralFail;
    SmallStrainQ2P1::SetGlobalShape();

    felem_val = 0.0;
    felem_mass = 0.0;

    const double* jac = fShapes->IPDets();;
    const double* weight = fShapes->IPWeights();
    fShapes->TopIP();
    while(fShapes->NextIP())
    {
	const dArrayT& internalstrains = fCurrSSMat->InternalStrainVars();
	//	cout << "\n internalstrains: "<< internalstrains<<endl;
        const double* pQbU = fShapes->IPShapeU();
	
	for (int i=0; i<NumElementNodes(); i++)
        {
            const double* pQaU = fShapes->IPShapeU();

            /*lumped element mass matrix*/
            for (int j = 0; j<NumElementNodes(); j++)
                felem_mass[i] += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  

            /*element forcing function*/
            for (int cnt = 0; cnt < fNumInternalVal; cnt++)
	            felem_val(i,cnt) += (*pQbU)*(*jac)*(*weight)*internalstrains[cnt];

	        pQbU++;
        }
        weight++;
        jac++;
    }    
    AssembleArray(felem_mass, fGlobalMass, CurrentElement().NodesX());
    AssembleArray2D(felem_val, fGlobalVal, CurrentElement().NodesX());
  }  
  for (int i = 0; i< fNumGroupNodes; i++)
    for (int j = 0; j< fNumInternalVal; j++)
        fGlobalVal(i,j) /= fGlobalMass[i];
  //  cout << "\nfGlobalVal: "<<fGlobalVal;
}

