/* $Id: SSMF.cpp,v 1.13 2005-03-02 17:41:47 paklein Exp $ */
#include "SSMF.h"

#include "OutputSetT.h"
#include "Traction_CardT.h"
#include "ScheduleT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"
#include "GraphT.h"
#include "GeometryT.h"
#include "ModelManagerT.h"
#include "CommunicatorT.h"
#include "SSMatSupportT.h"
#include "KBC_ControllerT.h"
#include "KBC_CardT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ParameterUtils.h"
#include "eIntegratorT.h"
#include "SolidMatListT.h"

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

using namespace Tahoe;

/* constructor */
SSMF::SSMF(const ElementSupportT& support):
  SmallStrainT(support),
  MFSupportT(support),
  fdynamic(false),
  ftraction(LocalArrayT::kUnspecified),
  fsurf_disp(LocalArrayT::kDisp),
  fsurf_coords(LocalArrayT::kInitCoords) 
{
	SetName("small_strain_material_force");
}

/* information about subordinate parameter lists */
void SSMF::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SmallStrainT::DefineSubs(sub_list);
	
	/* nodes where force will be calculated */
	sub_list.AddSub("mf_nodes_ID_list");

	/* nodes on the boundary */
	sub_list.AddSub("boundary_mf_nodes_ID_list", ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void SSMF::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SmallStrainT::TakeParameterList(list);

	/* dimension workspace */
	fEshelby.Dimension(NumSD());

	fBodyForce.Dimension(NumSD()*NumElementNodes());
	fip_body.Dimension(NumSD());
	fVel.Dimension(NumSD());
	fAcc.Dimension(NumSD());
	fGradVel.Dimension(NumSD());
    
	fGradU_List.Dimension(NumIP());
	for (int i = 0; i< NumIP(); i++)
		fGradU_List[i].Dimension(NumSD());

	b_dil.Dimension(NumSD()*NumElementNodes());

	/* nodes over which to calculate material forces */
	const ParameterListT& mf_nodes = list.GetList("mf_nodes_ID_list");
	StringListT::Extract(mf_nodes, fNID);

	/* boundary nodes */
	const ParameterListT* boundary_mf_nodes = list.List("boundary_mf_nodes_ID_list");
	if (boundary_mf_nodes)
		StringListT::Extract(*boundary_mf_nodes, fBoundID);

	fhas_dissipation = fMaterialList->HasHistoryMaterials();
}

void SSMF::SetGlobalShape(void)
{
	SmallStrainT::SetGlobalShape();

	for (int ip = 0; ip < NumIP(); ip++) 
		fShapes->GradU(fLocDisp, fGradU_List[ip], ip);

	if (fStrainDispOpt == kMeanDilBbar) {
		SetMeanGradient(fMeanGradient);
		for (int ip = 0; ip< NumIP(); ip++) {
	
			double mean_correction = (grad_dil(fShapes->Derivatives_U(ip), fMeanGradient))/3.0;
			fGradU_List[ip](0,0) -= mean_correction;
			fGradU_List[ip](1,1) -= mean_correction;
		}
	}
}

/***************************outputs managers***********************************/
/* register self for output */
void SSMF::RegisterOutput(void)
{
  /* inherited */
  SolidElementT::RegisterOutput();

  ArrayT<StringT> n_labels(4*NumSD());
  ArrayT<StringT> e_labels;
  
  StringT mf_label = "MF";
  StringT mfd_label = "DF";
  StringT mfdd_label = "KF";
  StringT disp_label = "D";
  const char* suffix[3] = {"_X", "_Y", "_Z"};
  int dex = 0;
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(disp_label, suffix[i]);
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(mf_label, suffix[i]);
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(mfd_label, suffix[i]);
  for (int i = 0; i < NumSD(); i++)
    n_labels[dex++].Append(mfdd_label, suffix[i]);
  
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
void SSMF::WriteOutput(void)
{
  /* inherited */
  SolidElementT::WriteOutput();
  
  /* calculate output values */
  dArray2DT n_values; 
  dArray2DT e_values;

  MapOutput();
  ComputeMatForce(n_values);
   /* send to output */
  const CommunicatorT& comm = ElementSupport().Communicator();
  if (comm.Size() == 1)
      WriteSummary(n_values);
 
  /* send to output */
  ElementSupport().WriteOutput(fMatForceOutputID, n_values, e_values);
}

void SSMF::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
			AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
  /*inherrited function*/
  ElementBaseT::ConnectsU(connects_1, connects_2);
  const CommunicatorT& comm = ElementSupport().Communicator();
  
  /*check for parallel execution*/
  /*make graph from element connectivities*/
  bool verbose = true;
  GraphT graph(true);
  
  for (int i = 0; i<fConnectivities.Length(); i++)
    graph.AddGroup(*fConnectivities[i]);
  graph.MakeGraph();
  
  int nnd = graph.NumNodes(); 
  iArrayT first_edges;
  AutoFill2DT<int> fill(nnd, 1, 5, 16);
  for (int i = 0; i<nnd; i++)
  {
    graph.GetEdges(i, first_edges);
    fill.Append(i,first_edges);
    for (int j = 0; j<first_edges.Length(); j++)
    {  
      iArrayT second_edges;
      graph.GetEdges(first_edges[j],second_edges);
      /*Do we need to reorder the neighbor list so that the node number appears in the first entry?*/
      fill.AppendUnique(i,second_edges);
    }
  }
  
  /* non-const this */
  SSMF* non_const_this = (SSMF*) this; 
  non_const_this->fXConnects.Copy(fill); 
  
  /*
    iArrayT tmp(fXConnects.Length(), fXConnects.Pointer());
    tmp++;
    ofstreamT& out = ElementSupport().Output();
    out << "\nextended connectivities: \n";
    fXConnects.WriteNumbered(ElementSupport().Output());
    tmp--;
  */
  
  connects_2.Append(&fXConnects);
}

/**********************Material Force Driver**************/
/*driver*/
void SSMF::ComputeMatForce(dArray2DT& output)
{
  const char caller[] = "SSMF::ComputeMatForce";

  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nen = NumElementNodes();
  int nmf = nnd*NumSD();

  /*dimension output array and workspace*/
  output.Dimension(nnd,4*NumSD());
  output = 0.0;

  const dArray2DT& global_disp = Field()[0];
  
  fMatForce.Dimension(nmf);
  fDissipForce.Dimension(nmf);
  fDynForce.Dimension(nmf);
  fGroupDisp.Dimension(nmf);

  felem_rhs.Dimension(NumSD()*nen);  

  fMatForce = 0.0;
  fDissipForce = 0.0;
  fDynForce = 0.0;
  fGroupDisp = 0.0;

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

  /*check for dynamic analysis*/
  fdynamic = (fIntegrator->Order() > 0);

  /*evaluate volume contributions to material and dissipation force*/
  Top();
  while (NextElement())
  {
    if (CurrentElement().Flag() != ElementCardT::kOFF)
    {  
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
    if (!fCurrSSMat) ExceptionT::GeneralFail(caller);
    
    /*Set Global Shape Functions for current element*/
    SetGlobalShape();
    SetLocalU(fLocDisp);
 
    if (fdynamic)
    {
      SetLocalU(fLocAcc);
      SetLocalU(fLocVel);
      MatForceDynamic(felem_rhs);
      AssembleArray(felem_rhs, fDynForce, CurrentElement().NodesX());
    }

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
    /*gather displacements*/
    GatherDisp(global_disp,fGroupDisp, CurrentElement().NodesX());   
  }  
  
  /*add surface contribution*/
  MatForceSurfMech(fMatForce);

  /*assemble material forces and displacements into output array*/
  double* pout_disp = output.Pointer();
  double* pout_force = output.Pointer(NumSD());
  double* pout_dissip = output.Pointer(2*NumSD());
  double* pout_dyn = output.Pointer(3*NumSD());

  double* pmat_force = fMatForce.Pointer();
  double* pmat_fdissip = fDissipForce.Pointer();
  double* pmat_fdyn = fDynForce.Pointer();
  double* pmat_fdisp = fGroupDisp.Pointer();

  ModelManagerT& model = ElementSupport().ModelManager();
 			
 const ArrayT<KBC_ControllerT*>& kbc_controllers = Field().KBC_Controllers();
 for (int i = 0; i < kbc_controllers.Length(); i++)
 {
	const ArrayT<KBC_CardT>& kbc_cards = kbc_controllers[i]->KBC_Cards();
	for (int j = 0; j < kbc_cards.Length(); j++) 
	{
		int mapped_node = fMap[kbc_cards[j].Node()];
		if (mapped_node > -1)
		{
			fExclude[mapped_node] = 1;
		}
	}
 }
 
 for (int i = 0; i<nnd; i++)
 {
	for (int j = 0; j<NumSD(); j++)
	{
		if (fExclude[i] == 1) {
			*pout_force++ = 0.0;
			pmat_force++;
		}
		else
			*pout_force++ = (*pmat_force++) - (*pmat_fdissip) - (*pmat_fdyn);

		*pout_dissip++ = (*pmat_fdissip++);
		*pout_dyn++ = (*pmat_fdyn++);
		*pout_disp++ = (*pmat_fdisp++);
	}
	pout_force += 3*NumSD();
	pout_dissip += 3*NumSD();
	pout_dyn += 3*NumSD();
	pout_disp += 3*NumSD();
  }
}

void SSMF::MatForceVolMech(dArrayT& elem_val)
{
  const char caller[] = "SSMF::MatForceVolMech";
  int nen = NumElementNodes();
  int elem = CurrElementNumber();

  elem_val = 0.0;
  
  /*get density*/
  double density = fCurrSSMat->Density();
  
  fBodyForce = 0.0;
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
 
    const dSymMatrixT& Cauchy = fCurrSSMat->s_ij();
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
      fEshelby(0,0) = gradU(0,0)*Cauchy(0,0) + gradU(1,0)*Cauchy(1,0) - energy;
      fEshelby(0,1) = gradU(0,0)*Cauchy(0,1) + gradU(1,0)*Cauchy(1,1);
      fEshelby(1,0) = gradU(0,1)*Cauchy(0,0) + gradU(1,1)*Cauchy(1,0);
      fEshelby(1,1) = gradU(0,1)*Cauchy(0,1) + gradU(1,1)*Cauchy(1,1) - energy;
      
      if (fdynamic)
      {
		fShapes->InterpolateU(fLocVel, fVel);
		fEshelby(0,0) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]);
		fEshelby(1,1) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]);
		if (elem ==6 && 0) { 
			cout << "\nstress: "<<Cauchy;
			cout << "\nenergy: "<<energy;
			cout << "\nfVel: "<<fVel; 
			cout << "\nfEshelby: "<<fEshelby;
		}
      }

      const double* pDQaX = DQa(0); 
      const double* pDQaY = DQa(1);
      
      for (int j = 0; j<nen; j++)
      {
	/*add nEshelby volume integral contribution*/

       	*(pforce++) += (fEshelby(0,0)*(*pDQaX) + fEshelby(0,1)*(*pDQaY)
          +(gradU(0,0)*fip_body[0]+gradU(1,0)*fip_body[1])*(*pQa))
          *(*jac)*(*weight);
	*(pforce++) += (fEshelby(1,0)*(*pDQaX++) + fEshelby(1,1)*(*pDQaY++)
	  +(gradU(0,1)*fip_body[0]+gradU(1,1)*fip_body[1])*(*pQa++))
          *(*jac)*(*weight); 
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
      
      fEshelby(0,0) = gradU(0,0)*Cauchy(0,0) + gradU(1,0)*Cauchy(1,0)
        + gradU(2,0)*Cauchy(2,0) - energy;
      fEshelby(0,1) = gradU(0,0)*Cauchy(0,1) + gradU(1,0)*Cauchy(1,1)
        + gradU(2,0)*Cauchy(2,1);
      fEshelby(0,2) = gradU(0,0)*Cauchy(0,2) + gradU(1,0)*Cauchy(1,2)
        + gradU(2,0)*Cauchy(2,2);      
      fEshelby(1,0) = gradU(0,1)*Cauchy(0,0) + gradU(1,1)*Cauchy(1,0)
        + gradU(2,1)*Cauchy(2,0);      
      fEshelby(1,1) = gradU(0,1)*Cauchy(0,1) + gradU(1,1)*Cauchy(1,1)
        + gradU(2,1)*Cauchy(2,1) - energy;      
      fEshelby(1,2) = gradU(0,1)*Cauchy(0,2) + gradU(1,1)*Cauchy(1,2)
        + gradU(2,1)*Cauchy(2,2);      
      fEshelby(2,0) = gradU(0,2)*Cauchy(0,0) + gradU(1,2)*Cauchy(1,0)
        + gradU(2,2)*Cauchy(2,0);     
      fEshelby(2,1) = gradU(0,2)*Cauchy(0,1) + gradU(1,2)*Cauchy(1,1)
        + gradU(2,2)*Cauchy(2,1);      
      fEshelby(2,2) = gradU(0,2)*Cauchy(0,2) + gradU(1,2)*Cauchy(1,2)
        + gradU(2,2)*Cauchy(2,2)-energy;

      if (fdynamic)
      {
	fShapes->InterpolateU(fLocVel, fVel);
	fEshelby(0,0) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]
				      +fVel[2]*fVel[2]);
	fEshelby(1,1) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]
				      +fVel[2]*fVel[2]);
	fEshelby(2,2) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]
				      +fVel[2]*fVel[2]);
      }

      const double* pDQaX = DQa(0); 
      const double* pDQaY = DQa(1);
      const double* pDQaZ = DQa(2);
      
      for (int j = 0; j<nen; j++)
      {
	/*add Eshelby volume integral contribution*/
	*(pforce++) += (fEshelby[0]*(*pDQaX) + fEshelby[3]*(*pDQaY) 
	 +fEshelby[6]*(*pDQaZ)+(gradU[0]*fip_body[0] + gradU[1]*fip_body[1] 
         +gradU[2]*fip_body[2])*(*pQa) )*(*jac)*(*weight);

    	*(pforce++) += (fEshelby[1]*(*pDQaX) + fEshelby[4]*(*pDQaY) 
	 +fEshelby[7]*(*pDQaZ)+(gradU[3]*fip_body[0] + gradU[4]*fip_body[1] 
	 +gradU[5]*fip_body[2])*(*pQa) )*(*jac)*(*weight);

	*(pforce++) += (fEshelby[2]*(*pDQaX++) + fEshelby[5]*(*pDQaY++) 
	 +fEshelby[8]*(*pDQaZ++) + (gradU[6]*fip_body[0]+gradU[7]*fip_body[1]
	 +gradU[8]*fip_body[2])*(*pQa++))*(*jac)*(*weight);
      }
    }
    weight++;
    jac++;
  }
  if (elem ==6 && 0) cout << "VolElem: "<<elem_val;
}

void SSMF::MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch)
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

      /*integrate material force*/
      const dArrayT& internalstress = fCurrSSMat->InternalStressVars();
      const double* pstress = internalstress.Pointer();
      double xval = -ScalarProduct(pstress, pGradX, fInternalDOF);
      double yval = -ScalarProduct(pstress, pGradY, fInternalDOF);

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
      double xval = -ScalarProduct(pstress, pGradX, fInternalDOF);
      double yval = -ScalarProduct(pstress, pGradY, fInternalDOF);
      double zval = -ScalarProduct(pstress, pGradZ, fInternalDOF);

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

void SSMF::MatForceDynamic(dArrayT& elem_val)
{
  const char caller[] = "SSMF::MatForceDynamic";
  int nen = NumElementNodes();
  int elem = CurrElementNumber();
  elem_val = 0.0;  

  double density = fCurrSSMat->Density();

  /*intialize shape function data*/
  const double* jac = fShapes->IPDets();
  const double* weight = fShapes->IPWeights();

  fShapes->TopIP();
  while(fShapes->NextIP())
  {
    /*get shape function and derivatives at integration point*/
    const double* pQa = fShapes->IPShapeX();
//    const double* pQaU = fShapes->IPShapeU();
//    const dArray2DT& DQa = fShapes->Derivatives_X();

    /*integration point values*/

    const dMatrixT& gradU = DisplacementGradient();
    fShapes->GradU(fLocVel,fGradVel); 
    fShapes->InterpolateU(fLocVel, fVel);
    fShapes->InterpolateU(fLocAcc, fAcc);

    double* pelem_val = elem_val.Pointer();
    if (NumSD() ==2)
    {
      for (int i = 0; i<nen; i++)
      {
	double xval = density*(fGradVel[0]*fVel[0]+fGradVel[1]*fVel[1]
			       -gradU[0]*fAcc[0]-gradU[1]*fAcc[1]);
	double yval = density*(fGradVel[2]*fVel[0]+fGradVel[3]*fVel[1]
			       -gradU[2]*fAcc[0]-gradU[3]*fAcc[1]);

    	*pelem_val++ += xval*(*pQa)*(*jac)*(*weight);
	*pelem_val++ += yval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    else if (NumSD() == 3)
    {
      for (int i = 0; i<nen; i++)
      {
	double xval = density*(fGradVel(0,0)*fVel[0]+fGradVel(1,0)*fVel[1]+fGradVel(2,0)*fVel[2]
			       -gradU(0,0)*fAcc[0]-gradU(1,0)*fAcc[1]-gradU(2,0)*fAcc[2]);
	double yval = density*(fGradVel(0,1)*fVel[0]+fGradVel(1,1)*fVel[1]+fGradVel(2,1)*fVel[2]
			       -gradU(0,1)*fAcc[0]-gradU(1,1)*fAcc[1]-gradU(2,1)*fAcc[2]);
	double zval = density*(fGradVel(0,2)*fVel[0]+fGradVel(1,2)*fVel[1]+fGradVel(2,2)*fVel[2]
			       -gradU(0,2)*fAcc[0]-gradU(1,2)*fAcc[1]-gradU(2,2)*fAcc[2]); 
	*pelem_val++ += xval*(*pQa)*(*jac)*(*weight);
	*pelem_val++ += yval*(*pQa)*(*jac)*(*weight);      
	*pelem_val++ += zval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    jac++;
    weight++;
  }
    if (elem ==6 && 0) cout<<"\nDynElem: "<<elem_val;
}

void SSMF::MatForceSurfMech(dArrayT& global_val)
{

  if (fTractionList.Length() > 0)
  {
    
    /*dimensions*/
    int nen = NumElementNodes();
    int nsd = NumSD();
    const Traction_CardT& BC_card = fTractionList[0];
    const iArrayT& surf_nodes = BC_card.Nodes();
    int nfn = surf_nodes.Length();

    /*dimension workspace for surface quantities*/
    ftraction.Dimension(nfn,nsd);                
    fsurf_disp.Dimension(nfn,nsd);               
    Field().RegisterLocal(fsurf_disp);
    fsurf_coords.Dimension(nfn, nsd);            
    ElementSupport().RegisterCoordinates(fsurf_coords);

    /*facet contribution to material force*/ 
    fsurf_val.Dimension(nfn*nsd);                            

    /*integration point values*/
    fip_tract.Dimension(nsd);                                
    fgradU.Dimension(nsd,nsd-1);                
    fjacobian.Dimension(nsd, nsd-1);             
    fjac_loc.Dimension(nsd-1);
    fQ.Dimension(nsd);                           

    /*loop through traction card*/
    for (int k = 0; k<fTractionList.Length(); k++)
    {
      /*retrieve information from traction card*/
      const Traction_CardT& BC_card = fTractionList[k];
      const iArrayT& surf_nodes = BC_card.Nodes();
      int elem, facet;

      /*retrieve traction information*/
      BC_card.Destination(elem, facet);
      BC_card.CurrentValue(ftraction);      
      fsurf_coords.SetLocal(surf_nodes);
      fsurf_disp.SetLocal(surf_nodes);

      //	cout << "\n Elem "<<elem<<" MatSurfVol";
      /*get surface shape function*/
      const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
      int nip = surf_shape.NumIP();
      fsurf_val = 0.0;
      const double* ip_w = surf_shape.Weight();
      for (int j = 0; j<nip; j++)
      {
	/*calculate surface jacobian and jacobian determinant*/
	surf_shape.DomainJacobian(fsurf_coords,j,fjacobian);
	double detj = surf_shape.SurfaceJacobian(fjacobian, fQ);
	if (nsd ==2)
	{ 
	  double* t = fQ(0);
	  double* n = fQ(1);
	  /*interpolate tractions to integration points*/
	  surf_shape.Interpolate(ftraction, fip_tract,j);
	  
	  /*rotate tractions from global to local coordinates*/
	  if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    double local_tract1 = t[0]*fip_tract[0]+t[1]*fip_tract[1];
	    double local_tract2 = n[0]*fip_tract[0]+n[1]*fip_tract[1];
	    fip_tract[0] = local_tract1;
	    fip_tract[1] = local_tract2;
	  }
	  //      	  cout <<"\nIP traction: "<<fip_tract;
	  //	  cout <<"\nfjacobian: "<<fjacobian;
	  //	  cout <<"\ntangent: ("<<t[0]<<", "<<t[1]<<")";
	  /*get displacement gradient*/
	  fgradU = 0;
	  double& du1 = fgradU[0];
	  double& du2 = fgradU[1];
	  
	  /*spatial derivative along the facet direction 1.0/(ds/dz)*/
	  double jac = (t[0]*fjacobian[0]+t[1]*fjacobian[1]);
	  if (jac <= 0.0)
	  {
	    //	      cout <<"\nUplagMF::MatForceSurf"<<endl;
	    throw ExceptionT::kBadJacobianDet;
	  }
	  fjac_loc[0] = 1.0/jac;
	  /*parent domain shape function derivatives*/
	  const double* pdz = surf_shape.DShape(j,0);
	  const double& j11 = fjac_loc[0];
	  const double* pu1 = fsurf_disp(0);
	  const double* pu2 = fsurf_disp(1);
	    
	  //       	  cout << "\n U "<< fsurf_disp;
	  //       	  cout << "\n fjac_loc "<<fjac_loc[0];
	  for (int i = 0; i<nfn; i++)
	  {
	    /*Rotate global displacement vector to local coordinates*/
	    double local_u1 = t[0]*(*pu1)+t[1]*(*pu2);
	    double local_u2 = n[0]*(*pu1)+n[1]*(*pu2);
	    //       	    cout << "\n local disp: ("<<local_u1<<", "<<local_u2<<")";
	    //       	    cout << "\n Dz: "<< *pdz;
	    /*calculate displacement gradient along local direction*/
	    du1 += (*pdz)*j11*(local_u1);
	    du2 += (*pdz)*j11*(local_u2);
	    pdz++; pu1++; pu2++;
	  }
	  //       	  cout << "\n GradU" <<fgradU;
	  
	  /*calculate material traction*/
	  double mat_tract = -(fgradU[0])*fip_tract[0]-fgradU[1]*fip_tract[1];
	  /*rotate material traction to global coordinates*/
	  //       	  cout << "\n mat_tract: "<< mat_tract;
	  
	  double mat_tractx = t[0]*mat_tract;
	  double mat_tracty = t[1]*mat_tract;
	  //       	  cout << "\n mat_tractx: "<< mat_tractx;
	  //       	  cout << "\n mat_tracty: "<< mat_tracty;
	  
	  /*integrate material force*/
	  const double* pQa = surf_shape.Shape(j);
	  double* psurf_val = fsurf_val.Pointer();
	  for (int i = 0; i<nfn; i++)
	  {
	    (*psurf_val++) += (*pQa)*(mat_tractx)*detj*(*ip_w);
	    (*psurf_val++) += (*pQa++)*(mat_tracty)*detj*(*ip_w);
	  }
	  ip_w++;
	}
	else if (NumSD() == 3)
	{
	  double* n1 = fQ(0);
	  double* n2 = fQ(1);
	  double* n3 = fQ(2);
	  /*interpolate tractions to integration points*/
	  surf_shape.Interpolate(ftraction, fip_tract,j);
	  
	  /*rotate tractions from global to local coordinates*/
	  if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    double local_tract1 = n1[0]*fip_tract[0]+n1[1]*fip_tract[1]+n1[2]*fip_tract[2];
	    double local_tract2 = n2[0]*fip_tract[0]+n2[1]*fip_tract[1]+n2[2]*fip_tract[2];
	    double local_tract3 = n3[0]*fip_tract[0]+n3[1]*fip_tract[1]+n2[2]*fip_tract[2];
	    fip_tract[0] = local_tract1;
	    fip_tract[1] = local_tract2;
	    fip_tract[2] = local_tract3;
	  }
	  /*get displacement gradient*/
	  fgradU = 0;
	  double& du11 = fgradU[0];
	  double& du21 = fgradU[1];
	  double& du31 = fgradU[2];
	  double& du12 = fgradU[3];
	  double& du22 = fgradU[4];
	  double& du32 = fgradU[5];
	  
	  /*spatial derivative along the facet direction 1.0/(ds/dz)*/
	  fjac_loc[0] = (n1[0]*fjacobian[0]+n1[1]*fjacobian[1]+n1[2]*fjacobian[2]);
	  fjac_loc[1] = (n2[0]*fjacobian[0]+n2[1]*fjacobian[1]+n2[2]*fjacobian[2]);
	  fjac_loc[2] = (n1[0]*fjacobian[4]+n1[1]*fjacobian[5]+n1[2]*fjacobian[6]);
	  fjac_loc[3] = (n2[0]*fjacobian[4]+n2[1]*fjacobian[5]+n2[2]*fjacobian[6]);
	  double jac = fjac_loc.Det();
	  if (jac <= 0.0)
	  {
	    cout <<"\nUplagMF::MatForceSurfMech";
	    throw ExceptionT::kBadJacobianDet;
	  }	  
	  fjac_loc.Inverse();
	  
	  /*parent domain shape function derivatives*/
	  const double* pdz = surf_shape.DShape(j,0);
	  const double* pdn = surf_shape.DShape(j,1);
	  
	  const double& j11 = fjac_loc[0];
	  const double& j21 = fjac_loc[1];
	  const double& j12 = fjac_loc[2];
	  const double& j22 = fjac_loc[3];
	  
	  const double* pu1 = fsurf_disp(0);
	  const double* pu2 = fsurf_disp(1);
	  const double* pu3 = fsurf_disp(2);
	  
	  for (int i = 0; i<nfn; i++)
	  {
	    /*Rotate global displacement vector to local coordinates*/
	    double local_u1 = n1[0]*(*pu1)+n1[1]*(*pu2)+n1[2]*(*pu3);
	    double local_u2 = n2[0]*(*pu1)+n2[1]*(*pu2)+n2[2]*(*pu3);
	    double local_u3 = n3[0]*(*pu1)+n3[1]*(*pu2)+n3[2]*(*pu3);
	    
	    /*calculate displacement gradient along local direction*/
	    du11 += (*pdz)*j11*(local_u1) + (*pdn)*j21*(local_u1);
	    du21 += (*pdz)*j11*(local_u2) + (*pdn)*j21*(local_u2);
	    du31 += (*pdz)*j11*(local_u3) + (*pdn)*j21*(local_u3);
	    
	    du12 += (*pdz)*j12*(local_u1) + (*pdn)*j22*(local_u1);
	    du22 += (*pdz)*j12*(local_u2) + (*pdn)*j22*(local_u2);
	    du32 += (*pdz)*j12*(local_u3) + (*pdn)*j22*(local_u3);
	    
	    pdz++; pdn++; pu1++; pu2++; pu3++;
	  }
	  /*calculate material traction*/
	  double mat_tract1 = -(fgradU[0])*fip_tract[0]-fgradU[1]*fip_tract[1]+fgradU[2]*fip_tract[2];
	  double mat_tract2 = -fgradU[3]*fip_tract[0]-(fgradU[4])*fip_tract[1]+fgradU[5]*fip_tract[2];
	  
	  /*rotate material traction to global coordinates*/
	  double mat_tractx = n1[0]*mat_tract1 + n2[0]*mat_tract2;
	  double mat_tracty = n1[1]*mat_tract1 + n2[1]*mat_tract2;
	  double mat_tractz = n1[2]*mat_tract1 + n2[2]*mat_tract2;
	  
	  /*integrate material force*/
	  const double* pQa = surf_shape.Shape(j);
	  double* psurf_val = fsurf_val.Pointer();
	  for (int i = 0; i<nfn; i++)
	  {
	    (*psurf_val++) += (*pQa)*(mat_tractx)*detj*(*ip_w);
	    (*psurf_val++) += (*pQa)*(mat_tracty)*detj*(*ip_w);
	    (*psurf_val++) += (*pQa++)*(mat_tractz)*detj*(*ip_w);
	  }
	  ip_w++;
	}
      }
      AssembleArray(fsurf_val, global_val, surf_nodes);
      //      cout << "\nsurf_val: "<<fsurf_val;
    }
  }
}

/****************utitlity functions******************************/
void SSMF::Extrapolate(void)
{
  const char caller[] = "SSMF::Extrapolate";   

  int nen = NumElementNodes();

  Top();
  while (NextElement())
  {
    /*initialize element matrices*/
    felem_mass = 0.0;

    if (CurrentElement().Flag() != ElementCardT::kOFF)
    {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
    if (!fCurrSSMat) throw ExceptionT::kGeneralFail;
    SmallStrainT::SetGlobalShape();

    felem_val = 0.0;

    const double* jac = fShapes->IPDets();;
    const double* weight = fShapes->IPWeights();
    fShapes->TopIP();
    while(fShapes->NextIP())
    {
	const dArrayT& internalstrains = fCurrSSMat->InternalStrainVars();
        const double* pQbU = fShapes->IPShapeU();
	
	    for (int i=0; i<nen; i++)
        {
            const double* pQaU = fShapes->IPShapeU();

            /*lumped element mass matrix*/
            for (int j = 0; j<NumElementNodes(); j++)
	    felem_mass[i] += (*pQaU++)*(*pQbU)*(*jac)*(*weight); 

            /*element forcing function*/
            for (int cnt = 0; cnt < fNumInternalVal; cnt++)
	      felem_val(i,cnt) += (*pQbU)*(*jac)*(*weight)
		*internalstrains[cnt];
	    pQbU++;
        }
        weight++;
        jac++;
    } 
    AssembleArray(felem_mass, fGlobalMass, CurrentElement().NodesX());
    AssembleArray2D(felem_val, fGlobalVal, CurrentElement().NodesX());
  }
  }  
  for (int i = 0; i< fNumGroupNodes; i++)
    for (int j = 0; j< fNumInternalVal; j++)
        fGlobalVal(i,j) /= fGlobalMass[i];
}

const double& SSMF::grad_dil(const dArray2DT& DNa, const dArray2DT& mean_gradient)
{

	int nnd = DNa.MinorDim();
	double* pB = b_dil.Pointer();

	/* 2D */
	if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
			
		const double* pBmx = mean_gradient(0);
		const double* pBmy = mean_gradient(1);
			
		for (int i = 0; i < nnd; i++)
		{
			*pB++ = ((*pBmx++) - (*pNax++))/3.0;
			*pB++ = ((*pBmy++) - (*pNay++))/3.0;
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);

		const double* pBmx = mean_gradient(0);
		const double* pBmy = mean_gradient(1);
		const double* pBmz = mean_gradient(2);
			
		for (int i = 0; i < nnd; i++)
		{
			*pB++ = ((*pBmx++) - (*pNax++))/3.0;
			*pB++ = ((*pBmy++) - (*pNay++))/3.0;
			*pB++ = ((*pBmz++) - (*pNaz++))/3.0;
		}
	}
	fLocDisp.ReturnTranspose(fLocDispTranspose);
	tmp = dArrayT::Dot(b_dil, fLocDispTranspose);

	return(tmp);
}
