/* $Id: UpLagMF.cpp,v 1.17 2004-07-15 08:28:31 paklein Exp $ */
#include <ctype.h>

#include "UpLagMF.h"
#include "OutputSetT.h"
#include "Traction_CardT.h"
#include "ScheduleT.h"
#include "ShapeFunctionT.h"
#include "FSSolidMatT.h"
#include "GeometryT.h"
#include "ModelManagerT.h"
#include "CommunicatorT.h"
#include "GraphT.h"
#include "FEManagerT.h"
#include "ifstreamT.h"

/* materials lists */
#include "MaterialListT.h"
#if 0
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"
#endif

using namespace Tahoe;

/* constructor */
UpLagMF::UpLagMF(const ElementSupportT& support, const FieldT& field):
  UpdatedLagrangianT(support),
  MFSupportT(support),
  LocalizeT(support),
  fdynamic(false),
  ftraction(LocalArrayT::kUnspecified),
  fsurf_disp(LocalArrayT::kDisp),
  fsurf_coords(LocalArrayT::kInitCoords){

  fMassType = kConsistentMass;
}

void UpLagMF::Initialize(void)
{
ExceptionT::GeneralFail("UpLagMF::Initialize", "out of date");
#if 0	
  UpdatedLagrangianT::Initialize();

  int nel = ElementBaseT::NumElements();
  ostream& out = ElementSupport().Output();

  /*localization*/
  fcheckflag.Dimension(nel);
  flocflag.Dimension(nel);
  flocflagtot.Dimension(nel);
  felem_centers.Dimension(nel,NumSD());
  fnormals.Dimension(nel, NumSD());

  fcheckflag = 0;
  flocflag = 0;
  flocflagtot = 0;
  felem_centers = 0.0;
  fnormals = 0.0;
  
  if (fCheck == 1)
  {
    /*set check localization flags*/
    Top();
    while(NextElement())
    {
      int elem = CurrElementNumber();
      int cnt = 0;
      bool match=false;
      while (cnt < fBlockList.Length() && !match)
      {
	if (fBlockList[cnt]+1 > fBlockData.Length()) {
	  cout<<"\nUpLagMF::Initialize:Block number exceeds group dimension.";
	  throw ExceptionT::kGeneralFail;
	}
	const StringT& A = fBlockData[fBlockList[cnt]].ID();
	const StringT& B = ElementBlockID(elem);
	if (strlen(A) == strlen(B) && strncmp(A,B,strlen(A)) == 0) {
	  match = true;
	}
	cnt ++;
      }
      if (match)
	fcheckflag[elem] = 1;
    }
  }
    
  /*dimension workspace for bulk quantities*/
  fEshelby.Dimension(NumSD());
  fC.Dimension(NumSD());
  fBodyForce.Dimension(NumSD()*NumElementNodes());
  fip_body.Dimension(NumSD());
  fVel.Dimension(NumSD());
  fAcc.Dimension(NumSD());
  fGradVel.Dimension(NumSD());
#endif
}

void UpLagMF::SetGlobalShape(void)
{
  UpdatedLagrangianT::SetGlobalShape();
}

/***************************outputs managers***********************************/
/* register self for output */
void UpLagMF::RegisterOutput(void)
{
  /* inherited */
  SolidElementT::RegisterOutput();

  ArrayT<StringT> n_labels(4*NumSD());
  ArrayT<StringT> e_labels;
  
  StringT mf_label = "mF";
  StringT mfd_label = "mF_dissip";
  StringT mfdd_label = "mF_dynamic";
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
void UpLagMF::WriteOutput(void)
{
  /* inherited */
  SolidElementT::WriteOutput();
  
  /* calculate output values */
  dArray2DT n_values; 
  dArray2DT e_values;
  
  if (fCheck==1)
    WriteLocalize(flocflagtot, felem_centers, fnormals);
  MapOutput();
  ComputeMatForce(n_values);

  /* send to output */
  const CommunicatorT& comm = ElementSupport().Communicator();
  if (comm.Size() == 1)
     WriteSummary(n_values);

  ElementSupport().WriteOutput(fMatForceOutputID, n_values, e_values);
}

void UpLagMF::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
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
  UpLagMF* non_const_this = (UpLagMF*) this; 
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

GlobalT::RelaxCodeT UpLagMF::RelaxSystem(void)
{
  const char caller[] = "UpLagMF::RelaxSystem";
 
  /*inherited function*/
  GlobalT::RelaxCodeT relax = UpdatedLagrangianT::RelaxSystem();
  ostream& out = ElementSupport().Output();
  out << "\nRelaxation: Localization Check\n";

  flocflag = 0;
  Top();
  while (NextElement() && fCheck == 1)
  {
    int elem = CurrElementNumber();
    if (fcheckflag[elem] == 1)
    {
      SetGlobalShape();
      SetLocalX(fLocInitCoords);
 
      fCurrShapes->TopIP();
      while (fCurrShapes->NextIP() && flocflag[elem] == 0)
      {
	const dSymMatrixT& stress = fCurrMaterial->s_ij();
	const dMatrixT& modulus = fCurrMaterial->c_ijkl();
	
	int loc = CheckLocalizeFS(stress, modulus,fLocInitCoords);
       	if (loc == 1)
       	{
	  out << "Localization detected in element " << elem 
	      << " IP " << CurrIP() << endl;
	  flocflag[elem] = 1;
	  flocflagtot[elem] = 1;
	  felem_centers.SetRow(elem, LocalizedElemCenter());
	  fnormals.SetRow(elem, LocalizedNormal());
	}
      }
    }
  }
  return(GlobalT::MaxPrecedence(relax, GlobalT::kNoRelax));
}
/**********************Material Force Driver**************/
/*driver*/
void UpLagMF::ComputeMatForce(dArray2DT& output)
{
  const char caller[] = "UpLagMF::ComputeMatForce";

  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nen = NumElementNodes();
  int nmf = nnd*NumSD();

  /*dimension output array and workspace*/
  output.Dimension(nnd,4*NumSD());
  output = 0.0;
  
  const dArray2DT& disp = Field()[0];
  if (disp.MajorDim() != output.MajorDim()) throw ExceptionT::kGeneralFail;

  fMatForce.Dimension(nmf);
  fDissipForce.Dimension(nmf);
  fDynForce.Dimension(nmf);
  felem_rhs.Dimension(NumSD()*nen);  

  fMatForce = 0.0;
  fDissipForce = 0.0;
  fDynForce = 0.0;

  /*if internal dissipation vars exists, extrapolate from element ip to nodes*/
  if (fhas_dissipation) {
    /*assume all materials within elem group have the same internal dissipation variables*/
    const MaterialListT& matlist = ContinuumElementT::MaterialsList();
    if (matlist.Length() > 1){
      cout << "\nUpLagMF::ComputeMatForce: Currently inelastic material force class requires that each element group has only one element block.";
      throw ExceptionT::kGeneralFail;
    }

    ContinuumMaterialT* pmat0 = (*fMaterialList)[ElementCard(0).MaterialNumber()];
    fCurrFSMat = dynamic_cast<FSSolidMatT*>(pmat0);
    if (!fCurrFSMat) throw ExceptionT::kGeneralFail;

    fInternalDOF = fCurrFSMat->InternalDOF();
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
  int analysiscode = ElementSupport().FEManager().Analysis();
  if (analysiscode ==  GlobalT::kLinExpDynamic  ||
      analysiscode == GlobalT::kNLExpDynamic    ||
      analysiscode == GlobalT::kVarNodeNLExpDyn ||
      analysiscode == GlobalT::kPML)
    fdynamic = true;

  /*evaluate volume contributions to material and dissipation force*/
  Top();
  while (NextElement())
  {
    int elem = CurrElementNumber();
    if (flocflag[elem] == 0)
    {
      ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
      fCurrFSMat = dynamic_cast<FSSolidMatT*>(pmat);
      if (!fCurrFSMat) ExceptionT::GeneralFail(caller);
      
      /*Set Global Shape Functions for current element*/
      SetGlobalShape();

      if (fdynamic)
      {
	SetLocalU(fLocAcc);
	SetLocalU(fLocVel);
	SetLocalU(fLocDisp);
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
    else cout << "\nElement "<<CurrElementNumber()<<" localized, excluded from material force calculations.\n";
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

  for (int i = 0; i<nnd; i++)
  {
    for (int j = 0; j<NumSD(); j++)
    {
      *pout_disp++ = disp[i*NumSD()+j];

      /*material force set to zero on external boundary*/
      if (fExclude[i] == 1)
      {
	    *pout_force++ = 0.0;
	    *pout_dissip++ = 0.0;
	    *pout_dyn++ = 0.0;
	    pmat_force++;
	    pmat_fdissip++;
	    pmat_fdyn++;
      }
      else
      {
	    *pout_force++ = (*pmat_force++) + (*pmat_fdissip) + (*pmat_fdyn);
	    *pout_dissip++ = (*pmat_fdissip++);
	    *pout_dyn++ = (*pmat_fdyn++);
      }
    }
    pout_disp += 3*NumSD();
    pout_force += 3*NumSD();
    pout_dissip += 3*NumSD();
    pout_dyn += 3*NumSD();
  }
}

void UpLagMF::MatForceVolMech(dArrayT& elem_val)
{
  const char caller[] = "UpLagMF::MatForceVolMech";
  int nen = NumElementNodes();
  int elem = CurrElementNumber();
  elem_val = 0.0;
  
  /*get density*/
  double density = fCurrFSMat->Density();
  
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
    double energy = fCurrFSMat->StrainEnergyDensity();
 
    const dSymMatrixT& S = fCurrFSMat->S_IJ();
    const dMatrixT& F = fCurrFSMat->F_mechanical();
    fCurrFSMat->Stretch(fC);

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
      fEshelby(0,0) = fC(0,0)*S(0,0) + fC(0,1)*S(1,0)- energy;
      fEshelby(0,1) = fC(0,0)*S(0,1) + fC(0,1)*S(1,1);
      fEshelby(1,0) = fC(1,0)*S(0,0) + fC(1,1)*S(1,0);
      fEshelby(1,1) = fC(1,0)*S(0,1) + fC(1,1)*S(1,1) - energy;


      if (fdynamic)
      {
	fShapes->InterpolateU(fLocVel, fVel);
	fEshelby(0,0) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]);
	fEshelby(1,1) -= 0.5*density*(fVel[0]*fVel[0]+fVel[1]*fVel[1]);
      }

      if (elem == 0 && 0)
      {
	cout << "\nstress: "<<S;
	cout << "\nenergy: "<<energy;
	cout << "\nfVel: "<<fVel;
	cout << "\nEshelby: "<< fEshelby;
      }    

      const double* pDQaX = DQa(0); 
      const double* pDQaY = DQa(1);
      
      for (int j = 0; j<nen; j++)
      {
	/*add negative of Eshelby stress and body force contribution*/
       	*(pforce++) += (fEshelby(0,0)*(*pDQaX) + fEshelby(0,1)*(*pDQaY)
	  +(F(0,0)*fip_body[0]+F(1,0)*fip_body[1])*(*pQa))*(*jac)*(*weight);
	*(pforce++) += (fEshelby(1,0)*(*pDQaX++) + fEshelby(1,1)*(*pDQaY++)
	  +(F(0,1)*fip_body[0]+F(1,1)*fip_body[1])*(*pQa++))*(*jac)*(*weight); 
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
      fEshelby(0,0) = fC[0]*S[0] + fC[5]*S[5] + fC[4]*S[4] - energy;
      fEshelby(1,1) = fC[5]*S[5] + fC[1]*S[1] + fC[3]*S[3] - energy;
      fEshelby(2,2) = fC[4]*S[4] + fC[3]*S[3] + fC[2]*S[2] - energy;

      fEshelby(1,0) = fC[5]*S[0] + fC[1]*S[5] + fC[3]*S[4];
      fEshelby(2,0) = fC[4]*S[0] + fC[3]*S[5] + fC[2]*S[4];
      fEshelby(0,1) = fC[0]*S[5] + fC[5]*S[1] + fC[4]*S[3];
      fEshelby(2,1) = fC[4]*S[5] + fC[3]*S[1] + fC[2]*S[3];
      fEshelby(0,2) = fC[0]*S[4] + fC[5]*S[3] + fC[4]*S[2];
      fEshelby(1,2) = fC[5]*S[4] + fC[1]*S[3] + fC[3]*S[2];

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
	*(pforce++) += (fEshelby(0,0)*(*pDQaX) + fEshelby(0,1)*(*pDQaY) 
	  + fEshelby(0,2)*(*pDQaZ) 
	  +(F(0,0)*fip_body[0] + F(1,0)*fip_body[1] + F(2,0)*fip_body[2])
	  *(*pQa) )*(*jac)*(*weight);
	
    	*(pforce++) += (fEshelby(1,0)*(*pDQaX) + fEshelby(1,1)*(*pDQaY) 
          + fEshelby(1,2)*(*pDQaZ)
	  +(F(0,1)*fip_body[0] + F(1,1)*fip_body[1] + F(2,1)*fip_body[2])
	  *(*pQa) )*(*jac)*(*weight);
	
	*(pforce++) += (fEshelby(2,0)*(*pDQaX++) + fEshelby(2,1)*(*pDQaY++) 
          + fEshelby(2,2)*(*pDQaZ++) 
	  +(F(0,2)*fip_body[0]+F(1,2)*fip_body[1]+F(2,2)*fip_body[2])
	  *(*pQa++))*(*jac)*(*weight);
      }
    }
    weight++;
    jac++;
  }
  if (elem == 0 && 0)      cout<< "\nVolElem: "<<elem_val;
}

void UpLagMF::MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch)
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
      const dArrayT& internalstress = fCurrFSMat->InternalStressVars();
      const double* pstress = internalstress.Pointer();
      double xval = ScalarProduct(pstress, pGradX, fInternalDOF);
      double yval = ScalarProduct(pstress, pGradY, fInternalDOF);

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
      const dArrayT& internalstress = fCurrFSMat->InternalStressVars();
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

void UpLagMF::MatForceDynamic(dArrayT& elem_val)
{
  const char caller[] = "UpLagMF::MatForceDynamic";
  int nen = NumElementNodes();
  int elem = CurrElementNumber();
  elem_val = 0.0;  

  double density = fCurrFSMat->Density();

  /*intialize shape function data*/
  const double* jac = fShapes->IPDets();
  const double* weight = fShapes->IPWeights();
  if (elem == 0 && 0) {
    cout << "\nAcc: "<<fLocAcc;
    cout << "\nVel: "<<fLocVel;
    cout <<"\n fDisp: "<<fLocDisp;
  }
  fShapes->TopIP();
  while(fShapes->NextIP())
  {
    /*get shape function and derivatives at integration point*/
    const double* pQa = fShapes->IPShapeX();
    const double* pQaU = fShapes->IPShapeU();
    const dArray2DT& DQa = fShapes->Derivatives_X();

    /*integration point values*/
    const dMatrixT& F = fCurrFSMat->F_mechanical();
    fShapes->GradU(fLocVel,fGradVel); 
    fShapes->InterpolateU(fLocVel, fVel);
    fShapes->InterpolateU(fLocAcc, fAcc);
    if (elem == 0 && 0)
    {
      cout << "\nelem "<<elem<<" ip "<<CurrIP()<<endl; 
      cout << "\nfGradVel: "<<fGradVel;
      cout <<"\n F: "<< F;
      cout <<"\n fVel: "<<fVel;
      cout <<"\n fAcc: "<<fAcc;
    }
    double* pelem_val = elem_val.Pointer();
    if (NumSD() ==2)
    {
      for (int i = 0; i<nen; i++)
      {
	double xval = density*(-fGradVel[0]*fVel[0]-fGradVel[1]*fVel[1]
			       +F[0]*fAcc[0]+F[1]*fAcc[1]);
	double yval = density*(-fGradVel[2]*fVel[0]-fGradVel[3]*fVel[1]
			       +F[2]*fAcc[0]+F[3]*fAcc[1]);
    	*pelem_val++ += xval*(*pQa)*(*jac)*(*weight);
	*pelem_val++ += yval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    else if (NumSD() == 3)
    {
      for (int i = 0; i<nen; i++)
      {
	double xval = density*(-fGradVel(0,0)*fVel[0]-fGradVel(1,0)*fVel[1]
			       -fGradVel(2,0)*fVel[2]+F(0,0)*fAcc[0]
			       +F(1,0)*fAcc[1]+F(2,0)*fAcc[2]);
	double yval = density*(-fGradVel(0,1)*fVel[0]-fGradVel(1,1)*fVel[1]
			       -fGradVel(2,1)*fVel[2]+F(0,1)*fAcc[0]
			       +F(1,1)*fAcc[1]+F(2,1)*fAcc[2]);
	double zval = density*(-fGradVel(0,2)*fVel[0]-fGradVel(1,2)*fVel[1]
			       -fGradVel(2,2)*fVel[2] + F(0,2)*fAcc[0]
			       +F(1,2)*fAcc[1]+F(2,2)*fAcc[2]); 
    	*pelem_val++ += xval*(*pQa)*(*jac)*(*weight);
	*pelem_val++ += yval*(*pQa)*(*jac)*(*weight);      
	*pelem_val++ += zval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    jac++;
    weight++;
  }
  if (elem == 0 && 0)  cout<<"\nDynnElem: "<<elem_val;
}

void UpLagMF::MatForceSurfMech(dArrayT& global_val)
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

      if (flocflag[elem]==0)
      {
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
	    double mat_tract = -(1+fgradU[0])*fip_tract[0]-fgradU[1]*fip_tract[1];
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
	    double mat_tract1 = -(1+fgradU[0])*fip_tract[0]-fgradU[1]*fip_tract[1]+fgradU[2]*fip_tract[2];
	    double mat_tract2 = -fgradU[3]*fip_tract[0]-(1+fgradU[4])*fip_tract[1]+fgradU[5]*fip_tract[2];
	    
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
}

/****************utitlity functions******************************/
void UpLagMF::Extrapolate(void)
{
  const char caller[] = "UpLagMF::Extrapolate";   
  
  Top();
  while (NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrFSMat = dynamic_cast<FSSolidMatT*>(pmat);
    if (!fCurrFSMat) throw ExceptionT::kGeneralFail;
    UpdatedLagrangianT::SetGlobalShape();

    felem_val = 0.0;
    felem_mass = 0.0;

    const double* jac = fShapes->IPDets();;
    const double* weight = fShapes->IPWeights();
    fShapes->TopIP();
    while(fShapes->NextIP())
    {
	    const dArrayT& internalstrains = fCurrFSMat->InternalStrainVars();
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
}

