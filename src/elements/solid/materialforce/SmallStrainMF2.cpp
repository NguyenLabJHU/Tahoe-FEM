/* $Id: SmallStrainMF2.cpp,v 1.1 2003-05-15 05:16:53 thao Exp $ */
#include "SmallStrainMF2.h"

#include "OutputSetT.h"
#include "Traction_CardT.h"
#include "ScheduleT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"
#include "GeometryT.h"
#include "ModelManagerT.h"
#include "SSMatSupportT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

/* materials lists */
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
using namespace Tahoe;

/* constructor */
SmallStrainMF2::SmallStrainMF2(const ElementSupportT& support, const FieldT& field):
  SmallStrainT(support, field),
  fNumGroupNodes(0),
  fMatForceOutputID(-1)
{
  ifstreamT& in = ElementSupport().Input();
  ostream&  out = ElementSupport().Output();

  in >> fhas_dissipation;
  if (fhas_dissipation > 1  && fhas_dissipation < 0)
  {
    cout << "\nSmallStrainMF2:: SmallStrainMF2 invalid input for dissipation flag: ";
    throw ExceptionT::kBadInputValue;
  }
  ModelManagerT& model = ElementSupport().Model();
  const ArrayT<StringT>& nsetIDs = model.NodeSetIDs();
  in >> fnumset;
  fNID.Dimension (fnumset);
  
  for (int i=0; i < fnumset; i++)
  {
    StringT name;
    in >> name;
    int index = model.NodeSetIndex(name);
    if (index < 0) 
    {
      cout << "\nSmallStrainMF22::SmallStrainMF2:  Node set " << name << " is undefined: ";
      throw ExceptionT::kDatabaseFail;
    }
    else
    {
      fNID[i] = nsetIDs[index];
    }
  }
  out << "\n Number of nodesets for summing material force: "<<fnumset;
  
  for (int j = 0; j<fnumset; j++) out << "\n\tNodeset: "<<fNID[j];
  out <<'\n';
  fopen = false;
}
	
SmallStrainMF2::~SmallStrainMF2(void)
{
  delete fOutputSet;
}

void SmallStrainMF2::Initialize(void)
{
  SmallStrainT::Initialize();
  
  fGradU_List.Dimension(NumIP());
  for (int i = 0; i< NumIP(); i++)
    fGradU_List[i].Dimension(NumSD());
}

void SmallStrainMF2::SetGlobalShape(void)
{
  SmallStrainT::SetGlobalShape();
  for (int i = 0; i < NumIP(); i++)
    fShapes->GradU(fLocDisp, fGradU_List[i], i);
}

/***************************outputs managers***********************************/
/* register self for output */
void SmallStrainMF2::RegisterOutput(void)
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
void SmallStrainMF2::WriteOutput(void)
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

void SmallStrainMF2::WriteSummary(dArray2DT& output)
{
  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nsd = NumSD();

  /*write summary of MF results to external file*/  
  ModelManagerT& model = ElementSupport().Model();
  int numnset = model.NumNodeSets();
 
  ifstreamT& in = ElementSupport().Input();
  const StringT& input_file = in.filename();
  fsummary_file.Root(input_file);
  fsummary_file.Append(".sum");
 
  double time = ElementSupport().Time();
  int precision = 6;
  int doublewidth = kPrecision+kDoubleExtra;
  int intwidth = kIntWidth;
  if (nsd == 2)
  { 
    double* pFx = output.Pointer();
    double* pFy = output.Pointer(1);
    double* pDFx = output.Pointer(2);
    double* pDFy = output.Pointer(3);
   
    /*sum components of material force over a given nodeset*/ 
    double MFx, MFy, DFx, DFy;  
    MFx = MFy = DFx = DFy = 0.0;
    for (int i = 0; i<fnumset; i++)
    {
      StringT& ID = fNID[i];
      const int nlength = model.NodeSetLength(ID);
      const iArrayT& nset = model.NodeSet(ID);
      for (int j = 0; j<nlength; j++)
      {
        int index = fMap[nset[j]];
        MFx += *(pFx+index*3*nsd);
        MFy += *(pFy+index*3*nsd);
        DFx += *(pDFx+index*3*nsd);
        DFy += *(pDFy+index*3*nsd);
      }
    }
    
    /*find maximum material force*/
    const iArrayT& nodes_used = fOutputSet->NodesUsed();
    double maxFx, maxFy;
    int nFx, nFy;
    maxFx = maxFy = 0.0;
    for (int i = 0; i<nnd; i++)
    {
      if (fabs(maxFx) < fabs(*(pFx+i*3*nsd))) 
        {maxFx = *(pFx+i*3*nsd); nFx = nodes_used[i];}
      if (fabs(maxFy) < fabs(*(pFy+i*3*nsd))) 
        {maxFy = *(pFy+i*3*nsd); nFy = nodes_used[i];}
    }
    
    /*write summary output file*/
    if (fopen)
    {
      fout.open_append(fsummary_file);
    }  
    else
    {
      fout.open(fsummary_file);
      fopen = true;
      
      fout <<"\nSummary of material force results:" 
	   <<"\n\t Summed x component . . . . . . . . . . . . . . . . F_X"
           <<"\n\t Summed y component . . . . . . . . . . . . . . . . F_Y"
           <<"\n\t Summed x component of dissipation contribution . . Fd_X"
           <<"\n\t Summed y component of dissipation contribution . . Fd_Y"   
           <<"\n\t Maximum x component . . . . . . . . . . . . . . . . maxF_X"
           <<"\n\t Maximum y component . . . . . . . . . . . . . . . . maxF_Y"
           <<endl<<endl<<endl;

      fout <<setw(intwidth) << "Time" 
           <<setw(doublewidth) << "F_X" 
           <<setw(doublewidth) << "F_Y" 
           <<setw(doublewidth) << "Fd_X" 
           <<setw(doublewidth) << "Fd_Y" 
           <<setw(doublewidth) << "Node" 
           <<setw(doublewidth) << "(F_X)max" 
           <<setw(doublewidth) << "Node" 
           <<setw(doublewidth) << "(F_Y)max"
           <<endl;        
    }

    fout <<setw(intwidth) << time 
         <<setw(doublewidth) << MFx 
         <<setw(doublewidth) << MFy 
         <<setw(doublewidth) << DFx 
         <<setw(doublewidth) << DFy 
         <<setw(intwidth+2) << nFx+1 
         <<setw(doublewidth+1) << maxFx 
         <<setw(intwidth+6) << nFy+1
         <<setw(doublewidth+2) << maxFy
         <<endl;        
  }
  else if (nsd == 3)
  { 
    double* pFx = output.Pointer();
    double* pFy = output.Pointer(1);
    double* pFz = output.Pointer(2);
    double* pDFx = output.Pointer(3);
    double* pDFy = output.Pointer(4);
    double* pDFz = output.Pointer(5);
     
    double MFx,MFy, MFz, DFx, DFy, DFz;  
    MFx = MFy = MFz = DFx = DFy = DFz = 0.0;
 
    for (int i = 0; i<fnumset; i++)
    {
      StringT& ID = fNID[i];
      const int nlength = model.NodeSetLength(ID);
      const iArrayT& nset = model.NodeSet(ID);
      for (int j = 0; j<nlength; j++)
      {
        int index = fMap[nset[j]];
        MFx += *(pFx+index*3*nsd);
        MFy += *(pFy+index*3*nsd);
        MFz += *(pFz+index*3*nsd);
        DFx += *(pDFx+index*3*nsd);
        DFy += *(pDFy+index*3*nsd);
        DFz += *(pDFz+index*3*nsd);
      }
    }
    /*find maximum material force*/
    const iArrayT& nodes_used = fOutputSet->NodesUsed();
    double maxFx, maxFy, maxFz;
    int nFx, nFy, nFz;
    maxFx = maxFy = maxFz = 0.0;
    for (int i = 0; i<nnd; i++)
    {
      if (maxFx < *(pFx+i*3*nsd)) {maxFx = *(pFx+i*3*nsd); nFx = nodes_used[i];}
      if (maxFy < *(pFy+i*3*nsd)) {maxFy = *(pFy+i*3*nsd); nFy = nodes_used[i];}
      if (maxFz < *(pFz+i*3*nsd)) {maxFz = *(pFz+i*3*nsd); nFz = nodes_used[i];}
    }
    /*write summary output file*/
    if (fopen)
    {
      fout.open_append(fsummary_file);
    }  
    else
    {
      fout.open(fsummary_file);
      fopen = true;
      
      fout <<"\nSummary of material force results:"
           <<"\n\t Summed x component . . . . . . . . . . . . . . . . F_X"
           <<"\n\t Summed y component . . . . . . . . . . . . . . . . F_Y"
           <<"\n\t Summed z component . . . . . . . . . . . . . . . . F_Z"
           <<"\n\t Summed x component of dissipation contribution . . Fd_X"
           <<"\n\t Summed y component of dissipation contribution . . Fd_Y"   
           <<"\n\t Summed z component of dissipation contribution . . Fd_Z"   
           <<"\n\t Maximum x component . . . . . . . . . . . . . . . . maxF_X"
           <<"\n\t Maximum y component . . . . . . . . . . . . . . . . maxF_Y"
           <<"\n\t Maximum z component . . . . . . . . . . . . . . . . maxF_Z"
           <<endl<<endl<<endl;

      fout <<setw(intwidth) << "Time" 
           <<setw(doublewidth) << "F_X" 
           <<setw(doublewidth) << "F_Y" 
           <<setw(doublewidth) << "F_Z" 
           <<setw(doublewidth) << "Fd_X" 
           <<setw(doublewidth) << "Fd_Y" 
           <<setw(doublewidth) << "Fd_Z" 
           <<setw(doublewidth) << "Node" 
           <<setw(doublewidth) << "(F_X)max" 
           <<setw(doublewidth) << "Node" 
           <<setw(doublewidth) << "(F_Y)max"
           <<setw(doublewidth) << "Node" 
           <<setw(doublewidth) << "(F_Z)max"
           <<endl;        
    }

    fout <<setw(intwidth) << time 
         <<setw(doublewidth) << MFx 
         <<setw(doublewidth) << MFy 
         <<setw(doublewidth) << MFz 
         <<setw(doublewidth) << DFx 
         <<setw(doublewidth) << DFy 
         <<setw(doublewidth) << DFz 
         <<setw(intwidth+2) << nFx+1 
         <<setw(doublewidth+1) << maxFx 
         <<setw(intwidth+6) << nFy+1
         <<setw(doublewidth+2) << maxFy
         <<setw(intwidth+6) << nFz+1
         <<setw(doublewidth+2) << maxFz
         <<endl;        
  }
}

/**********************evaluates material force**************/
/*driver*/
void SmallStrainMF2::ComputeMatForce(dArray2DT& output)
{
  const char caller[] = "SmallStrainT::ComputeMatForce";

  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nen = NumElementNodes();
  int nsd = NumSD();
  int nmf = nnd*nsd;

  /*dimension output array and workspace*/
  output.Dimension(nnd,3*nsd);

  const dArray2DT& disp = Field()[0];
  if (disp.MajorDim() != output.MajorDim()) throw ExceptionT::kGeneralFail;

  dArrayT mat_force(nmf);
  dArrayT mat_fdissip(nmf);
  mat_force = 0;
  mat_fdissip = 0;
  
  /*if internal dissipation vars exists, extrapolate from element ip to nodes*/
  if (fhas_dissipation) Extrapolate();

  /*evaluate volume contributions to material force*/
  Top();
  dArrayT elem_data(nsd*nen);  
  while (NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
    if (!fCurrSSMat) ExceptionT::GeneralFail(caller);
    
    /*Set Global Shape Functions for current element*/
    SetGlobalShape();
    MatForceVolMech(elem_data);
    AssembleArray(elem_data, mat_force, CurrentElement().NodesX());
    if (fhas_dissipation) 
    {
      felem_val.Free();
      ExtractArray2D(fglobal_val, felem_val, CurrentElement().NodesX());
      MatForceDissip(elem_data, felem_val);
      AssembleArray(elem_data, mat_fdissip, CurrentElement().NodesX());
    }
  }

  /*add surface contribution*/
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
      /*material force set to zero for kinematically constrained nodes*/
      if(eqno[i*nsd+j] < 1)
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
  WriteSummary(output);
}

void SmallStrainMF2::MatForceVolMech(dArrayT& elem_val)
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
  double density = fCurrSSMat->Density();
  
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
    /*gather material data*/
    double energy = fCurrSSMat->StrainEnergyDensity();
    const dSymMatrixT& sig = fCurrSSMat->s_ij();
    const dMatrixT& gradU = DisplacementGradient();

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

      /*form negative of Eshelby stress -SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
      nEshelby(0,0) = gradU(0,0)*sig(0,0)+gradU(1,0)*sig(1,0) - energy;
      nEshelby(0,1) = gradU(0,0)*sig(0,1)+gradU(1,0)*sig(1,1);
      nEshelby(1,0) = gradU(0,1)*sig(0,0)+gradU(1,1)*sig(1,0);
      nEshelby(1,1) = gradU(0,1)*sig(0,1)+gradU(1,1)*sig(1,1) - energy;

      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      
      for (int j = 0; j<nen; j++)
      {
	/*add nEshelby volume integral contribution*/
       	*(pforce++) += (nEshelby[0]*(*pDQaX)+nEshelby[2]*(*pDQaY)
	       +(gradU[0]*ip_body[0]+gradU[1]*ip_body[1])*(*pQa))*(*jac)*(*weight);
	*(pforce++) += (nEshelby[1]*(*pDQaX++)+nEshelby[3]*(*pDQaY++)
	       +(gradU[2]*ip_body[0]+gradU[3]*ip_body[1])*(*pQa++))*(*jac)*(*weight); 

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
	        
      /*form negative of Eshelby stress -SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
      
      nEshelby(0,0)=gradU(0,0)*sig(0,0)+gradU(1,0)*sig(1,0)+gradU(2,0)*sig(2,0)-energy;
      nEshelby(0,1)=gradU(0,0)*sig(0,1)+gradU(1,0)*sig(1,1)+gradU(2,0)*sig(2,1);
      nEshelby(0,2)=gradU(0,0)*sig(0,2)+gradU(1,0)*sig(1,2)+gradU(2,0)*sig(2,2);      
      nEshelby(1,0)=gradU(0,1)*sig(0,0)+gradU(1,1)*sig(1,0)+gradU(2,1)*sig(2,0);      
      nEshelby(1,1)=gradU(0,1)*sig(0,1)+gradU(1,1)*sig(1,1)+gradU(2,1)*sig(2,1)-energy;      
      nEshelby(1,2)=gradU(0,1)*sig(0,2)+gradU(1,1)*sig(1,2)+gradU(2,1)*sig(2,2);      
      nEshelby(2,0)=gradU(0,2)*sig(0,0)+gradU(1,2)*sig(1,0)+gradU(2,2)*sig(2,0);     
      nEshelby(2,1)=gradU(0,2)*sig(0,1)+gradU(1,2)*sig(1,1)+gradU(2,2)*sig(2,1);      
      nEshelby(2,2)=gradU(0,2)*sig(0,2)+gradU(1,2)*sig(1,2)+gradU(2,2)*sig(2,2)-energy;

      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      double* pDQaZ = DQa(2);
      
      for (int j = 0; j<nen; j++)
      {
	/*add Eshelby volume integral contribution*/
	*(pforce++)+=(nEshelby[0]*(*pDQaX)+nEshelby[3]*(*pDQaY)+nEshelby[6]*(*pDQaZ) 
	 +(gradU[0]*ip_body[0]+gradU[1]*ip_body[1]+gradU[2]*ip_body[2])*(*pQa) )*(*jac)*(*weight);

	*(pforce++) += (nEshelby[1]*(*pDQaX)+nEshelby[4]*(*pDQaY)+nEshelby[7]*(*pDQaZ)
	 +(gradU[3]*ip_body[0]+gradU[4]*ip_body[1]+gradU[5]*ip_body[2])*(*pQa) )*(*jac)*(*weight);

	*(pforce++) += (nEshelby[2]*(*pDQaX++)+nEshelby[5]*(*pDQaY++)+nEshelby[8]*(*pDQaZ++)
	 +(gradU[6]*ip_body[0]+gradU[7]*ip_body[1]+gradU[8]*ip_body[2])*(*pQa++))*(*jac)*(*weight);
      }
    }
    weight++;
    jac++;
  }
}

void SmallStrainMF2::MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch)
{
  bool print = false;
  int pos = fElementCards.Position();
  if (pos == 1 && 0) 
    print = true;
  
  /*obtain dimensions*/
  int nen = NumElementNodes();
  int nsd = NumSD();
  int nip = fShapes->NumIP();
  int varsets = finternaldof.Length();

  fgrad_intstrain.Dimension(nsd,fnumval);
  elem_val = 0;

  const double* jac = fShapes->IPDets();
  const double* weight = fShapes->IPWeights();
  fShapes->TopIP();

  while(fShapes->NextIP())
  { 
    const double* pQa = fShapes->IPShapeX();
    const dArray2DT& DQa = fShapes->Derivatives_X();
    if (nsd ==2)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);

      /*Interpolate grad of inverse inelastic stretch to ip*/
      fgrad_intstrain = 0;
      double* pGradX = fgrad_intstrain(0);
      double* pGradY = fgrad_intstrain(1);
      for (int i = 0; i<nen; i++)
      {
        for (int cnt = 0; cnt < fnumval; cnt++)
        {
          pGradX[cnt] += (*pDQaX)*internalstretch(i,cnt);
          pGradY[cnt] += (*pDQaY)*internalstretch(i,cnt);
	}
    	pDQaX++;
    	pDQaY++;
      }

      /*integrate material force*/
      const dArrayT& internalstress = fCurrSSMat->InternalStressVars();
      double* pstress = internalstress.Pointer();
      double xval = ScalarProduct(pstress, pGradX, finternaldof);
      double yval = ScalarProduct(pstress, pGradY, finternaldof);

      double* pelem_val = elem_val.Pointer();
      for (int i = 0; i<nen; i++)
      {
	(*pelem_val++) += xval*(*pQa)*(*jac)*(*weight);
	(*pelem_val++) += yval*(*pQa++)*(*jac)*(*weight);      
      }
    }
    else if (nsd ==3)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);
      const double* pDQaZ = DQa(2);

      /*Interpolate grad of iverse inelastic stretch to ip*/
      fgrad_intstrain = 0;
      double* pGradX = fgrad_intstrain(0);
      double* pGradY = fgrad_intstrain(1);
      double* pGradZ = fgrad_intstrain(2);
      for (int i = 0; i<nen; i++)
      {	
        for (int cnt = 0; cnt < fnumval; cnt++)
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
      double* pstress = internalstress.Pointer();
      double xval = ScalarProduct(pstress, pGradX, finternaldof);
      double yval = ScalarProduct(pstress, pGradY, finternaldof);
      double zval = ScalarProduct(pstress, pGradZ, finternaldof);

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

void SmallStrainMF2::MatForceSurfMech(dArrayT& global_val)
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
      ContinuumMaterialT* pmat = 
(*fMaterialList)[elem_card.MaterialNumber()];
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
        pval.Allocate(nsd);
        pval = 0;
      }
      dArrayT ExtrapEnergy(nen);
      ExtrapEnergy = 0;
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
        nval.Allocate(nsd);
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
      dArrayT ip_tract(nsd);
      dArrayT global_tract(nsd);
      dArrayT ip_eshelby(nsd);     
      dMatrixT ip_gradU(nsd);
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
          ip_gradU = 0;
          ip_tract = 0;
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
          ip_gradU = 0;
          ip_tract = 0;
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

/****************utitlity functions******************************/
void SmallStrainMF2::Extrapolate(void)
{
  const char caller[] = "SmallStrainT::Extrapolate";

  /*dimensions*/
  int nen = NumElementNodes();
  int nnd = fNumGroupNodes;
  int nip = fShapes->NumIP();
    
  Top();
  NextElement();
  /*assume all materials within elemblock have the same internal dissipation variables*/
  ContinuumMaterialT* pmat0 = (*fMaterialList)[CurrentElement().MaterialNumber()];
  fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat0);
  if (!fCurrSSMat) ExceptionT::GeneralFail(caller);

  finternaldof = fCurrSSMat->InternalDOF();
  int varsets = finternaldof.Length();
  fnumval = 0;
  for (int i = 0; i<varsets; i++)
    fnumval += finternaldof[i];

  felem_mass.Dimension(nen);
  fglobal_mass.Dimension(nnd);
  fglobal_val.Dimension(nnd, fnumval);
  felem_val.Dimension(nen, fnumval);

  fglobal_val = 0.0;
  fglobal_mass = 0.0;
   
  Top(); 
  while (NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    fCurrSSMat = dynamic_cast<SSSolidMatT*>(pmat);
    if (!fCurrSSMat) ExceptionT::GeneralFail(caller);

    felem_val = 0.0;
    felem_mass = 0.0;

    SetGlobalShape();
    const double* jac = fShapes->IPDets();;
    const double* weight = fShapes->IPWeights();
    fShapes->TopIP();
    while(fShapes->NextIP())
    {
	const dArrayT& internalstrains = fCurrSSMat->InternalStrainVars();
	if (internalstrains.Length()!=fnumval) ExceptionT::GeneralFail(caller);
        const double* pQbU = fShapes->IPShapeU();
	
	for (int i=0; i<nen; i++)
        {
          const double* pQaU = fShapes->IPShapeU();
          for (int j = 0; j<nen; j++)
            felem_mass[i] += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
          for (int cnt = 0; cnt < fnumval; cnt++)
	    felem_val(i,cnt) += (*pQbU)*(*jac)*(*weight)*internalstrains[cnt];
          pQbU++;
        }
        weight++;
        jac++;
    }    
    AssembleArray(felem_mass, fglobal_mass, CurrentElement().NodesX());
    AssembleArray2D(felem_val, fglobal_val, CurrentElement().NodesX());
  }  
  for (int i = 0; i< nnd; i++)
    for (int j = 0; j< fnumval; j++)
        fglobal_val(i,j) /= fglobal_mass[i];
}

void SmallStrainMF2::AssembleArray(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& elem_nodes)
{
  int nen = elem_nodes.Length();
  int nsd = NumSD();
  int index;
  int numval = elem_val.Length();
  
  if (numval != nen*nsd && numval != nen) throw ExceptionT::kGeneralFail;
  
  for (int i = 0; i< nen; i++)
  {  
    if (numval == nen*nsd)
    {
      index=fMap[elem_nodes[i]]*nsd;    
      for (int j = 0; j<nsd; j++)
      	 global_val[index+j] += elem_val[i*nsd+j];
    }
    else
    {
        index = fMap[elem_nodes[i]];
        global_val[index] += elem_val[i];
    }
  }
}

void SmallStrainMF2::AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val, const iArrayT& elem_nodes)
{
  int nen = elem_nodes.Length();
  int numval = elem_val.MinorDim();
  
  if (global_val.MinorDim() != numval) throw ExceptionT::kGeneralFail;
  if (elem_val.MajorDim() != nen) throw ExceptionT::kGeneralFail;

  int index;
  for(int i = 0; i < nen; i++)
  {
    index=fMap[elem_nodes[i]];
    for (int j = 0; j < numval; j++)
      	global_val(index,j) += elem_val(i,j);
  }
}

void SmallStrainMF2::ExtractArray2D(const dArray2DT& global_val, dArray2DT& elem_val, const iArrayT& elem_nodes)
{
  int nen = elem_nodes.Length();
  int numval = global_val.MinorDim();
  elem_val.Dimension(nen, numval);

  int index;
  for (int i = 0; i<nen; i++)
  {
    index = fMap[elem_nodes[i]];
    for (int j = 0; j< numval; j++)
        elem_val(i,j) = global_val(index,j);
  }
}

void SmallStrainMF2::MapOutput(void)
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

double SmallStrainMF2::ScalarProduct(double* pa, double* pb, const iArrayT& dims)
{
  bool print = false;
  int pos = fElementCards.Position();
  if (pos == 1 && 0) 
    print = true;
  
  int varsets = dims.Length();
  double val = 0;
  for (int i = 0; i<varsets; i++)
  {
    int numval = dims[i];
    switch(numval)
    {
      case 1:{ 
	val += pa[0]*pb[0];
	break;}
      case 4:{
	val +=pa[0]*pb[0]+pa[1]*pb[1]+2.0*pa[2]*pb[2];
	break;}
      case 6:{
       	val +=pa[0]*pb[0] + pa[1]*pb[1] + pa[2]*pb[2] +
	 2.0*(pa[3]*pb[3] + pa[4]*pb[4]+ pa[5]*pb[5]);
	break;}
      default:
	throw ExceptionT::kGeneralFail;
    }
    pa+=numval;
    pb+=numval;
  }
  return(val);
} 

