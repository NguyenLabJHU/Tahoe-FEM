/* $Id: SSMFBaseT.cpp,v 1.1 2003-08-10 23:29:12 thao Exp $ */
#include "MFBaseT.h"

#include "OutputSetT.h"
#include "ScheduleT.h"
#include "ShapeFunctionT.h"
#include "GeometryT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
using namespace Tahoe;

/* constructor */
MFBaseT::MFBaseT(const ElementSupportT& support, const FieldT& field):
  fNumGroupNodes(0),
  fMatForceOutputID(-1)
{
  ifstreamT& in = ElementSupport().Input();
  ostream&  out = ElementSupport().Output();

  /*does constitutive model have dissipation variables*/
  in >> fhas_dissipation;
  if (fhas_dissipation > 1  && fhas_dissipation < 0)
  {
    cout << "\nMFBaseT:: MFBaseT invalid input for dissipation flag: ";
    throw ExceptionT::kBadInputValue;
  }

  /*read in node sets over which material forces are summed*/
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
      cout << "\nMFBaseT2::MFBaseT:  Node set " << name << " is undefined: ";
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
  
  /*initialize fio boolean*/
  fopen = false;
  
  /*map output*/
  
}
	
MFBaseT::~MFBaseT(void)
{
  delete fOutputSet;
}

void MFBaseT::MapOutput(void)
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

void MFBaseT::SetGlobalShape(void)
{
  SmallStrainT::SetGlobalShape();
  for (int i = 0; i < NumIP(); i++)
    fShapes->GradU(fLocDisp, fGradU_List[i], i);
}

/***************************outputs managers***********************************/
/* register self for output */
void MFBaseT::RegisterOutput(void)
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
void MFBaseT::WriteOutput(void)
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

void MFBaseT::WriteSummary(dArray2DT& output)
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
void MFBaseT::ComputeMatForce(dArray2DT& output)
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


/****************utitlity functions******************************/
void MFBaseT::Extrapolate(void)
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
  if (!fCurrSSMat) throw ExceptionT::kGeneralFail;

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
    if (!fCurrSSMat) throw ExceptionT::kGeneralFail;

    bool print = false;
    int pos = fElementCards.Position();
    if (pos == 1&&0) 
      print = true;

    felem_val = 0.0;
    felem_mass = 0.0;

    SetGlobalShape();
    const double* jac = fShapes->IPDets();;
    const double* weight = fShapes->IPWeights();
    fShapes->TopIP();
    while(fShapes->NextIP())
    {
	const dArrayT& internalstrains = fCurrSSMat->InternalStrainVars();
	
	if (print)
	{
	  cout<<"\nIP: "<<fShapes->CurrIP();
	  cout<<"\ninternalstrains: "<<internalstrains;
	}

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

void MFBaseT::AssembleArray(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& elem_nodes)
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

void MFBaseT::AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val, const iArrayT& elem_nodes)
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

void MFBaseT::ExtractArray2D(const dArray2DT& global_val, dArray2DT& elem_val, const iArrayT& elem_nodes)
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


double MFBaseT::ScalarProduct(double* pa, double* pb, const iArrayT& dims)
{
  bool print = false;
  int pos = fElementCards.Position();
  if (pos == 1&&0) 
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
        if (print)
	  cout<<"\ncase1: "
	      <<"\npa0: "<<pa[0]
	      <<"\npb0: "<<pb[0]
	      <<"\nval: "<<val;
	break;}
      case 4:{
	val +=pa[0]*pb[0]+pa[1]*pb[1]+2.0*pa[2]*pb[2];
        if (print)
	  cout<<"\ncase2: "
	      <<"\npa0: "<<pa[0]
	      <<"\npb0: "<<pb[0]
	      <<"\nval: "<<val;
	break;}
      case 6:{
       	val +=pa[0]*pb[0] + pa[1]*pb[1] + pa[2]*pb[2] +
	 2.0*(pa[3]*pb[3] + pa[4]*pb[4]+ pa[5]*pb[5]);
        if (print)
	  cout<<"\ncase3: "
	      <<"\npa0: "<<pa[0]
	      <<"\npb0: "<<pb[0]
	      <<"\nval: "<<val;
	break;}
      default:
	throw ExceptionT::kGeneralFail;
    }
    pa+=numval;
    pb+=numval;
  }
  return(val);
} 

