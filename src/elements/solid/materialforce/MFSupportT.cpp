/* $Id: MFSupportT.cpp,v 1.11 2005-03-02 17:41:47 paklein Exp $ */
#include "MFSupportT.h"

#include "dArrayT.h"
#include "dArray2DT.h"
#include "ElementSupportT.h"
#include "OutputSetT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
using namespace Tahoe;

/* constructor */
MFSupportT::MFSupportT(const ElementSupportT& support):
	fSupport(support),
	fNumGroupNodes(0),
	fMatForceOutputID(-1),
	fOutputSet(NULL),
	fhas_dissipation(false),
	fopen(false)
{
#if 0
  ifstreamT& in = fSupport.Input();
  ostream&  out = fSupport.Output();

  /*does constitutive model have dissipation variables*/
  in >> fhas_dissipation;
  out << "\nMaterial has dissipation: "<<fhas_dissipation;
  if (fhas_dissipation > 1  && fhas_dissipation < 0)
  {
    cout << "\nMFSupportT:: MFSupportT invalid input for dissipation flag: ";
    throw ExceptionT::kBadInputValue;
  }

  /*read in node sets over which material forces are summed*/
  ModelManagerT& model = fSupport.Model();
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
      cout << "\nMFSupportT::MFSupportT:  Node set " << name << " is undefined: ";
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

  /*boundary node set*/
  StringT name;
  in >> name;

  int index = model.NodeSetIndex(name);
  if (index < 0)
  {
    cout << "\nMFSupport::MFSupportT: Node set "<< name << " is undefinded: ";
    fBoundID = "0";
  }
  else
    fBoundID = nsetIDs[index];

  /*initialize fio boolean*/
  fopen = false;
#endif
}
	
MFSupportT::~MFSupportT(void)
{
  delete fOutputSet;
}

/* set nodes over which material force is calculated and nodes on the boundary */
void MFSupportT::SetNodes(const ArrayT<StringT>& mf_nodes, const ArrayT<StringT>& boundary_nodes)
{
	fNID = mf_nodes;
	fBoundID = boundary_nodes;
}

/***************************I/0 functions***********************************/
/* register self for output */
void MFSupportT::MapOutput(void)
{
	/* check */
	if (!fOutputSet) ExceptionT::GeneralFail("MFSupportT::MapOutput",
		"no output set");

  ModelManagerT& model = fSupport.ModelManager();  
  int nnd = model.NumNodes();
  
  const iArrayT& nodes_used = fOutputSet->NodesUsed();

  fNumGroupNodes = nodes_used.Length();
  
  /*map ordering*/
  fMap.Dimension(nnd);
  fMap = -1;
  for (int i = 0; i<fNumGroupNodes; i++)
  {
    fMap[nodes_used[i]]=i;
  }

  fExclude.Dimension(fNumGroupNodes);
  fExclude = 0;

	if (fBoundID.Length() > 0) {
		iArrayT bound_nodes;
		model.ManyNodeSets(fBoundID, bound_nodes);
		for (int i = 0; i< bound_nodes.Length(); i++)
			fExclude[fMap[bound_nodes[i]]] = 1;
	}
}

void MFSupportT::WriteSummary(dArray2DT& output)
{
	/* check */
	if (!fOutputSet) ExceptionT::GeneralFail("MFSupportT::WriteSummary",
		"no output set");

  /*obtain dimensions*/
  int nnd = fNumGroupNodes;
  int nsd = fSupport.NumSD();

  /*write summary of MF results to external file*/  
  ModelManagerT& model = fSupport.ModelManager();
  int numnset = model.NumNodeSets();

  const StringT& input_file = fSupport.InputFile();
  fsummary_file.Root(input_file);
  fsummary_file.Append(".sum");
 
  double time = fSupport.Time();
  int precision = 6;
  int doublewidth = kPrecision+kDoubleExtra;
  int intwidth = kIntWidth;
  if (nsd == 2)
  { 
    double* pFx = output.Pointer(2);
    double* pFy = output.Pointer(3);
    double* pDFx = output.Pointer(4);
    double* pDFy = output.Pointer(5);
    double* pKFx = output.Pointer(6);
    double* pKFy = output.Pointer(7);
 
    /*sum components of material force over a given nodeset*/ 
    double MFx, MFy, DFx, DFy, KFx, KFy;  
    MFx = MFy = DFx = DFy = KFx = KFy = 0.0;
    for (int i = 0; i < fNID.Length(); i++)
    {
      StringT& ID = fNID[i];
      const int nlength = model.NodeSetLength(ID);
      const iArrayT& nset = model.NodeSet(ID);
      for (int j = 0; j<nlength; j++)
      {
        int index = fMap[nset[j]];
        MFx += *(pFx+index*4*nsd);
        MFy += *(pFy+index*4*nsd);
        DFx += *(pDFx+index*4*nsd);
        DFy += *(pDFy+index*4*nsd);
	KFx += *(pKFx+index*4*nsd);
	KFy += *(pKFy+index*4*nsd);
      }
    }
    
    /*find maximum material force*/
    const iArrayT& nodes_used = fOutputSet->NodesUsed();
    double maxFx, maxFy;
    int nFx, nFy;
    nFx = nFy = 0;
    maxFx = maxFy = 0.0;
    for (int i = 0; i<nnd; i++)
    {
      if (fabs(maxFx) < fabs(*(pFx+i*4*nsd))) 
        {maxFx = *(pFx+i*4*nsd); nFx = nodes_used[i];}
      if (fabs(maxFy) < fabs(*(pFy+i*4*nsd))) 
        {maxFy = *(pFy+i*4*nsd); nFy = nodes_used[i];}
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
           <<"\n\t Summed x component of dissipation contribution . . Fk_X"
           <<"\n\t Summed y component of dissipation contribution . . Fk_Y"   
           <<endl<<endl<<endl;

      fout <<setw(intwidth) << "Time" 
           <<setw(doublewidth) << "F_X" 
           <<setw(doublewidth) << "F_Y" 
           <<setw(doublewidth) << "Fd_X" 
           <<setw(doublewidth) << "Fd_Y" 
           <<setw(doublewidth) << "Fk_X" 
           <<setw(doublewidth) << "Fk_Y" 
           <<endl;        
    }

    fout <<setw(intwidth) << time 
         <<setw(doublewidth) << MFx 
         <<setw(doublewidth) << MFy 
         <<setw(doublewidth) << DFx 
         <<setw(doublewidth) << DFy 
         <<setw(doublewidth) << KFx 
         <<setw(doublewidth) << KFy 
         <<endl; 
       
    fout.close();
  }
  else if (nsd == 3)
  { 
    double* pFx = output.Pointer(3);
    double* pFy = output.Pointer(4);
    double* pFz = output.Pointer(5);
    double* pDFx = output.Pointer(6);
    double* pDFy = output.Pointer(7);
    double* pDFz = output.Pointer(8);
     
    double MFx,MFy, MFz, DFx, DFy, DFz;  
    MFx = MFy = MFz = DFx = DFy = DFz = 0.0;
 
    for (int i = 0; i < fNID.Length(); i++)
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

    fout.close();
  }
}

/****************utitlity functions******************************/
void MFSupportT::GatherDisp(const dArray2DT& global_disp, dArrayT& group_disp, const iArrayT& elem_nodes)
{
  int nen = elem_nodes.Length();
  int nsd = fSupport.NumSD();
  int group_index, node;
  
  for (int i = 0; i< nen; i++)
  {  
    node = elem_nodes[i]; 
   group_index=fMap[node]*nsd;  
    for (int j = 0; j<nsd; j++)
      group_disp[group_index+j] = global_disp(node,j);
  }
}

void MFSupportT::AssembleArray(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& elem_nodes)
{
  int nen = elem_nodes.Length();
  int nsd = fSupport.NumSD();
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

void MFSupportT::AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val, const iArrayT& elem_nodes)
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

void MFSupportT::ExtractArray2D(const dArray2DT& global_val, dArray2DT& elem_val, const iArrayT& elem_nodes)
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


double MFSupportT::ScalarProduct(const double* pa, const double* pb, const iArrayT& dims)
{  
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
	    val +=pa[0]*pb[0]+pa[1]*pb[1]+pa[3]*pb[3]+2.0*pa[2]*pb[2];
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

