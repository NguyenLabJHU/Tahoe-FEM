/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.cpp,v 1.4 2006-07-12 17:48:40 a-kopacz Exp $ */
/* created: a-kopacz (07/04/2006) */
#include "FluidElementT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "eIntegratorT.h"
#include "ShapeFunctionT.h"
#include "XDOF_ManagerT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* initialize static data */
const int FluidElementT::NumNodalOutputCodes = 4;
static const char* NodalOutputNames[] = {
  "coordinates",
  "velocities",
  "accelerations",
  "pressures"};

const int FluidElementT::NumElementOutputCodes = 0;
static const char* ElementOutputNames[] = {
  "NONE"};

const int FluidElementT::NumStabParamCodes = 1;
static const char* StabParamNames[] = {
  "tau_m_is_tau_c"}; /* T.E. Tezduyar, Stabilized Finite Element Formulations for Incompressible Flow Computations, Adv. Appl. Mech. 28(1991) 1-44. */
  
/* parameters */
const int FluidElementT::kPressureNDOF = 1;

/* constructor */
FluidElementT::FluidElementT(const ElementSupportT& support):
  ContinuumElementT(support),
  /* velocity */
  fLocCurVel(LocalArrayT::kDisp),  
  fLocOldVel(LocalArrayT::kLastDisp),
  /* accelaration */
  fLocCurAcc(LocalArrayT::kVel),
  /* pressure */   
  fLocPrs(LocalArrayT::kUnspecified)
{
  SetName("incompressible_newtonian_fluid_element");
}

/* destructor */
FluidElementT::~FluidElementT(void)
{
  delete fCurrShapes;
  fCurrShapes = NULL;
}

/* compute nodal force */
void FluidElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
  WriteCallLocation("AddNodalForce: not implemented"); //DEBUG  
  //not implemented
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* returns the energy as defined by the derived class types */
double FluidElementT::InternalEnergy(void)
{
  WriteCallLocation("InternalEnergy: not implemented"); //DEBUG
  //not implemented
  double energy = 0.0;
  return energy;
}

/** compute specified output parameter and send for smoothing */
void FluidElementT::SendOutput(int kincode)
{
  WriteCallLocation("SendOutput"); //DEBUG 
  /* output flags */
  iArrayT flags(fNodalOutputCodes.Length());

  /* set flags to get desired output */
  flags = IOBaseT::kAtNever;
  
  switch (kincode)
  {
    case iNodalCrd:
      flags[iNodalCrd] = 1;
    break;
    case iNodalVel:
      flags[iNodalVel] = 1;
    break;
    case iNodalAcc:
      flags[iNodalAcc] = 1;
    break;
    case iNodalPrs:
      flags[iNodalPrs] = 1;
    break;
    default:
      cout << "\n FluidElementT::SendOutput: invalid output code: "
        << kincode << endl;
    }

  /* number of output values */
  iArrayT n_counts;
  SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

  /* reset averaging workspace */
  ElementSupport().ResetAverage(n_counts.Sum());

  /* no element output */
  iArrayT e_counts(fElementOutputCodes.Length());
  e_counts = 0;

  /* generate output */
  dArray2DT n_values, e_values;
  ComputeOutput(n_counts, n_values, e_counts, e_values);
}
/* driver for calculating output values */
void FluidElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
  const iArrayT& e_codes, dArray2DT& e_values)
{
  WriteCallLocation("ComputeOutput"); //DEBUG 
  /* number of output values */
  int n_out = n_codes.Sum();
  int e_out = e_codes.Sum();

  /* nothing to output */
  if (n_out == 0 && e_out == 0) return;

  // not implemented
  return;
}  
/***********************************************************************
 * Protected
 ***********************************************************************/

/** initialize local arrays */
void FluidElementT::SetLocalArrays(void)
{
  WriteCallLocation("SetLocalArrays"); //DEBUG
  /* inherited */
  ContinuumElementT::SetLocalArrays();

  /* allocate */
  int nen = NumElementNodes();
  fLocCurVel.Dimension(nen, NumDOF());
  fLocOldVel.Dimension(nen, NumDOF());
  fLocCurAcc.Dimension(nen, NumDOF());

  /* set source */
  Field().RegisterLocal(fLocCurVel);
  Field().RegisterLocal(fLocOldVel);

  if (fIntegrator->Order() > 0)
    Field().RegisterLocal(fLocCurAcc);

  /* allocate */
  fLocPrs.Dimension(NumElementNodes(), 1);

  /* no need to set source for fLocPrs; handled by XDOF_ManagerT */
}

void FluidElementT::SetShape(void)
{
  WriteCallLocation("SetShape"); //DEBUG
  /* link shape functions to current velocities */
  fCurrShapes = new ShapeFunctionT(GeometryCode(), NumIP(), fLocCurVel);
  if (!fCurrShapes ) throw ExceptionT::kOutOfMemory;

  fCurrShapes->Initialize();
}

/** form shape functions and derivatives */
void FluidElementT::SetGlobalShape(void)
{
  WriteCallLocation("SetGlobalShape"); //DEBUG
  /* shape function wrt current velocities; Eularian formulation */
  SetLocalX(fLocCurVel);

  /* compute shape function derivatives */
  fCurrShapes->SetDerivatives();
}

void FluidElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
  WriteCallLocation("LHSDriver"); //DEBUG
}
void FluidElementT::RHSDriver(void)
{
  WriteCallLocation("RHSDriver"); //DEBUG  
}

/** appends group connectivities to the array (X -> geometry, U -> field) */
void FluidElementT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
  AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
  WriteCallLocation("ConnectsU"); //DEBUG
  /* inherited */
  ContinuumElementT::ConnectsU(connects_1, connects_2);
  
  bool found = false;
  for (int i = connects_1.Length() - 1; i > -1 && !found; i--)
    if (connects_1[i] == fConnectivities[0])
    {
      connects_1[i] = &fXDOFConnectivities;
      found = true;
    }
  /* check */
  if (!found) connects_1.AppendUnique(&fXDOFConnectivities);  
}

/** collecting element group equation numbers */
void FluidElementT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
  AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
  WriteCallLocation("Equations"); //DEBUG
#pragma unused(eq_2)
  /* collect using method allowing mixed node/tag numbers */
  ElementSupport().XDOF_Manager().XDOF_SetLocalEqnos(Group(), fXDOFConnectivities, fXDOFEqnos);
  /* add to equation list */
  eq_1.Append(&fXDOFEqnos);
}
    
/** determine number of tags needed */
void FluidElementT::SetDOFTags(void)
{
  WriteCallLocation("SetDOFTags"); //DEBUG
  /* need a single tag w/ 1 DOF for pressure node
   * based on NodesUsed in elment group
   */
  ElementBaseT::NodesUsed(fNodesUsed);   
  fPressureDOFtags.Dimension(fNodesUsed.Length());
}

/** return an array of tag numbers */
iArrayT& FluidElementT::DOFTags(int tag_set)
{
  WriteCallLocation("DOFTags"); //DEBUG
  /* check */
#if __option(extended_errorcheck)
  if (tag_set != 0)
  {
    cout << "\n FluidElementT::DOFTags: expecting tag set 0:"
      << tag_set << endl;
    throw ExceptionT::kOutOfRange;
  }
#endif  
  return fPressureDOFtags;
}

/** generate nodal connectivities */
void FluidElementT::GenerateElementData(void)
{
  WriteCallLocation("GenerateElementData"); //DEBUG
  iArrayT& DOFTagSet=DOFTags(0);

  /* consider elmenet within some group
   *
   *  08---01---04
   *  |     |    |
   *  09---11---05
   *  |     |    |
   *  12---02---03
   *
   * fNodesUsed [ 12 02 03 09 11 05 08 01 04 ]
   * fPressureDOFtags [ 01 02 03 04 05 06 07 08 09 ]
   * fNodesUsedInverse [ 08 02 03 09 06 -1 -1 07 04 -1 05 01 ]
   
   * AK!? this might cause memory problems, lets say node 1 is '89'
   * meaning the size of fNodesUsedInverse is 89, with 80 values of -1

   *
   */
    
  fNodesUsedInverse.Dimension(fNodesUsed.Max()-1);
  fNodesUsedInverse=-1;
  for(int i=0; i<fNodesUsed.Length(); i++)
  {
     fNodesUsedInverse[fNodesUsed[i]]=DOFTagSet[i];
  }

  fXDOFConnectivities.Dimension(fConnectivities[0]->MajorDim(),2*NumElementNodes());
     
  const int *pelem = fConnectivities[0]->Pointer();
  int rowlength = fConnectivities[0]->MinorDim();
 
  for(int i=0; i<fConnectivities[0]->MajorDim(); i++, pelem+=rowlength)
  {
    int* pxelem = fXDOFConnectivities(i);
    for(int j=0; j<2*NumElementNodes(); j++)
    {
      if(j<NumElementNodes())
      {
        pxelem[j]=  pelem[j];
      }
      else
      {
        pxelem[j]= fNodesUsedInverse[pelem[j-NumElementNodes()]];
      }
    }
  }
}

/** return the connectivities associated with the node */
const iArray2DT& FluidElementT::DOFConnects(int tag_set) const
{
  WriteCallLocation("DOFConnects"); //DEBUG
  /* check */
#if __option(extended_errorcheck)  
  if (tag_set != 0)
  {
    cout << "\n FluidElementT::DOFTags: expecting tag set 0:"
      << tag_set << endl;
    throw ExceptionT::kOutOfRange;
  }
#endif  
  return fXDOFConnectivities;
}

/** describe the parameters needed by the interface */
void FluidElementT::DefineParameters(ParameterListT& list) const
{
  WriteCallLocation("DefineParameters"); //DEBUG
  /* inherited */
  ElementBaseT::DefineParameters(list);
}

/** information about subordinate parameter lists */
void FluidElementT::DefineSubs(SubListT& sub_list) const
{
  WriteCallLocation("DefineSubs"); //DEBUG
  /* inherited */
  ContinuumElementT::DefineSubs(sub_list);

  sub_list.AddSub("fluid_element_nodal_output", ParameterListT::ZeroOrOnce);
  sub_list.AddSub("fluid_element_element_output", ParameterListT::ZeroOrOnce);
  sub_list.AddSub("fluid_element_stab_param", ParameterListT::ZeroOrOnce);
}

/** a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FluidElementT::NewSub(const StringT& name) const
{
  WriteCallLocation("NewSub"); //DEBUG
  if (name == "fluid_element_nodal_output")
  {
    ParameterContainerT* node_output = new ParameterContainerT(name);
    /* all false by default */
    for (int i = 0; i < NumNodalOutputCodes; i++)
    {
      ParameterT output(ParameterT::Integer, NodalOutputNames[i]);
      output.SetDefault(1);
      node_output->AddParameter(output, ParameterListT::ZeroOrOnce);
    }
    return node_output;
  }
  else if (name == "fluid_element_element_output")
  {
    ParameterContainerT* element_output = new ParameterContainerT(name);
    /* all false by default */
    for (int i = 0; i < NumElementOutputCodes; i++)
    {
      ParameterT output(ParameterT::Integer, ElementOutputNames[i]);
      output.SetDefault(1);
      element_output->AddParameter(output, ParameterListT::ZeroOrOnce);
    }
    return element_output;
  }
  else if (name == "fluid_element_stab_param")
  {
    ParameterContainerT* stab_param = new ParameterContainerT(name);
    /* all false by default */
    for (int i = 0; i < NumStabParamCodes; i++)
    {
      ParameterT output(ParameterT::Integer, StabParamNames[i]);
      output.SetDefault(1);
      stab_param->AddParameter(output, ParameterListT::ZeroOrOnce);
    }
    return stab_param;
  }  
  else /* inherited */
    return ContinuumElementT::NewSub(name);  
}

/** accept parameter list */
void FluidElementT::TakeParameterList(const ParameterListT& list)
{
  WriteCallLocation("TakeParameterList"); //DEBUG
  const char caller[] = "FluidElementT::TakeParameterList";

  /* inherited */
  ContinuumElementT::TakeParameterList(list);

  /* nodal output codes */
  fNodalOutputCodes.Dimension(NumNodalOutputCodes);
  fNodalOutputCodes = IOBaseT::kAtNever;
  const ParameterListT* node_output = list.List("fluid_element_nodal_output");
  if (node_output)
  {
    /* set flags */
    for (int i = 0; i < NumNodalOutputCodes; i++)
    {
      /* look for entry */
      const ParameterT* nodal_value = node_output->Parameter(NodalOutputNames[i]);
      if (nodal_value)
      {
        int do_write = *nodal_value;        
        if (do_write == 0)
        {
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: coordinates"); //DEBUG
        }
        else if (do_write == 1)
        {
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: velocities"); //DEBUG
        }
        else if (do_write == 2)
        {
          /* check order of the time integrator */
          if (fIntegrator->Order() < 1)
            ExceptionT::GeneralFail(caller, "expecting time integrator order of 1 or higher");          
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: accelerations"); //DEBUG
        }
        else if (do_write == 3)
        {
          fNodalOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fNodalOutputCodes: pressures"); //DEBUG
        }
      }
    }
  }
  
  /* element output codes */
  fElementOutputCodes.Dimension(NumElementOutputCodes);
  fElementOutputCodes = IOBaseT::kAtNever;
  const ParameterListT* element_output = list.List("fluid_element_element_output");
  if (element_output)
  {
    /* set flags */
    for (int i = 0; i < NumElementOutputCodes; i++)
    {
      /* look for entry */
      const ParameterT* element_value = element_output->Parameter(ElementOutputNames[i]);
      if (element_value)
      {
        int do_write = *element_value;
        if (do_write == 0)
        {
          fElementOutputCodes[i] = IOBaseT::kAtInc;
          WriteCallLocation("Picked: fElementOutputCodes: NONE"); //DEBUG
        }
      }
    }
  }

  /* stabilization parameter codes */
  const ParameterListT* stab_param = list.List("fluid_element_stab_param");
  if (stab_param)
  {
    /* set flags */
    for (int i = 0; i < NumStabParamCodes; i++)
    {
      /* look for entry */
      const ParameterT* stab_param_value = stab_param->Parameter(StabParamNames[i]);
      if (stab_param_value)
      {
        int param_value = *stab_param_value;
        if (param_value == 0)
        {
          /* implement \tau_m = \tau_c = \tau_PSPG = \tau_SUPG \tau
           *
           * \tau = [ (2/(del t))^2 + (2/h * ||v||) + (4\mu / h^2) ]^(-1/2)
           * ** check reference
           */
           WriteCallLocation("Picked: \tau_m = \tau_c"); //DEBUG
        }
      }
    }
  }  
    
  /* only 1 tag set for the group */
  int NumTagSet = 1;
  iArrayT numDOF(NumTagSet);
  numDOF[0] = FluidElementT::kPressureNDOF;

  /* register with node manager - sets initial fPressureDOFtags */
  ElementSupport().XDOF_Manager().XDOF_Register(this, numDOF); 
}

/* construct output labels array */
void FluidElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
  iArrayT& counts) const
{
  WriteCallLocation("SetNodalOutputCodes"); //DEBUG
  /* initialize */
  counts.Dimension(flags.Length());
  counts = 0;

  /* set output flags */
  if (flags[iNodalCrd] == mode) counts[iNodalCrd] = NumSD();
  if (flags[iNodalVel] == mode) counts[iNodalVel] = NumSD();
  if (flags[iNodalAcc] == mode) counts[iNodalAcc] = NumSD();
  if (flags[iNodalPrs] == mode) counts[iNodalPrs] = FluidElementT::kPressureNDOF;
}

void FluidElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
  iArrayT& counts) const
{
  WriteCallLocation("SetElementOutputCodes"); //DEBUG
  /* initialize */
  counts.Dimension(flags.Length());
  counts = 0;

  /* set output flags */
  if (fElementOutputCodes[iNONE] == mode) counts[iNONE] = 0;
}

void FluidElementT::GenerateOutputLabels(const iArrayT& n_codes,
    ArrayT<StringT>& n_labels, const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
  WriteCallLocation("GenerateOutputLabels"); //DEBUG
  const char caller[] = "FluidElementT::GenerateOutputLabels";

  /* allocate node labels */
  n_labels.Dimension(n_codes.Sum());

  int count = 0;
  if (n_codes[iNodalCrd])
  {
    const char* xlabels[] = {"x1", "x2", "x3"};
    for (int i = 0; i < NumSD(); i++)
      n_labels[count++] = xlabels[i];
  }

  if (n_codes[iNodalVel])
  {
    const char* xlabels[] = {"v1", "v2", "v3"};
    for (int i = 0; i < NumSD(); i++)
      n_labels[count++] = xlabels[i];
  }

  if (n_codes[iNodalAcc])
  {
    const char* xlabels[] = {"a1", "a2", "a3"};
    for (int i = 0; i < NumSD(); i++)
      n_labels[count++] = xlabels[i];
  }

  if (n_codes[iNodalPrs])
  {
    const char* xlabels[] = {"p"};
    for (int i = 0; i < FluidElementT::kPressureNDOF; i++)
      n_labels[count++] = xlabels[i];
  }

  if (e_codes.Sum() != 0)
    ExceptionT::GeneralFail("FluidElementT::GenerateOutputLabels",
      "not expecting any element output codes");  
}

/***********************************************************************
 * Private
 ***********************************************************************/


/** FOR DEBUGGING PURPOSES ONLY */
void FluidElementT::WriteCallLocation( char* loc ) const {
cout << "Inside of FluidElementT::" << loc << endl;
}