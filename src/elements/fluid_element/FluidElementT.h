/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.h,v 1.4 2006-07-12 17:48:41 a-kopacz Exp $ */
/* created: a-kopacz (07/04/2006) */
#ifndef _FLUID_ELEMENT_H_
#define _FLUID_ELEMENT_H_

/* base class */
#include "ContinuumElementT.h"
#include "DOFElementT.h"

/* direct members
 * NONE
 */

namespace Tahoe {

/* forward declarations
 * NONE
 */
 
/* Fluid element; 4 DOF's per node. Fourth degree of freedom
 * is supported by the inheritance of the DOFElementT interface.
 */
class FluidElementT: public ContinuumElementT, public DOFElementT
{
public:

  /** list/index of nodal outputs */
  enum NodalOutputCodeT {
    iNodalCrd = 0, /**< (reference) coordinates */
    iNodalVel = 1, /**< velocitiets */
    iNodalAcc = 2, /**< accelerations */
    iNodalPrs = 3 /**< pressures */
  };

  /** list/index of element outputs */
  enum ElementOutputCodeT {
    iNONE = 0 /**< NONE */
  };

  /** list/index of stabilization parameters */
  enum StabParamCodeT {
    iStabParamOne = 0 /**< \tau_m = \tau_c = \tau_PSPG = \tau_SUPG */
  };

  /** constructor */
  FluidElementT(const ElementSupportT& support);

  /** destructor */
  virtual ~FluidElementT(void);
  
  /** \name access to nodal values */
  /*@{*/
  const LocalArrayT& OldVelocities(void) const;
  const LocalArrayT& Velocities(void) const;
  const LocalArrayT& Accelerations(void) const;
  const LocalArrayT& Pressures(void) const;
  /*@}*/

  /** compute nodal force */
  virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

  /** returns the stored energy */
  virtual double InternalEnergy(void);

  /** compute specified output parameter and send for smoothing */
  virtual void SendOutput(int kincode);

  /** driver for calculating output values */
  virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
    const iArrayT& e_codes, dArray2DT& e_values);  

protected:

  /* parameters */
  static const int NumNodalOutputCodes;
  static const int NumElementOutputCodes;
  static const int NumStabParamCodes;
  static const int kPressureNDOF;

  /** pointer to shape functions wrt current velocities */
  ShapeFunctionT* fCurrShapes;

  /** initialization functions */
  virtual void SetLocalArrays(void);
  virtual void SetShape(void);

  /** form shape functions and derivatives */
  virtual void SetGlobalShape(void);

  /** drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
  virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
  virtual void RHSDriver(void);

  /** tag for the pressure DOF (1 DOF/tag; 1 tag/node) */
  iArrayT fPressureDOFtags;

  /** map vectors */
  iArrayT fNodesUsed;
  iArrayT fNodesUsedInverse;

  /** extended interaction data */
  iArray2DT fXDOFConnectivities;
  iArray2DT fXDOFEqnos;

  /** appends group connectivities to the array (X -> geometry, U -> field) */  
  virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const { /** ConnectX is called if geometry changes */ }
  virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
    AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
  
  /** collecting element group equation numbers */
  virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
    AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
      
  /**************************************************************/
  /************ implementing interface of DOFElementT ************/
  /**************************************************************/
  /* virtual int Group(void) const;
   * virtual void SetDOFTags(void);
   * virtual iArrayT& DOFTags(int tag_set);
   * virtual void GenerateElementData(void);
   * virtual const iArray2DT& DOFConnects(int tag_set) const;
   * virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;
   * virtual int Reconfigure(void);
   * virtual void ResetState(void);
  /**************************************************************/
  /**************************************************************/

  /** \name implementation of the DOFElementT interface */
  /*@{*/  
  /** return the equation group to which the generate degrees of freedom belong */
  virtual int Group(void) const { return ElementBaseT::Group(); };

  /** determine number of tags needed */
  virtual void SetDOFTags(void);

  /** return an array of tag numbers */
  virtual iArrayT& DOFTags(int tag_set);

  /** generate nodal connectivities */
  virtual void GenerateElementData(void);  

  /** return the connectivities associated with the node */
  virtual const iArray2DT& DOFConnects(int tag_set) const;

  /** restore/initialize the values of the element DOF */
  virtual void ResetDOF(dArray2DT& DOF, int) const { /* DOF[0] = */ };

  /** check element group for tag reconfiguration */
  virtual int Reconfigure(void) { return 0; };

  /** restore any state data to the previous converged state */
  virtual void ResetState(void) {};
  /*@}*/
  
  /**************************************************************/
  /******* implementing interface of ParameterInterfaceT ********/
  /**************************************************************/
  /* virtual void DefineParameters(ParameterListT& list) const;
   * virtual void DefineSubs(SubListT& sub_list) const;
   * virtual ParameterInterfaceT* NewSub(const StringT& name) const;
   * virtual virtual void TakeParameterList(const ParameterListT& list);
  /**************************************************************/
  /**************************************************************/

  /** \name implementation of the ParameterInterfaceT interface */
  /*@{*/
  /** describe the parameters needed by the interface */
  virtual void DefineParameters(ParameterListT& list) const;

  /** information about subordinate parameter lists */
  virtual void DefineSubs(SubListT& sub_list) const;

  /** a pointer to the ParameterInterfaceT of the given subordinate */
  virtual ParameterInterfaceT* NewSub(const StringT& name) const;

  /** accept parameter list */
  virtual void TakeParameterList(const ParameterListT& list);
  /*@}*/
   
private:

  /** reference to current element shape functions, over undeformed configuration*/
  const ShapeFunctionT& CurrShapeFunction(void) const;

  /** nodal pressure values with local ordering */
  LocalArrayT fLocPrs;

  /** nodal current/old velocities with local ordering */
  LocalArrayT fLocCurVel;
  LocalArrayT fLocOldVel;

  /** nodal current accelerations with local ordering */
  LocalArrayT fLocCurAcc; // post-processing only

  /** \name construct output labels array */
  /*@{*/
  virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
    iArrayT& counts) const;
  virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
    iArrayT& counts) const;
  virtual void GenerateOutputLabels(const iArrayT& n_counts,
    ArrayT<StringT>& n_labels, const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;
  /*@}*/
      
  /** FOR DEBUGGING PURPOSES ONLY */
  void WriteCallLocation( char* loc ) const;
};

/* accessors */
inline const ShapeFunctionT& FluidElementT::CurrShapeFunction(void) const
{
#if __option(extended_errorcheck)
if (!fCurrShapes)
  ExceptionT::GeneralFail("FluidElementT::CurrShapeFunctionT", "no shape functions");
#endif
  return *fCurrShapes;
}
inline const LocalArrayT& FluidElementT::OldVelocities(void) const { return fLocOldVel; }
inline const LocalArrayT& FluidElementT::Velocities(void) const { return fLocCurVel; }
inline const LocalArrayT& FluidElementT::Accelerations(void) const { return fLocCurAcc; }
inline const LocalArrayT& FluidElementT::Pressures(void) const { return fLocPrs; }

} // namespace Tahoe
#endif /* _FLUID_ELEMENT_H_ */
