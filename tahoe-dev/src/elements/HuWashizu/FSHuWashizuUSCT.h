//
// $Id: FSHuWashizuUSCT.h,v 1.1 2009-05-05 23:33:42 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2009/04/02 00:49:52  amota
// Changes in handling of internal variables.
//
// Revision 1.1  2008/12/12 18:56:58  amota
// Moved to new localtion.
//
// Revision 1.1  2008/07/14 16:13:25  lxmota
// Initial sources (disabled for now)
//
//

#if !defined(_FSHuWashizuUSCT_)
#define _FSHuWashizuUSCT_

#include <cassert>

#include "FiniteStrainT.h"
#include "FSSolidMatT.h"
#include "MaterialListT.h"
#include "MixedShapeFunctionT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSMatSupportT;

  //
  //
  //
  class FSHuWashizuUSCT: public FiniteStrainT {

  public:

    //
    // constructor
    //
    FSHuWashizuUSCT(const ElementSupportT& support);

    //
    // destructor
    //
    virtual ~FSHuWashizuUSCT();

    //
    // parameters needed by interface
    //
    virtual void DefineParameters(ParameterListT& list) const;

    //
    // accept parameter list
    //
    virtual void TakeParameterList(const ParameterListT& list);

    //
    // increment current element
    //
    virtual bool NextElement();

    //
    // strain-displacement operator
    //
    virtual void Set_B(const dArray2DT& DNaX, dMatrixT& B);

    //
    // element stiffness matrix
    //
    virtual void FormStiffness(double constK);

    //
    // internal force
    //
    virtual void FormKd(double constK);

    //
    // @}
    //

    virtual const ShapeFunctionT& MixedShapeFunction() const;
    virtual int NumberMixedNodes() const;
    virtual int NumberMixedIP() const;

    virtual void InitStep();
    virtual void CloseStep();

    // restore last converged state
    virtual GlobalT::RelaxCodeT ResetStep();

  protected:

    //
    // initialize shape functions
    //
    virtual void SetShape();

    //
    // form shape functions and derivatives
    //
    virtual void SetGlobalShape();

    //
    // write all current element information to the stream. used to
    // generate debugging information after runtime errors
    //
    virtual void CurrElementInfo(ostream& out) const;

    //
    // Initialize local arrays
    //
    virtual void SetLocalArrays();

    //
    // Compute deformation measures at an arbitrary location
    // given in parametric coordinates.
    //
    virtual dMatrixT DeformationGradient(const dArrayT& paramPt) const;

    virtual dSymMatrixT
        RightCauchyGreenDeformation(const dArrayT& paramPt) const;

    //
    // Compute 2nd PK stress at an arbitrary location
    // given in parametric coordinates.
    //
    virtual dSymMatrixT Stress2PK(const dArrayT& paramPt);

  private:

    //
    // Internal variables
    //
    void Initialize();

    //
    //
    //
    void Workspace();

    void ComputeInverseWeightedVolume();

    //
    // Given deformation gradient and interpolation functions
    // derivatives at arbitrary parametric point, compute
    // strain displacement operator.
    //
    dMatrixT StrainDisplacementOperator(const dMatrixT& F,
        const dArray2DT& DNaX) const;

    //
    // Strain-displacement operator at arbitrary parametric point.
    //
    dMatrixT StrainDisplacementOperator(const dArrayT& paramPt) const;

    //
    // Mixed strain-displacement operator at arbitrary parametric point.
    //
    dMatrixT MixedStrainDisplacementOperator(const dArrayT& paramPt) const;

    //
    // Strain measure at arbitrary parametric point.
    //
    dSymMatrixT Strain(const dArrayT& paramPt);

    //
    // Stress at arbitrary parametric point.
    //
    dSymMatrixT Stress(const dArrayT& paramPt);

    //
    // Mixed strain measure at arbitrary parametric point.
    //
    dSymMatrixT MixedStrain(const dArrayT& paramPt);

    //
    // Mixed stress at arbitrary parametric point.
    //
    dSymMatrixT MixedStress(const dArrayT& paramPt);

    //
    // Interpolate internal variables to arbitrary parametric point.
    //
    dArrayT InterpolateInternal(const dArrayT& paramPt);

    //
    // Internal variables held by the material.
    //
    const dArrayT& InternalVariables();

    //
    // Compute nodal values of internal variables by extrapolating
    // from integration points.
    //
    void NodalInternalVariables();

    //
    //
    //
    void AccumulateGeometricStiffness(dMatrixT& Kg, const dArray2DT& DNaX,
        dSymMatrixT& S);

    LocalArrayT GetMixedNodalPositions();

  protected:

    //
    // \name work space
    //

    //
    // @{
    //
    dMatrixT fC_nodal;

    dMatrixT fS_nodal;

    //
    // @}
    //

    //
    // The material support used to construct materials lists. This
    // pointer is only set the first time
    // FSHuWashizuUSCT::NewMaterialList is called.
    //
    FSMatSupportT* fFSMatSupport;

  private:

    int fNumberMixedNodes;
    int fNumberMixedIP;

    bool fUseMixed;

    LocalArrayT fLocC;
    LocalArrayT fLocS;
    LocalArrayT fLocMixedInitCoords;
    LocalArrayT fNodalInternalVariables;
    FSSolidMatT* fCurrMaterial;

    //
    // Stiffness storage
    //
    dMatrixT fMaterialTangent;
    dMatrixT fGeometricTangent;

    //
    // Mixed interpolation functions
    //
    MixedShapeFunctionT* fShapesMixed;

    //
    // inverse weighted volume (H inverse)
    //
    dMatrixT fInvH;

  };

} // namespace Tahoe

#include "FSHuWashizuUSCT.i.h"

#endif // _FSHuWashizuUSCT_
