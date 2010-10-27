
#if !defined(_FSDielectricElastomerT_)
#define _FSDielectricElastomerT_

#include <cassert>

#include "FiniteStrainT.h"
#include "FSDEMatT.h"

namespace Tahoe {

  /* Forward declarations */
  class FSDEMatSupportT;

  //
  // interface for finite deformation dielectric elastomers 
  // based on 2008 JMPS paper of Suo et al.
  
  class FSDielectricElastomerT: public FiniteStrainT {

  public:

    //
    // constructor
    //
    FSDielectricElastomerT(const ElementSupportT& support);

    //
    // destructor
    //
    virtual ~FSDielectricElastomerT();

    //
    // specify parameters needed by the interface
    //
    virtual void DefineParameters(ParameterListT& list) const;

    //
    // accept parameter list
    //
    virtual void TakeParameterList(const ParameterListT& list);

    //
    // information about subordinate parameter lists
    //
    virtual void DefineSubs(SubListT& sub_list) const;

    //
    //
    //
    virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
        AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

    //
    // \name Electric displacements
    // @{
    //

    //
    // Electric displacement at current integration point
    //
    const dArrayT& ElectricDisplacement() const;

    //
    // Electric displacement at given integration point
    //
    const dArrayT& ElectricDisplacement(int ip) const;

    //
    // increment current element
    //
    virtual bool NextElement();

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

  protected:

    //
    // Construct a new material support and return a
    // pointer. Recipient is responsible for freeing the pointer.
    //


    //
    // \param p an existing MaterialSupportT to be initialized. If
    // 0, allocate a new MaterialSupportT and initialize it.
    //
    virtual MaterialSupportT*
    NewMaterialSupport(MaterialSupportT* p = 0) const;

    //
    // Return a pointer to a new material list. Recipient is
    // responsible for freeing the pointer.
    //
    // \param name list identifier
    // \param size length of the list
    //
    virtual MaterialListT* NewMaterialList(const StringT& name, int size);

    //
    // form shape functions and derivatives
    //
    virtual void SetGlobalShape(void);

    //
    // write all current element information to the stream. used to
    // generate debugging information after runtime errors
    //
    virtual void CurrElementInfo(ostream& out) const;

    //
    // driver for calculating output values
    //
    virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
        const iArrayT& e_codes, dArray2DT& e_values);

  private:

    //
    //
    //
    void Workspace();

  protected:

    //
    // \name work space
    //

    //
    // @{
    //

    // electric displacement
    ArrayT<dArrayT> fD_List;
    dArrayT fD_all;

    //
    // @}
    //

    //
    // The material support used to construct materials lists. This
    // pointer is only set the first time
    // FSDielectricElastomerT::NewMaterialList is called.
    //
    FSDEMatSupportT* fFSDEMatSupport;

  private:

    FSDEMatT* fCurrMaterial;
    //
    // Stiffness storage
    //
    dMatrixT fAmm;	// purely mechanical part of Hessian matrix
    dMatrixT fAme;	// mechanical-electrical coupling part of Hessian matrix
    dMatrixT fAem;	// electrical-mechanical coupling part of Hessian matrix
    dMatrixT fAee;	// electrical-electrical coupling part of Hessian matrix
    dMatrixT fDefGrad;	// deformation gradient
    dMatrixT fStretch;	// stretch (F^T * F)
    /* these stiffnesses may not be needed */
    dMatrixT fCmm;	// mechanical stiffness
    dMatrixT fCme;	// mechanical-electrical stiffness
    dMatrixT fCem;	// electrical-mechanical stiffness
    dMatrixT fCee;	// electrical stiffness
    dMatrixT fGradNa;	// shape function gradients matrix
    
  };

} // namespace Tahoe

#include "FSDielectricElastomerT.i.h"
#endif // _FSDielectricElastomerT_
