//
// $Id: FSHuWashizuUSCT.h,v 1.1 2008-09-03 18:40:50 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/07/14 16:13:25  lxmota
// Initial sources (disabled for now)
//
//

#if !defined(_FSHuWashizuUSCT_)
#define _FSHuWashizuUSCT_

#include <cassert>

#include "FiniteStrainT.h"
#include "FSIsotropicMatT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSMatSupportT;

  //
  // interface for finite deformation piezolectric and field gradients
  //
  class FSHuWashizuUSCT : public FiniteStrainT
  {
    
  public:
      
    //
    // constructor
    //
    FSHuWashizuUSCT(const ElementSupportT& support);

    //
    // destructor
    //
    virtual ~FSHuWashizuUSCT();

    // \name implementation of the ParameterInterfaceT interface
    // @{
    // information about subordinate parameter lists
    virtual void DefineSubs(SubListT& sub_list) const;

    //
    // return the description of the given inline subordinate
    // parameter list
    //
    virtual void DefineInlineSub(const StringT& name,
				 ParameterListT::ListOrderT& order, 
				 SubListT& sub_lists) const;

    //
    // return the description of the given inline subordinate
    // parameter list
    //
    virtual ParameterInterfaceT* NewSub(const StringT& name) const;

    //
    // accept parameter list
    //
    virtual void TakeParameterList(const ParameterListT& list);

    //
    // strain-displacement operator
    //
    virtual void Set_B(const dArray2DT& DNaX, dMatrixT& B);

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

    //
    // extract the list of material parameters
    //
    virtual void CollectMaterialInfo(const ParameterListT& all_params,
				     ParameterListT& mat_params) const;
	
    //
    //
    //
    const ShapeFunctionT& CurrShapeFunction(void) const;

  protected:

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
    // Initialize local arrays
    //
    virtual void SetLocalArrays();

  private:

    //
    //
    //
    void Initialize();

    void Set_B_bar(const dArray2DT& DNaX, dMatrixT& B_bar);

    void AccumulateGeometricStiffness(dMatrixT& Kg,
				      const dArray2DT& DNaX,
				      dSymMatrixT& S);

  protected:

    //
    // \name work space
    //

    //
    // @{
    //
    ArrayT<dArrayT> fC_List;       
    dArrayT         fC_all;        
				   
    ArrayT<dArrayT> fS_List;  
    dArrayT         fS_all;   
				   
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

    LocalArrayT fLocC;
    LocalArrayT fLocS;
    FSIsotropicMatT* fCurrMaterial;

    //
    // Stiffness storage
    //
    dMatrixT fMaterialTangent;
    dMatrixT fGeometricTangent;

    //
    // Mixed interpolation functions
    //
    ShapeFunctionT* fShapesHW;

  };
  
} // namespace Tahoe

#include "FSHuWashizuUSCT.i.h"

#endif // _FSHuWashizuUSCT_
