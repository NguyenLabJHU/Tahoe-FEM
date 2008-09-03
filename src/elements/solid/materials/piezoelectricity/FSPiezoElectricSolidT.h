//
// $Id: FSPiezoElectricSolidT.h,v 1.1 2008-09-03 18:40:50 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2008/07/14 17:37:23  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:15:10  lxmota
// Piezoelectric solid. Initial source.
//
//


#if !defined(_FSPiezoElectricSolidT_)
#define _FSPiezoElectricSolidT_

#include <cassert>

#include "FiniteStrainT.h"
#include "FSNeoHookePZLinT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSPZMatSupportT;

  //
  // interface for finite deformation piezolectric and field gradients
  //
  class FSPiezoElectricSolidT : public FiniteStrainT
  {

  public:

    //
    // constructor
    //
    FSPiezoElectricSolidT(const ElementSupportT& support);

    //
    // destructor
    //
    virtual ~FSPiezoElectricSolidT();

    //
    // specify parameters needed by the interface
    //
    virtual void DefineParameters(ParameterListT& list) const;

    //
    //
    //
    virtual int TotalNumDOF() const;

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
    // Electric displacement from end of previous time step at current
    // integration point
    //
    const dArrayT& ElectricDisplacement_last() const;

    //
    // Electric displacement from end of previous time step at given
    // integration point
    //
    const dArrayT& ElectricDisplacement_last(int ip) const;

#if 0
    //
    // @}
    //

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
    // extract the list of material parameters
    //
    virtual void CollectMaterialInfo(const ParameterListT& all_params,
             ParameterListT& mat_params) const;

#endif

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
    //
    //
    const ShapeFunctionT& CurrShapeFunction(void) const;

#if 0
    virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
       AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
#endif

#if 0
    //
    //
    //
    static const int kNeedVectorPotential     = KNeedLastDisp + 1;
    static const int kNeedLastVectorPotential = KNeedLastDisp + 2;

    static const int kD      = kF_last + 1;
    static const int kD_last = kF_last + 2;

    bool Needs_D(int materialNumber) const;
    bool Needs_D_last(int materialNumber) const;
#endif

    const int ManifoldDim() const;
    const int StrainDim() const;
    const int ElectricalDim() const;

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
    // Initialize local arrays
    //
    virtual void SetLocalArrays();

  private:

    //
    //
    //
    void Initialize();

    void Set_B_C(const dArray2DT& DNaX, dMatrixT& B_C);
    void Set_B_D(const dArray2DT& DNaX, dMatrixT& B_D);
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
    ArrayT<dArrayT> fD_List;       // < deformation gradient
    dArrayT         fD_all;        // < grouped memory for all
				                           // < electric displacements
    ArrayT<dArrayT> fD_last_List;  // < last deformation gradient
    dArrayT         fD_last_all;   // < grouped memory for all last
				                           // < electric displacements
    //
    // @}
    //

    //
    // The material support used to construct materials lists. This
    // pointer is only set the first time
    // FSPiezoElectricSolidT::NewMaterialList is called.
    //
    FSPZMatSupportT* fFSPZMatSupport;

  private:

#if 0
    int fNeedsOffset;
#endif
    LocalArrayT fLocVectorPotential;
    LocalArrayT fLocLastVectorPotential;
    FSNeoHookePZLinT* fCurrMaterial;

    //
    // Stiffness storage
    //
    dMatrixT fMaterialTangent;
    dMatrixT fGeometricTangent;
    dMatrixT fMechanical2ElectricTangent;
    dMatrixT fElectric2MechanicalTangent;
    dMatrixT fElectricTangent;

    //
    // Electric field
    //
    const FieldT* fElectricVectorPotentialField;

  };

} // namespace Tahoe

#include "FSPiezoElectricSolidT.i.h"

#endif // _FSPiezoElectricSolidT_
