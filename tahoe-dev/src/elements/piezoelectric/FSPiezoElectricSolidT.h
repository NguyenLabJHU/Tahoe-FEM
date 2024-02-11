//
// $Id: FSPiezoElectricSolidT.h,v 1.1 2009-05-05 23:34:54 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/12/12 18:56:28  amota
// Moved to new localtion.
//
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
  // interface for finite deformation piezoelectric and field gradients
  //
  class FSPiezoElectricSolidT: public FiniteStrainT {

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
    // accept parameter list
    //
    virtual void TakeParameterList(const ParameterListT& list);

    //
    //
    //
    virtual int TotalNumDOF() const;

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
    // accessors for divergence of vector potential
    //
    const double DivergenceVectorPotential() const;
    const double DivergenceVectorPotential(int ip) const;
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

    void Set_B_C(const dArray2DT& DNaX, dMatrixT& B_C);
    void Set_B_D(const dArray2DT& DNaX, dMatrixT& B_D);
    void Set_B_K(const dArray2DT& DNaX, dMatrixT& B_K);
    void AccumulateGeometricStiffness(dMatrixT& Kg, const dArray2DT& DNaX,
        dSymMatrixT& S);

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

    // divergence of vector potential
    dArrayT fDivPhi_List;
    dArrayT fDivPhi_all;

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

    LocalArrayT fLocVectorPotential;
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
