#if !defined(_FSDEMatSupportQ1P0ViscoT_)
#define _FSDEMatSupportQ1P0ViscoT_

#include <cassert>

#include "FSMatSupportT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSDielectricElastomerQ1P0ViscoT;

  //
  // support for dielectric elastomer, finite strain Tahoe materials classes
  //
  class FSDEMatSupportQ1P0ViscoT: public FSMatSupportT
  {

  public:

    //
    // constructors
    //
    FSDEMatSupportQ1P0ViscoT(int ndof, int nip);

    //
    // \name host code information
    // @{
    // return a pointer to the host element. Returns 0 if no
    // element information in available.
    //
    const FSDielectricElastomerQ1P0ViscoT* FSDielectricElastomerQ1P0Visco() const;

    //
    // set the element group pointer
    //
    virtual void SetContinuumElement(const ContinuumElementT* p);

    //
    // \name Electric fields
    // @{
    //

    //
    // Electric field at current integration point
    //
    const dArrayT& ElectricField() const;

    //
    // Electric field at given integration point
    //
    const dArrayT& ElectricField(int ip) const;

    //
    // Set source for electric field
    //
    void SetElectricField(const ArrayT<dArrayT>* E_List);

    //
    // Return pointer to specified local array
    //
    virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;

    //
    // Set pointer to local array
    //
    virtual void SetLocalArray(const LocalArrayT& array);

    //
    // Nodal electric scalar potentials
    //
    const LocalArrayT* ScalarPotentials() const;

    //
    //
    //
    void SetScalarPotentials(const LocalArrayT& scalarPotentials);

  private:

    //
    // \name Sources for electric field and pointers to local
    // arrays
    //
    // @{
    //

    //
    // Electric field
    //
   const ArrayT<dArrayT>* fE_List;

    //
    // Pointers to local arrays
    //
    const LocalArrayT* fESP;

    // pointer to the host element
    const FSDielectricElastomerQ1P0ViscoT* fFSDielectricElastomerQ1P0Visco;

  };

} // namespace Tahoe

#include "FSDEMatSupportQ1P0ViscoT.i.h"
#endif // _FSDEMatSupportQ1P0ViscoT_
