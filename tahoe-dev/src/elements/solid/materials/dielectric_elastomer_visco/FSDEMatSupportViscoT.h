#if !defined(_FSDEMatSupportViscoT_)
#define _FSDEMatSupportViscoT_

#include <cassert>

#include "FSMatSupportT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSDielectricElastomerViscoT;

  //
  // support for dielectric elastomer, finite strain Tahoe materials classes
  //
  class FSDEMatSupportViscoT: public FSMatSupportT
  {

  public:

    //
    // constructors
    //
    FSDEMatSupportViscoT(int ndof, int nip);

    //
    // \name host code information
    // @{
    // return a pointer to the host element. Returns 0 if no
    // element information in available.
    //
    const FSDielectricElastomerViscoT* FSDielectricElastomerVisco() const;

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
    const FSDielectricElastomerViscoT* fFSDielectricElastomerVisco;

  };

} // namespace Tahoe

#include "FSDEMatSupportViscoT.i.h"
#endif // _FSDEMatSupportViscoT_
