#if !defined(_FSDEMatSupportT_)
#define _FSDEMatSupportT_

#include <cassert>

#include "FSMatSupportT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSDielectricElastomerT;

  //
  // support for dielectric elastomer, finite strain Tahoe materials classes
  //
  class FSDEMatSupportT: public FSMatSupportT
  {

  public:

    //
    // constructors
    //
    FSDEMatSupportT(int ndof, int nip);

    //
    // \name host code information
    // @{
    // return a pointer to the host element. Returns 0 if no
    // element information in available.
    //
    const FSDielectricElastomerT* FSDielectricElastomer() const;

    //
    // set the element group pointer
    //
    virtual void SetContinuumElement(const ContinuumElementT* p);

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
    // Set source for electric displacement
    //
    void SetElectricDisplacement(const ArrayT<dArrayT>* D_List);

    // @}
    //

    static const int ManifoldDim() { return 3; };
    static const int StrainDim() { return 6; };
    static const int ElectricalDim() { return 3; };

  private:

    //
    // \name Sources for electric displacements and pointers to local
    // arrays
    //
    // @{
    //

    //
    // Electric displacement
    //
    const ArrayT<dArrayT>* fD_List;

    //
    // @}
    //

    // pointer to the host element
    const FSDielectricElastomerT* fFSDielectricElastomer;

  };

} // namespace Tahoe

#include "FSDEMatSupportT.i.h"
#endif // _FSDEMatSupportT_
