//
// $Id: FSPZMatSupportT.h,v 1.2 2008-12-12 18:58:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2008/07/14 17:38:53  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:21:41  lxmota
// Piezoelectric material support. Initial sources.
//
//

#if !defined(_FSPZMatSupportT_)
#define _FSPZMatSupportT_

#include <cassert>

#include "FSMatSupportT.h"

namespace Tahoe {

  //
  // forward declarations
  //
  class FSPiezoElectricSolidT;

  //
  // support for piezoelectric, finite strain Tahoe materials classes
  //
  class FSPZMatSupportT: public FSMatSupportT
  {

  public:

    //
    // constructors
    //
    FSPZMatSupportT(int ndof, int nip);

    //
    // \name host code information
    // @{
    // return a pointer to the host element. Returns 0 if no
    // element information in available.
    //
    const FSPiezoElectricSolidT* FSPiezoElectricSolid() const;

    //
    // set the element group pointer
    //
    virtual void SetContinuumElement(const ContinuumElementT* p);

    //
    // @}
    //

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

    //
    // accessors for divergence of vector potential
    //
    double DivergenceVectorPotential() const;
    double DivergenceVectorPotential(int ip) const;

    //
    // mutators for divergence of vector potential
    //
    void SetDivergenceVectorPotential(const dArrayT* DivPhi_List);

    //
    // Return pointer to specified local array
    //
    virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;

    //
    // Set pointer to local array
    //
    virtual void SetLocalArray(const LocalArrayT& array);

    //
    // Nodal electric vector potentials
    //
    const LocalArrayT* VectorPotentials() const;

    //
    //
    //
    void SetVectorPotentials(const LocalArrayT& vectorPotentials);

    //
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

    // divergence of vector potential
    const dArrayT* fDivPhi_List;

    //
    // Pointers to local arrays
    //
    const LocalArrayT* fEVP;

    //
    // @}
    //

    // pointer to the host element
    const FSPiezoElectricSolidT* fFSPiezoElectricSolid;

  };

} // namespace Tahoe

#include "FSPZMatSupportT.i.h"

#endif // _FSPZMatSupportT_
