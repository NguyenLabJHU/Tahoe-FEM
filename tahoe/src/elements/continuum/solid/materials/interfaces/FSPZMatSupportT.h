//
// $Id: FSPZMatSupportT.h,v 1.1 2008-06-16 18:21:41 lxmota Exp $
//
// $Log: not supported by cvs2svn $
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
    // Electric displacement from end of previous time step at current
    // integration point
    //
    const dArrayT& ElectricDisplacement_last() const;

    //
    // Electric displacement from end of previous time step at given
    // integration point
    //
    const dArrayT& ElectricDisplacement_last(int ip) const;

    //
    // Set source for electric displacement
    //
    void SetElectricDisplacement(const ArrayT<dArrayT>* D_List);

    //
    // Set source for electric displacement from end of previous time
    // step
    //
    void SetElectricDisplacement_last(const ArrayT<dArrayT>* D_last_List);

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
    // Last nodal electric vector potentials
    //
    const LocalArrayT* LastVectorPotentials() const;

    //
    //
    //
    void SetVectorPotentials(const LocalArrayT& vectorPotentials);

    //
    //
    //
    void SetLastVectorPotentials(const LocalArrayT& last_vectorPotentials);

    //
    // @}
    //
    
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
    // Last electric displacement
    //
    const ArrayT<dArrayT>* fD_last_List;

    //
    // Pointers to local arrays
    //
    const LocalArrayT* fLastEVP;
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
