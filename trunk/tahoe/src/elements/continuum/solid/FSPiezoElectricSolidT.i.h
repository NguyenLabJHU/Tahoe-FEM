//
// $Id: FSPiezoElectricSolidT.i.h,v 1.1 2008-06-16 18:15:10 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

namespace Tahoe {

  //
  //
  //
  inline
  FSPiezoElectricSolidT::FSPiezoElectricSolidT(const ElementSupportT& support):
    FiniteStrainT(support),
    fFSPZMatSupport(0),
    fLocVectorPotential(LocalArrayT::kEVP),
    fLocLastVectorPotential(LocalArrayT::kLastEVP)
  {

    SetName("Finite-deformation piezoelectric");
    Initialize();

  }

  //
  //
  //
  inline const dArrayT&
  FSPiezoElectricSolidT::ElectricDisplacement() const
  {

    return fD_List[CurrIP()];

  }

  //
  //
  //
  inline const dArrayT&
  FSPiezoElectricSolidT::ElectricDisplacement(int ip) const
  {

    return fD_List[ip];

  }

  //
  //
  //
  inline const dArrayT&
  FSPiezoElectricSolidT::ElectricDisplacement_last() const
  {

    return fD_last_List[CurrIP()];

  }

  //
  //
  //
  inline const dArrayT&
  FSPiezoElectricSolidT::ElectricDisplacement_last(int ip) const
  {

    return fD_last_List[ip];

  }

  //
  //
  //
  inline bool
  FSPiezoElectricSolidT::Needs_D(int materialNumber) const
  {

    const ArrayT<bool>& needs = fMaterialNeeds[materialNumber];

    return needs[fNeedsOffset + kD];

  }

  //
  //
  //
  inline bool
  FSPiezoElectricSolidT::Needs_D_last(int materialNumber) const
  {

    const ArrayT<bool>& needs = fMaterialNeeds[materialNumber];

    return needs[fNeedsOffset + kD_last];

  }

} // namespace Tahoe
