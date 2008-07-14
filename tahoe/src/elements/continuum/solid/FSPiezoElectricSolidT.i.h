//
// $Id: FSPiezoElectricSolidT.i.h,v 1.2 2008-07-14 17:37:23 lxmota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/06/16 18:15:10  lxmota
// Piezoelectric solid. Initial source.
//
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
    fLocLastVectorPotential(LocalArrayT::kLastEVP),
    fCurrMaterial(0),
    fElectricVectorPotentialField(0)
  {

    SetName("piezoelectric");

  }

  //
  //
  //
  inline int
  FSPiezoElectricSolidT::TotalNumDOF() const
  {
    return ManifoldDim() + ElectricalDim();
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

#if 0
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

#endif


  //
  //
  //
  inline const int
  FSPiezoElectricSolidT::ManifoldDim() const
  {
    return FSPZMatSupportT::ManifoldDim();
  }


  //
  //
  //
  inline const int
  FSPiezoElectricSolidT::StrainDim() const
  {
    return FSPZMatSupportT::StrainDim();
  }


  //
  //
  //
  inline const int
  FSPiezoElectricSolidT::ElectricalDim() const
  {
    return FSPZMatSupportT::ElectricalDim();
  }

} // namespace Tahoe
