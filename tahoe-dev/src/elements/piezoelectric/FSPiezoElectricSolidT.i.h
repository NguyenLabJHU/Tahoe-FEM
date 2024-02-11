//
// $Id: FSPiezoElectricSolidT.i.h,v 1.1 2009-05-05 23:34:54 beichuan Exp $
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

namespace Tahoe {

  //
  //
  //
  inline FSPiezoElectricSolidT::FSPiezoElectricSolidT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSPZMatSupport(0), fLocVectorPotential(
        LocalArrayT::kEVP), fCurrMaterial(0), fElectricVectorPotentialField(0)
  {
    SetName("piezoelectric");
  }

  //
  //
  //
  inline void FSPiezoElectricSolidT::Workspace()
  {

    const int nen = NumElementNodes();
    const int nsd = NumSD();
    const int nel = nen * nsd;
    const int nme = nen * nsd;
    const int dof = ManifoldDim() + ElectricalDim();
    const int neq = nen * dof;

    fMaterialTangent.Dimension(nme, nme);
    fGeometricTangent.Dimension(nme, nme);
    fMechanical2ElectricTangent.Dimension(nel, nme);
    fElectric2MechanicalTangent.Dimension(nme, nel);
    fElectricTangent.Dimension(nel, nel);

    fLHS.Dimension(neq);
    fRHS.Dimension(neq);

  }

  //
  //
  //
  inline int FSPiezoElectricSolidT::TotalNumDOF() const
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
  inline const double FSPiezoElectricSolidT::DivergenceVectorPotential() const
  {
    return fDivPhi_List[CurrIP()];
  }

  //
  //
  //
  inline const double FSPiezoElectricSolidT::DivergenceVectorPotential(int ip) const
  {
    return fDivPhi_List[ip];
  }

  //
  //
  //
  inline const int FSPiezoElectricSolidT::ManifoldDim() const
  {
    return FSPZMatSupportT::ManifoldDim();
  }

  //
  //
  //
  inline const int FSPiezoElectricSolidT::StrainDim() const
  {
    return FSPZMatSupportT::StrainDim();
  }

  //
  //
  //
  inline const int FSPiezoElectricSolidT::ElectricalDim() const
  {
    return FSPZMatSupportT::ElectricalDim();
  }

} // namespace Tahoe
