
namespace Tahoe {

  //
  //
  //
  inline FSDielectricElastomerT::FSDielectricElastomerT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSDEMatSupport(0), fCurrMaterial(0)
  {
    SetName("dielectric_elastomer");
  }

  //
  //
  //
  inline const dArrayT&
  FSDielectricElastomerT::ElectricDisplacement() const
  {
    return fD_List[CurrIP()];
  }

  //
  //
  //
  inline const dArrayT&
  FSDielectricElastomerT::ElectricDisplacement(int ip) const
  {
    return fD_List[ip];
  }


} // namespace Tahoe
