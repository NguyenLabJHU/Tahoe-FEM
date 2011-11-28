
namespace Tahoe {

  //
  inline FSDielectricElastomer2DViscoT::FSDielectricElastomer2DViscoT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSDEMatSupport2DVisco(0), fCurrMaterial(0),
    fLocScalarPotential(LocalArrayT::kESP), fElectricScalarPotentialField(0)
  {
    SetName("dielectric_elastomer_2D_visco");
  }

  //
  inline int FSDielectricElastomer2DViscoT::TotalNumDOF() const
  {
 	int mechdof = 2;
 	int elecdof = 1;
    return (mechdof+elecdof);
  }

  inline const dArrayT&
  FSDielectricElastomer2DViscoT::ElectricField() const
  {
  	cout << "FSDielectricElastomer2DViscoT::ElectricField()" << endl;
    return fE_List[CurrIP()];
  }

  //
  inline const dArrayT&
  FSDielectricElastomer2DViscoT::ElectricField(int ip) const
  {
  	cout << "FSDielectricElastomer2DViscoT::ElectricField(int ip)" << endl;
    return fE_List[ip];
  }


} // namespace Tahoe
