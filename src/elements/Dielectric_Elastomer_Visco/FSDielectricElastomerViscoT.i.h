
namespace Tahoe {

  //
  //
  //
  inline FSDielectricElastomerViscoT::FSDielectricElastomerViscoT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSDEMatSupport(0), fCurrMaterial(0),
    fLocScalarPotential(LocalArrayT::kESP), fElectricScalarPotentialField(0)
  {
    SetName("dielectric_elastomer");
  }

  //
  //
  //
  inline int FSDielectricElastomerViscoT::TotalNumDOF() const
  {
 	int mechdof = 3;
 	int elecdof = 1;
    return (mechdof+elecdof);
  }

  inline const dArrayT&
  FSDielectricElastomerViscoT::ElectricField() const
  {
  	cout << "FSDielectricElastomerViscoT::ElectricField()" << endl;
    return fE_List[CurrIP()];
  }

  //
  //
  //
  inline const dArrayT&
  FSDielectricElastomerViscoT::ElectricField(int ip) const
  {
  	cout << "FSDielectricElastomerViscoT::ElectricField(int ip)" << endl;
    return fE_List[ip];
  }


} // namespace Tahoe
