
namespace Tahoe {

  //
  //
  //
  inline FSDielectricElastomerT::FSDielectricElastomerT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSDEMatSupport(0), fCurrMaterial(0),
    fLocScalarPotential(LocalArrayT::kESP), fElectricScalarPotentialField(0)
  {
    SetName("dielectric_elastomer");
  }

  //
  //
  //
  inline int FSDielectricElastomerT::TotalNumDOF() const
  {
 	int mechdof = 3;
 	int elecdof = 1;
    return (mechdof+elecdof);
  }

  inline const dArrayT&
  FSDielectricElastomerT::ElectricField() const
  {
  	cout << "FSDielectricElastomerT::ElectricField()" << endl;
    return fE_List[CurrIP()];
  }

  //
  //
  //
  inline const dArrayT&
  FSDielectricElastomerT::ElectricField(int ip) const
  {
  	cout << "FSDielectricElastomerT::ElectricField(int ip)" << endl;
    return fE_List[ip];
  }


} // namespace Tahoe
