
namespace Tahoe {

  //
  //
  //
  inline FSDielectricElastomer2DT::FSDielectricElastomer2DT(
      const ElementSupportT& support) :
    FiniteStrainT(support), fFSDEMatSupport2D(0), fCurrMaterial(0),
    fLocScalarPotential(LocalArrayT::kESP), fElectricScalarPotentialField(0)
  {
    SetName("dielectric_elastomer_2D");
  }

  //
  //
  //
  inline int FSDielectricElastomer2DT::TotalNumDOF() const
  {
 	int mechdof = 2;
 	int elecdof = 1;
    return (mechdof+elecdof);
  }

  inline const dArrayT&
  FSDielectricElastomer2DT::ElectricField() const
  {
  	cout << "FSDielectricElastomer2DT::ElectricField()" << endl;
    return fE_List[CurrIP()];
  }

  //
  //
  //
  inline const dArrayT&
  FSDielectricElastomer2DT::ElectricField(int ip) const
  {
  	cout << "FSDielectricElastomer2DT::ElectricField(int ip)" << endl;
    return fE_List[ip];
  }


} // namespace Tahoe
