
namespace Tahoe {

  //
  //
  //
  inline FSDEMatSupportViscoT::FSDEMatSupportViscoT(int ndof, int nip) :
    FSMatSupportT(ndof, nip), fE_List(0), fFSDielectricElastomerVisco(0)
  {

  }

  //
  //
  //
  inline const FSDielectricElastomerViscoT*
  FSDEMatSupportViscoT::FSDielectricElastomerVisco() const
  {
    return fFSDielectricElastomerVisco;
  }

  inline const dArrayT&
  FSDEMatSupportViscoT::ElectricField() const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[CurrIP()];
  }

  inline const dArrayT&
  FSDEMatSupportViscoT::ElectricField(int ip) const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[ip];
  }

  inline void FSDEMatSupportViscoT::SetElectricField(
      const ArrayT<dArrayT>* E_List)
  {
    fE_List = E_List;
  }

  //
  // Return pointer to specified local array
  //
  inline const LocalArrayT*
  FSDEMatSupportViscoT::LocalArray(LocalArrayT::TypeT t) const
  {
    const LocalArrayT* pla = 0;

    switch (t) {

    case LocalArrayT::kESP:
      pla = fESP;
      break;

    default:
      //
      // Inherited
      //
      pla = FSMatSupportT::LocalArray(t);

    }

    return pla;

  }

  //
  // Set pointer to local array
  //
  inline void FSDEMatSupportViscoT::SetLocalArray(const LocalArrayT& array)
  {
    switch (array.Type()) {

    case LocalArrayT::kESP:
      fESP = &array;
      break;

    default:
      //
      // Inherited
      //
      FSMatSupportT::SetLocalArray(array);

    }

  }

  //
  // Nodal electric scalar potentials
  //
  inline const LocalArrayT*
  FSDEMatSupportViscoT::ScalarPotentials() const
  {
    return fESP;
  }

  //
  //
  //
  inline void FSDEMatSupportViscoT::SetScalarPotentials(
      const LocalArrayT& scalarPotentials)
  {
    fESP = &scalarPotentials;
  }


} //namespace Tahoe
