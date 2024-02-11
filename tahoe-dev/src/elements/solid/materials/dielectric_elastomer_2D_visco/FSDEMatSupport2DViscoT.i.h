
namespace Tahoe {

  //
  //
  //
  inline FSDEMatSupport2DViscoT::FSDEMatSupport2DViscoT(int ndof, int nip) :
    FSMatSupportT(ndof, nip), fE_List(0), fFSDielectricElastomer2DVisco(0)
  {

  }

  //
  //
  //
  inline const FSDielectricElastomer2DViscoT*
  FSDEMatSupport2DViscoT::FSDielectricElastomer2DVisco() const
  {
    return fFSDielectricElastomer2DVisco;
  }

  inline const dArrayT&
  FSDEMatSupport2DViscoT::ElectricField() const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[CurrIP()];
  }

  inline const dArrayT&
  FSDEMatSupport2DViscoT::ElectricField(int ip) const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[ip];
  }

  inline void FSDEMatSupport2DViscoT::SetElectricField(
      const ArrayT<dArrayT>* E_List)
  {
    fE_List = E_List;
  }

  //
  // Return pointer to specified local array
  //
  inline const LocalArrayT*
  FSDEMatSupport2DViscoT::LocalArray(LocalArrayT::TypeT t) const
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
  inline void FSDEMatSupport2DViscoT::SetLocalArray(const LocalArrayT& array)
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
  FSDEMatSupport2DViscoT::ScalarPotentials() const
  {
    return fESP;
  }

  //
  //
  //
  inline void FSDEMatSupport2DViscoT::SetScalarPotentials(
      const LocalArrayT& scalarPotentials)
  {
    fESP = &scalarPotentials;
  }


} //namespace Tahoe
