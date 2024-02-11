
namespace Tahoe {

  //
  //
  //
  inline FSDEMatSupport2DT::FSDEMatSupport2DT(int ndof, int nip) :
    FSMatSupportT(ndof, nip), fE_List(0), fFSDielectricElastomer2D(0)
  {

  }

  //
  //
  //
  inline const FSDielectricElastomer2DT*
  FSDEMatSupport2DT::FSDielectricElastomer2D() const
  {
    return fFSDielectricElastomer2D;
  }

  inline const dArrayT&
  FSDEMatSupport2DT::ElectricField() const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[CurrIP()];
  }

  inline const dArrayT&
  FSDEMatSupport2DT::ElectricField(int ip) const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[ip];
  }

  inline void FSDEMatSupport2DT::SetElectricField(
      const ArrayT<dArrayT>* E_List)
  {
    fE_List = E_List;
  }

  //
  // Return pointer to specified local array
  //
  inline const LocalArrayT*
  FSDEMatSupport2DT::LocalArray(LocalArrayT::TypeT t) const
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
  inline void FSDEMatSupport2DT::SetLocalArray(const LocalArrayT& array)
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
  FSDEMatSupport2DT::ScalarPotentials() const
  {
    return fESP;
  }

  //
  //
  //
  inline void FSDEMatSupport2DT::SetScalarPotentials(
      const LocalArrayT& scalarPotentials)
  {
    fESP = &scalarPotentials;
  }


} //namespace Tahoe
