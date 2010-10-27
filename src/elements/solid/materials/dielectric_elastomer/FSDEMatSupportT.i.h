
namespace Tahoe {

  //
  //
  //
  inline FSDEMatSupportT::FSDEMatSupportT(int ndof, int nip) :
    FSMatSupportT(ndof, nip), fD_List(0), 
        fFSDielectricElastomer(0)
  {
  }

  //
  //
  //
  inline const FSDielectricElastomerT*
  FSDEMatSupportT::FSDielectricElastomer() const
  {
    return fFSDielectricElastomer;
  }

  //
  //
  //
  inline const dArrayT&
  FSDEMatSupportT::ElectricDisplacement() const
  {
    if (fD_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fD_List)[CurrIP()];
  }

  //
  //
  //
  inline const dArrayT&
  FSDEMatSupportT::ElectricDisplacement(int ip) const
  {
    if (fD_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fD_List)[ip];
  }

  //
  //
  //
  inline void FSDEMatSupportT::SetElectricDisplacement(
      const ArrayT<dArrayT>* D_List)
  {
    fD_List = D_List;
  }


} //namespace Tahoe
