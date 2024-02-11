//
// $Id: FSPZMatSupportT.i.h,v 1.2 2008-12-12 18:58:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2008/07/14 17:38:53  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:21:41  lxmota
// Piezoelectric material support. Initial sources.
//
//

namespace Tahoe {

  //
  //
  //
  inline FSPZMatSupportT::FSPZMatSupportT(int ndof, int nip) :
    FSMatSupportT(ndof, nip), fD_List(0), fDivPhi_List(0),
        fFSPiezoElectricSolid(0)
  {
  }

  //
  //
  //
  inline const FSPiezoElectricSolidT*
  FSPZMatSupportT::FSPiezoElectricSolid() const
  {
    return fFSPiezoElectricSolid;
  }

  //
  //
  //
  inline const dArrayT&
  FSPZMatSupportT::ElectricDisplacement() const
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
  FSPZMatSupportT::ElectricDisplacement(int ip) const
  {
    if (fD_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fD_List)[ip];
  }

  //
  //
  //
  inline void FSPZMatSupportT::SetElectricDisplacement(
      const ArrayT<dArrayT>* D_List)
  {
    fD_List = D_List;
  }

  //
  //
  //
  inline double FSPZMatSupportT::DivergenceVectorPotential() const
  {
    if (fDivPhi_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fDivPhi_List)[CurrIP()];
  }

  //
  //
  //
  inline double FSPZMatSupportT::DivergenceVectorPotential(int ip) const
  {
    if (fDivPhi_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fDivPhi_List)[ip];
  }

  //
  //
  //
  inline void FSPZMatSupportT::SetDivergenceVectorPotential(
      const dArrayT* DivPhi_List)
  {
    fDivPhi_List = DivPhi_List;
  }

  //
  // Return pointer to specified local array
  //
  inline const LocalArrayT*
  FSPZMatSupportT::LocalArray(LocalArrayT::TypeT t) const
  {

    const LocalArrayT* pla = 0;

    switch (t) {

    case LocalArrayT::kEVP:
      pla = fEVP;
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
  inline void FSPZMatSupportT::SetLocalArray(const LocalArrayT& array)
  {

    switch (array.Type()) {

    case LocalArrayT::kEVP:
      fEVP = &array;
      break;

    default:
      //
      // Inherited
      //
      FSMatSupportT::SetLocalArray(array);

    }

  }

  //
  // Nodal electric vector potentials
  //
  inline const LocalArrayT*
  FSPZMatSupportT::VectorPotentials() const
  {

    return fEVP;

  }

  //
  //
  //
  inline void FSPZMatSupportT::SetVectorPotentials(
      const LocalArrayT& vectorPotentials)
  {

    fEVP = &vectorPotentials;

  }

} //namespace Tahoe
