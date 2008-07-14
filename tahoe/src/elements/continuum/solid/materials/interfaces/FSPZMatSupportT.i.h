//
// $Id: FSPZMatSupportT.i.h,v 1.2 2008-07-14 17:38:53 lxmota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/06/16 18:21:41  lxmota
// Piezoelectric material support. Initial sources.
//
//

namespace Tahoe{
  
  //
  //
  //
  inline
  FSPZMatSupportT::FSPZMatSupportT(int ndof, int nip) :
    FSMatSupportT(ndof, nip),
    fD_List(0),
    fD_last_List(0),
    fFSPiezoElectricSolid(0)
  {
  }

  //
  //
  //
  inline const FSPiezoElectricSolidT*
  FSPZMatSupportT::FSPiezoElectricSolid() const {

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
  inline const dArrayT&
  FSPZMatSupportT::ElectricDisplacement_last() const
  {

    if (fD_List == 0) {

      throw ExceptionT::kGeneralFail;

    }

    return (*fD_last_List)[CurrIP()];

  }
    
  //
  //
  //
  inline const dArrayT&
  FSPZMatSupportT::ElectricDisplacement_last(int ip) const
  {

    if (fD_List == 0) {

      throw ExceptionT::kGeneralFail;

    }

    return (*fD_last_List)[ip];

  }

  //
  //
  //
  inline void
  FSPZMatSupportT::SetElectricDisplacement(const ArrayT<dArrayT>* D_List)
  {

    fD_List = D_List;

  }
    
  //
  //
  //
  inline void
  FSPZMatSupportT::SetElectricDisplacement_last(const ArrayT<dArrayT>*
						D_last_List)
  {

    fD_last_List = D_last_List;

  }
  
  //
  // Return pointer to specified local array
  //
  inline const LocalArrayT*
  FSPZMatSupportT::LocalArray(LocalArrayT::TypeT t) const
  {

    const LocalArrayT* pla = 0;

    switch (t) {

    case LocalArrayT::kLastEVP:
      pla = fLastEVP;
      break;

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
  inline void
  FSPZMatSupportT::SetLocalArray(const LocalArrayT& array)
  {

    switch (array.Type()) {

    case LocalArrayT::kLastEVP:
      fLastEVP = &array;
      break;

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
  // Last nodal electric vector potentials
  //
  inline const LocalArrayT*
  FSPZMatSupportT::LastVectorPotentials() const
  {

    return fLastEVP;

  }

  //
  //
  //
  inline void
  FSPZMatSupportT::SetVectorPotentials(const LocalArrayT& vectorPotentials)
  {

    fEVP = &vectorPotentials;

  }
  
  //
  //
  //
  inline void
  FSPZMatSupportT::SetLastVectorPotentials(const LocalArrayT&
					 last_vectorPotentials)
  {

    fLastEVP = &last_vectorPotentials;

  }
  
} //namespace Tahoe
