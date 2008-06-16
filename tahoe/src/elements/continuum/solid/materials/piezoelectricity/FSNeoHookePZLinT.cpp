//
// $Id: FSNeoHookePZLinT.cpp,v 1.1 2008-06-16 18:10:49 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

#include "ExceptionT.h"
#include "FSNeoHookePZLinT.h"

namespace Tahoe {

  //
  //
  //
  static const char NPZ[] = "Neohookean-piezoelectric";
  const char* FSNeoHookePZLinT::Name = NPZ;

  //
  //
  //
  void
  FSNeoHookePZLinT::initialize()
  {

    const int nsd = NumSD();
    const int ElectricalDim = nsd;
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    fShearModulus = 0.0;
    fBulkModulus = 0.0;
    fElectricPermittivity = 0.0;
    fPiezoelectricTensor = dMatrixT(ElectricalDim, StrainDim);
    fPiezoelectricTensor = 0.0;

  }

  //
  //
  //
  void
  FSNeoHookePZLinT::DefineParameters(ParameterListT& list) const
  {
    
    FSIsotropicMatT::DefineParameters(list);

    list.AddParameter(fShearModulus, "mu");
    list.AddParameter(fBulkModulus, "kappa");
    list.AddParameter(fElectricPermittivity, "epsilon");
  
    const int nsd = NumSD();
    const int ElectricalDim = nsd;
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    char p[4] = "g00";

    for (int i = 0; i < ElectricalDim; ++i) {

      p[1] = '0' + i;

      for (int j = 0; j < StrainDim; ++j) {

        p[2] = '0' + j;

        double gij = fPiezoelectricTensor(i,j);
        list.AddParameter(gij, p);

      }

    }

        
    //
    // set the description
    //
    list.SetDescription("Psi(C)=0.5*mu*(I1bar-3)+0.25*kappa*(J^2-1-2*log(J))");

  }

  //
  //
  //
  void
  FSNeoHookePZLinT::TakeParameterList(const ParameterListT& list)
  {
    
    FSIsotropicMatT::TakeParameterList(list);
  
    fShearModulus = list.GetParameter("mu");
    fBulkModulus  = list.GetParameter("kappa");
    fElectricPermittivity = list.GetParameter("epsilon");

    const int nsd = NumSD();
    const int ElectricalDim = nsd;
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    char p[4] = "g00";

    for (int i = 0; i < ElectricalDim; ++i) {

      p[1] = '0' + i;

      for (int j = 0; j < StrainDim; ++j) {

        p[2] = '0' + j;

        fPiezoelectricTensor(i,j) = list.GetParameter(p);

      }

    }

    //
    // check
    //
    if (fShearModulus < -kSmall) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::TakeParameterList",
				"expecting a non-negative value mu: %e",
				fShearModulus);
    }
  
    if (fBulkModulus < -kSmall) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::TakeParameterList",
				"expecting a non-negative value kappa: %e",
				fBulkModulus);
    }

    if (fElectricPermittivity < -kSmall) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::TakeParameterList",
				"expecting a non-negative value epsilon: %e",
				fElectricPermittivity);
    }

  }

  //
  // information about subordinate parameter lists
  //
  void
  FSNeoHookePZLinT::DefineSubs(SubListT& sub_list) const
  {
    return;
  }

  //
  // Set piezoelectric constants
  //
  void
  FSNeoHookePZLinT::setPiezoelectricConstant(int i, int j, double gij)
  {

    const int nsd = NumSD();
    const int ElectricalDim = nsd;
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    bool validI = (0 <= i && i < StrainDim) == true;

    bool validJ = (0 <= j && j < ElectricalDim) == true;

    if (validI == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
				"invalid piezoelectric index i: %d",
				i);
    }
  
    if (validJ == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
				"invalid piezoelectric index j: %d",
				j);
    }

    fPiezoelectricTensor(i,j) = gij;

  }

  //
  // Get piezoelectric constants
  //
  double
  FSNeoHookePZLinT::getPiezoelectricConstant(int i, int j) const
  {

    const int nsd = NumSD();
    const int ElectricalDim = nsd;
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    bool validI = (0 <= i && i < StrainDim) == true;

    bool validJ = (0 <= j && j < ElectricalDim) == true;

    if (validI == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
				"invalid piezoelectric index i: %d",
				i);
    }
  
    if (validJ == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
				"invalid piezoelectric index j: %d",
				j);
    }

    return fPiezoelectricTensor(i,j);

  }

} //namespace Tahoe
