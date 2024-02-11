//
// $Id: FSNeoHookePZLinT.cpp,v 1.2 2008-12-12 18:58:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2008/07/14 17:37:44  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:10:49  lxmota
// Piezoelectric material. Initial sources.
//
//

#include "ExceptionT.h"
#include "FSNeoHookePZLinT.h"

namespace Tahoe {

  //
  //
  //
  static const char NPZ[] = "Neohookean-piezoelectric";
  const bool FSNeoHookePZLinT::dependsOnC = true;
  const char* FSNeoHookePZLinT::Name = NPZ;

  //
  //
  //
  void FSNeoHookePZLinT::Initialize()
  {

    fElectricPermittivity = 0.0;
    fPiezoelectricTensor = dMatrixT(ElectricalDim(), StrainDim());
    fPiezoelectricTensor = 0.0;
    fPenaltyCoefficient = 0.0;

  }

  //
  //
  //
  void FSNeoHookePZLinT::DefineParameters(ParameterListT& list) const
  {

    FSIsotropicMatT::DefineParameters(list);

    list.AddParameter(fElectricPermittivity, "epsilon");
    list.AddParameter(fPenaltyCoefficient, "penalty");

    char p[4] = "h11";

    for (int i = 0; i < ElectricalDim(); ++i) {

      p[1] = '1' + i;

      for (int j = 0; j < StrainDim(); ++j) {

        p[2] = '1' + j;

        const double & gij = fPiezoelectricTensor(i, j);
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
  void FSNeoHookePZLinT::TakeParameterList(const ParameterListT& list)
  {

    FSIsotropicMatT::TakeParameterList(list);

    fElectricPermittivity = list.GetParameter("epsilon");
    fPenaltyCoefficient = list.GetParameter("penalty");

    char p[4] = "h11";

    for (int i = 0; i < ElectricalDim(); ++i) {

      p[1] = '1' + i;

      for (int j = 0; j < StrainDim(); ++j) {

        p[2] = '1' + j;

        fPiezoelectricTensor(i, j) = list.GetParameter(p);

      }

    }

    //
    // check
    //
    if (fElectricPermittivity < -kSmall) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::TakeParameterList",
          "expecting a non-negative epsilon: %e", fElectricPermittivity);
    }

    if (fPenaltyCoefficient < -kSmall) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::TakeParameterList",
          "expecting a non-negative penalty coefficient: %e",
          fPenaltyCoefficient);
    }

  }

  //
  // information about subordinate parameter lists
  //
  void FSNeoHookePZLinT::DefineSubs(SubListT& sub_list) const
  {
    FSIsotropicMatT::DefineSubs(sub_list);
    return;
  }

  //
  // Set piezoelectric constants
  //
  void FSNeoHookePZLinT::SetPiezoelectricConstant(int i, int j, double gij)
  {

    bool validI = (0 <= i && i < StrainDim()) == true;
    bool validJ = (0 <= j && j < ElectricalDim()) == true;

    if (validI == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
          "invalid piezoelectric index i: %d", i);
    }

    if (validJ == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
          "invalid piezoelectric index j: %d", j);
    }

    fPiezoelectricTensor(i, j) = gij;

  }

  //
  // Get piezoelectric constants
  //
  double FSNeoHookePZLinT::GetPiezoelectricConstant(int i, int j) const
  {

    bool validI = (0 <= i && i < StrainDim()) == true;
    bool validJ = (0 <= j && j < ElectricalDim()) == true;

    if (validI == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
          "invalid piezoelectric index i: %d", i);
    }

    if (validJ == false) {
      ExceptionT::BadInputValue("FSNeoHookePZLinT::setPiezoelectricConstant",
          "invalid piezoelectric index j: %d", j);
    }

    return fPiezoelectricTensor(i, j);

  }

  //
  //
  //
  const dSymMatrixT FSNeoHookePZLinT::StressNum(const dSymMatrixT& C,
      const dArrayT& D) const
  {
    const double eps = 1e-6;
    dSymMatrixT Cm(ManifoldDim());
    dSymMatrixT Cp(ManifoldDim());
    dSymMatrixT S(ManifoldDim());

    for (int i = 0; i < StrainDim(); ++i) {
      Cm = C;
      Cm[i] -= i < ManifoldDim() ? (0.5 * eps) : (0.25 * eps);
      const double Wm = EnergyDensity(Cm, D);

      Cp = C;
      Cp[i] += i < ManifoldDim() ? (0.5 * eps) : (0.25 * eps);
      const double Wp = EnergyDensity(Cp, D);

      S[i] = 2.0 * (Wp - Wm) / eps;
    }

    return S;

  }

  //
  //
  //
  const dMatrixT FSNeoHookePZLinT::TangentMechanicalNum(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    const double eps = 1e-6;
    dSymMatrixT Cm(ManifoldDim());
    dSymMatrixT Cp(ManifoldDim());
    dMatrixT CC(StrainDim());

    for (int i = 0; i < StrainDim(); ++i) {
      Cm = C;
      Cm[i] -= i < ManifoldDim() ? (0.5 * eps) : (0.25 * eps);
      const dSymMatrixT Sm = Stress(Cm, D);

      Cp = C;
      Cp[i] += i < ManifoldDim() ? (0.5 * eps) : (0.25 * eps);
      const dSymMatrixT Sp = Stress(Cp, D);

      for (int j = 0; j < StrainDim(); ++j) {
        CC(i, j) = 2.0 * (Sp[j] - Sm[j]) / eps;
      }

    }

    return CC;

  }

} //namespace Tahoe
