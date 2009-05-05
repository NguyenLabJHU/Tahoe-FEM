//
// $Id: FSHuWashizuUSCT.cpp,v 1.1 2009-05-05 23:33:42 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/12/12 18:56:58  amota
// Moved to new localtion.
//
// Revision 1.1  2008/07/14 16:13:25  lxmota
// Initial sources (disabled for now)
//
//

#include "FSSolidMatT.h"
#include "FSHuWashizuUSCT.h"
#include "FSMatSupportT.h"
#include "ParameterContainerT.h"
#include "OutputSetT.h"

//
// materials lists
//
#include "FSSolidMatList1DT.h"
#include "FSSolidMatList2DT.h"
#include "FSSolidMatList3DT.h"

namespace Tahoe {

  //
  //
  //
  FSHuWashizuUSCT::~FSHuWashizuUSCT()
  {

    if (0 != fFSMatSupport) delete fFSMatSupport;
    if (0 != fShapesMixed) delete fShapesMixed;

  }

  //
  // specify parameters needed by the interface
  //
  void FSHuWashizuUSCT::DefineParameters(ParameterListT& list) const
  {
    // inherited
    FiniteStrainT::DefineParameters(list);

    // additional fields
    list.AddParameter(ParameterT::Integer, "number_mixed_nodes");
    list.AddParameter(ParameterT::Integer, "number_mixed_int_pts");
    list.AddParameter(ParameterT::Boolean, "use_mixed_fields");

  }

  //
  // accept parameter list
  //
  void FSHuWashizuUSCT::TakeParameterList(const ParameterListT& list)
  {

    //
    //
    //
    fNumberMixedNodes = list.GetParameter("number_mixed_nodes");
    fNumberMixedIP = list.GetParameter("number_mixed_int_pts");
    fUseMixed = list.GetParameter("use_mixed_fields");

    //
    // inherited
    //
    FiniteStrainT::TakeParameterList(list);

    // Allocate workspace
    Initialize();
    Workspace();

  }

  //
  // initialize shape functions
  //
  void FSHuWashizuUSCT::SetShape()
  {

    // inherited
    if (fShapes == 0) {
      FiniteStrainT::SetShape();
    }

    // construct mixed shape functions
    if (fShapesMixed == 0) {

      const int nmn = NumberMixedNodes();
      const int nsd = NumSD();

      fLocMixedInitCoords.Dimension(nmn, nsd);
      fLocMixedInitCoords.SetType(LocalArrayT::kInitCoords);
      fShapesMixed
          = new MixedShapeFunctionT(*fShapes, GeometryCode(), NumberMixedIP(),
              fLocMixedInitCoords, false);
      if (0 == fShapesMixed) throw ExceptionT::kOutOfMemory;

      // initialize
      fShapesMixed->Initialize();
    }

  }

  //
  // form shape functions and derivatives
  //
  void FSHuWashizuUSCT::SetGlobalShape()
  {

    //
    // inherited
    //
    FiniteStrainT::SetGlobalShape();

    fShapes->SetIPParamCoords();

    fLocMixedInitCoords = GetMixedNodalPositions();
    fLocMixedInitCoords.SetType(LocalArrayT::kInitCoords);

    fShapesMixed->SetDerivatives();
    fShapesMixed->SetIPParamCoords();

    ComputeInverseWeightedVolume();

  }

  //
  // Determine the position of the mixed nodes in terms of the
  // positions of the primary nodes and primary and
  // mixed interpolation schemes.
  //
  LocalArrayT FSHuWashizuUSCT::GetMixedNodalPositions()
  {

    const int nen = NumElementNodes();
    const int nmn = NumberMixedNodes();
    const int nsd = NumSD();

    assert(nen >= nmn);

    LocalArrayT mixedPositions;
    mixedPositions.Dimension(nmn, nsd);
    mixedPositions.SetType(LocalArrayT::kInitCoords);

    // If only one mixed node, calculate position by averaging all primary
    // positions, otherwise take the position of the first nmn primary
    // nodes.
    if (1 == nmn) {
      dArrayT center;
      fLocInitCoords.Average(center);
      for (int j = 0; j < nsd; ++j) {
        mixedPositions(0, j) = center[j];
      }
    } else {
      for (int i = 0; i < nmn; ++i) {
        for (int j = 0; j < nsd; ++j) {
          mixedPositions(i, j) = fLocInitCoords(i, j);
        }
      }
    }

    return mixedPositions;

  }

  //
  // write all current element information to the stream
  //
  void FSHuWashizuUSCT::CurrElementInfo(ostream& out) const
  {

    //
    // inherited
    //
    FiniteStrainT::CurrElementInfo(out);

    out << std::endl << "Mixed Cauchy-Green deformation at node:" << std::endl;

    for (int i = 0; i < fC_nodal.Rows(); ++i) {

      out << " node: " << i + 1 << std::endl << fC_nodal[i] << std::endl;

    }

    out << std::endl;

    out << std::endl << "Mixed 2PK stress at node:" << std::endl;

    for (int i = 0; i < fS_nodal.Rows(); ++i) {

      out << " node: " << i + 1 << std::endl << fS_nodal[i] << std::endl;

    }

    out << std::endl;

  }

  //
  // Initialize local arrays
  //
  void FSHuWashizuUSCT::SetLocalArrays()
  {

    //
    // Inherited
    //
    FiniteStrainT::SetLocalArrays();

    //
    // Allocate storage space
    //
    const int nmn = NumberMixedNodes();
    const int ndof = dSymMatrixT::NumValues(NumDOF());

    fLocC.Dimension(nmn, ndof);
    fLocS.Dimension(nmn, ndof);

  }

  //
  // increment current element
  //
  bool FSHuWashizuUSCT::NextElement()
  {

    bool isThereNext = FiniteStrainT::NextElement();

    if (isThereNext == true) {

      const int index = CurrentElement().MaterialNumber();

      ContinuumMaterialT* pMaterial = (*fMaterialList)[index];

      fCurrMaterial = dynamic_cast<FSSolidMatT*> (pMaterial);

    }

    return isThereNext;

  }

  //
  //
  //
  void FSHuWashizuUSCT::Initialize()
  {

    bool nothingToDo = (fMaterialList == 0)
        || (fMaterialList->HasHistoryMaterials() == false);

    if (nothingToDo) {
      return;
    }

    // initialize internal variables for each integration point
    Top();

    while (NextElement() == true) {

      ElementCardT& element = CurrentElement();

      SetShape();

      ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];
      iArrayT iidPoint = pmat->InitialIntegerData();
      dArrayT iddPoint = pmat->InitialDoubleData();

      const int nip = fShapes->NumIP();

      iArrayT iidElement;
      const int nid = iidPoint.Length();
      if (nid > 0) {

        iidElement.Dimension(nid * nip);

        for (int ip = 0; ip < nip; ++ip) {
          iidElement.CopyIn(nid * ip, iidPoint);
        }

      }

      dArrayT iddElement;
      const int ndd = iddPoint.Length();
      if (ndd > 0) {

        iddElement.Dimension(ndd * nip);

        for (int ip = 0; ip < nip; ++ip) {
          iddElement.CopyIn(ndd * ip, iddPoint);
        }

      }

      element.Dimension(nid * nip, ndd * nip);
      element.IntegerData() = iidElement;
      element.DoubleData() = iddElement;

    }

  }

  //
  //
  //
  void FSHuWashizuUSCT::InitStep()
  {

    FiniteStrainT::InitStep();

    bool nothingToDo = (fMaterialList == 0)
        || (fMaterialList->HasHistoryMaterials() == false);

    if (nothingToDo) {
      return;
    }

    // initialize internal variables for each integration point
    // (differs from implementation elsewhere in Tahoe)
    Top();
    while (NextElement() == true) {

      ElementCardT& element = CurrentElement();

      SetShape();

      if (element.IsAllocated() == 0) {
        continue;
      }

      ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];
      dArrayT& elementDoubleData = element.DoubleData();
      iArrayT& elementIntegerData = element.IntegerData();

      const int nip = NumIP();
      for (int ip = 0; ip < nip; ++ip) {
        pmat->BeginStep(elementIntegerData, elementDoubleData, nip, ip);
      }

    }

  }

  //
  //
  //
  void FSHuWashizuUSCT::CloseStep()
  {

    FiniteStrainT::CloseStep();

    bool nothingToDo = (fMaterialList == 0)
        || (fMaterialList->HasHistoryMaterials() == false);

    if (nothingToDo) {
      return;
    }

    // update internal variables for each integration point
    // (differs from implementation elsewhere in Tahoe)
    Top();

    while (NextElement() == true) {

      ElementCardT& element = CurrentElement();

      if (element.IsAllocated() == 0) {
        continue;
      }

      ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];
      dArrayT& elementDoubleData = element.DoubleData();
      iArrayT& elementIntegerData = element.IntegerData();
      const int nip = fShapes->NumIP();
      fShapes->TopIP();

      while (fShapes->NextIP() != 0) {

        const int ip = fShapes->CurrIP();

        pmat->UpdateHistory(elementIntegerData, elementDoubleData, nip, ip);
        pmat->EndStep(elementIntegerData, elementDoubleData, nip, ip);

      }

    }

  }

  //
  // resets to the last converged solution
  //
  GlobalT::RelaxCodeT FSHuWashizuUSCT::ResetStep()
  {

    GlobalT::RelaxCodeT relax = FiniteStrainT::ResetStep();

    bool nothingToDo = (fMaterialList == 0)
        || (fMaterialList->HasHistoryMaterials() == false);

    if (nothingToDo) {
      return relax;
    }

    // update material internal variables for each integration point
    // (differs from implementation elsewhere in Tahoe)
    Top();

    while (NextElement() == true) {

      ElementCardT& element = CurrentElement();

      if (element.IsAllocated() == 0) {
        continue;
      }

      ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];
      dArrayT& elementDoubleData = element.DoubleData();
      iArrayT& elementIntegerData = element.IntegerData();
      const int nip = fShapes->NumIP();
      fShapes->TopIP();

      while (fShapes->NextIP() != 0) {

        const int ip = fShapes->CurrIP();

        pmat->ResetHistory(elementIntegerData, elementDoubleData, nip, ip);

      }

    }

    return relax;

  }

  //
  // Geometric stiffness
  //
  void FSHuWashizuUSCT::AccumulateGeometricStiffness(dMatrixT& Kg,
      const dArray2DT& DNaX, dSymMatrixT& S)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();

    dMatrixT GradNa(nsd, nen);
    fShapes->GradNa(DNaX, GradNa);

    dMatrixT DNaS(nen);

    dMatrixT Sm(nsd);
    S.ToMatrix(Sm);

    DNaS.MultQTBQ(GradNa, Sm);

    dMatrixT SI(nsd);

    for (int i = 0; i < nen; ++i) {

      for (int j = 0; j < nen; ++j) {

        SI.Identity(DNaS(i, j));

        Kg.AddBlock(nsd * i, nsd * j, SI);

      }

    }

  }

  //
  // strain-displacement operator at current point
  // of primary integration scheme (interpolation
  // function derivatives are assumed to be evaluated there)
  //
  void FSHuWashizuUSCT::Set_B(const dArray2DT& DNaX, dMatrixT& B)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    const dMatrixT F = FiniteStrainT::DeformationGradient();

    B = StrainDisplacementOperator(F, DNaX);

  }

  //
  // element stiffness matrix
  //
  void FSHuWashizuUSCT::FormStiffness(double constK)
  {

    //
    // matrix format
    //
    dMatrixT::SymmetryFlagT format = (fLHS.Format()
        == ElementMatrixT::kNonSymmetric)
        ? dMatrixT::kWhole
        : dMatrixT::kUpperOnly;

    //
    // integrate over element
    //
    const int nsd = NumSD();
    const int nen = NumElementNodes();

    fGeometricTangent = 0.0;
    fMaterialTangent = 0.0;

    fShapes->TopIP();

    while (fShapes->NextIP() != 0) {

      //
      // scale/weighting factor for integration
      //
      const double w = constK * fShapes->IPDet() * fShapes->IPWeight();

      //
      // get material tangent moduli
      //
      dMatrixT C = fCurrMaterial->C_IJKL();

      dArrayT primaryIP(nsd);
      fShapes->IPParamCoords(primaryIP);

      dSymMatrixT S_bar = fUseMixed == true
          ? MixedStress(primaryIP)
          : Stress(primaryIP);

      C *= (0.25 * w);
      S_bar *= w;

      //
      // prepare derivatives of interpolation functions
      //
      const dArray2DT & DNaX = fShapes->Derivatives_X();

      //
      // Primary variables -> gradient operators
      //

      dMatrixT B_bar = fUseMixed == true
          ? MixedStrainDisplacementOperator(primaryIP)
          : StrainDisplacementOperator(primaryIP);

      //
      // Material stiffness
      //
      fMaterialTangent.MultQTBQ(B_bar, C, format, dMatrixT::kAccumulate);

      //
      // Geometric stiffness
      //
      AccumulateGeometricStiffness(fGeometricTangent, DNaX, S_bar);

    }

    //
    // Add geometric stiffness
    //
    fMaterialTangent.Expand(fGeometricTangent, 1, dMatrixT::kAccumulate);

    //
    // Assemble into element stiffness matrix
    //
    fLHS.AddBlock(0, 0, fMaterialTangent);

  }

  //
  // internal force
  //
  void FSHuWashizuUSCT::FormKd(double constK)
  {
    //
    //
    //
    const int neq = NumElementNodes() * NumDOF();
    const int nsd = NumSD();

    fShapes->TopIP();

    while (fShapes->NextIP() != 0) {

      //
      // integration weight
      //
      const double w = constK * fShapes->IPDet() * fShapes->IPWeight();

      //
      // Deformation power P_D := \int S : 0.5 \dot(C) dV
      //
      dSymMatrixT S = fCurrMaterial->S_IJ();
      S *= (0.5 * w);

      dArrayT primaryIP(nsd);
      fShapes->IPParamCoords(primaryIP);

      dMatrixT B_bar = fUseMixed == true
          ? MixedStrainDisplacementOperator(primaryIP)
          : StrainDisplacementOperator(primaryIP);

      B_bar.MultTx(S, fRHS, 1.0, dMatrixT::kAccumulate);

    }

  }

} // namespace Tahoe

