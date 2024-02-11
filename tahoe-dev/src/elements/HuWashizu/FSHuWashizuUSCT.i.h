//
// $Id: FSHuWashizuUSCT.i.h,v 1.1 2009-05-05 23:33:42 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2009/04/02 00:49:52  amota
// Changes in handling of internal variables.
//
// Revision 1.1  2008/12/12 18:56:58  amota
// Moved to new localtion.
//
// Revision 1.1  2008/07/14 16:13:25  lxmota
// Initial sources (disabled for now)
//
//

namespace Tahoe {

  //
  //
  //
  inline FSHuWashizuUSCT::FSHuWashizuUSCT(const ElementSupportT& support) :
    FiniteStrainT(support), fFSMatSupport(0), fCurrMaterial(0), fShapesMixed(0)
  {

    SetName("Hu_Washizu_USC");

  }

  //
  //
  //
  inline void FSHuWashizuUSCT::Workspace()
  {

    const int nen = NumElementNodes();
    const int nmn = NumberMixedNodes();
    const int nsd = NumSD();
    const int nel = nen * nsd;
    const int nme = nen * nsd;

    fMaterialTangent.Dimension(nme, nme);
    fGeometricTangent.Dimension(nme, nme);

    fC_nodal.Dimension(nmn, dSymMatrixT::NumValues(nsd));
    fS_nodal.Dimension(nmn, dSymMatrixT::NumValues(nsd));

  }

  //
  //
  //
  inline const ShapeFunctionT&
  FSHuWashizuUSCT::MixedShapeFunction() const
  {
    return *fShapesMixed;
  }

  //
  //
  //
  inline int FSHuWashizuUSCT::NumberMixedNodes() const
  {
    return fNumberMixedNodes;
  }

  //
  //
  //
  inline int FSHuWashizuUSCT::NumberMixedIP() const
  {
    return fNumberMixedIP;
  }

  //
  // Compute deformation gradient at an arbitrary location
  // given in parametric coordinates.
  //
  inline dMatrixT FSHuWashizuUSCT::DeformationGradient(const dArrayT& paramPt) const
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();

    dMatrixT F(nsd);
    dMatrixT DX_xi(nsd);
    dMatrixT Dx_xi(nsd);
    dArrayT Na;
    dArray2DT DNa;

    // get values of interpolation functions and derivatives at
    // given parametric point and compute F there
    fShapes->GradU(fLocDisp, Dx_xi, paramPt, Na, DNa);
    fShapes->GradU(fLocInitCoords, DX_xi, paramPt, Na, DNa);
    Dx_xi += DX_xi;

    F.MultAB(Dx_xi, DX_xi.Inverse());

    return F;

  }

  //
  //
  //
  inline dSymMatrixT FSHuWashizuUSCT::RightCauchyGreenDeformation(
      const dArrayT& paramPt) const
  {
    const dMatrixT F = DeformationGradient(paramPt);
    const int nsd = NumSD();
    dMatrixT FTF(nsd);
    FTF.MultATB(F, F);
    dSymMatrixT C(nsd);
    C.Symmetrize(FTF);

    return C;
  }

  //
  // Compute 2nd PK stress at an arbitrary location
  // given in parametric coordinates.
  //
  inline dSymMatrixT FSHuWashizuUSCT::Stress2PK(const dArrayT& paramPt)
  {
    const dSymMatrixT C = RightCauchyGreenDeformation(paramPt);
    fStress = fCurrMaterial->S_IJ(C);

    return fStress;
  }

  //
  //
  //
  inline void FSHuWashizuUSCT::ComputeInverseWeightedVolume()
  {

    const int nmn = NumberMixedNodes();

    //
    // integration scheme
    //
    const double* Det = fShapesMixed->IPDets();
    const double* Weight = fShapesMixed->IPWeights();

    fShapesMixed->TopIP();

    fInvH.Dimension(nmn);
    fInvH = 0.0;

    while (fShapesMixed->NextIP() != 0) {

      //
      // integration weight
      //
      const double w = (*Weight++) * (*Det++);

      dMatrixT PsiI_x_PsiJ(nmn);
      PsiI_x_PsiJ.Outer(fShapesMixed->IPShapeU(), fShapesMixed->IPShapeU());

      PsiI_x_PsiJ *= w;

      fInvH += PsiI_x_PsiJ;

    }

    fInvH.Inverse();

  }

  //
  // Given deformation gradient and interpolation functions
  // derivatives at arbitrary parametric point, compute
  // strain displacement operator.
  //
  inline dMatrixT FSHuWashizuUSCT::StrainDisplacementOperator(
      const dMatrixT& F, const dArray2DT& DNaX) const
  {
    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int StrainDim = dSymMatrixT::NumValues(nsd);

    dMatrixT B(StrainDim, nsd * nen);

    for (int ij = 0; ij < StrainDim; ++ij) {

      int i;
      int j;

      dSymMatrixT::ExpandIndex(nsd, ij, i, j);

      for (int a = 0; a < nen; ++a) {

        for (int k = 0; k < nsd; ++k) {

          const int ak = a * nsd + k;

          B(ij, ak) = DNaX(i, a) * F(k, j) + DNaX(j, a) * F(k, i);

          // Shear components doubled to conform with Voigt
          // convention
          if (ij >= nsd) {
            B(ij, ak) *= 2.0;
          }


        }

      }

    }

    return B;

  }

  //
  // Strain-displacement operator at an arbitrary parametric point.
  //
  inline dMatrixT FSHuWashizuUSCT::StrainDisplacementOperator(
      const dArrayT& paramPt) const
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();

    dMatrixT DX_xi(nsd);
    dArrayT Na;
    dArray2DT DNa;

    fShapes->GradU(fLocInitCoords, DX_xi, paramPt, Na, DNa);

    dArray2DT DNaX(nsd, nen);
    DX_xi.Inverse().Transpose();

    DNaX = 0.0;
    for (int i = 0; i < nsd; ++i) {
      for (int j = 0; j < nen; ++j) {
        for (int k = 0; k < nsd; ++k) {
          DNaX(i,j) += DX_xi(i,k) * DNa(k,j);
        }
      }
    }

    const dMatrixT F = DeformationGradient(paramPt);

    dMatrixT B = StrainDisplacementOperator(F, DNaX);

    return B;

  }

  //
  // Mixed strain-displacement operator at arbitrary parametric point.
  //
  inline dMatrixT FSHuWashizuUSCT::MixedStrainDisplacementOperator(
      const dArrayT& paramPt) const
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int StrainDim = dSymMatrixT::NumValues(nsd);
    const int nmn = NumberMixedNodes();

    //
    // mixed integration scheme
    //
    const double* Det = fShapesMixed->IPDets();
    const double* Weight = fShapesMixed->IPWeights();

    dMatrixT Beta(nmn, StrainDim * nsd * nen);

    Beta = 0.0;

    fShapesMixed->TopIP();

    while (fShapesMixed->NextIP() != 0) {

      //
      // integration weight
      //
      const double w = (*Weight++) * (*Det++);

      dArrayT mixedIP(nsd);
      fShapesMixed->IPParamCoords(mixedIP);

      dMatrixT B = StrainDisplacementOperator(mixedIP);

      dArrayT Psi;
      Psi.Alias(nmn, fShapesMixed->IPShapeU());

      dMatrixT Psi_x_B(nmn, StrainDim * nsd * nen);

      Psi_x_B.Outer(Psi, B);
      Psi_x_B *= w;
      Beta += Psi_x_B;

    }

    dMatrixT invHBeta(nmn, StrainDim * nsd * nen);
    invHBeta.MultAB(fInvH, Beta);

    dArrayT Psi(nmn);

    fShapesMixed->EvaluateShapeFunctions(paramPt, Psi);

    dMatrixT B_bar(StrainDim, nsd * nen);

    invHBeta.MultTx(Psi.Pointer(), B_bar.Pointer());

    return B_bar;

  }

  //
  // Strain measure at arbitrary parametric point.
  //
  inline dSymMatrixT FSHuWashizuUSCT::Strain(const dArrayT& paramPt)
  {
    return RightCauchyGreenDeformation(paramPt);
  }

  //
  // Mixed strain measure at arbitrary parametric point.
  //
  inline dSymMatrixT FSHuWashizuUSCT::MixedStrain(const dArrayT& paramPt)
  {

    const int nsd = NumSD();
    const int nmn = NumberMixedNodes();

    //
    // integration scheme
    //
    const double* Det = fShapesMixed->IPDets();
    const double* Weight = fShapesMixed->IPWeights();

    dMatrixT Xi(nmn, dSymMatrixT::NumValues(nsd));
    Xi = 0.0;

    fShapesMixed->TopIP();

    while (fShapesMixed->NextIP() != 0) {

      //
      // integration weight
      //
      const double w = (*Weight++) * (*Det++);

      dArrayT mixedIP(nsd);
      fShapesMixed->IPParamCoords(mixedIP);

      dSymMatrixT C = RightCauchyGreenDeformation(mixedIP);

      dArrayT Psi;
      Psi.Alias(nmn, fShapesMixed->IPShapeU());

      dMatrixT Psi_x_C(nmn, dSymMatrixT::NumValues(nsd));

      Psi_x_C.Outer(Psi, C);
      Psi_x_C *= w;
      Xi += Psi_x_C;

    }

    fC_nodal.MultAB(fInvH, Xi);

    dSymMatrixT C_bar(nsd);

    dArrayT Psi(nmn);
    fShapesMixed->EvaluateShapeFunctions(paramPt, Psi);

    fC_nodal.MultTx(Psi, C_bar);

    return C_bar;

  }

  //
  // Stress at arbitrary parametric point.
  //
  inline dSymMatrixT FSHuWashizuUSCT::Stress(const dArrayT& paramPt)
  {
    return fCurrMaterial->S_IJ(Strain(paramPt));
  }

  //
  // Mixed stress at arbitrary parametric point.
  //
  inline dSymMatrixT FSHuWashizuUSCT::MixedStress(const dArrayT& paramPt)
  {

    const int nsd = NumSD();
    const int nmn = NumberMixedNodes();

    //
    // integration scheme
    //
    const double* Det = fShapesMixed->IPDets();
    const double* Weight = fShapesMixed->IPWeights();

    dMatrixT Sigma(nmn, dSymMatrixT::NumValues(nsd));
    Sigma = 0.0;

    fShapesMixed->TopIP();

    while (fShapesMixed->NextIP() != 0) {

      //
      // integration weight
      //
      const double w = (*Weight++) * (*Det++);

      dArrayT mixedIP(nsd);
      fShapesMixed->IPParamCoords(mixedIP);

      const int savedIP = fShapesMixed->CurrIP();
      dSymMatrixT C_bar = MixedStrain(mixedIP);
      dArrayT iv = InterpolateInternal(mixedIP);
      fShapesMixed->SetIP(savedIP);

      dSymMatrixT S = fCurrMaterial->S_IJ(C_bar);

      dArrayT Psi;
      Psi.Alias(nmn, fShapesMixed->IPShapeU());

      dMatrixT Psi_x_S(nmn, dSymMatrixT::NumValues(nsd));

      Psi_x_S.Outer(Psi, S);
      Psi_x_S *= w;
      Sigma += Psi_x_S;

    }

    fS_nodal.MultAB(fInvH, Sigma);

    dSymMatrixT S_bar(nsd);

    dArrayT Psi(nmn);
    fShapesMixed->EvaluateShapeFunctions(paramPt, Psi);

    fS_nodal.MultTx(Psi, S_bar);

    return S_bar;

  }

  //
  // Internal variables held by the material.
  //
  inline const dArrayT& FSHuWashizuUSCT::InternalVariables()
  {
    return CurrentElement().DoubleData();
  }

  //
  // Compute nodal values of internal variables by extrapolating
  // from integration points.
  //
  inline void FSHuWashizuUSCT::NodalInternalVariables()
  {

    const int nip = NumIP();
    const int nen = NumElementNodes();

    // Get internal variables at integration points.
    const dArrayT& aiv = InternalVariables();
    const int niv = aiv.Length() / nip;

    fNodalInternalVariables.Dimension(nen, niv);

    //Extrapolate one by one.
    for (int iv = 0; iv < niv; ++iv) {

      dArrayT ipiv(nip);
      dArrayT ndiv(nen);

      for (int ip = 0; ip < nip; ++ip) {
        ipiv[ip] = aiv[ip * niv + iv];
      }

      fShapes->ExtrapolateAll(ipiv, ndiv);

      for (int en = 0; en < nen; ++en) {
        fNodalInternalVariables(en,iv) = ndiv[en];
      }

    }

  }

  //
  // Interpolate internal variables to arbitrary parametric point.
  //
  inline dArrayT FSHuWashizuUSCT::InterpolateInternal(const dArrayT& paramPt) {

    const int nip = NumIP();
    const dArrayT& aiv = InternalVariables();
    const int niv = aiv.Length() / nip;
    NodalInternalVariables();

    dArrayT iv(niv);

    fShapes->InterpolateU(fNodalInternalVariables, iv, paramPt);

    return iv;

  }


} // namespace Tahoe
