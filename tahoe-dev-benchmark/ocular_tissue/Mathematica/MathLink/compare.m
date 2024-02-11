Off[General::spell1]; Off[General::spell]; 
errchk = False; 
compare[alpha_, beta_, mu_, taulo_, taudif_, c_] := 
   Module[{}, TimeConstrained[If[errchk, Export["parameters.dat", 
        parameters]]; StressScale = 6.8947; nu = 0.49; gamma = nu/(1 - 2*nu); 
      tau1 = 10^taulo; tau2 = 10^(taulo + taudif); 
      Goo = 1/(1 + c*Log[tau2/tau1]); Wf[l_] = 
       alpha*(Exp[beta*(l - 1)] - beta*l); Wfp[l_] = D[Wf[l], l]; 
      Wfpp[l_] = D[Wfp[l], l]; Theta = {Cos[t], Sin[t]}; 
      MatrixForm[M = Outer[Times, Theta, Theta]]; 
      I4 = l1^2*Cos[t]^2 + l2^2*Sin[t]^2; sf1[l1_, l2_, l3_] = 
       2*Wfp[I4]*M[[1,1]]; sf2[l1_, l2_, l3_] = 2*Wfp[I4]*M[[2,2]]; 
      kf11[l1_, l2_, l3_] = 4*Wfpp[I4]*M[[1,1]]*M[[1,1]]; 
      kf12[l1_, l2_, l3_] = 4*Wfpp[I4]*M[[1,1]]*M[[2,2]]; 
      kf22[l1_, l2_, l3_] = 4*Wfpp[I4]*M[[2,2]]*M[[2,2]]; 
      H = (1./(2.*Pi))*NIntegrate[#1, {t, -Pi, Pi}, AccuracyGoal -> 10] & ; 
      J = l1*l2*l3; Um[l1_, l2_, l3_] = (mu/2)*(l1^2 + l2^2 + l3^2 - 3) + 
        (mu/(2*gamma))*(J^(-2*gamma) - 1); sm1[l1_, l2_, l3_] = 
       D[Um[l1, l2, l3], l1]/l1; sm2[l1_, l2_, l3_] = D[Um[l1, l2, l3], l2]/
        l2; sm3[l1_, l2_, l3_] = D[Um[l1, l2, l3], l3]/l3; 
      Cm11[l1_, l2_, l3_] = D[sm1[l1, l2, l3], l1]/l1; 
      Cm12[l1_, l2_, l3_] = D[sm1[l1, l2, l3], l2]/l2; 
      Cm13[l1_, l2_, l3_] = D[sm1[l1, l2, l3], l3]/l3; 
      Cm22[l1_, l2_, l3_] = D[sm2[l1, l2, l3], l2]/l2; 
      Cm23[l1_, l2_, l3_] = D[sm2[l1, l2, l3], l3]/l3; 
      Cm33[l1_, l2_, l3_] = D[sm3[l1, l2, l3], l3]/l3; 
      stress[l1m_, l2m_, l3m_] := {H[sf1[l1m, l2m, l3m]], 
         H[sf2[l1m, l2m, l3m]], 0} + {sm1[l1m, l2m, l3m], sm2[l1m, l2m, l3m], 
         sm3[l1m, l2m, l3m]}; stiffness[l1m_, l2m_, l3m_] := 
       Module[{CF11, CF12, CF22, CM11, CM12, CM13, CM22, CM23, CM33, C}, 
        CF11 = H[kf11[l1m, l2m, l3m]]; CF12 = H[kf12[l1m, l2m, l3m]]; 
         CF22 = H[kf22[l1m, l2m, l3m]]; CM11 = Cm11[l1m, l2m, l3m]; 
         CM12 = Cm12[l1m, l2m, l3m]; CM13 = Cm13[l1m, l2m, l3m]; 
         CM22 = Cm22[l1m, l2m, l3m]; CM23 = Cm23[l1m, l2m, l3m]; 
         CM33 = Cm33[l1m, l2m, l3m]; C = {{CF11 + CM11, CF12 + CM12, CM13}, 
           {CF12 + CM12, CF22 + CM22, CM23}, {CM13, CM23, CM33}}; C]; 
      E1[x_] = ExpIntegralE[1, x]; 
      G[t_] = Goo*(1 + c*(E1[t/tau2] - E1[t/tau1])); 
      G0 = Goo*(1 + c*Limit[E1[t/tau2] - E1[t/tau1], t -> 0]); 
      Goo = Limit[G[t], t -> Infinity]; dG[t_] = D[G[t], t]; 
      dG0 = Limit[dG[t], t -> 0]; tolG = 10^(-8); 
      tmax = t /. FindRoot[dG[t] == -tolG, {t, 2}]; 
      strainHistory[dt_, sh_] := Module[{maxClockTime = 60*5, ntsteps = 0, 
         ts = {}, tolR = 10^(-8), sigma, sigma0, sigint, errorR, Gscale, 
         Gs = {}, dGs = {}, l1m = 1, l2m = 1, l3m = 1, outInt = 100}, 
        ntsteps = First[Dimensions[sh]]; ts = Range[0, dt*(ntsteps - 1), dt]; 
         sigint = Table[{0, 0, 0}, {k, ntsteps}]; stretches = sigint; 
         If[Last[ts] > tmax, If[errchk, Print["can truncate dG"]]]; 
         dGs = dG /@ Drop[ts, 1]; Gs = G /@ Drop[ts, 1]; 
         Gscale = G0 + dt*dG0; For[m = 1, m <= ntsteps, m++, 
          time0 = TimeUsed[]; sigma = stress[l1m, l2m, l3m]; 
           sigma0 = dt*Take[dGs, m - 1] . Reverse[Take[sigint, m - 1]]; 
           sigma = Gscale*sigma + sigma0; Rext = sh[[m]]/{l1m, l2m, l3m}; 
           Res0 = Rext - sigma; errorR = Sqrt[Res0 . Res0]; errors = {}; 
           While[errorR > tolR, K = stiffness[l1m, l2m, l3m] . DiagonalMatrix[
                {l1m, l2m, l3m}]; K *= Gscale; K -= DiagonalMatrix[-sh[[m]]/
                {l1m*l1m, l2m*l2m, l3m*l3m}]; {l1m, l2m, l3m} += 
              Inverse[K] . Res0; If[l1m <= 0 || l2m <= 0 || l3m <= 0, 
              If[errchk, Print["negative stretch : ", {l1m, l2m, l3m}]]; 
               Abort[]]; sigma = sigint[[m]] = stress[l1m, l2m, l3m]; 
             sigma = Gscale*sigma + sigma0; Rext = sh[[m]]/{l1m, l2m, l3m}; 
             Res0 = Rext - sigma; errorR = Sqrt[Res0 . Res0]; 
             AppendTo[errors, errorR]; ]; If[errchk && Mod[m, outInt] == 0, 
            Print[m, " / ", ntsteps]; ]; sigint[[m]] = sigma; 
           stretches[[m]] = {l1m, l2m, l3m}; ]; stretches]; error = 0.; 
      datasets = {"is01b-creep"}; For[i = 1, i <= Length[datasets], i++, 
       If[errchk, Print[datasets[[i]]]]; stream = OpenRead[datasets[[i]]]; 
        header = Read[stream, String]; data = ReadList[stream, 
          {Number, Number, Number, Number, Number, Number}]; Close[stream]; 
        data = Drop[data, 4]; time0 = First[Transpose[data][[2]]]; 
        expTime = Transpose[data][[2]] - time0; 
        edts = Drop[RotateLeft[expTime] - expTime, -1]; 
        expStress = Transpose[data][[5]]; expStress *= StressScale; 
        expStrain = 1. + Transpose[data][[6]]; 
        Interpolate[th_, thB_, yhB_] := Module[{t1, t2, y1, y2, j = 1, t, y, 
           yh = {}}, t1 = thB[[j]]; t2 = thB[[j + 1]]; y1 = yhB[[j]]; 
           y2 = yhB[[j + 1]]; For[i = 1, i <= Length[th], i++, 
            t = th[[i]]; While[t > t2, j++; t1 = thB[[j]]; t2 = thB[[j + 1]]; 
               y1 = yhB[[j]]; y2 = yhB[[j + 1]]]; AppendTo[yh, 
              y = y1 + (y2 - y1)*((t - t1)/(t2 - t1))]; ]; yh]; dt = 1.; 
        th = Range[First[expTime], Last[expTime], dt]; 
        eshA = Interpolate[th, expTime, expStress]; 
        Th = Transpose[{0*th, eshA, 0*th}]; lh = strainHistory[dt, Th]; 
        data = Transpose[{th, Transpose[lh][[1]], Transpose[lh][[2]], 
           Transpose[lh][[2]]}]; Export["strain.dat", data]; 
        intStrain = Interpolate[th, expTime, expStrain]; 
        l2h = Transpose[lh][[2]]; diff = l2h - intStrain; 
        error += dt*Sqrt[diff . diff]; If[errchk, 
         RGB = {RGBColor[1, 0, 0], RGBColor[0, 1, 0], RGBColor[0, 0, 1]}; 
          Show[pls = Table[ListPlot[Transpose[{th, Transpose[lh][[i]]}], 
              PlotJoined -> True, PlotLabel -> "strain history", 
              PlotRange -> All, PlotStyle -> RGB[[i]]], {i, 3}], 
           ListPlot[Transpose[{expTime, expStrain}]]]; ]; ]; 
      If[errchk, Export["results.dat", {error}]]; Return[error], 60*5, 
     Return[1.]]]; 
