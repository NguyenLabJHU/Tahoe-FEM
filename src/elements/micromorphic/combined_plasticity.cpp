                              PlasticityCondition=4;
                                //initialize before iteration
                                fYield_function=fYield_function_tr;
                                fMicroYield_function=fMicroYield_function_tr;
                                fFe=fFe_tr;
                                fChie=fChie_tr;
                                fCchie=fCchie_tr;
                                //initial values for Fp is assumed the same with previous step
                                fFp=fFp_n;
                                //fTemp_matrix_nsd_x_nsd.Inverse(fFe);
                                //fFp.MultAB(fTemp_matrix_nsd_x_nsd,fDeformation_Gradient);

                                SPK=fSPK_tr;
                                devSPK=fdevSPK_tr;
                                devfSPKinv=devfSPKinv_tr;

                                SIGMA_S=SIGMA_S_tr;
                                devSIGMA_S_inv=devSIGMA_S_inv_tr;
                                devSIGMA_S=devSIGMA_S_tr;
                                PSIe=PSIe_tr;
                                Stress_Norm=Stress_Norm_tr;


                                fdelDelgamma = 0.0;
                                fDelgamma = 0.0;

                                //iterate using Newton-Raphson to solve for fDelgamma
                                iter_count = 0;
                                fs_micromorph3D_out<< "Gauss Point = "<< IP <<endl;
                                fs_micromorph3D_out << "Current  Macro Yield function = " << fYield_function << endl;
                                while (fabs(fYield_function) > dAbsTol && fabs(fYield_function/fYield_function_tr) > dRelTol && iter_count < iIterationMax)
                                {
                                	iter_count += 1;
                                	//Form  dFe/dDgamma
                                	fFp_inverse.Inverse(fFp);
                                	dFedDelgamma=0.0;
                                	fdGdS_n_transpose.Transpose(fdGdS_n);
                                	fCe_n_inverse.Inverse(fCe_n);
                                	fTemp_matrix_nsd_x_nsd.MultATBC(fdGdS_n,fFp_n,fFp_inverse);
                                	dFedDelgamma.MultABC(fFe,fCe_n_inverse,fTemp_matrix_nsd_x_nsd);
                                	dFedDelgamma*=-1;

                                	//Forming  dE^e/dDgamma  E^e: Elas. Lag. stn tensor
                                	dEedDelgamma.MultATB(dFedDelgamma,fFe);
                                	fTemp_matrix_nsd_x_nsd.MultATB(fFe,dFedDelgamma);
                                	dEedDelgamma+=fTemp_matrix_nsd_x_nsd;
                                	dEedDelgamma*=0.5;

                                	//Inverse of fChip
                                    fChip_inverse.Inverse(fChip);


                                    /* Form dChip/dDelgamma */
                                    fTemp_matrix_nsd_x_nsd.MultATBC(PSIe_n_inverse,fCchie_n,fChip_n);
                                    dChipdDelgamma.MultABC(PSIe_n_inverse,fdGchidSIGMA_S_n_transpose,fTemp_matrix_nsd_x_nsd);

                                    /* Forming dChie/dDgamma*/
                                    dChiedDelgamma.MultABC(fChie,dChipdDelgamma,fChip_inverse);
                                    dChiedDelgamma*=-1;


                                    /* Forming  dEpsilon^e/dDgamma  Epsilone^e: Elastic micro strain tensor */
                                    dEpsilonedDelgamma.MultATB(dFedDelgamma,fChie);
                                    fTemp_matrix_nsd_x_nsd.MultATB(fFe,dChiedDelgamma);
                                    dEpsilonedDelgamma+=fTemp_matrix_nsd_x_nsd;


                                	//Forming  dS/dDgamma  S= SPK tensor
                                	Temp_inv=dEedDelgamma.Trace();
                                	dSdDelgamma.SetToScaled((fMaterial_Params[kLambda]+fMaterial_Params[kTau])*Temp_inv,fIdentity_matrix);

                                	fTemp_matrix_nsd_x_nsd.SetToScaled(2*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]),dEedDelgamma);
                                	dSdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                	Temp_inv=dEpsilonedDelgamma.Trace();
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kEta]*Temp_inv,fIdentity_matrix);
                                	dSdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                	fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kKappa],dEpsilonedDelgamma);
                                	dSdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                	fTemp_matrix_nsd_x_nsd2.Transpose(dEpsilonedDelgamma);
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kNu],fTemp_matrix_nsd_x_nsd2);
                                	dSdDelgamma+=fTemp_matrix_nsd_x_nsd;


                                	//Forming  dP/dDgamma (scalar) P: pressure  dP/dDgamma= (1/3)1:dS/dDgamma
                                	//One error was here; dPdDelgamma*=1/3 makes dPdDelgamma=0
                                	dPdDelgamma=dSdDelgamma.Trace()/3;

                                	//Forming  dc/dDgamma  c: cohesion
                                	//fState_variables_n_IPs(IP,khc) =Aphi;
                                	dcdDelgamma=fState_variables_n_IPs(IP,khc)*fMaterial_Params[kHc];

                                	//Forming  d(devS)/dDgamma  devS: dev. part of SPK tensor
                                	ddevSdDelgamma.SetToScaled(dPdDelgamma,fIdentity_matrix);
                                	ddevSdDelgamma*=-1;
                                	ddevSdDelgamma+=dSdDelgamma;

                                    /* Forming  d(SIGMA-S)/dDgamma tensor*/
                                    Temp_inv=dEedDelgamma.Trace();
                                    dSIGMA_SdDelgamma.SetToScaled(fMaterial_Params[kTau]*Temp_inv,fIdentity_matrix);

                                    fTemp_matrix_nsd_x_nsd.SetToScaled(2*fMaterial_Params[kSigma_const],dEedDelgamma);
                                    dSIGMA_SdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                    Temp_inv=dEpsilonedDelgamma.Trace();
                                    fTemp_matrix_nsd_x_nsd.SetToScaled((fMaterial_Params[kEta]-fMaterial_Params[kTau])*Temp_inv,fIdentity_matrix);
                                    dSIGMA_SdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                    fTemp_matrix_nsd_x_nsd.SetToScaled((fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]),dEpsilonedDelgamma);
                                    dSIGMA_SdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                    fTemp_matrix_nsd_x_nsd2.Transpose(dEpsilonedDelgamma);
                                    fTemp_matrix_nsd_x_nsd.SetToScaled((fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]),fTemp_matrix_nsd_x_nsd2);
                                    dSIGMA_SdDelgamma+=fTemp_matrix_nsd_x_nsd;

                                    /*Forming  dPchi/dDgammachi (scalar) Pchi: pressure for micro-scale  dPchi/dDgammachi= (1/3)1:dSIGMA_S/
                                    dDgammachi*/
                                    dPchidDelgamma=dSIGMA_SdDelgamma.Trace()/3;

                                    /* Forming  dc/dDgamma  c: cohesion*/
                                    dcdDelgamma=fState_variables_n_IPs(IP,khc)*fMaterial_Params[kHc];

                                    /* Forming  d(dev(SIGMA_S))/dDgammachi  dev(SIGMA_S): dev. part of SIGMA-S (relative stress) tensor*/
                                    ddevSIGMA_SdDelgamma.SetToScaled(dPchidDelgamma,fIdentity_matrix);
                                    ddevSIGMA_SdDelgamma*=-1;
                                    ddevSIGMA_SdDelgamma+=dSIGMA_SdDelgamma;



                                    //CHECK STRESS NORM
                                	//Forming  d(||devS||)/dDgamma  devS: dev. part of SPK tensor
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(1/Stress_Norm,devSPK);
                                	InvddevSdDelgamma=dMatrixT::Dot(ddevSdDelgamma,fTemp_matrix_nsd_x_nsd);
                                    fTemp_matrix_nsd_x_nsd.SetToScaled(1/Stress_Norm,devSIGMA_S);
                                    Temp_inv=dMatrixT::Dot(ddevSIGMA_SdDelgamma,fTemp_matrix_nsd_x_nsd);
                                    InvddevSdDelgamma+=Temp_inv;

									// assemble the consistent tangent
                                	//dFYdDelgamma=InvddevSdDelgamma-(Aphi*dcdDelgamma-Bphi*dPdDelgamma);
                                	dFYdDelgamma=InvddevSdDelgamma-(Aphi*dcdDelgamma-Bphi*dPdDelgamma+Aphi_chi*dcchidDelgamma-Bphi_chi*dPchidDelgamma);


                                	//solve for fdelDelgamma
                                	if (dFYdDelgamma != 0.0) fdelDelgamma = -fYield_function/dFYdDelgamma;
                                	else fdelDelgamma = 0.0;

                                	//update fDelgamma
                                	fDelgamma += fdelDelgamma;


                                	if (fDelgamma < 0.0) fDelgamma = 0.0;
                                	fState_variables_IPs(IP,kDelgamma) = fDelgamma;


                                	//update c ISVs

                                	fState_variables_IPs(IP,kc)= fState_variables_n_IPs(IP,kc)
                                	+ fDelgamma*fState_variables_n_IPs(IP,khc)*fMaterial_Params[kHc];
                                	if (fState_variables_IPs(IP,kc) < 0.0) fState_variables_IPs(IP,kc) = 0.0;


                                    // update cx (c_chi) ISVs
                                    fState_variables_IPs(IP,kc_chi)= fState_variables_n_IPs(IP,kc_chi)
                                    + fDelgamma*fState_variables_n_IPs(IP,khc_chi)*fMaterial_Params[kHc_chi];
                                    if (fState_variables_IPs(IP,kc_chi) < 0.0) fState_variables_IPs(IP,kc_chi) = 0.0;


                                	//  	                              update fFp
                                	fdGdS_n_transpose.Transpose(fdGdS_n);
                                	fCe_n_inverse.Inverse(fCe_n);
                                	//fTemp_matrix_nsd_x_nsd.MultAB(fCe_n_inverse,fdGdS_n_transpose);
                                	fTemp_matrix_nsd_x_nsd.MultABT(fCe_n_inverse,fdGdS_n);
                                	fTemp_matrix_nsd_x_nsd*=fDelgamma;

                                	fTemp_matrix_nsd_x_nsd += fIdentity_matrix;
                                	fFp.MultAB(fTemp_matrix_nsd_x_nsd,fFp_n);

                                	//calculate fFp_Inverse
                                	fFp_inverse.Inverse(fFp);

                                	//calculate Fe
                                	fFe.MultAB(fDeformation_Gradient,fFp_inverse);

                                	//[fElastic_Right_Cauchy_Green_tensor] will be formed
                                	fRight_Elastic_Cauchy_Green_tensor.MultATB(fFe,fFe);
                                	if (fRight_Elastic_Cauchy_Green_tensor.Det()==0)
                                		fRight_Elastic_Cauchy_Green_tensor = fIdentity_matrix;

                                	//[fMicroElastic_Right_Cauchy_Green_tensor] will be formed
                                	fMicroRight_Elastic_Cauchy_Green_tensor.MultATB(fChie,fChie);
                                	if (fMicroRight_Elastic_Cauchy_Green_tensor.Det()==0)
                                		fMicroRight_Elastic_Cauchy_Green_tensor = fIdentity_matrix;


                                    /* update fChip */
                                    fTemp_matrix_nsd_x_nsd.MultATBC(PSIe_n_inverse,fCchie_n,fChip_n);
                                    fTemp_matrix_nsd_x_nsd2.MultABC(PSIe_n_inverse,fdGchidSIGMA_S_n_transpose,fTemp_matrix_nsd_x_nsd);
                                    fChip.SetToScaled(fDelgamma,fTemp_matrix_nsd_x_nsd2);
                                    fChip+=fChip_n;


                                    /* Form inverse of Chi^p*/
                                    fChip_inverse.Inverse(fChip);

                                    /* Calculate Chie */
                                    fChie.MultAB(ChiM,fChip_inverse);


                                    /* [fMicroElastic_Right_Cauchy_Green_tensor] will be formed */
                                    fMicroRight_Elastic_Cauchy_Green_tensor.MultATB(fChie,fChie);
                                    if (fMicroRight_Elastic_Cauchy_Green_tensor.Det()==0)
                                            fMicroRight_Elastic_Cauchy_Green_tensor = fIdentity_matrix;


                                    //Update fCchie
                                    fCchie.MultATB(fChie,fChie);

                                    /* Update PSIe */
                                    PSIe.MultATB(fFe,fChie);



                                	//Update Elastic Lagrangian strain tensor E
                                	Elastic_LagrangianStn=fIdentity_matrix;
                                	Elastic_LagrangianStn*=-1;
                                	Elastic_LagrangianStn+=fRight_Elastic_Cauchy_Green_tensor;
                                	Elastic_LagrangianStn*=0.5;

                                	//Update Elastic micro strain tensor will be formed in Bbar
                                	Elastic_MicroStnTensor = fIdentity_matrix;
                                	Elastic_MicroStnTensor *= -1;
                                	Elastic_MicroStnTensor += PSIe;

                                	//update S stress
                                	Temp_inv=0.0;
                                	Temp_inv=Elastic_LagrangianStn.Trace();//Calculating the tr(E) and keep in temp. var.
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_inv*(fMaterial_Params[kLambda]+fMaterial_Params[kTau]),fIdentity_matrix);

                                	SPK.SetToScaled(2*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]),Elastic_LagrangianStn);
                                	SPK+=fTemp_matrix_nsd_x_nsd;

                                	Temp_inv=0.0;
                                	Temp_inv=Elastic_MicroStnTensor.Trace();
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_inv*fMaterial_Params[kEta],fIdentity_matrix);
                                	SPK+=fTemp_matrix_nsd_x_nsd;

                                	fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kKappa],Elastic_MicroStnTensor);
                                	SPK+=fTemp_matrix_nsd_x_nsd;

                                	fTemp_matrix_nsd_x_nsd2.Transpose(Elastic_MicroStnTensor);
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kNu],fTemp_matrix_nsd_x_nsd2);
                                	SPK+=fTemp_matrix_nsd_x_nsd;

                                	//Update Relative stress SIGMA_S
                                	Temp_inv=Elastic_LagrangianStn.Trace();
                                	SIGMA_S.SetToScaled(Temp_inv*fMaterial_Params[kTau],fIdentity_matrix);
                                	// 2sigmaE
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(2*fMaterial_Params[kSigma_const],Elastic_LagrangianStn);
                                	SIGMA_S+=fTemp_matrix_nsd_x_nsd;
                                	//(eta-Tau)trEpsilon.1
                                	Temp_inv=Elastic_MicroStnTensor.Trace();
                                	fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_inv*(fMaterial_Params[kEta]-fMaterial_Params[kTau]),fIdentity_matrix);
                                	SIGMA_S+=fTemp_matrix_nsd_x_nsd;
                                	//(nu-sigma)*Epsilon
                                	fTemp_matrix_nsd_x_nsd.SetToScaled((fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]),Elastic_MicroStnTensor);
                                	SIGMA_S+=fTemp_matrix_nsd_x_nsd;
                                	//(kappa-sigma)*Epsilon^T
                                	fTemp_matrix_nsd_x_nsd2.Transpose(Elastic_MicroStnTensor);
                                	fTemp_matrix_nsd_x_nsd.SetToScaled((fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]),fTemp_matrix_nsd_x_nsd2);
                                	SIGMA_S+=fTemp_matrix_nsd_x_nsd;

                                	cout<<"Elastic_MicroStnTensor.det="<<Elastic_MicroStnTensor.Det()<<endl;

                                	//calculate  devS stress
                                	Pbar=SPK.Trace()/3;//Calculating the pressure term
                                	devSPK.SetToScaled(Pbar,fIdentity_matrix);
                                	devSPK*=-1;
                                	devSPK+=SPK;

                                    /* Form the deviatoric SIGMA-S */
                                    Pchibar=SIGMA_S.Trace()/3;//Calculating the pressure term
                                    devSIGMA_S.SetToScaled(Pchibar,fIdentity_matrix);
                                    devSIGMA_S*=-1;
                                    devSIGMA_S+=SIGMA_S;




                                    /* Calculate ||devS:devS+devR:devR||  */
                                    Temp_inv= devSPK.ScalarProduct();
                                    Stress_Norm=Temp_inv;
                                    Temp_inv= devSIGMA_S.ScalarProduct();
                                    Stress_Norm+=Temp_inv;
                                    Temp_inv=Stress_Norm;
                                    Stress_Norm=sqrt(Temp_inv);



                                	// Calculate yield function with updated parameters
                                	//fYield_function=Stress_Norm-(Aphi*(fState_variables_IPs(IP,kc))-Bphi*Pbar);
                                	fYield_function=Stress_Norm-(Aphi*fState_variables_IPs(IP,kc)-Bphi*Pbar+Aphi_chi*fState_variables_IPs(IP,kc_chi)-Bphi_chi*Pchibar);

                                	fs_micromorph3D_out  << "Current relative residual = " << fabs(fYield_function/fYield_function_tr) << endl;
