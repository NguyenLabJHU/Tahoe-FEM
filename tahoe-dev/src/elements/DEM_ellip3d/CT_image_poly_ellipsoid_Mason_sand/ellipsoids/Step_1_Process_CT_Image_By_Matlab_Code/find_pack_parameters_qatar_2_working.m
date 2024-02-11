	
	%===============================================================
    % Developed by Dr. Riyadh Al-Raoush, 12/31/2011
    % Southern University and A&M College
    % Departement of Civil and Environmental Engineering
    % Baton Rouge, LA 70813
    % email: riyadh_alraoush@subr.edu
    %
    %  find_pack_parameters_fast
    %===============================================================
	clc
    clear all;
    load objects_image_MasonSand.mat;
    %objects_image = objects_image(1:60,1:60,1:60);
    
    connect=8;   %%%  Connectivity
	
    nrows=size(objects_image,1);
	ncolumns=size(objects_image,2);
	ndepth=size(objects_image,3);
    n_objects=max(objects_image(:))
    
         
        
    objects_surface_area=[];
    objects_centroids=[];
    objects_sizes=[];
    objects_coord_numbers=[];
    all_coord_numbers=[];
    objects_angles=[];
    objects_lengths=[];
    objects_diameter_sphere=[];
    objects_sphereness=[];
    objects_roundness=[];
    objects_normal_vector=[];
    objects_tangent_vector=[];
    objects_center_center_vector=[];
    collect_i_zone=[];
    collect_object_voxels_final_index=[];
    pointer_collect_object_voxels_final_index=[];
    collect_coord_voxels=[];
    pointer_collect_coord_voxels=[];
    collect_coord_number=[];
    pointer_collect_coord_number=[];
    contact_points_coordinats=[];
    objects_contacts_center=[];
    
    object_poly_lengths = [];  % the x+, x-, y+, y-, z+ and z- of the equivalent polyellipsoids
    objects_new_center = [];    % the new center of the poly-ellipsoid when we use the shifted center
    objects_principal_vectors = []; % the coordinates of the principal vectors, 
                                    % unlike objects_angles, this
                                    % principla_vectors will give us 9
                                    % components of principals
    
     i_zone=1;
     
    % use: if we uncomment the shif-center part, then it will get
    % shift-center poly-ellipsoids & ellipsoids, comment, just same center
    % poly-ellipsoids and elliopsoids.
    % ellipsoids: objects_centroids, object_lengths, objects_principal_vectors
    % (the same both for poly and ellipsoids)
    % shift-center poly: objects_new_center, object_poy_lengths,
    % objects_principal_vectors
    % same center poly: objects_centroids, object_pol_lengths,
    % objects_principal_vectors
        
	while i_zone <= n_objects
	
	i_zone
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Centers of particels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    object_voxels_final_index=find(objects_image(:)==i_zone);
    
   
    d_z_object_final=floor((object_voxels_final_index-1)/(nrows*ncolumns))+1;
    d_y_object_final=ceil(object_voxels_final_index/nrows)-(ncolumns*(d_z_object_final-1));
    d_x_object_final=object_voxels_final_index-nrows*(d_y_object_final-1)-nrows*ncolumns*(d_z_object_final-1);
    
   
    object_voxels_final=[d_x_object_final,d_y_object_final,d_z_object_final];
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%      start get connect       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
         
       
    n_1=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_2=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_3=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_4=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_5=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_6=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_7=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_8=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];   
    n_9=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)-1];
    n_10=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_11=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_12=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_13=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_14=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_15=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_16=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_17=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)];
    n_18=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_19=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_20=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1]; 
    n_21=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_22=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_23=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1),object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_24=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)-1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_25=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2),object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
    n_26=[object_voxels_final(1:length(object_voxels_final(:,1,:)),1)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),2)+1,object_voxels_final(1:length(object_voxels_final(:,1,:)),3)+1];
              
    
       
    n_1_fake=n_1;
    n_2_fake=n_2;
    n_3_fake=n_3;
    n_4_fake=n_4;
    n_5_fake=n_5;
    n_6_fake=n_6;
    n_7_fake=n_7;
    n_8_fake=n_8;
    n_9_fake=n_9;
    n_10_fake=n_10;
    n_11_fake=n_11;
    n_12_fake=n_12;
    n_13_fake=n_13;
    n_14_fake=n_14;
    n_15_fake=n_15;
    n_16_fake=n_16;
    n_17_fake=n_17;
    n_18_fake=n_18;
    n_19_fake=n_19;
    n_20_fake=n_20;
    n_21_fake=n_21;
    n_22_fake=n_22;
    n_23_fake=n_23;
    n_24_fake=n_24; 
    n_25_fake=n_25;
    n_26_fake=n_26;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_1_fake(find(n_1(:,1)<=0),1)=n_1_fake(find(n_1(:,1)<=0),1)+1;         
    n_1_fake(find(n_1(:,2)<=0),2)=n_1_fake(find(n_1(:,2)<=0),2)+1;
	n_1_fake(find(n_1(:,3)<=0),3)=n_1_fake(find(n_1(:,3)<=0),3)+1;
  	n_2_fake(find(n_2(:,1)<=0),1)=n_2_fake(find(n_2(:,1)<=0),1)+1;
 	n_2_fake(find(n_2(:,2)<=0),2)=n_2_fake(find(n_2(:,2)<=0),2)+1;
 	n_2_fake(find(n_2(:,3)<=0),3)=n_2_fake(find(n_2(:,3)<=0),3)+1;
    n_3_fake(find(n_3(:,3)<=0),3)=n_3_fake(find(n_3(:,3)<=0),3)+1;
    n_3_fake(find(n_3(:,2)>ncolumns),2)=n_3_fake(find(n_3(:,2)>ncolumns),2)-1;
   	n_3_fake(find(n_3(:,1)<=0),1)=n_3_fake(find(n_3(:,1)<=0),1)+1;
  	n_4_fake(find(n_4(:,3)<=0),3)=n_4_fake(find(n_4(:,3)<=0),3)+1;
 	n_4_fake(find(n_4(:,2)<=0),2)=n_4_fake(find(n_4(:,2)<=0),2)+1;
    n_5_fake(find(n_5(:,3)<=0),3)=n_5_fake(find(n_5(:,3)<=0),3)+1;
    n_6_fake(find(n_6(:,3)<=0),3)=n_6_fake(find(n_6(:,3)<=0),3)+1;
   	n_6_fake(find(n_6(:,2)>ncolumns),2)=n_6_fake(find(n_6(:,2)>ncolumns),2)-1;
    n_7_fake(find(n_7(:,3)<=0),3)=n_7_fake(find(n_7(:,3)<=0),3)+1;
   	n_7_fake(find(n_7(:,2)<=0),2)=n_7_fake(find(n_7(:,2)<=0),2)+1;
 	n_7_fake(find(n_7(:,1)>nrows),1)=n_7_fake(find(n_7(:,1)>nrows),1)-1;
    n_8_fake(find(n_8(:,3)<=0),3)=n_8_fake(find(n_8(:,3)<=0),3)+1;
   	n_8_fake(find(n_8(:,1)>nrows),1)=n_8_fake(find(n_8(:,1)>nrows),1)-1;
    n_9_fake(find(n_9(:,3)<=0),3)=n_9_fake(find(n_9(:,3)<=0),3)+1;
  	n_9_fake(find(n_9(:,2)>ncolumns),2)=n_9_fake(find(n_9(:,2)>ncolumns),2)-1;
   	n_9_fake(find(n_9(:,1)>nrows),1)=n_9_fake(find(n_9(:,1)>nrows),1)-1;
   	n_10_fake(find(n_10(:,2)<=0),2)=n_10_fake(find(n_10(:,2)<=0),2)+1;
	n_10_fake(find(n_10(:,1)<=0),1)=n_10_fake(find(n_10(:,1)<=0),1)+1;     
    n_11_fake(find(n_11(:,1)<=0),1)=n_11_fake(find(n_11(:,1)<=0),1)+1;
   	n_12_fake(find(n_12(:,2)>ncolumns),2)=n_12_fake(find(n_12(:,2)>ncolumns),2)-1;
  	n_12_fake(find(n_12(:,1)<=0),1)=n_12_fake(find(n_12(:,1)<=0),1)+1;
    n_13_fake(find(n_13(:,2)<=0),2)=n_13_fake(find(n_13(:,2)<=0),2)+1;
   	n_14_fake(find(n_14(:,2)>ncolumns),2)=n_14_fake(find(n_14(:,2)>ncolumns),2)-1;
    n_15_fake(find(n_15(:,2)<=0),2)=n_15_fake(find(n_15(:,2)<=0),2)+1;
 	n_15_fake(find(n_15(:,1)>nrows),1)=n_15_fake(find(n_15(:,1)>nrows),1)-1;
    n_16_fake(find(n_16(:,1)>nrows),1)=n_16_fake(find(n_16(:,1)>nrows),1)-1;
    n_17_fake(find(n_17(:,2)>ncolumns),2)=n_17_fake(find(n_17(:,2)>ncolumns),2)-1;
   	n_17_fake(find(n_17(:,1)>nrows),1)=n_17_fake(find(n_17(:,1)>nrows),1)-1;
   	n_18_fake(find(n_18(:,3)>ndepth),3)=n_18_fake(find(n_18(:,3)>ndepth),3)-1;
   	n_18_fake(find(n_18(:,2)<=0),2)=n_18_fake(find(n_18(:,2)<=0),2)+1;
   	n_18_fake(find(n_18(:,1)<=0),1)=n_18_fake(find(n_18(:,1)<=0),1)+1;
    n_19_fake(find(n_19(:,3)>ndepth),3)=n_19_fake(find(n_19(:,3)>ndepth),3)-1;
   	n_19_fake(find(n_19(:,1)<=0),1)=n_19_fake(find(n_19(:,1)<=0),1)+1;
   	n_20_fake(find(n_20(:,3)>ndepth),3)=n_20_fake(find(n_20(:,3)>ndepth),3)-1;
  	n_20_fake(find(n_20(:,2)>ncolumns),2)=n_20_fake(find(n_20(:,2)>ncolumns),2)-1;
  	n_20_fake(find(n_20(:,1)<=0),1)=n_20_fake(find(n_20(:,1)<=0),1)+1;
    n_21_fake(find(n_21(:,3)>ndepth),3)=n_21_fake(find(n_21(:,3)>ndepth),3)-1;
	n_21_fake(find(n_21(:,2)<=0),2)=n_21_fake(find(n_21(:,2)<=0),2)+1;
   	n_22_fake(find(n_22(:,3)>ndepth),3)=n_22_fake(find(n_22(:,3)>ndepth),3)-1;
    n_23_fake(find(n_23(:,3)>ndepth),3)=n_23_fake(find(n_23(:,3)>ndepth),3)-1;
   	n_23_fake(find(n_23(:,2)>ncolumns),2)=n_23_fake(find(n_23(:,2)>ncolumns),2)-1;
    n_24_fake(find(n_24(:,3)>ndepth),3)=n_24_fake(find(n_24(:,3)>ndepth),3)-1;
 	n_24_fake(find(n_24(:,2)<=0),2)=n_24_fake(find(n_24(:,2)<=0),2)+1;
   	n_24_fake(find(n_24(:,1)>nrows),1)=n_24_fake(find(n_24(:,1)>nrows),1)-1;
    n_25_fake(find(n_25(:,3)>ndepth),3)=n_25_fake(find(n_25(:,3)>ndepth),3)-1;
  	n_25_fake(find(n_25(:,1)>nrows),1)=n_25_fake(find(n_25(:,1)>nrows),1)-1;
    n_26_fake(find(n_26(:,3)>ndepth),3)=n_26_fake(find(n_26(:,3)>ndepth),3)-1;
   	n_26_fake(find(n_26(:,2)>ncolumns),2)=n_26_fake(find(n_26(:,2)>ncolumns),2)-1;
	n_26_fake(find(n_26(:,1)>nrows),1)=n_26_fake(find(n_26(:,1)>nrows),1)-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    
    %index=x+r(y-1)+r*c(z-1)
    
   
        n_1=n_1_fake(:,1)+nrows*(n_1_fake(:,2)-1)+nrows*ncolumns*(n_1_fake(:,3)-1);
        n_2=n_2_fake(:,1)+nrows*(n_2_fake(:,2)-1)+nrows*ncolumns*(n_2_fake(:,3)-1);
        n_3=n_3_fake(:,1)+nrows*(n_3_fake(:,2)-1)+nrows*ncolumns*(n_3_fake(:,3)-1);
        n_4=n_4_fake(:,1)+nrows*(n_4_fake(:,2)-1)+nrows*ncolumns*(n_4_fake(:,3)-1);
        n_5=n_5_fake(:,1)+nrows*(n_5_fake(:,2)-1)+nrows*ncolumns*(n_5_fake(:,3)-1);
        n_6=n_6_fake(:,1)+nrows*(n_6_fake(:,2)-1)+nrows*ncolumns*(n_6_fake(:,3)-1);
        n_7=n_7_fake(:,1)+nrows*(n_7_fake(:,2)-1)+nrows*ncolumns*(n_7_fake(:,3)-1);
        n_8=n_8_fake(:,1)+nrows*(n_8_fake(:,2)-1)+nrows*ncolumns*(n_8_fake(:,3)-1);
        n_9=n_9_fake(:,1)+nrows*(n_9_fake(:,2)-1)+nrows*ncolumns*(n_9_fake(:,3)-1);
        
        n_10=n_10_fake(:,1)+nrows*(n_10_fake(:,2)-1)+nrows*ncolumns*(n_10_fake(:,3)-1);
        n_11=n_11_fake(:,1)+nrows*(n_11_fake(:,2)-1)+nrows*ncolumns*(n_11_fake(:,3)-1);
        n_12=n_12_fake(:,1)+nrows*(n_12_fake(:,2)-1)+nrows*ncolumns*(n_12_fake(:,3)-1);
        n_13=n_13_fake(:,1)+nrows*(n_13_fake(:,2)-1)+nrows*ncolumns*(n_13_fake(:,3)-1);
        n_14=n_14_fake(:,1)+nrows*(n_14_fake(:,2)-1)+nrows*ncolumns*(n_14_fake(:,3)-1);
        n_15=n_15_fake(:,1)+nrows*(n_15_fake(:,2)-1)+nrows*ncolumns*(n_15_fake(:,3)-1);
        n_16=n_16_fake(:,1)+nrows*(n_16_fake(:,2)-1)+nrows*ncolumns*(n_16_fake(:,3)-1);
        n_17=n_17_fake(:,1)+nrows*(n_17_fake(:,2)-1)+nrows*ncolumns*(n_17_fake(:,3)-1);
        
        n_18=n_18_fake(:,1)+nrows*(n_18_fake(:,2)-1)+nrows*ncolumns*(n_18_fake(:,3)-1);
        n_19=n_19_fake(:,1)+nrows*(n_19_fake(:,2)-1)+nrows*ncolumns*(n_19_fake(:,3)-1);
        n_20=n_20_fake(:,1)+nrows*(n_20_fake(:,2)-1)+nrows*ncolumns*(n_20_fake(:,3)-1);
        n_21=n_21_fake(:,1)+nrows*(n_21_fake(:,2)-1)+nrows*ncolumns*(n_21_fake(:,3)-1);
        n_22=n_22_fake(:,1)+nrows*(n_22_fake(:,2)-1)+nrows*ncolumns*(n_22_fake(:,3)-1);
        n_23=n_23_fake(:,1)+nrows*(n_23_fake(:,2)-1)+nrows*ncolumns*(n_23_fake(:,3)-1);
        n_24=n_24_fake(:,1)+nrows*(n_24_fake(:,2)-1)+nrows*ncolumns*(n_24_fake(:,3)-1);
        n_25=n_25_fake(:,1)+nrows*(n_25_fake(:,2)-1)+nrows*ncolumns*(n_25_fake(:,3)-1);
        n_26=n_26_fake(:,1)+nrows*(n_26_fake(:,2)-1)+nrows*ncolumns*(n_26_fake(:,3)-1);
      
    if connect == 4    %%% 2D
        
        n_1=[];
        n_2=[];
        n_3=[];
        n_4=[];
        n_5=[];
        n_6=[];
        n_7=[];
        n_8=[];
        n_9=[];
        n_10=[];
        n_12=[];
        n_15=[];
        n_17=[];
        n_18=[];
        n_19=[];
        n_20=[];
        n_21=[];
        n_22=[];
        n_23=[];
        n_24=[];
        n_25=[];
        n_26=[];
        
    end
    
        
    if connect == 8   %%%% 2D
        
        n_1=[];
        n_2=[];
        n_3=[];
        n_4=[];
        n_5=[];
        n_6=[];
        n_7=[];
        n_8=[];
        n_9=[];
        n_18=[];
        n_19=[];
        n_20=[];
        n_21=[];
        n_22=[];
        n_23=[];
        n_24=[];
        n_25=[];
        n_26=[];
        
    end
    
        
      if connect == 6   %%%3D
        
        n_1=[];
        n_2=[];
        n_3=[];
        n_4=[];
        n_6=[];
        n_7=[];
        n_8=[];
        n_9=[];
        n_10=[];
        n_12=[];
        n_15=[];
        n_17=[];
        n_18=[];
        n_19=[];
        n_20=[];
        n_21=[];
        n_23=[];
        n_24=[];
        n_25=[];
        n_26=[];
        
    end
        
           
    connect_voxels_index=[n_1;n_2;n_3;n_4;n_5;n_6;n_7;n_8;n_9;n_10;n_11;n_12;
        n_13;n_14;n_15;n_16;n_17;n_18;n_19;n_20;n_21;n_22;n_23;n_24;n_25;n_26];     
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     end get connect
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    outside_voxels_index=connect_voxels_index(find(objects_image(connect_voxels_index)~=i_zone));
    
    d_z_outside_voxels=floor((outside_voxels_index-1)/(nrows*ncolumns))+1;
    d_y_outside_voxels=ceil(outside_voxels_index/nrows)-(ncolumns*(d_z_outside_voxels-1));
    d_x_outside_voxels=outside_voxels_index-nrows*(d_y_outside_voxels-1)-nrows*ncolumns*(d_z_outside_voxels-1);
    
    outside_voxels=[d_x_outside_voxels,d_y_outside_voxels,d_z_outside_voxels];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     start get connect index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      
    n_1=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_2=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_3=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_4=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_5=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_6=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_7=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_8=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];   
    n_9=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)-1];
    n_10=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_11=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_12=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_13=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_14=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_15=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_16=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_17=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)];
    n_18=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_19=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_20=[outside_voxels(1:length(outside_voxels(:,1,:)),1)-1,outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)+1]; 
    n_21=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_22=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_23=[outside_voxels(1:length(outside_voxels(:,1,:)),1),outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_24=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2)-1,outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_25=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2),outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
    n_26=[outside_voxels(1:length(outside_voxels(:,1,:)),1)+1,outside_voxels(1:length(outside_voxels(:,1,:)),2)+1,outside_voxels(1:length(outside_voxels(:,1,:)),3)+1];
              
           
    n_1_fake=n_1;
    n_2_fake=n_2;
    n_3_fake=n_3;
    n_4_fake=n_4;
    n_5_fake=n_5;
    n_6_fake=n_6;
    n_7_fake=n_7;
    n_8_fake=n_8;
    n_9_fake=n_9;
    n_10_fake=n_10;
    n_11_fake=n_11;
    n_12_fake=n_12;
    n_13_fake=n_13;
    n_14_fake=n_14;
    n_15_fake=n_15;
    n_16_fake=n_16;
    n_17_fake=n_17;
    n_18_fake=n_18;
    n_19_fake=n_19;
    n_20_fake=n_20;
    n_21_fake=n_21;
    n_22_fake=n_22;
    n_23_fake=n_23;
    n_24_fake=n_24; 
    n_25_fake=n_25;
    n_26_fake=n_26;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_1_fake(find(n_1(:,1)<=0),1)=n_1_fake(find(n_1(:,1)<=0),1)+1;         
    n_1_fake(find(n_1(:,2)<=0),2)=n_1_fake(find(n_1(:,2)<=0),2)+1;
	n_1_fake(find(n_1(:,3)<=0),3)=n_1_fake(find(n_1(:,3)<=0),3)+1;
  	n_2_fake(find(n_2(:,1)<=0),1)=n_2_fake(find(n_2(:,1)<=0),1)+1;
 	n_2_fake(find(n_2(:,2)<=0),2)=n_2_fake(find(n_2(:,2)<=0),2)+1;
 	n_2_fake(find(n_2(:,3)<=0),3)=n_2_fake(find(n_2(:,3)<=0),3)+1;
    n_3_fake(find(n_3(:,3)<=0),3)=n_3_fake(find(n_3(:,3)<=0),3)+1;
    n_3_fake(find(n_3(:,2)>ncolumns),2)=n_3_fake(find(n_3(:,2)>ncolumns),2)-1;
   	n_3_fake(find(n_3(:,1)<=0),1)=n_3_fake(find(n_3(:,1)<=0),1)+1;
  	n_4_fake(find(n_4(:,3)<=0),3)=n_4_fake(find(n_4(:,3)<=0),3)+1;
 	n_4_fake(find(n_4(:,2)<=0),2)=n_4_fake(find(n_4(:,2)<=0),2)+1;
    n_5_fake(find(n_5(:,3)<=0),3)=n_5_fake(find(n_5(:,3)<=0),3)+1;
    n_6_fake(find(n_6(:,3)<=0),3)=n_6_fake(find(n_6(:,3)<=0),3)+1;
   	n_6_fake(find(n_6(:,2)>ncolumns),2)=n_6_fake(find(n_6(:,2)>ncolumns),2)-1;
    n_7_fake(find(n_7(:,3)<=0),3)=n_7_fake(find(n_7(:,3)<=0),3)+1;
   	n_7_fake(find(n_7(:,2)<=0),2)=n_7_fake(find(n_7(:,2)<=0),2)+1;
 	n_7_fake(find(n_7(:,1)>nrows),1)=n_7_fake(find(n_7(:,1)>nrows),1)-1;
    n_8_fake(find(n_8(:,3)<=0),3)=n_8_fake(find(n_8(:,3)<=0),3)+1;
   	n_8_fake(find(n_8(:,1)>nrows),1)=n_8_fake(find(n_8(:,1)>nrows),1)-1;
    n_9_fake(find(n_9(:,3)<=0),3)=n_9_fake(find(n_9(:,3)<=0),3)+1;
  	n_9_fake(find(n_9(:,2)>ncolumns),2)=n_9_fake(find(n_9(:,2)>ncolumns),2)-1;
   	n_9_fake(find(n_9(:,1)>nrows),1)=n_9_fake(find(n_9(:,1)>nrows),1)-1;
   	n_10_fake(find(n_10(:,2)<=0),2)=n_10_fake(find(n_10(:,2)<=0),2)+1;
	n_10_fake(find(n_10(:,1)<=0),1)=n_10_fake(find(n_10(:,1)<=0),1)+1;     
    n_11_fake(find(n_11(:,1)<=0),1)=n_11_fake(find(n_11(:,1)<=0),1)+1;
   	n_12_fake(find(n_12(:,2)>ncolumns),2)=n_12_fake(find(n_12(:,2)>ncolumns),2)-1;
  	n_12_fake(find(n_12(:,1)<=0),1)=n_12_fake(find(n_12(:,1)<=0),1)+1;
    n_13_fake(find(n_13(:,2)<=0),2)=n_13_fake(find(n_13(:,2)<=0),2)+1;
   	n_14_fake(find(n_14(:,2)>ncolumns),2)=n_14_fake(find(n_14(:,2)>ncolumns),2)-1;
    n_15_fake(find(n_15(:,2)<=0),2)=n_15_fake(find(n_15(:,2)<=0),2)+1;
 	n_15_fake(find(n_15(:,1)>nrows),1)=n_15_fake(find(n_15(:,1)>nrows),1)-1;
    n_16_fake(find(n_16(:,1)>nrows),1)=n_16_fake(find(n_16(:,1)>nrows),1)-1;
    n_17_fake(find(n_17(:,2)>ncolumns),2)=n_17_fake(find(n_17(:,2)>ncolumns),2)-1;
   	n_17_fake(find(n_17(:,1)>nrows),1)=n_17_fake(find(n_17(:,1)>nrows),1)-1;
   	n_18_fake(find(n_18(:,3)>ndepth),3)=n_18_fake(find(n_18(:,3)>ndepth),3)-1;
   	n_18_fake(find(n_18(:,2)<=0),2)=n_18_fake(find(n_18(:,2)<=0),2)+1;
   	n_18_fake(find(n_18(:,1)<=0),1)=n_18_fake(find(n_18(:,1)<=0),1)+1;
    n_19_fake(find(n_19(:,3)>ndepth),3)=n_19_fake(find(n_19(:,3)>ndepth),3)-1;
   	n_19_fake(find(n_19(:,1)<=0),1)=n_19_fake(find(n_19(:,1)<=0),1)+1;
   	n_20_fake(find(n_20(:,3)>ndepth),3)=n_20_fake(find(n_20(:,3)>ndepth),3)-1;
  	n_20_fake(find(n_20(:,2)>ncolumns),2)=n_20_fake(find(n_20(:,2)>ncolumns),2)-1;
  	n_20_fake(find(n_20(:,1)<=0),1)=n_20_fake(find(n_20(:,1)<=0),1)+1;
    n_21_fake(find(n_21(:,3)>ndepth),3)=n_21_fake(find(n_21(:,3)>ndepth),3)-1;
	n_21_fake(find(n_21(:,2)<=0),2)=n_21_fake(find(n_21(:,2)<=0),2)+1;
   	n_22_fake(find(n_22(:,3)>ndepth),3)=n_22_fake(find(n_22(:,3)>ndepth),3)-1;
    n_23_fake(find(n_23(:,3)>ndepth),3)=n_23_fake(find(n_23(:,3)>ndepth),3)-1;
   	n_23_fake(find(n_23(:,2)>ncolumns),2)=n_23_fake(find(n_23(:,2)>ncolumns),2)-1;
    n_24_fake(find(n_24(:,3)>ndepth),3)=n_24_fake(find(n_24(:,3)>ndepth),3)-1;
 	n_24_fake(find(n_24(:,2)<=0),2)=n_24_fake(find(n_24(:,2)<=0),2)+1;
   	n_24_fake(find(n_24(:,1)>nrows),1)=n_24_fake(find(n_24(:,1)>nrows),1)-1;
    n_25_fake(find(n_25(:,3)>ndepth),3)=n_25_fake(find(n_25(:,3)>ndepth),3)-1;
  	n_25_fake(find(n_25(:,1)>nrows),1)=n_25_fake(find(n_25(:,1)>nrows),1)-1;
    n_26_fake(find(n_26(:,3)>ndepth),3)=n_26_fake(find(n_26(:,3)>ndepth),3)-1;
   	n_26_fake(find(n_26(:,2)>ncolumns),2)=n_26_fake(find(n_26(:,2)>ncolumns),2)-1;
	n_26_fake(find(n_26(:,1)>nrows),1)=n_26_fake(find(n_26(:,1)>nrows),1)-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
     
    
    %index=x+r(y-1)+r*c(z-1)
    
   
        n_1=n_1_fake(:,1)+nrows*(n_1_fake(:,2)-1)+nrows*ncolumns*(n_1_fake(:,3)-1);
        n_2=n_2_fake(:,1)+nrows*(n_2_fake(:,2)-1)+nrows*ncolumns*(n_2_fake(:,3)-1);
        n_3=n_3_fake(:,1)+nrows*(n_3_fake(:,2)-1)+nrows*ncolumns*(n_3_fake(:,3)-1);
        n_4=n_4_fake(:,1)+nrows*(n_4_fake(:,2)-1)+nrows*ncolumns*(n_4_fake(:,3)-1);
        n_5=n_5_fake(:,1)+nrows*(n_5_fake(:,2)-1)+nrows*ncolumns*(n_5_fake(:,3)-1);
        n_6=n_6_fake(:,1)+nrows*(n_6_fake(:,2)-1)+nrows*ncolumns*(n_6_fake(:,3)-1);
        n_7=n_7_fake(:,1)+nrows*(n_7_fake(:,2)-1)+nrows*ncolumns*(n_7_fake(:,3)-1);
        n_8=n_8_fake(:,1)+nrows*(n_8_fake(:,2)-1)+nrows*ncolumns*(n_8_fake(:,3)-1);
        n_9=n_9_fake(:,1)+nrows*(n_9_fake(:,2)-1)+nrows*ncolumns*(n_9_fake(:,3)-1);
        n_10=n_10_fake(:,1)+nrows*(n_10_fake(:,2)-1)+nrows*ncolumns*(n_10_fake(:,3)-1);
        n_11=n_11_fake(:,1)+nrows*(n_11_fake(:,2)-1)+nrows*ncolumns*(n_11_fake(:,3)-1);
        n_12=n_12_fake(:,1)+nrows*(n_12_fake(:,2)-1)+nrows*ncolumns*(n_12_fake(:,3)-1);
        n_13=n_13_fake(:,1)+nrows*(n_13_fake(:,2)-1)+nrows*ncolumns*(n_13_fake(:,3)-1);
        n_14=n_14_fake(:,1)+nrows*(n_14_fake(:,2)-1)+nrows*ncolumns*(n_14_fake(:,3)-1);
        n_15=n_15_fake(:,1)+nrows*(n_15_fake(:,2)-1)+nrows*ncolumns*(n_15_fake(:,3)-1);
        n_16=n_16_fake(:,1)+nrows*(n_16_fake(:,2)-1)+nrows*ncolumns*(n_16_fake(:,3)-1);
        n_17=n_17_fake(:,1)+nrows*(n_17_fake(:,2)-1)+nrows*ncolumns*(n_17_fake(:,3)-1);
        n_18=n_18_fake(:,1)+nrows*(n_18_fake(:,2)-1)+nrows*ncolumns*(n_18_fake(:,3)-1);
        n_19=n_19_fake(:,1)+nrows*(n_19_fake(:,2)-1)+nrows*ncolumns*(n_19_fake(:,3)-1);
        n_20=n_20_fake(:,1)+nrows*(n_20_fake(:,2)-1)+nrows*ncolumns*(n_20_fake(:,3)-1);
        n_21=n_21_fake(:,1)+nrows*(n_21_fake(:,2)-1)+nrows*ncolumns*(n_21_fake(:,3)-1);
        n_22=n_22_fake(:,1)+nrows*(n_22_fake(:,2)-1)+nrows*ncolumns*(n_22_fake(:,3)-1);
        n_23=n_23_fake(:,1)+nrows*(n_23_fake(:,2)-1)+nrows*ncolumns*(n_23_fake(:,3)-1);
        n_24=n_24_fake(:,1)+nrows*(n_24_fake(:,2)-1)+nrows*ncolumns*(n_24_fake(:,3)-1);
        n_25=n_25_fake(:,1)+nrows*(n_25_fake(:,2)-1)+nrows*ncolumns*(n_25_fake(:,3)-1);
        n_26=n_26_fake(:,1)+nrows*(n_26_fake(:,2)-1)+nrows*ncolumns*(n_26_fake(:,3)-1);
      
    if connect == 4    %%% 2D
        
        n_1=[];
        n_2=[];
        n_3=[];
        n_4=[];
        n_5=[];
        n_6=[];
        n_7=[];
        n_8=[];
        n_9=[];
        n_10=[];
        n_12=[];
        n_15=[];
        n_17=[];
        n_18=[];
        n_19=[];
        n_20=[];
        n_21=[];
        n_22=[];
        n_23=[];
        n_24=[];
        n_25=[];
        n_26=[];
        
    end
    
        
    if connect == 8   %%%% 2D
        
        n_1=[];
        n_2=[];
        n_3=[];
        n_4=[];
        n_5=[];
        n_6=[];
        n_7=[];
        n_8=[];
        n_9=[];
        n_18=[];
        n_19=[];
        n_20=[];
        n_21=[];
        n_22=[];
        n_23=[];
        n_24=[];
        n_25=[];
        n_26=[];
        
    end
    
        
      if connect == 6   %%%3D
        
        n_1=[];
        n_2=[];
        n_3=[];
        n_4=[];
        n_6=[];
        n_7=[];
        n_8=[];
        n_9=[];
        n_10=[];
        n_12=[];
        n_15=[];
        n_17=[];
        n_18=[];
        n_19=[];
        n_20=[];
        n_21=[];
        n_23=[];
        n_24=[];
        n_25=[];
        n_26=[];
        
    end
          
       
    connect_voxels_index=[n_1;n_2;n_3;n_4;n_5;n_6;n_7;n_8;n_9;n_10;n_11;n_12;
        n_13;n_14;n_15;n_16;n_17;n_18;n_19;n_20;n_21;n_22;n_23;n_24;n_25;n_26];     
             
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     end get connect index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    surface_area_index=find(objects_image(connect_voxels_index)==i_zone);
    
    temp_image=objects_image;
    temp_image(:)=0;
    temp_image(connect_voxels_index(surface_area_index))=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Boundaries of particels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    object_boundary_index=find(temp_image(:)==1);
    d_z_boundary=floor((object_boundary_index-1)/(nrows*ncolumns))+1;
    d_y_boundary=ceil(object_boundary_index/nrows)-(ncolumns*(d_z_boundary-1));
    d_x_boundary=object_boundary_index-nrows*(d_y_boundary-1)-nrows*ncolumns*(d_z_boundary-1);
    
    object_boundary=[d_x_boundary,d_y_boundary,d_z_boundary];    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
              
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%
    %%%%%        Centers, sizes, orintations, lengths, diameter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    temp_objects_sizes=length(object_voxels_final_index);
    surface_area=length(surface_area_index);
        
    if isempty(surface_area_index)==0
    
      
    centroid_x=(sum(d_x_object_final)/temp_objects_sizes);
	centroid_y=(sum(d_y_object_final)/temp_objects_sizes);
	centroid_z=(sum(d_z_object_final)/temp_objects_sizes);
    
    temp_objects_centroids=([centroid_x,centroid_y,centroid_z]);
   	objects_centroids=[objects_centroids;temp_objects_centroids];
    objects_sizes=[objects_sizes;temp_objects_sizes];
    objects_surface_area=[objects_surface_area;length(surface_area_index)];
    
    min_x_boundary=min(object_voxels_final(:,1));
	max_x_boundary=max(object_voxels_final(:,1));
	min_y_boundary=min(object_voxels_final(:,2));
	max_y_boundary=max(object_voxels_final(:,2));
	min_z_boundary=min(object_voxels_final(:,3));
	max_z_boundary=max(object_voxels_final(:,3));
	
    clear aa* n_*_*  min_*_* max_*_* 
   
	diameter_sphere_matrix=[(centroid_x-object_boundary(:,1)).^2 ,...
     (centroid_y-object_boundary(:,2)).^2,(centroid_z-object_boundary(:,3)).^2];
	
     
    diameter_sphere_matrix=sum(diameter_sphere_matrix,2);
	   
    diameter_sphere=(round(sqrt(diameter_sphere_matrix)));
   
    clear d_*_* dd_*_*
    
    Irr=sum((object_voxels_final(:,2)-centroid_y).^2)+sum((object_voxels_final(:,3)-centroid_z).^2);
	Icc=sum((object_voxels_final(:,1)-centroid_x).^2)+sum((object_voxels_final(:,3)-centroid_z).^2);
	Idd=sum((object_voxels_final(:,1)-centroid_x).^2)+sum((object_voxels_final(:,2)-centroid_y).^2);
	Irc=sum(((object_voxels_final(:,1)-centroid_x)).*((object_voxels_final(:,2)-centroid_y)));
	Ird=sum(((object_voxels_final(:,1)-centroid_x)).*((object_voxels_final(:,3)-centroid_z)));
	Icd=sum(((object_voxels_final(:,2)-centroid_y)).*((object_voxels_final(:,3)-centroid_z)));
	
	moment_matrix=[Irr -Irc -Ird; -Irc Icc -Icd; -Ird -Icd Idd];
	
	[v,d]=eig(moment_matrix);
	
	length_points=[];
	
	length_points=v*object_voxels_final(1:length(object_voxels_final(:,1)),:)';
	length_points=length_points';
    % get the centroid in the pricinpla coordinates
    length_centroid = v*temp_objects_centroids';
    length_centroid = length_centroid';
    
	
    temp_objects_angles=[acos(abs(v(3,:)))*360/(2*pi)];
    objects_angles=[objects_angles;temp_objects_angles];
    
    temp_principal_vectors = [v(:,1)' v(:,2)' v(:,3)'];
    objects_principal_vectors = [objects_principal_vectors; temp_principal_vectors];
  
	%temp_objects_lengths=[range(length_points(:,1))+1,range(length_points(:,2))+1,...
    %range(length_points(:,3))+1];
    
	range_x=(max(length_points(:,1))-min(length_points(:,1)))+1;
    range_y=(max(length_points(:,2))-min(length_points(:,2)))+1;
    range_z=(max(length_points(:,3))-min(length_points(:,3)))+1;
    temp_objects_lengths=[range_x,range_y,range_z];

	objects_lengths=[objects_lengths;temp_objects_lengths];
    
    % translate the length_points into the principla coordinates with
    % lenght_centroid as the origin
    length_points(:,1) = length_points(:,1)-length_centroid(1);
    length_points(:,2) = length_points(:,2)-length_centroid(2);
    length_points(:,3) = length_points(:,3)-length_centroid(3);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % include below to have the shift center
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % find the shift x, y and z
    
    % find the 6 boundary points
    x_max = max(length_points(:,1));
    x_min = min(length_points(:,1));
    
    y_max = max(length_points(:,2));
    y_min = min(length_points(:,2));
    
    z_max = max(length_points(:,3));
    z_min = min(length_points(:,3));
    % these 6 points are used to construct the eight planes
    
    % find the volume_above and volume_below of that plane in each octant
    % if the volume of the present voxel and the other three points is
    % positive, then above volume; if negative, then below volume; if zero,
    % 0.5 for above and 0.5 for below
    above_v1 = 0;   % the first octant
    below_v1 = 0;
    
    above_v2 = 0;
    below_v2 = 0;
    
    above_v3 = 0;
    below_v3 = 0;
    
    above_v4 = 0;
    below_v4 = 0;
    
    above_v5 = 0;
    below_v5 = 0;
    
    above_v6 = 0;
    below_v6 = 0;
    
    above_v7 = 0;
    below_v7 = 0;
    
    above_v8 = 0;
    below_v8 = 0;
    
    for i_p=1:length(length_points(:,1))    % loop over all the voxel points
        x_p = length_points(i_p,1);
        y_p = length_points(i_p,2);
        z_p = length_points(i_p,3);

        if(x_p>=0 && y_p>=0 && z_p>=0)  % 1st quadrant
            temp_matrix = [1 1 1 1; x_max 0 0 x_p; 0 y_max 0 y_p; 0 0 z_max z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0 % means this voxel point is above the plane
                above_v1 = above_v1+1;
            end
            if temp_v<0 % means this voxel point is below the plane
                below_v1 = below_v1+1;
            end
            if temp_v==0    % means this voxel point is in the plane
                above_v1 = above_v1+0.5;
                below_v1 = below_v1+0.5;
            end
        end
        
        if(x_p<=0 && y_p>=0 && z_p>=0)  % 2nd quadrant
            temp_matrix = [1 1 1 1; x_min 0 0 x_p; 0 0 y_max y_p; 0 z_max 0 z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v2 = above_v2+1;
            end
            if temp_v<0
                below_v2 = below_v2+1;
            end
            if temp_v==0
                above_v2 = above_v2+0.5;
                below_v2 = below_v2+0.5;
            end
        end
        
        if(x_p<=0 && y_p<=0 && z_p>=0)  % 3rd quadrant
            temp_matrix = [1 1 1 1; x_min 0 0 x_p; 0 y_min 0 y_p; 0 0 z_max z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v3 = above_v3+1;
            end
            if temp_v<0
                below_v3 = below_v3+1;
            end
            if temp_v==0
                above_v3 = above_v3+0.5;
                below_v3 = below_v3+0.5;
            end
        end
        
        if(x_p>=0 && y_p<=0 && z_p>=0)  % 4th quadrant
            temp_matrix = [1 1 1 1; 0 x_max 0 x_p; y_min 0 0 y_p; 0 0 z_max z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v4 = above_v4+1;
            end
            if temp_v<0
                below_v4 = below_v4+1;
            end
            if temp_v==0
                above_v4 = above_v4+0.5;
                below_v4 = below_v4+0.5;
            end
        end
        
        if(x_p>=0 && y_p>=0 && z_p<=0)  % 5th quadrant
            temp_matrix = [1 1 1 1; x_max 0 0 x_p; 0 0 y_max y_p; 0 z_min 0 z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v5 = above_v5+1;
            end
            if temp_v<0
                below_v5 = below_v5+1;
            end
            if temp_v==0
                above_v5 = above_v5+0.5;
                below_v5 = below_v5+0.5;
            end
        end
        
        if(x_p<=0 && y_p>=0 && z_p<=0)  % 6th quadrant
            temp_matrix = [1 1 1 1; x_min 0 0 x_p; 0 y_max 0 y_p; 0 0 z_min z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v6 = above_v6+1;
            end
            if temp_v<0
                below_v6 = below_v6+1;
            end
            if temp_v==0
                above_v6 = above_v6+0.5;
                below_v6 = below_v6+0.5;
            end
        end
        
        if(x_p<=0 && y_p<=0 && z_p<=0)  % 7th quadrant
            temp_matrix = [1 1 1 1; x_min 0 0 x_p; 0 0 y_min y_p; 0 z_min 0 z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v7 = above_v7+1;
            end
            if temp_v<0
                below_v7 = below_v7+1;
            end
            if temp_v==0
                above_v7 = above_v7+0.5;
                below_v7 = below_v7+0.5;
            end
        end
        
        if(x_p>=0 && y_p<=0 && z_p<=0)  % 8th quadrant
            temp_matrix = [1 1 1 1; x_max 0 0 x_p; 0 y_min 0 y_p; 0 0 z_min z_p];
            temp_v = det(temp_matrix)/6;
            if temp_v>0
                above_v8 = above_v8+1;
            end
            if temp_v<0
                below_v8 = below_v8+1;
            end
            if temp_v==0
                above_v8 = above_v8+0.5;
                below_v8 = below_v8+0.5;
            end
        end
    end
    
    % calculate the volume of the tet that constrained by the the three
    % points and the origin in each octant
    plane_v1 = det([1 1 1 1; 0 x_max 0 0; 0 0 y_max 0; 0 0 0 z_max])/6;
    plane_v2 = det([1 1 1 1; x_min 0 0 0; 0 0 y_max 0; 0 0 0 z_max])/6;
    plane_v3 = det([1 1 1 1; x_min 0 0 0; 0 y_min 0 0; 0 0 0 z_max])/6;
    plane_v4 = det([1 1 1 1; 0 x_max 0 0; y_min 0 0 0; 0 0 0 z_max])/6;
    
    plane_v5 = det([1 1 1 1; 0 0 x_max 0; 0 y_max 0 0; 0 0 0 z_min])/6;
    plane_v6 = det([1 1 1 1; 0 x_min 0 0; 0 0 y_max 0; 0 0 0 z_min])/6;
    plane_v7 = det([1 1 1 1; 0 0 x_min 0; 0 y_min 0 0; 0 0 0 z_min])/6;
    plane_v8 = det([1 1 1 1; 0 x_max 0 0; 0 0 y_min 0; 0 0 0 z_min])/6;
    
    % find volume_1 to volume_8
    volume_1 = above_v1-(plane_v1-below_v1);
    volume_2 = above_v2-(plane_v2-below_v2);
    volume_3 = above_v3-(plane_v3-below_v3);
    volume_4 = above_v4-(plane_v4-below_v4);
    volume_5 = above_v5-(plane_v5-below_v5);
    volume_6 = above_v6-(plane_v6-below_v6);
    volume_7 = above_v7-(plane_v7-below_v7);
    volume_8 = above_v8-(plane_v8-below_v8);
    
    % find shift x, y and z
    shift_x = (volume_1+volume_4+volume_8+volume_5...
              -volume_2-volume_3-volume_7-volume_6)/temp_objects_sizes*range_x;
    shift_y = (volume_1+volume_2+volume_6+volume_5...
              -volume_4-volume_3-volume_7-volume_8)/temp_objects_sizes*range_y;
    shift_z = (volume_1+volume_2+volume_3+volume_4...
              -volume_5-volume_6-volume_7-volume_8)/temp_objects_sizes*range_z;
    
    % translate the origin to the new shifted x, y and z
    length_points(:,1) = length_points(:,1)-shift_x;
    length_points(:,2) = length_points(:,2)-shift_y;
    length_points(:,3) = length_points(:,3)-shift_z;
    
    %%%%%%find the new center in the original x, y and z axises
    %%%%%%temp_objects_new_center = temp_objects_centroids+(v\[shift_x; shift_y; shift_z])';
    %
    %%%%%%objects_new_center = [objects_new_center; temp_objects_new_center];   
    
while(1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % include above to have the shift center
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % find the parts in eight quadrants of this particle
    length_points_1 = [];
    length_points_2 = [];
    length_points_3 = [];
    length_points_4 = [];
    length_points_5 = [];
    length_points_6 = [];
    length_points_7 = [];
    length_points_8 = [];   % the eight quadrant parts
    for i_p=1:length(length_points(:,1))    % loop over all the voxel points
        x_p = length_points(i_p,1);
        y_p = length_points(i_p,2);
        z_p = length_points(i_p,3);
        
        temp_p = [x_p,y_p,z_p];
        if(x_p>=0 && y_p>=0 && z_p>=0)  % 1st quadrant
            length_points_1 = [length_points_1;temp_p];
        end
        if(x_p<=0 && y_p>=0 && z_p>=0)  % 2nd quadrant
            length_points_2 = [length_points_2;temp_p];
        end
        if(x_p<=0 && y_p<=0 && z_p>=0)  % 3rd quadrant
            length_points_3 = [length_points_3;temp_p];
        end
        if(x_p>=0 && y_p<=0 && z_p>=0)  % 4th quadrant
            length_points_4 = [length_points_4;temp_p];
        end
        if(x_p>=0 && y_p>=0 && z_p<=0)  % 5th quadrant
            length_points_5 = [length_points_5;temp_p];
        end
        if(x_p<=0 && y_p>=0 && z_p<=0)  % 6th quadrant
            length_points_6 = [length_points_6;temp_p];
        end
        if(x_p<=0 && y_p<=0 && z_p<=0)  % 7th quadrant
            length_points_7 = [length_points_7;temp_p];
        end
        if(x_p>=0 && y_p<=0 && z_p<=0)  % 8th quadrant
            length_points_8 = [length_points_8;temp_p];
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % include below to have the shift center
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if(length(length_points_1) > 0 & length(length_points_2) > 0 ...
     & length(length_points_3) > 0 & length(length_points_4) > 0 ...
     & length(length_points_5) > 0 & length(length_points_6) > 0 ...
     & length(length_points_7) > 0 & length(length_points_8) > 0)
        break;  % no empty octant, then need not to change shift_x, go out while
    end
    
    disp('shift too much at particle '), disp(i_zone)
    length_points(:,1) = length_points(:,1)+0.1*shift_x;
    length_points(:,2) = length_points(:,2)+0.1*shift_y;
    length_points(:,3) = length_points(:,3)+0.1*shift_z;
    shift_x = 0.9*shift_x;
    shift_y = 0.9*shift_y;
    shift_z = 0.9*shift_z;
end
    % find the new center in the original x, y and z axises
    temp_objects_new_center = temp_objects_centroids+(v\[shift_x; shift_y; shift_z])';
    
    objects_new_center = [objects_new_center; temp_objects_new_center];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % include above to have the shift center
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % the lengths of the equivalent polyellipsoid
    poly_x1 = max(length_points_1(:,1))+0.5; % the 1st quadrant
    poly_y1 = max(length_points_1(:,2))+0.5;
    poly_z1 = max(length_points_1(:,3))+0.5;
    
    poly_x2 = -min(length_points_2(:,1))-0.5; % the 2nd quadrant
    poly_y2 = max(length_points_2(:,2))+0.5;
    poly_z2 = max(length_points_2(:,3))+0.5;
    
    poly_x3 = -min(length_points_3(:,1))-0.5; % the 3rd quadrant
    poly_y3 = -min(length_points_3(:,2))-0.5;
    poly_z3 = max(length_points_3(:,3))+0.5;
    
    poly_x4 = max(length_points_4(:,1))+0.5; % the 4th quadrant
    poly_y4 = -min(length_points_4(:,2))-0.5;
    poly_z4 = max(length_points_4(:,3))+0.5;
    
    poly_x5 = max(length_points_5(:,1))+0.5; % the 5th quadrant
    poly_y5 = max(length_points_5(:,2))+0.5;
    poly_z5 = -min(length_points_5(:,3))-0.5;
    
    poly_x6 = -min(length_points_6(:,1))-0.5; % the 6th quadrant
    poly_y6 = max(length_points_6(:,2))+0.5;
    poly_z6 = -min(length_points_6(:,3))-0.5;
    
    poly_x7 = -min(length_points_7(:,1))-0.5; % the 7th quadrant
    poly_y7 = -min(length_points_7(:,2))-0.5;
    poly_z7 = -min(length_points_7(:,3))-0.5;
    
    poly_x8 = max(length_points_8(:,1))+0.5; % the 8th quadrant
    poly_y8 = -min(length_points_8(:,2))-0.5;
    poly_z8 = -min(length_points_8(:,3))-0.5;
    
    % the pricinpal lengths of the equivalent poly-ellipsoid
    x_plus = (poly_x1+poly_x4+poly_x8+poly_x5)/4;
    x_minus = (poly_x2+poly_x3+poly_x7+poly_x6)/4;
    
    y_plus = (poly_y1+poly_y2+poly_y6+poly_y5)/4;
    y_minus = (poly_y4+poly_y3+poly_y7+poly_y8)/4;
    
    z_plus = (poly_z1+poly_z2+poly_z3+poly_z4)/4;
    z_minus = (poly_z5+poly_z6+poly_z7+poly_z8)/4;
    
    temp_poly_lengths = [x_plus, x_minus, y_plus, y_minus, z_plus, z_minus];
    object_poly_lengths = [object_poly_lengths; temp_poly_lengths];
    
    % end of pricinpal lengths of poly-ellipsoid
    
    
    clear I*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% End centers, sizes, orintations, lengths, diameter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    sphere_index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D=diameter_sphere;
    D_max=.5*max(temp_objects_lengths);
    D_min=.5*min(temp_objects_lengths);
    
    temp_sphere_index=mean(abs((D/D_max)-(D/D_min)));
    
   % temp_sphere_index=(sum(abs((diameter_sphere/.5*min(temp_objects_lengths))-(diameter_sphere/.5*max(temp_objects_lengths)))))/length(boundary);
    
    objects_sphereness=[objects_sphereness;temp_sphere_index];
    temp_round_index=(length(object_boundary)/(4*pi*(.25*min(temp_objects_lengths)+.25*max(temp_objects_lengths))^2));
    objects_roundness=[objects_roundness;temp_round_index];
	
    clear D D_*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     End sphere_index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%  Start Coordination Number
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
    coord_voxels=outside_voxels_index;
    coord_voxels_null=find(objects_image(coord_voxels)==0);
    coord_voxels(coord_voxels_null)=[];
    coord_voxels_data=sort(objects_image(coord_voxels));
    match_matrix=[double(min(coord_voxels_data)):double(max(coord_voxels_data))];
    [coord_voxels_data_index,match_matrix_index]=wcommon(coord_voxels_data,match_matrix);
    coord_numbers_index=find(match_matrix_index==1);
    coord_number=match_matrix(coord_numbers_index);  % The contact numbers of spheres
    temp_coord_number=[ i_zone length(coord_numbers_index)];   %coord number of the current i_zone (i.e., shpere)
    
    objects_coord_numbers=[objects_coord_numbers;temp_coord_number];
    all_coord_numbers=[all_coord_numbers;coord_number'];
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%  End Coordination Number
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%  Start Collecting data
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   collect_i_zone=[collect_i_zone;i_zone];
   collect_object_voxels_final_index=[collect_object_voxels_final_index;object_voxels_final_index];
   pointer_collect_object_voxels_final_index=[pointer_collect_object_voxels_final_index;length(object_voxels_final_index)];
   collect_coord_voxels=[collect_coord_voxels;coord_voxels];
   pointer_collect_coord_voxels=[pointer_collect_coord_voxels;length(coord_voxels)];
   collect_coord_number=[collect_coord_number;coord_number'];
   pointer_collect_coord_number=[pointer_collect_coord_number;length(coord_number)];
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%  End Collecting data
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   
    else
        
    objects_image(object_voxels_final_index)=0;
  
    higher_values_2=find(object_voxels_final_index> i_zone);
    objects_image(higher_values_2)=double(objects_image(higher_values_2))-1;
    
    
    end
    
    
    i_zone=i_zone+1;
    
    end
   
      
    cumsum_object_voxels_final_index=[0;cumsum(pointer_collect_object_voxels_final_index)];
    cumsum_coord_voxels=[0;cumsum(pointer_collect_coord_voxels)];
    cumsum_coord_number=[0;cumsum(pointer_collect_coord_number)];
    
    
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%  Find contact information
   %%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%  Coordination Number
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i_zone=1:length(collect_i_zone)
%         
%     fake_i_zone=collect_i_zone(i_zone);
%         
%               
%     coord_voxels=collect_coord_voxels(cumsum_coord_voxels(i_zone)+1:cumsum_coord_voxels(i_zone+1));
%    
%     coord_number=collect_coord_number(cumsum_coord_number(i_zone)+1:cumsum_coord_number(i_zone+1));
%     
%     for i=1:length(coord_number)
%         
%         xxxx=find(collect_i_zone==coord_number(i));
%         
%         all_contact_points_index=collect_object_voxels_final_index((cumsum_object_voxels_final_index(xxxx)+1):cumsum_object_voxels_final_index(xxxx+1));
%         
%         some_contact_points_index=coord_voxels(find(objects_image(coord_voxels)==coord_number(i)));   
%         
%         %%%%%%%%%%%%%%%  Controid of contacs of points %%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         z_some_contact_points=floor((some_contact_points_index-1)/(nrows*ncolumns))+1;
%         y_some_contact_points=ceil(some_contact_points_index/nrows)-(ncolumns*(z_some_contact_points-1));
%         x_some_contact_points=some_contact_points_index-nrows*(y_some_contact_points-1)-nrows*ncolumns*(z_some_contact_points-1);
%         
%         center_some_contact_coordinats=[x_some_contact_points,y_some_contact_points,z_some_contact_points];
%         
%         center_some_contact=[sum(center_some_contact_coordinats(:,1))/length(center_some_contact_coordinats(:,1)),...
%                           sum(center_some_contact_coordinats(:,2))/length(center_some_contact_coordinats(:,1)),...
%                           sum(center_some_contact_coordinats(:,3))/length(center_some_contact_coordinats(:,1))];
%         
%         
%         temp_objects_contacts_center=[fake_i_zone, coord_number(i),center_some_contact(1),center_some_contact(2),center_some_contact(3)];
%         
%                                                                     
%          objects_contacts_center=[ objects_contacts_center; temp_objects_contacts_center];
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        
%         %%%%%%%%%%%%%%%%%%%
%         z_all_contact_points=floor((all_contact_points_index-1)/(nrows*ncolumns))+1;
%         y_all_contact_points=ceil(all_contact_points_index/nrows)-(ncolumns*(z_all_contact_points-1));
%         x_all_contact_points=all_contact_points_index-nrows*(y_all_contact_points-1)-nrows*ncolumns*(z_all_contact_points-1);
%        
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         %%%%%%%  Center of  contacts and vectors  %%%%%%%%%%%%%%%%%%%%%
%         %center_contact_coordinats=[x_all_contact_points,...
%         %y_all_contact_points-ncolumns*(floor(all_contact_points_index/(nrows*ncolumns))),...
%         %floor(all_contact_points_index/(nrows*ncolumns))+1];
%         
%         center_contact_coordinats=[x_all_contact_points,y_all_contact_points,z_all_contact_points];
%         
%         centroid_contact=[sum(center_contact_coordinats(:,1))/length(center_contact_coordinats(:,1)),...
%                           sum(center_contact_coordinats(:,2))/length(center_contact_coordinats(:,1)),...
%                           sum(center_contact_coordinats(:,3))/length(center_contact_coordinats(:,1))];
%         
%                           
%                           
%        temp_center_center_vector=[fake_i_zone, coord_number(i), centroid_contact(:,1)-objects_centroids(i_zone,1),...
%                                      centroid_contact(:,2)-objects_centroids(i_zone,2),...
%                                     centroid_contact(:,3)-objects_centroids(i_zone,3)];
%                                  
%          objects_center_center_vector=[ objects_center_center_vector; temp_center_center_vector];
%          
%          
%          
%         
%         
%         %%%%  objects_center_center_vector=[ particle #, Contact Label, X
%         %%%% component of vector,Y component of vector,Z component of vector]
%                                      
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         [any_value, some_contact_points_index_match]=wcommon(some_contact_points_index,all_contact_points_index);
%        
%         contact_points_matrix=[x_all_contact_points(find(some_contact_points_index_match==1)),...
%         y_all_contact_points(find(some_contact_points_index_match==1)),...
%         all_contact_points_index(find(some_contact_points_index_match==1))];
%         
%         temp_contact_points_coordinats=[contact_points_matrix(:,1),...
%         contact_points_matrix(:,2),...
%         floor(contact_points_matrix(:,3)/(nrows*ncolumns))+1];
%     
%         %%%%%%%% coordinates of contact points
%         
%         %contact_points_coordinats=[contact_points_coordinats;temp_contact_points_coordinats];
%         contact_points_coordinats=temp_contact_points_coordinats;
%         
%           
%         [coeff]=princomp(contact_points_coordinats);
%           
%         current_normal_vector=coeff(:,3)';
%              
%         temp_normal_vector=[fake_i_zone, coord_number(i), current_normal_vector];
%         objects_normal_vector=[ objects_normal_vector; temp_normal_vector];
%          
%         %%%%  objects_normal_vector=[ particle #, Contact Label, X, component of vector,Y component of vector,Z component of
%         %%%% vector]
%         
%         current_tangent_vector=coeff(:,1)';
%               
%         temp_tangent_vector=[fake_i_zone,coord_number(i), current_tangent_vector];
%         objects_tangent_vector=[ objects_tangent_vector; temp_tangent_vector];
%         
%         
%         %%%%  objects_tangent_vector=[ particle #, Contact Label, X, component of vector,Y component of vector,Z component of
%         %%%% vector]
%               
%     end
%     end
%     
    

    
    clear small_* object_index object_voxels_final objects_diameter_sphere row_* range_* v
    clear column_* depth_* copy_image
    
   
    save('all_result_10_11.mat')
    
    
    
    
    
    
    
    
    
    
