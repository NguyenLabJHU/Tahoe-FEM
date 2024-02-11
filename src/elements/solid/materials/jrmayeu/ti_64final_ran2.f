c===================================================================
c
c  This UMAT subroutine implements the following model in ABAQUS.
c-------------------------------------------------------------------

      subroutine ti_umat (	stress,  statev,  ddsdde,  sse,     spd,
     &			scd,     rpl,     ddsddt,  drplde,  drpldt,
     &			strain,  dstrain, time,    dtime,   temp,
     &			dtemp,   predef,  dpred,   cmname,  ndi,
     &			nshr,    ntens,   nstatv,  props,   nprops,
     &			coords,  drot,    pnewdt,  celent,  dfgrd0,
     &			dfgrd1,  noel,    npt,     layer,   kspt,
     &			kstep,   kinc )

      include 'ABA_PARAM.INC'

      parameter(num_slip_sys	= 24,
     & 		max_grains	= 1,
     &		max_loops	= 10,
     &		Tolerance	= 0.00001,
     &        noel1=9)

c-------------------------------------------------------------------
c  Common Block variables are used for writing out to the text files,   
c  using the ABAQUS user subroutine uexternaldb()
c-------------------------------------------------------------------

      common/flag1/storsdv(9,4),jnoel(9),xcord(9,3),kincp,
     &icount
       
      character*7 cmname

      logical Converged, Improved

c-------------------------------------------------------------------
c  Dimension arrays passed into the UMAT sub
c-------------------------------------------------------------------

      dimension
     &	coords(3),	! Coordinates of Gauss pt. being evaluated
     &	ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     &	ddsddt(ntens),	! Change in stress per change in temperature
     &	dfgrd0(3,3),	! Deformation gradient at beginning of step
     &	dfgrd1(3,3),	! Deformation gradient at end of step
     &	dpred(1),	! Change in predefined state variables
     &	drplde(ntens),	! Change in heat generation per change in strain
     &	drot(3,3),	! Rotation matrix
     &	dstrain(ntens),	! Strain increment tensor stored in vector form
     &	predef(1),	! Predefined state vars dependent on field variables
     &	props(nprops),	! Material properties passed in
     &	statev(nstatv),	! State Variables
     &	strain(ntens),	! Strain tensor stored in vector form
     &	stress(ntens),	! Cauchy stress tensor stored in vector form
     &	time(2)		! Step Time and Total Time

c-------------------------------------------------------------------
    
c  Dimension other arrays used in this UMAT
c-------------------------------------------------------------------
    

      dimension
     &	array1	(3,3),		! Dummy array
     &	array2	(3,3),		! Dummy array
     &	array3	(num_slip_sys,num_slip_sys), ! A_matrix of Newton Raphson
     &	array4	(6,6),		! Dummy array used in Voigt notation
     &	array5	(6,6),		! Inverse of array4().
     &	array6	(3,3,3,3),	! 4th rank dummy array
     &	array7	(3,3,3,3),	! Another 4th rank dummy array
     &	a0	(num_slip_sys), ! Kinematic stress at beginning of step
     &	a	(num_slip_sys), ! Kinematic stress at end of step
     &    C_IV(6,6),      ! Local elastic stiffness tensor in Voight notation
     &	C0	(3,3,3,3),	! Local  4th rank elastic stiffness tensor
     &	C	(3,3,3,3),	! Global 4th rank elastic stiffness tensor
     &	C_avg	(3,3,3,3),	! Average over all grains
     &	del	(3,3),		! Kronecker delta tensor
     &	ddpdsig	(3,3,3,3),	! deriv of D_p wrt sig * dt
     &	ddsdde_4th(3,3,3,3),	! 4th rank tangent stiffness tensor.
     &	daadga	(num_slip_sys),	! deriv of a_alpha wrt del_gamma_alpha
     &	ddadgb	(num_slip_sys,num_slip_sys),! deriv of d_alpha wrt del_gamma_beta
     &	dtadgb	(num_slip_sys,num_slip_sys),! deriv of t_alpha wrt del_gamma_beta
     &	dt1adgb	(num_slip_sys,num_slip_sys),! deriv of t_alpha wrt del_gamma_beta
     &	d_gamma_dot(num_slip_sys),! Gamma_dot step from Newt-Raphson
     &	dlpdsd	(3,3,3,3),	! deriv of xL_p wrt sigma_dot
     &	dir_cos	(3,3),! Direction cosines 
     &	dtaudgd	(num_slip_sys,num_slip_sys), ! Deriv of Tau wrt Gamma_dot
     &	E_el	(3,3),		! Elastic Green strain tensor
     &	E_p	(3,3),		! Plastic strain tensor
     &	E_tot	(3,3),		! Green strain tensor E_el + E_p
     &    E_dev   (3,3),  ! Deviatoric part of Green strain tensor
     &    E_p_dev   (3,3),  ! Deviatoric part of Plastic Green strain tensor
     &	F0	(3,3),		! F at beginning of sub increment
     &	F1	(3,3),		! F at end of sub increment
     &	F_el	(3,3),		! Elastic part of F
     &	F_el_inv(3,3),		! Inverse of elastic part of F
     &	F_dot	(3,3),		! F-dot
     &	F_inv	(3,3),		! Inverse of deformation gradient
     &	F_p_inv_0(3,3),! Inverse of plastic part of F at beginning
     &	F_p_inv(3,3),! Inverse of plastic part of F at end of step
     &	F_p(3,3),!  plastic part of F at end of step
     &	F_p_t(3,3),!  Transpose of plastic part of F at end of step
     &	func	(num_slip_sys),	! Function that is solved to get gamma_dot(k)
     &	d0	(num_slip_sys), ! drag stress at beginning of step
     &	d	(num_slip_sys), ! drag stress at end of step
     &	t0	(num_slip_sys), ! threshold stress at beginning of step
     &	t	(num_slip_sys), ! threshold stress at end of step 
     &    t1  (num_slip_sys), ! threshold stress at end of step
     &	gamma_dot(num_slip_sys),! shear strain rate on system
     &	gamma_try(num_slip_sys),! Candidate gamma_dots for line search
     &	grad(num_slip_sys),	! gradient of sum-sqares-error
     &	psi	(3),	! Euler angles
     &	sig0	(3,3),		! Stress tensor at beginning of step.
     &	sig	(3,3),		! Stress tensor
     &	sig_avg	(3,3),		! Rate of stress change for a grain
     &	Spk2	(3,3),		! 2nd Piola Kirkhoff stress
     &	tau	(num_slip_sys),! Resolved shear stress for slip dir.
     &	xL	(3,3),		! Eulerian velocity gradient
     &	xL_p	(3,3),		! Plastic vel grad in current configuration
     &	xs0	(3,num_slip_sys),! Inter config slip directions in global coords
     &	xs	(3,num_slip_sys),! Current config slip directions in global coords
     &	xm0	(3,num_slip_sys),! Inter config plane normals in global coords
     &	xm	(3,num_slip_sys),! Current config plane normals in global coords
     &	y	(3,num_slip_sys),! Miller indices of slip plane normals 
     &	z	(3,num_slip_sys), ! Miller indices of slip directions
     &    dummy(num_slip_sys,num_slip_sys),
c  Other arrays
     &    delta_E_p(3,3),!increment of plastic strain tensor
     &	gamma(num_slip_sys),! Accumulated shear strain
     &    delta_gamma(num_slip_sys), !increment of shear strain on each slip sys
     &    D_p(3,3),
     &    euler(3),
     &    bor(12,3,3),
     &	psi_b(3,12),
     &    text(3,3),
     &    xs1(3,num_slip_sys),
     &    xm1(3,num_slip_sys),
     &    ps(3),
     &    an(3,3),
     &    gamma_p(6)

      Pi = 4. * atan(1.)

c-------------------------------------------------------------------
c   Determine number of grains.  Divide the total number of ISVs
c   by the number of ISVs per grain.
c   3 - Euler angles
c   9 - F_p(i,j)
c  24 - g(i) for 3 slip systems
c  24 - a(i) for 3 slip systems
c ----
c  60
c-------------------------------------------------------------------

      num_grains = 1
     
c-------------------------------------------------------------------
c  Assign props() array to logical variable names
c-------------------------------------------------------------------
 
       C11 = props(1)     
       C12 = props(2)
	 C13 = props(3)
	 C33 = props(4)
       C44 = props(5)
       gamma_dot_zero = props(6)
       flow_exp	= props(7)
       d_zero_p	= props(8)
       Hdir	= props(9)
       Hdyn	= props(10)
       xLatent	= props(11)
       a_zero    = props(12)
       Adir      = props(13)*props(14)
       Adyn	= props(13)
       psi_ang = props(15)
	 theta_ang = props(16)
	 phi_ang = props(17)
	 t_zero = props(18)
       a1 = props(19)
       phase = props(20)

       c_a = 1.588
c-------------------------------------------------------------------
c  Initialize Kronnecker delta tensor
c-------------------------------------------------------------------

      do i = 1,3
         do j = 1,3
            del(i,j) = 0.0
         end do
         del(i,i) = 1.0
      end do

c-------------------------------------------------------------------
c  Assign slip system normals and slip directions for an HCP material.
c  y(i,j) - plane normal / z(i,j) - slip direction
c-------------------------------------------------------------------

c  Basal Slip Systems (0001)<1120>
c  (0001)[-1 2 -1 0]  
	 y(1,1) = 0.
	 y(2,1) = 0.
	 y(3,1) = 1.
       z(1,1) = -0.5
	 z(2,1) = sqrt(3.)/2.
	 z(3,1) = 0.
c  (0001)[-1 2 -1 0]
	 y(1,2) = 0.
	 y(2,2) = 0.
	 y(3,2) = 1.
	 z(1,2) = -0.5
	 z(2,2) = -sqrt(3.)/2.
	 z(3,2) = 0.
c  (0001)[2 -1 -1 0]
	 y(1,3) = 0.
	 y(2,3) = 0.
	 y(3,3) = 1.
       z(1,3) = 1.
	 z(2,3) = 0.
	 z(3,3) = 0.

c  Prismatic Slip Systems {1010}<1120>
c  (1 0 -1 0)[-1 2 -1 0]
	 y(1,4) = sqrt(3.)/2.
	 y(2,4) = 0.5
	 y(3,4) = 0.
       z(1,4) = -0.5
	 z(2,4) = sqrt(3.)/2.
	 z(3,4) = 0.
c  (-1 1 0 0)[-1 -1 2 0]
	 y(1,5) = -sqrt(3.)/2.
	 y(2,5) = 0.5
	 y(3,5) = 0.
       z(1,5) = -0.5
	 z(2,5) = -sqrt(3.)/2.
	 z(3,5) = 0.
c  (0 -1 1 0)[2 -1 -1 0]
	 y(1,6) = 0.
	 y(2,6) = -1.
	 y(3,6) = 0.
       z(1,6) = 1.
	 z(2,6) = 0.
	 z(3,6) = 0.

c  Pyramidal Slip Systems {10-11}<1120>
c  (1 0 -1 1)[-1 2 -1 0]
	 y(1,7) = sqrt(3.)/2.*c_a
	 y(2,7) = 0.5*c_a
	 y(3,7) = sqrt(3.)/2.
	 z(1,7) = -0.5
	 z(2,7) = sqrt(3.)/2.
	 z(3,7) = 0.
c  (-1 0 1 1)[-1 2 -1 0]  
	 y(1,8) = -sqrt(3.)/2.*c_a
	 y(2,8) = -0.5*c_a
	 y(3,8) = sqrt(3.)/2.
       z(1,8) = -0.5
	 z(2,8) = sqrt(3.)/2
	 z(3,8) = 0.
c  (-1 1 0 1)[-1 -1 2 0]  
	 y(1,9) = -sqrt(3.)/2.*c_a
	 y(2,9) = 0.5*c_a
	 y(3,9) = sqrt(3.)/2.
       z(1,9) = -0.5
	 z(2,9) = -sqrt(3.)/2.
	 z(3,9) = 0.
c  (1 -1 0 1)[-1 -1 2 0]  
	 y(1,10) = sqrt(3.)/2.*c_a
	 y(2,10) = -0.5*c_a
	 y(3,10) = sqrt(3.)/2.
       z(1,10) = -0.5
	 z(2,10) = -sqrt(3.)/2.
	 z(3,10) = 0.
c  (0 -1 1 1)[2 -1 -1 0] 
	 y(1,11) = 0.
	 y(2,11) = -c_a
	 y(3,11) = sqrt(3.)/2.
       z(1,11) = 1.
	 z(2,11) = 0.
	 z(3,11) = 0.
c  (0 1 -1 1)[2 -1 -1 0] 
 	 y(1,12) = 0. 
	 y(2,12) = c_a
	 y(3,12) = sqrt(3.)/2.
	 z(1,12) = 1.
	 z(2,12) = 0.
	 z(3,12) = 0.

      if (phase.eq.1) then  
c  Pyramidal <c+a> Slip Systems {1011}<1123>
c  (1 0 -1 1)[-2 1 1 3]  
	 y(1,13) = y(1,7)
	 y(2,13) = y(2,7)
	 y(3,13) = y(3,7)
	 z(1,13) = -z(1,6)
	 z(2,13) = -z(2,6)
	 z(3,13) = c_a - z(3,6)
c  (1 0 -1 1)[-1 -1 2 3] 
	 y(1,14) = y(1,13)
	 y(2,14) = y(2,13)
	 y(3,14) = y(3,13)
       z(1,14) = z(1,5)
	 z(2,14) = z(2,5)
	 z(3,14) = c_a + z(3,5)
c  (-1 1 0 1)[2 -1 -1 3]   
	 y(1,15) = y(1,8)
	 y(2,15) = y(2,8)
	 y(3,15) = y(3,8)
       z(1,15) = z(1,6)
	 z(2,15) = z(2,6)
	 z(3,15) = c_a + z(3,6)
c  (-1 1 0 1)[1 1 -2 3]     
	 y(1,16) = y(1,15)
	 y(2,16) = y(2,15)
	 y(3,16) = y(3,15)
       z(1,16) = -z(1,5)
	 z(2,16) = -z(2,5)
	 z(3,16) = c_a - z(3,5)
c  (-1 1 0 1)[2 -1 -1 3]  
	 y(1,17) = y(1,9)
	 y(2,17) = y(2,9)
	 y(3,17) = y(3,9)
	 z(1,17) = z(1,6)
	 z(2,17) = z(2,6)
	 z(3,17) = c_a + z(3,6)
c  (-1 1 0 1)[1 -2 1 3]   
	 y(1,18) = y(1,17)
	 y(2,18) = y(2,17)
	 y(3,18) = y(3,17)
       z(1,18) = -z(1,4)
	 z(2,18) = -z(2,4)
	 z(3,18) = c_a - z(3,4)
c  (1 -1 0 1)[-2 1 1 3] 
	 y(1,19) = y(1,10)
	 y(2,19) = y(2,10)
	 y(3,19) = y(3,10)
       z(1,19) = -z(1,6)
	 z(2,19) = -z(2,6)
	 z(3,19) = c_a - z(3,6)
c  (1 -1 0 1)[-1 2 -1 3] 
	 y(1,20) = y(1,19)
	 y(2,20) = y(2,19)
	 y(3,20) = y(3,19)
	 z(1,20) = z(1,4)
	 z(2,20) = z(2,4)
	 z(3,20) = c_a + z(3,4)
c  (0 -1 1 1)[-1 2 -1 3] 
	 y(1,21) = y(1,11)
	 y(2,21) = y(2,11)
	 y(3,21) = y(3,11)
       z(1,21) = z(1,4)
	 z(2,21) = z(2,4)
	 z(3,21) = c_a + z(3,4)
c  (0 -1 1 1)[1 1 -2 3]   
	 y(1,22) = y(1,21)
	 y(2,22) = y(2,21)
	 y(3,22) = y(3,21)
       z(1,22) = -z(1,5)
	 z(2,22) = -z(2,5)
	 z(3,22) = c_a - z(3,5)
c  (0 1 -1 1)[1 -2 1 3]  
	 y(1,23) = y(1,12)
	 y(2,23) = y(2,12)
	 y(3,23) = y(3,12)
       z(1,23) = -z(1,4)
	 z(2,23) = -z(2,4)
	 z(3,23) = c_a - z(3,4)
c  (0 1 -1 1)[-1 -1 2 3]   
	 y(1,24) = y(1,23)
	 y(2,24) = y(2,23)
	 y(3,24) = y(3,23)
	 z(1,24) = z(1,5)
	 z(2,24) = z(2,5)
	 z(3,24) = c_a + z(3,5)

	else

c  {110}<111> BCC Slip Systems
c  (101)[-1 -1 1]
	 y(1,13) = 1.
	 y(2,13) = 0.
	 y(3,13) = 1. 
	 z(1,13) = -1.
	 z(2,13) = -1.
	 z(3,13) = 1.
c  (101)[-1 1 1]
	 y(1,14) = y(1,13)
	 y(2,14) = y(2,13)
	 y(3,14) = y(3,13)
       z(1,14) = -1.
	 z(2,14) = 1.
	 z(3,14) = 1.
c  (-1 0 1)[1 -1 1]
	 y(1,15) = -1.
	 y(2,15) = 0.
	 y(3,15) = 1.
       z(1,15) = 1.
	 z(2,15) = -1.
	 z(3,15) = 1.
c  (-1 0 1)[1 1 1] 
	 y(1,16) = y(1,15)
	 y(2,16) = y(2,15)
	 y(3,16) = y(3,15)
       z(1,16) = 1.
	 z(2,16) = 1.
	 z(3,16) = 1.
c  (110)[-1 1 1]
	 y(1,17) = 1.
	 y(2,17) = 1.
	 y(3,17) = 0.
       z(1,17) = -1.
	 z(2,17) = 1.
	 z(3,17) = 1.
c  (110)[1 -1 1] 
	 y(1,18) = y(1,17)
	 y(2,18) = y(2,17)
	 y(3,18) = y(3,17)
	 z(1,18) = 1.
	 z(2,18) = -1.
	 z(3,18) = 1.
c  (-1 1 0)[-1 -1 1]
	 y(1,19) = -1.
	 y(2,19) = 1.
	 y(3,19) = 0.
	 z(1,19) = -1.
	 z(2,19) = -1.
	 z(3,19) = 1.
c  (-1 1 0)[1 1 1]
	 y(1,20) = y(1,19)
	 y(2,20) = y(2,19)
	 y(3,20) = y(3,19)
       z(1,20) = 1.
	 z(2,20) = 1.
	 z(3,20) = 1.
c  (011)[1 -1 1]
	 y(1,21) = 0.
	 y(2,21) = 1.
	 y(3,21) = 1.
       z(1,21) = 1.
	 z(2,21) = -1.
	 z(3,21) = 1.
c  (011)[-1 -1 1]
	 y(1,22) = y(1,21)
	 y(2,22) = y(2,21)
	 y(3,22) = y(3,21)
       z(1,22) = -1.
	 z(2,22) = -1.
	 z(3,22) = 1.
c  (0 1 -1)[1 1 1]
	 y(1,23) = 0.
	 y(2,23) = 1.
	 y(3,23) = -1.
       z(1,23) = 1.
	 z(2,23) = 1.
	 z(3,23) = 1.
c  (0 1 -1)[-1 1 1]
	 y(1,24) = y(1,23)
	 y(2,24) = y(2,23)
	 y(3,24) = y(3,23)
	 z(1,24) = -1.
	 z(2,24) = 1.
	 z(3,24) = 1.

c (112)<111>
c BCC slip systems not used in the current model.
c	 y(1,13) = 1.
c	 y(2,13) = 1.
c	 y(3,13) = 2. 
c	 z(1,13) = -1.
c	 z(2,13) = -1.
c	 z(3,13) = 1.
c  
c	 y(1,14) = 1.
c	 y(2,14) = 1.
c	 y(3,14) = -2.
c       z(1,14) = 1.
c	 z(2,14) = 1.
c	 z(3,14) = 1.
cc 
c	 y(1,15) = 1.
c	 y(2,15) = -1.
c	 y(3,15) = 2.
c       z(1,15) = -1.
c	 z(2,15) = 1.
c	 z(3,15) = 1.
cc   
c	 y(1,16) = -1.
c	 y(2,16) = 1.
c	 y(3,16) = 2.
c       z(1,16) = 1.
c	 z(2,16) = -1.
c	 z(3,16) = 1.
cc  
c	 y(1,17) = 1.
c	 y(2,17) = 2.
c	 y(3,17) = 1.
c       z(1,17) = 1.
c	 z(2,17) = -1.
c	 z(3,17) = 1.
cc   
c	 y(1,18) = 1.
c	 y(2,18) = 2.
c	 y(3,18) = -1.
c	 z(1,18) = -1.
c	 z(2,18) = 1.
c	 z(3,18) = 1.
cc
c	 y(1,19) = 1.
c	 y(2,19) = -2.
c	 y(3,19) = 1.
c	 z(1,19) = 1.
c	 z(2,19) = 1.
c	 z(3,19) = 1.
cc  
c	 y(1,20) = -1.
c	 y(2,20) = 2.
c	 y(3,20) = 1.
c       z(1,20) = -1.
c	 z(2,20) = -1.
c	 z(3,20) = 1.
cc 
c	 y(1,21) = 2.
c	 y(2,21) = 1.
c	 y(3,21) = 1.
c       z(1,21) = -1.
c	 z(2,21) = 1.
c	 z(3,21) = 1.
cc   
c	 y(1,22) = 2.
c	 y(2,22) = 1.
c	 y(3,22) = -1.
c       z(1,22) = 1.
c	 z(2,22) = -1.
c	 z(3,22) = 1.
cc  
c	 y(1,23) = 2.
c	 y(2,23) = -1.
c	 y(3,23) = 1.
c       z(1,23) = -1.
c	 z(2,23) = -1.
c	 z(3,23) = 1.
cc  
c	 y(1,24) = -2.
c	 y(2,24) = 1.
c	 y(3,24) = 1.
c	 z(1,24) = 1.
c	 z(2,24) = 1.
c	 z(3,24) = 1.
c
	end if
c----- --------------------------------------------------------------
c  Normalize the slip plane normals and slip directions to length one.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
         call normalize_vector( y(1,i), y(2,i), y(3,i) )
         call normalize_vector( z(1,i), z(2,i), z(3,i) )
      end do

c-------------------------------------------------------------------
c  Check for normality of slip plane normals and slip directions.
c-------------------------------------------------------------------

      do k = 1,num_slip_sys
         sum = 0.0
         do i = 1,3
            sum = sum + y(i,k) * z(i,k)
         end do
         if (abs(sum) .gt. tolerance) then
          print*,'The Miller indices are WRONG!!!'
          print*,'on slip system # ',k
          STOP
         end if
      end do

c-------------------------------------------------------------------
c  Initialize internal variables for initial time step
c-------------------------------------------------------------------

      if (time(2) .eq. 0.0) then

        psi(1) = psi_ang*pi/180
        psi(2) = theta_ang*pi/180
        psi(3) = phi_ang*pi/180


c--- state variables 1-3
c-------------------------------------------------------------------
c  Initialize drag stress for each slip system of
c  each grain.
c-------------------------------------------------------------------
        if (phase.eq.1) then

         do i = 1,3
          d0(i) = d_zero_p
         end do

         do i = 4,6
          d0(i) = d_zero_p
         end do
 
         do i = 7,12
          d0(i) = 400.
         end do
 
         do i = 13,num_slip_sys
          d0(i) = 650.
         end do

	  else

	   do i = 1,6
	    d0(i) = d_zero_p
	   end do
	   do i = 7,12
	    d0(i) = 400
	   end do
	   do i = 13,24
	    d0(i) = 0.9*d_zero_p
	   end do

	  end if
          
c--- state variables 4-27
c-------------------------------------------------------------------
c  Initialize back stress for each slip system of
c  each grain.
c-------------------------------------------------------------------

        do i = 1,num_slip_sys
         a0(i) = a_zero
        end do

            
c--- state variables 28-51
c-------------------------------------------------------------------
c  Initialize threshold stress for each slip system of
c  each grain.
c-------------------------------------------------------------------
        if (phase.eq.1) then

          do i = 1,num_slip_sys
           t0(i) = t_zero
          end do
           
        else

          do i = 1,12
           t0(i) = 450.
          end do
 
          do i = 13,24
           t0(i) = 750.
          end do

c  The easy glide systems will change depending on which variant of the
c  BOR is selected and must be reflected here.  These are the easy glide
c  systems for the first variant as indicated later in the code.         
          t0(1) = 85.
          t0(2) = 85.
          t0(3) = 85.
          t0(6) = 100.
          t0(13) = 100.
          t0(14) = 100.

        end if
  
c--- state variables 52-75
c-------------------------------------------------------------------
c  Initialize F_p_inv_0
c-------------------------------------------------------------------

        do i = 1,3
          do j = 1,3
           F_p_inv_0(i,j) = 0.0
          end do
          F_p_inv_0(i,i) = 1.0
        end do
            
c--- state variables 76-84
c-------------------------------------------------------------------
c  Initialize E_p
c-------------------------------------------------------------------

        do i = 1,3
	    do j = 1,3
           E_p(i,j) = 0.
          end do
	  end do

c--- state variables 85-93
c-------------------------------------------------------------------
c  Initialize E_eff
c-------------------------------------------------------------------

        E_eff = 0.

c--- state variables 94
c-------------------------------------------------------------------
c  Initialize E_p_eff
c-------------------------------------------------------------------

        E_p_eff = 0.

c--- state variables 95
c-------------------------------------------------------------------
c  Initialize E_p_eff_cum
c-------------------------------------------------------------------

        E_p_eff_cum = 0.
c--- state variables 96
c-------------------------------------------------------------------
c  Initialize gamma_p
c-------------------------------------------------------------------
c        do i = 1,6
c         gamma_p(i) = 0.
c	  end do
c--- state variables 97-102

        tempvar1 = 0.
        tempvar2 = 0.

c-------------------------------------------------------------------
c  End of initializations.  Read in internal variables.
c-------------------------------------------------------------------

      else  ! time<>0

        n = 0
            
c-------------------------------------------------------------------
c  Read in Euler Angles
c-------------------------------------------------------------------

        do i = 1,3
         n = n + 1
         psi(i) = statev(n)
        end do
          
c--state variables 1-3
c-------------------------------------------------------------------
c  Read drag stress values
c-------------------------------------------------------------------

        do i = 1,num_slip_sys
         n = n + 1
         d0(i) = statev(n)
        end do

c--state variables 4-27
c-------------------------------------------------------------------
c  Read kinematic stress values
c-------------------------------------------------------------------

        do i = 1,num_slip_sys
         n = n + 1
         a0(i) = statev(n)
        end do

c--state variables 28-51
c-------------------------------------------------------------------
c  Read threshold stress values
c-------------------------------------------------------------------

        do i = 1,num_slip_sys
         n = n + 1
         t0(i) = statev(n)
        end do

c--state variables 52-75
c-------------------------------------------------------------------
c  Read inverse of the plastic part of F
c-------------------------------------------------------------------

        do i = 1,3
         do j = 1,3
          n = n + 1
          F_p_inv_0(i,j) = statev(n)
         end do
        end do
          
c--state variables 76-84
c-------------------------------------------------------------------
c  Read E_p
c-------------------------------------------------------------------

        do i = 1,3
         do j = 1,3
          n = n + 1
          E_p(i,j) = statev(n)
         end do
        end do
          
c--state variables 85-93
c-------------------------------------------------------------------
c  Read E_eff
c-------------------------------------------------------------------
        n=n+1
	  E_eff = statev(n)

c--- state variables 94
c-------------------------------------------------------------------
c  Read E_p_eff
c-------------------------------------------------------------------
        n=n+1
	  E_p_eff = statev(n)

c--- state variables 95
c-------------------------------------------------------------------
c  Read E_p_eff_cum
c-------------------------------------------------------------------
        n = n+1
        E_p_eff_cum = statev(n)    
c--- state variables 96
c-------------------------------------------------------------------
c  Read gamma_p
c-------------------------------------------------------------------
c        do i = 1,6
c         n = n+1
c         gamma_p(i) = statev(n) 
c	  end do   
c--- state variables 96
       
        n=n+1
        tempvar1=statev(n)
        n=n+1
        tempvar2=statev(n)
c-------------------------------------------------------------------
c  End of initializations
c-------------------------------------------------------------------
         
      end if ! (time = 0)

c-------------------------------------------------------------------
c  Calculate orientation matrix between bcc and hcp based on BOR.  
c  Euler angles are given for alignment of the indicated BCC slip 
c  system with the (0001)[2 -1 -1 0] system in the alpha phase.
c  The Euler angles are initially given Bunge notation
c  (cf. Randle and Engler p. 26) and are later converted to the 
c  Bunge-Roe convention as used in the UMAT orienation matrix.
c-------------------------------------------------------------------
c  Euler Angles for variant (101)[-1 -1 1]
	psi_b(1,1) = pi/2. 
	psi_b(2,1) = pi/4.
	psi_b(3,1) = pi/2. + acos(sqrt(2./3.))

c  Euler Angles for variant (101)[-1 1 1]
	psi_b(1,2) = pi/2. 
	psi_b(2,2) = pi/4.
	psi_b(3,2) = pi/2. - acos(sqrt(2./3.))

c  Euler Angles for variant
	psi_b(1,3) = pi/2.
	psi_b(2,3) = 7.*pi/4.
	psi_b(3,3) = 3.*pi/2. - acos(sqrt(2./3.))

c  Euler Angles for variant
	psi_b(1,4) = pi/2.
	psi_b(2,4) = 7.*pi/4.
	psi_b(3,4) = 3.*pi/2. + acos(sqrt(2./3.))

c  Euler Angles for variant
	psi_b(1,5) = 3.*pi/4.
	psi_b(2,5) = pi/2.
	psi_b(3,5) = pi/2. - acos(sqrt(1./3.))

c  Euler Angles for variant
	psi_b(1,6) = 3.*pi/4
	psi_b(2,6) = pi/2.
	psi_b(3,6) = pi/2. + acos(sqrt(1./3.))

c  Euler Angles for variant
	psi_b(1,7) = pi/4.
	psi_b(2,7) = 3.*pi/2.
	psi_b(3,7) = 3.*pi/2. - acos(sqrt(1./3.))

c  Euler Angles for variant
	psi_b(1,8) = pi/4.
	psi_b(2,8) = 3.*pi/2.
	psi_b(3,8) = 3.*pi/2. + acos(sqrt(1./3.))

c  Euler Angles for variant
	psi_b(1,9) = 0.
	psi_b(2,9) = 7.*pi/4.
	psi_b(3,9) = 3.*pi/2. + acos(sqrt(2./3.))

c  Euler Angles for variant
	psi_b(1,10) = 0.
	psi_b(2,10) = 7.*pi/4.
	psi_b(3,10) = 3.*pi/2. - acos(sqrt(2./3.))

c  Euler Angles for variant
	psi_b(1,11) = 0.
	psi_b(2,11) = 5.*pi/4.
	psi_b(3,11) = 3.*pi/2. + acos(sqrt(2./3.))

c  Euler Angles for variant
	psi_b(1,12) = 0.
	psi_b(2,12) = 5.*pi/4.
	psi_b(3,12) = 3.*pi/2. - acos(sqrt(2./3.))

c  Convert to Bunge-Roe convention used by bob
	do i = 1,12
	  psi_b(1,i) = psi_b(1,i) + 3.*pi/2. 
	  psi_b(2,i) = psi_b(2,i)
	  psi_b(3,i) = psi_b(3,i) + pi/2.
	end do

c  Initialize the orientation matrices related to the BOR
      do i = 1,12
        s1 = sin(psi_b(1,i))
        c1 = cos(psi_b(1,i))
        s2 = sin(psi_b(2,i))
        c2 = cos(psi_b(2,i))
        s3 = sin(psi_b(3,i))
        c3 = cos(psi_b(3,i))
            
        bor(i,1,1) = c1*c2*c3-s1*s3
        bor(i,1,2) = c3*c2*s1+s3*c1
        bor(i,1,3) = -c3*s2
        bor(i,2,1) = -s3*c2*c1-c3*s1
        bor(i,2,2) = -s3*c2*s1+c3*c1
        bor(i,2,3) = s3*s2
        bor(i,3,1) = s2*c1
        bor(i,3,2) = s2*s1
        bor(i,3,3) = c2
      end do

c-------------------------------------------------------------------
c  Calculate direction cosines based on Euler angles
c-------------------------------------------------------------------

        s1 = sin(psi(1))
        c1 = cos(psi(1))
        s2 = sin(psi(2))
        c2 = cos(psi(2))
        s3 = sin(psi(3))
        c3 = cos(psi(3))
            
        dir_cos(1,1) = c1*c2*c3-s1*s3
        dir_cos(2,1) = c3*c2*s1+s3*c1
        dir_cos(3,1) = -c3*s2
        dir_cos(1,2) = -s3*c2*c1-c3*s1
        dir_cos(2,2) = -s3*c2*s1+c3*c1
        dir_cos(3,2) = s3*s2
        dir_cos(1,3) = s2*c1
        dir_cos(2,3) = s2*s1
        dir_cos(3,3) = c2

c-------------------------------------------------------------------
c  Initialize ANISOTROPIC elastic stiffness tensor in Voight notation
c-------------------------------------------------------------------

      do i = 1,6
	  do j = 1,6
	    C_IV(i,j) = 0.0
	  end do
	end do

	C_IV(1,1) = C11
	C_IV(1,2) = C12
	C_IV(1,3) = C13
	C_IV(2,1) = C12
	C_IV(2,2) = C11
	C_IV(2,3) = C13
	C_IV(3,1) = C13
	C_IV(3,2) = C13
	C_IV(3,3) = C33
	C_IV(4,4) = C44
	C_IV(5,5) = C44
	C_IV(6,6) = 0.5 * (C11-C12)
      

c-------------------------------------------------------------------
c  Convert initial elastic stiffness tensor from Voigt notation to
c  the standard 4th rank stiffness tensor
c-------------------------------------------------------------------
 
      do i = 1,3
       do j = 1,3
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = 1,3
          ib = k
          if (k.ne.l) ib=9-k-l
          C0(i,j,k,l) = C_IV(ia,ib)
         end do
        end do
       end do
      end do

c--------------------------------------------------------------------
c  Rotate local isotropic elasticity tensor to
c  global coordinates.
c--------------------------------------------------------------------

      call rotate_4th(dir_cos(1,1),C0,C)

c--------------------------------------------------------------------
c  Rotate slip plane normals and slip directions to global coordinates.
c--------------------------------------------------------------------

      do n = 1,12
         do i = 1,3
            xs0(i,n) = 0.0
            xm0(i,n) = 0.0
            do j = 1,3
               xs0(i,n) = xs0(i,n) + dir_cos(i,j) * z(j,n)
               xm0(i,n) = xm0(i,n) + dir_cos(i,j) * y(j,n)
            end do
         end do
      end do
      
	if (phase.eq.1) then

       do n = 13,num_slip_sys
         do i = 1,3
            xs0(i,n) = 0.0
            xm0(i,n) = 0.0
            do j = 1,3
               xs0(i,n) = xs0(i,n) + dir_cos(i,j) * z(j,n)
               xm0(i,n) = xm0(i,n) + dir_cos(i,j) * y(j,n)
            end do
         end do
       end do

      else

c  Additional rotation required for the BCC systems as they must first be 
c  rotated into the local coordinate system attached to the hcp crystal 
c  prior to specification of "grain" orientation.

       do n = 13,num_slip_sys
	  do i = 1,3
         xs0(i,n) = 0.0
         xm0(i,n) = 0.0
         do j = 1,3
          do k = 1,3
		 xs0(i,n) = xs0(i,n) + dir_cos(i,j)*bor(1,j,k)*z(k,n)
		 xm0(i,n) = xm0(i,n) + dir_cos(i,j)*bor(1,j,k)*y(k,n)
	    end do
	   end do
	  end do
       end do

      end if
      
c  Junk used for debuggin purposes
c      if (kcnt.lt.24) then
c      do i = 1,24
c        print*,kcnt,time(2)
c        print 134,i,xs0(1,i),xs0(2,i),xs0(3,i)
c        print 134,i,xm0(1,i),xm0(2,i),xm0(3,i)
c        kcnt = kcnt + 1
c      end do
c      end if
c

134   format(i2,2x,f9.4,2x,f9.4,2x,f9.4)

c--------------------------------------------------------------------
c  Initialize number of subincrements.  Note that the
c  subincrement initially equals the total increment.  This
c  remains the case unless the process starts to diverge.
c--------------------------------------------------------------------

      N_incr       = 1
      N_incr_total = 1

c====================================================================
c  Top of Subincrement Time Step Integration Loop.
c  Calculate subincrement time step size.
c====================================================================

  100 dt_incr = dtime / N_incr_total

c-------------------------------------------------------------------
c  Initialize drag stress (isotropic hardening)
c-------------------------------------------------------------------

      do n = 1,num_slip_sys
        d(n) = d0(n)
      end do

c-------------------------------------------------------------------
c  Initialize back stress (kinematic hardening)
c-------------------------------------------------------------------

      do n = 1,num_slip_sys
        a(n) = a0(n)
      end do
c-------------------------------------------------------------------
c  Initialize threshold stress 
c-------------------------------------------------------------------

      do n = 1,num_slip_sys
        t(n) = t0(n)
      end do


c--------------------------------------------------------------------
c  Initialize deformation gradients for beginning and
c  end of subincrement.
c--------------------------------------------------------------------

      do i = 1,3
         do j = 1,3
            F0(i,j) = dfgrd0(i,j) + (dfgrd1(i,j) - dfgrd0(i,j)) *
     &	  (N_incr - 1) / N_incr_total
            F1(i,j) = dfgrd0(i,j) + (dfgrd1(i,j) - dfgrd0(i,j)) *
     &	  N_incr / N_incr_total
         end do
      end do

c--------------------------------------------------------------------
c  Multiply F() by F_p_inv_0() to get F_el() and take the inverse of 
c  F_el() to get F_el_inv
c--------------------------------------------------------------------

      call aa_dot_bb(3,F0,F_p_inv_0,F_el)            
c      call inverse_3x3(F_el,F_el_inv)

c--------------------------------------------------------------------
c  Rotate xs0 and xm0 to current coordinates, called xs and xm.
c--------------------------------------------------------------------

c      do n = 1,num_slip_sys
c         do i = 1,3
c            xs(i,n) = 0.0
c            xm(i,n) = 0.0
c            do j = 1,3
c               xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
c               xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
c            end do
c         end do
c      end do

c--------------------------------------------------------------------
c  Calculate elastic Green Strain E=(1/2)*(C-I)
c--------------------------------------------------------------------

      call transpose(3,F_el,array1)
      call aa_dot_bb(3,array1,F_el,E_el)
      do i = 1,3
         E_el(i,i) = E_el(i,i) - 1
         do j = 1,3
            E_el(i,j) = E_el(i,j) / 2
         end do
      end do 

c--------------------------------------------------------------------
c  Multiply the anisotropic stiffness tensor by the Green strain 
c  to get the 2nd Piola Kirkhhoff stress
c--------------------------------------------------------------------

      call aaaa_dot_dot_bb(3,C,E_el,Spk2)

c--------------------------------------------------------------------
c  Convert from Spk2 stress to Cauchy stress
c--------------------------------------------------------------------

c      det = determinant(F_el)
c      call transpose(3,F_el,array2)
c      call aa_dot_bb(3,F_el,Spk2,array1)
c      call aa_dot_bb(3,array1,array2,sig)
      
c	if (det.eq.0) then
c      do i = 1,3
c         do j = 1,3
c            sig(i,j) =0.
c         end do
c      end do
c      endif

c      if (det.ne.0) then
c      do i = 1,3
c         do j = 1,3
c            sig(i,j) = sig(i,j) / det
c         end do
c      end do
c      endif

c--------------------------------------------------------------------
c  Calculate resolved shear stress for each slip system.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
         tau(k) = 0.0
         do j = 1,3
            do i = 1,3
	        tau(k) = tau(k) + xs0(i,k) * xm0(j,k) * Spk2(i,j)
            end do
         end do
      end do

c--------------------------------------------------------------------
c   Modify threshold stress to include non-Schmid terms
c--------------------------------------------------------------------
      
	do i = 1,num_slip_sys
	   t1(i) = t(i)
	end do

	t1(4) = t(4) + a1*tau(7) - a1*tau(8)
	t1(5) = t(5) + a1*tau(9) - a1*tau(10)
	t1(6) = t(6) + a1*tau(11) - a1*tau(12)

c        print 15, t1(4),t1(5),t1(6)
c15      format(f10.4,2x,f10.4,2x,f10.4)
c        write(7,'(f10.4,2x,f10.4,2x,f10.4)'),t1(4),t1(5),t1(6)
c--------------------------------------------------------------------
c  Calculate 1st estimate of gamma_dot for each slip system.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
	  taud = abs(tau(k)-a(k))-t1(k)
	  if (taud.gt.0.) then
	    gamma_dot(k) = gamma_dot_zero*power((taud/d(k)),flow_exp)*
     &    sgn(tau(k)-a(k))
	  else 
	    gamma_dot(k) = 0.
	  end if
	end do

c	write(7,*) 'gamma_dot1',gamma_dot(7)
c--------------------------------------------------------------------
c  Calculate d(Tau)/d(Gamma_dot)
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
       do ib = 1,num_slip_sys
        dtaudgd(ia,ib) = 0.0

        do i = 1,3
         do j = 1,3
          array1(i,j) = xs0(i,ib) * xm0(j,ib)
         end do
        end do
        call aaaa_dot_dot_bb(3,C,array1,array2)

        do i = 1,3
         do j = 1,3
          array1(i,j) = xs0(i,ia) * xm0(j,ia)
         end do
        end do
        call aa_dot_dot_bb(3,array1,array2,dtaudgd(ia,ib))

        dtaudgd(ia,ib) = dtaudgd(ia,ib) * dt_incr * (-1.)
       end do !ib
      end do !ia

c====================================================================
c  Begin Newton-Raphson iterative loops.
c====================================================================

      converged = .false.

      do while (.not.converged)

      converged = .true.

c--------------------------------------------------------------------
c  Calculate g, F_p_inv, F_el, sig, tau, func, sse.
c--------------------------------------------------------------------
c        print 15, 'before first eval_func',noel,time(2)
c15      format(a,2x,i7,2x,f7.4)

      call eval_func(      xs0,	dt_incr,	gamma_dot,	
     &			    xm0,	F1,		num_slip_sys,
     &			    C,		F_p_inv,	F_p_inv_0,
     &			    d0,		Hdir,		Hdyn,
     &			    a0,		Adir,		Adyn,
     &			    g_sat,	F_el,		flow_exp,
     &			    sig,	tau,		gamma_dot_zero,
     &			    d,		func,		xL_p,
     &			    xLatent,	sse,	a, t, t0, t1, a1,
     &                Spk2,gamma,text)

      sse_ref = sse

c--------------------------------------------------------------------
c  Begin calculation of the partial derivatives needed for 
c  the Newton-Raphson step!!!
c--------------------------------------------------------------------
c  Calculate derivative of the hardening variable, d-alpha,
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------

      sum = 0
      do ia = 1,num_slip_sys      
        sum = sum + abs(gamma_dot(ia))
      end do
      do ia = 1,num_slip_sys
       do ib = 1,num_slip_sys
         temp = xLatent
         if (ia .eq. ib) temp = 1.0
         ddadgb(ia,ib) = (Hdir * temp - Hdyn*d(ia)) * dt_incr / 
     &	(1 + Hdyn * sum * dt_incr)
         if(gamma_dot(ib).lt.0.0)  ddadgb(ia,ib)=-ddadgb(ia,ib)
       end do
      end do        !If g=const > ddadgb = 0
      
c	do ia = 1,num_slip_sys
c	  do ib = 1,num_slip_sys
c	    ddadgb(ia,ib) = 0.
c	  end do
c	end do
c--------------------------------------------------------------------
c  Calculate derivative of kinematic stress, a-alpha
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
         daadga(ia) = (Adir - Adyn * a(ia) * sgn(gamma_dot(ia)))
     &   * dt_incr / (1 + Adyn * dt_incr * abs(gamma_dot(ia)))
      end do

c--------------------------------------------------------------------
c  Calculate derivative of threshold stress, t-alpha
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------

c      do ia = 1,num_slip_sys
c	  do ib = 1,num_slip_sys
c         dtadgb(ia,ib) = 0.
c	  end do
c      end do

      do ia = 1,num_slip_sys
	  do ib = 1,num_slip_sys
         dt1adgb(ia,ib) = 0.
	  end do
      end do

c--------------------------------------------------------------------
c  Form "A-matrix" of derivatives wrt d_gamma_beta.  
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
	  do ib = 1,num_slip_sys
	    array3(ia,ib) = 0.
	  end do
	end do
 
      do ia = 1,num_slip_sys
	  do ib = 1,num_slip_sys
	    taud = abs(tau(ia)-a(ia))-t1(ia)
	    if (taud.gt.0) then
	      if (ia.eq.ib) then
	        dtaudadgb = (dtaudgd(ia,ib)-daadga(ia))*sgn(tau(ia)-a(ia))
     &        - dt1adgb(ia,ib)
	      else
	        dtaudadgb = dtaudgd(ia,ib)*sgn(tau(ia)-a(ia))
     &        - dt1adgb(ia,ib)
            end if       
	      array3(ia,ib) = flow_exp*power((taud/d(ia)),flow_exp-1.)*
     &      (dtaudadgb/d(ia)-taud*ddadgb(ia,ib)/(d(ia)**2))
	      if (ia.eq.ib) then
	        array3(ia,ib) = array3(ia,ib) - sgn(tau(ia)-a(ia))/
     &        gamma_dot_zero
	      end if
	    end if	  
	  end do
	end do
	      
c--------------------------------------------------------------------
c  Calculate the gradient of sse wrt gamma_dot().  Will be used
c  later to ensure that line search is in the correct direction.
c--------------------------------------------------------------------

      do j = 1,num_slip_sys
         grad(j) = 0.0
         do i = 1,num_slip_sys
            grad(j) = grad(j) + func(i) * array3(i,j)
         end do
         grad(j) = 2. * grad(j)
      end do

c--------------------------------------------------------------------
c  Solve for increment of gamma_dot.  Solution is returned 
c  in the func() array.
c--------------------------------------------------------------------

c      call simeq(num_slip_sys,array3,func)

	    do ia=1,num_slip_sys
	      do ib=1,num_slip_sys
	      dummy(ia,ib) = array3(ia,ib)
	      end do
	    end do

	    do ia=1,num_slip_sys
	      taud = abs(tau(ia)-a(ia)) - t1(ia)
	      if (taud/d(ia).le.0.0) dummy(ia,ia)=1.0
	    end do

          call simeq(num_slip_sys, dummy, func)
c--------------------------------------------------------------------
c  Store offset in d_gamma_dot(k) 
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
         d_gamma_dot(k) = func(k)
      end do

c--------------------------------------------------------------------
c  Check to make sure that N-R step leads 'down hill' the 
c  sse surface.
c--------------------------------------------------------------------

      sum = 0.0
      do i = 1,num_slip_sys
         sum = sum - grad(i) * d_gamma_dot(i)
      end do

      if (sum .gt. 0.0) then
        do i = 1,num_slip_sys
           d_gamma_dot(i) = -d_gamma_dot(i)
        end do
      end if

c--------------------------------------------------------------------
c  Multiply step size by two because next loop will divide it by 2.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
         d_gamma_dot(k) = d_gamma_dot(k) * 2
      end do
      
c====================================================================
c  Begin line search.
c====================================================================

      improved = .false.

      do N_ctr = 1,max_loops

      sse_old = sse

c--------------------------------------------------------------------
c  Divide step size by two.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
         d_gamma_dot(k) = d_gamma_dot(k) / 2
         gamma_try(k)   = gamma_dot(k) + d_gamma_dot(k)
      end do
      
c--------------------------------------------------------------------
c  Calculate g, F_p_inv, F_el, sig, tau, func, and sse based
c  on gamma_try(k)
c--------------------------------------------------------------------
c        print 15, 'before second eval_func',noel,time(2)

      call eval_func(      xs0,	dt_incr,	gamma_try,	
     &			    xm0,	F1,		num_slip_sys,
     &			    C,		F_p_inv,	F_p_inv_0,
     &			    d0,		Hdir,		Hdyn,
     &			    a0,		Adir,		Adyn,
     &			    g_sat,	F_el,		flow_exp,
     &			    sig,	tau,		gamma_dot_zero,
     &			    d,		func,		xL_p,
     &			    xLatent,	sse,	a, t, t0, t1, a1,
     &                Spk2,gamma,text)

c        print 15, 'after second eval_func',noel,time(2)
c--------------------------------------------------------------------
c  Check for 'convergence' of the line search.  Note that the line
c  search doesn't waste time converging closely.  This is because
c  the N-R step does so much better.
c--------------------------------------------------------------------

      if ((sse_old.le.sse_ref).and.(sse.ge.sse_old).and.
     &	(N_ctr.gt.1)) improved=.true.

      if (improved) go to 200

      end do ! Line Search

c--------------------------------------------------------------------
c  Add "d_gamma_dot" to gamma_dot to get new values for
c  this iteration. 
c--------------------------------------------------------------------

  200 do k = 1,num_slip_sys
        gamma_dot(k) = gamma_dot(k) + d_gamma_dot(k) * 2.0
      end do

c	write(7,*) 'gamma_dot2',gamma_dot(7)

c--------------------------------------------------------------------
c  If (sse_old > tolerance) then this step has not converged.
c--------------------------------------------------------------------

      if (sse_old .gt. tolerance) converged = .false.

c--------------------------------------------------------------------
c  If (sse_old > sse_ref/2) then convergence is too slow and
c  increment is divided into two subincrements.
c--------------------------------------------------------------------
c      write(7,*) sse_old,sse_ref

      if ((sse_old.gt.sse_ref/2.0).and.(.not.converged)) then
        N_incr = 2 * N_incr - 1
        N_incr_total = 2 * N_incr_total
        go to 100
      end if

c--------------------------------------------------------------------
c  End iterative loop.
c--------------------------------------------------------------------

      end do ! 'Newton Raphson Iterative Loop'

c--------------------------------------------------------------------
c  If another subincrement remains to be done, then update counters.
c--------------------------------------------------------------------

      if (N_incr .lt. N_incr_total) then
        if (N_incr .eq. (N_incr/2)*2) then	! N_incr is 'even'
          N_incr = N_incr / 2 + 1
          N_incr_total = N_incr_total / 2
        else					! N_incr is 'odd'
          N_incr = N_incr + 1
        end if

c--------------------------------------------------------------------
c  Update F_p_inv_0 for next time sub-step.
c--------------------------------------------------------------------

        do i = 1,3
          do j = 1,3
            F_p_inv_0(i,j) = F_p_inv(i,j)
          end do
        end do

c--------------------------------------------------------------------
c  Update g() and a() for next time sub-step.
c--------------------------------------------------------------------

        do i = 1,num_slip_sys
          d0(i) = d(i)
          a0(i) = a(i)
          t0(i) = t(i)
        end do

        go to 100
      end if

c--------------------------------------------------------------------
c  Output stress, strain, etc. for post processing.
c--------------------------------------------------------------------
c      call aa_dot_bb(3,F_el,dir_cos,array1)
c      call kocks_angles(npt,time(2),array1)

c--------------------------------------------------------------------
c  Convert from Spk2 stress to Cauchy stress
c--------------------------------------------------------------------

      det = determinant(F_el)
      call transpose(3,F_el,array2)
      call aa_dot_bb(3,F_el,Spk2,array1)
      call aa_dot_bb(3,array1,array2,sig)
      
	if (det.eq.0) then
       do i = 1,3
        do j = 1,3
         sig(i,j) =0.
        end do
       end do
      endif

      if (det.ne.0) then
       do i = 1,3
        do j = 1,3
         sig(i,j) = sig(i,j) / det
        end do
       end do
      endif
c--------------------------------------------------------------------
c  Calculate Green Strain
c--------------------------------------------------------------------
      call transpose(3,F1,array1)
	call aa_dot_bb(3,array1,F1,E_tot)

      do i = 1,3
        E_tot(i,i) = E_tot(i,i) - 1
        do j = 1,3
          E_tot(i,j) = E_tot(i,j) / 2
        end do 
      end do 
c====================================================================
c  Begin calculation of the Jacobian (the tangent stiffness matrix).
c===================================================================
c--------------------------------------------------------------------
c  Rotate xs0 and xm0 to current coordinates, called xs and xm.
c--------------------------------------------------------------------

      call inverse_3x3(F_el,F_el_inv)
      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate the derivative of the plastic part of the rate of
c  deformation tensor in the current configuration wrt sigma.
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
           ddpdsig(i,j,k,l) = 0.0
         end do ! l
        end do ! k
       end do ! j
      end do ! i

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          do n = 1,num_slip_sys
	      taud = abs(tau(n)-a(n)) - t1(n)
	      if (taud.gt.0.) then
             ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l) + (xs(i,n)*xm(j,n) +
     &       xm(i,n) * xs(j,n)) * (xs(k,n) * xm(l,n) + xm(k,n) *
     &       xs(l,n)) * power((taud/d(n)),
     &	   flow_exp-1.0) / d(n)
	      end if
          end do ! num_slip_sys
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Calculate the inverse of F_el
c--------------------------------------------------------------------

c      call inverse_3x3(F_el,F_el_inv)

c--------------------------------------------------------------------
c  Scale by appropriate constants and divide by num_grains to
c  get the average value.  ALSO multiply by 'dtime' which is
c  d(sig)/d(sig_dot).
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l) * dtime * flow_exp *
     &	gamma_dot_zero / 4. / num_grains
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Multiply the 4th rank elastic stiffness tensor by the derivative
c  of the plastic part of the rate of deformation tensor wrt sig_dot.
c--------------------------------------------------------------------

      call aaaa_dot_dot_bbbb(3,C,ddpdsig,array6)

c--------------------------------------------------------------------
c  Add 4th rank identity tensor to array6()
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          array6(i,j,k,l) = array6(i,j,k,l) + 0.5 * 
     &		(del(i,k) * del(j,l) + del(i,l) * del(j,k))
         end do
        end do
       end do
      end do

c--------------------------------------------------------------------
c  Need to take the inverse of Array4.  Since it relates two 2nd
c  rank tensors that are both symmetric, Array4 can be xformed to 
c  Voigt notation style in order to do the inverse, then xformed back.
c--------------------------------------------------------------------

      call forth_to_Voigt (Array6,Array4)
      call inverse(6,Array4,Array5)
      call Voigt_to_forth (Array5,Array6)

c--------------------------------------------------------------------
c  Multiply Array6 by C, the elastic stiffness matrix to
c  finally get the Jacobian, but as a 4th rank tensor.
c--------------------------------------------------------------------

      call aaaa_dot_dot_bbbb (3,Array6,C,ddsdde_4th)

c--------------------------------------------------------------------
c  Store the stress tensor in the ABAQUS stress 'vector'
c--------------------------------------------------------------------

      do i = 1,ndi
         stress(i) = sig(i,i)
      end do
      if (nshr .eq. 1) stress(ndi+1) = sig(1,2)
      if (nshr .eq. 3) then
         stress(4) = sig(1,2)
         stress(5) = sig(1,3)
         stress(6) = sig(2,3)
      end if

c--------------------------------------------------------------------
c  Store the Jacobian in Voigt notation form.
c--------------------------------------------------------------------

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=i+j+1
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=k+l+1
          array4(ia,ib) = ddsdde_4th(i,j,k,l)
c          IF(IB.GE.4) ARRAY4(IA,IB) = 2 * ARRAY4(IA,IB)
         end do
        end do
       end do
      end do

      do i =1,6
        do j = 1,6
          ddsdde(i,j) = 0.0
        end do
      end do

      if (ndi .eq. 1) then			! 1-D
         ddsdde(1,1) = array4(1,1)
      else if (ndi .eq. 2) then			! 2-D plane stress & axi
         do i = 1,2
            do j = 1,2
               ddsdde(i,j) = array4(i,j)
            end do
         end do
         ddsdde(1,3) = array4(1,4)
         ddsdde(2,3) = array4(2,4)
         ddsdde(3,1) = array4(4,1)
         ddsdde(3,2) = array4(4,2)
         ddsdde(3,3) = array4(4,4)
      else if (ndi .eq. 3 .and. nshr .eq. 1) then ! plane strain
         do i = 1,4
            do j = 1,4
               ddsdde(i,j) = array4(i,j)
            end do 
           end do 
         else					! Full 3-D
         do i = 1,6
            do j = 1,6
               ddsdde(i,j) = array4(i,j)
            end do
         end do
      end if

     
c-------------------------------------------------------------------
c  Store the internal variables in the statev() array
c-------------------------------------------------------------------

      n = 0
	               
c-------------------------------------------------------------------
c  Store the Euler Angles
c-------------------------------------------------------------------

      do i = 1,3
        n = n + 1
        statev(n) = psi(i)
      end do
            
c   state variables 1-3
c-------------------------------------------------------------------
c  Store the drag stresses.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
       n = n + 1
       statev(n) = d(i)
      end do

c   state variables 4-27
c-------------------------------------------------------------------
c  Store the back stresses.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
       n = n + 1
       statev(n) = a(i)
      end do

c   state variables 28-51
c-------------------------------------------------------------------
c  Store the threshold stresses.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
       n = n + 1
       statev(n) = t(i)
      end do

c   state variables 52-75
c-------------------------------------------------------------------
c  Store F_p_inv
c-------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        n=n+1
        statev(n)= F_p_inv(i,j)
       end do
      end do

c   state variables 76-84
c-------------------------------------------------------------------
c  Plastic Strain Calculations
c-------------------------------------------------------------------
c  Increment of plastic shear strain accumulated on each slip system
c  over this time step.
      do i = 1,num_slip_sys
	  delta_gamma(i) = gamma_dot(i)*dtime
	end do

	do j = 1,3
	  do k = 1,3
         delta_E_p(j,k) =  0.0
	  end do
	end do

c  Increment of the plastic strain tensor
	do j = 1,3
	  do k = 1,3
	    do l = 1,num_slip_sys
	      delta_E_p(j,k) = delta_E_p(j,k) + 0.5*delta_gamma(l)*
     &     (xs0(j,l)*xm0(k,l)+xs0(k,l)*xm0(j,l))
          end do
	  end do
	end do 

c  Plastic strain tensor
      do i = 1,3
	  do j = 1,3
	    E_p(i,j) = E_p(i,j)+delta_E_p(i,j)
	  end do
	end do 
    
	do i = 1,3
	 do j = 1,3
	   n = n+1
	   statev(n) = E_p(i,j)
	 end do
	end do		      
c  state variables 85-93
 
c  Effective total strain
      call aa_dot_dot_bb(3,E_tot,E_tot,sum)
	E_eff = sqrt(2./3.*sum)
      
	n = n+1
      statev(n)= E_eff
c   state variables 94

c  Effective Plastic Strain
      call aa_dot_dot_bb(3,E_p,E_p,sum1)
	E_p_eff = sqrt(2./3.*sum1)

      if (E_p_eff.lt.0) E_p_eff = 0.
		     
	n = n+1
      statev(n)= E_p_eff
c   state variables 95

c  Cumulative Effective plastic strain 
      sum2 = 0.
      call aa_dot_dot_bb(3,delta_E_p,delta_E_p,sum2)
      sum2=sqrt(2./3.*sum2)
      E_p_eff_cum = sum2 + E_p_eff_cum

      n = n+1       
      statev(n) = E_p_eff_cum
c   state variable 96

c  Maximum plastic shear strain calculated for every 15 degrees
c      do i= 1,6
c	 theta = 15.*pi/180.*i
c	 xt1 = cos(theta)
c	 xt2 = sin(theta)
c	 xn1 = -sin(theta)
c	 xn2 = cos(theta)
c	 gamma_p(i) = xn1*E_p(1,1)*xt1+xn1*E_p(1,2)*xt2+xn2*E_p(2,1)*xt1+
c     & xn2*E_p(2,2)*xt2
c	 gamma_p(i) = gamma_p(i)/2.
c      end do
  
c      do i = 1,6
c	 n=n+1 
c	 statev(n) = gamma_p(i)
c	end do

      n=n+1
      statev(n) = xrot
      n=n+1
      statev(n) = F_p(2,2)



c  Steps for passing items into uexternaldb() for output to text files
      if (kincp.eq.0) icount = 1 
c         print 117,icount,noel,kinc,kincp,statev(85),statev(86),
c     &   statev(89),statev(93)
c         print 117,icount,noel,kinc,kincp,xcord(icount,1),
c     &   xcord(icount,2),xcord(icount,3),time(2)
      jnoel(icount) = noel
      xcord(icount,1) = coords(1)
      xcord(icount,2) = coords(2)
      xcord(icount,3) = coords(3)
      storsdv(icount,1) = statev(85)
      storsdv(icount,2) = statev(86)
      storsdv(icount,3) = statev(89)
      storsdv(icount,4) = statev(93)
c      print 121,noel,icount,kinc,time(2),dtime
c      print*,'*******************************************'
      kincp = kinc
      icount = icount+1
      if (icount.eq.(noel1+1)) icount = 1


117   format(i2,2x,i2,2x,i2,2x,i2,2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)
121   format(i4,2x,i4,2x,i2,2x,f12.6,2x,f12.6)
c-------------------------------------------------------------------

      return
      end
c====================================================================
c====================================================================
c====================== S U B R O U T I N E S =======================
c====================================================================
c====================================================================
c
c  Calculate a vector cross product.
c
c  c = a cross b
c
c--------------------------------------------------------------------

      subroutine kcross_product(a1,a2,a3,b1,b2,b3,c1,c2,c3)
      
      include 'ABA_PARAM.INC'
      
      c1 = a2 * b3 - a3 * b2
      c2 = a3 * b1 - a1 * b3
      c3 = a1 * b2 - a2 * b1

      return
      end

c====================================================================
c====================================================================
c
c  Normalize the length of a vector to one.
c
c--------------------------------------------------------------------

      subroutine normalize_vector(x,y,z)

      include 'ABA_PARAM.INC'

      xlength = sqrt(x*x+y*y+z*z)
      x = x / xlength
      y = y / xlength
      z = z / xlength

      return
      end

c====================================================================
c====================================================================
c
c  Transpose an ( n x n ) tensor.
c
c--------------------------------------------------------------------

      subroutine transpose(n,a,b)

      include 'ABA_PARAM.INC'

      dimension a(n,n), b(n,n)

      do i = 1,n
         do j = 1,n
            b(i,j) = a(j,i)
         end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Calculate the dot product of two 2nd rank tensors.
c  Result is stored in cc(i,j)
c
c--------------------------------------------------------------------

      subroutine aa_dot_bb(n,a,b,c)

      include 'ABA_PARAM.INC'

      dimension a(n,n), b(n,n), c(n,n)

      do i = 1,n
         do j = 1,n
            c(i,j) = 0
            do k = 1,n
               c(i,j) = c(i,j) + a(i,k) * b(k,j)
            end do
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 2nd rank tensors.
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bb(n,a,b,sum)

      include 'ABA_PARAM.INC'

      dimension a(n,n), b(n,n)

      sum = 0.0
      do i = 1,n
         do j = 1,n
            sum = sum + a(i,j) * b(i,j)
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 4th rank tensors.
c  Result is stored in c(i,j,k,l)
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bbbb(n,a,b,c)

      include 'ABA_PARAM.INC'

      dimension a(n,n,n,n), b(n,n,n,n), c(n,n,n,n)

      do i = 1,n
       do j = 1,n
        do k = 1,n
         do l = 1,n
          c(i,j,k,l) = 0
          do m1 = 1,n
           do m2 = 1,n
            c(i,j,k,l) = c(i,j,k,l) + a(i,j,m1,m2) * b(m1,m2,k,l)
           end do !m2
          end do !m1
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of a 4th rank tensor and
c  a 2nd rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bb(n,a,b,c)

      include 'ABA_PARAM.INC'

      dimension a(n,n,n,n), b(n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(i,j,k,l) * b(k,l)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of a 2nd rank tensor and
c  a 4th rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bbbb(n,a,b,c)

      include 'ABA_PARAM.INC'

      dimension a(n,n), b(n,n,n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(k,l) * b(k,l,i,j)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Rotates any 3x3x3x3 tensor by a rotation matrix.
c
c  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
c
c--------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

      do m = 1,3
       do n = 1,3
        do k = 1,3
         do l = 1,3
          d(m,n,k,l) = a(k,1) * (a(l,1) * b(m,n,1,1) + 
     &		a(l,2) * b(m,n,1,2) + a(l,3) * b(m,n,1,3)) +
     &		a(k,2) * (a(l,1) * b(m,n,2,1) + 
     &		a(l,2) * b(m,n,2,2) + a(l,3) * b(m,n,2,3)) +
     &		a(k,3) * (a(l,1) * b(m,n,3,1) + 
     &		a(l,2) * b(m,n,3,2) + a(l,3) * b(m,n,3,3))
         end do
        end do
       end do
      end do

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          c(i,j,k,l) = a(i,1) * (a(j,1) * d(1,1,k,l) + 
     &		a(j,2) * d(1,2,k,l) + a(j,3) * d(1,3,k,l)) +
     &		a(i,2) * (a(j,1) * d(2,1,k,l) + 
     &		a(j,2) * d(2,2,k,l) + a(j,3) * d(2,3,k,l)) +
     &		a(i,3) * (a(j,1) * d(3,1,k,l) + 
     &		a(j,2) * d(3,2,k,l) + a(j,3) * d(3,3,k,l))
         end do
        end do
       end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      subroutine inverse_3x3(a,b)

      include 'ABA_PARAM.INC'

      dimension a(3,3), b(3,3)

      b(1,1) = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b(1,2) = a(3,2) * a(1,3) - a(1,2) * a(3,3)
      b(1,3) = a(1,2) * a(2,3) - a(2,2) * a(1,3)
      b(2,1) = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b(2,2) = a(1,1) * a(3,3) - a(3,1) * a(1,3)
      b(2,3) = a(2,1) * a(1,3) - a(1,1) * a(2,3)
      b(3,1) = a(2,1) * a(3,2) - a(3,1) * a(2,2)
      b(3,2) = a(3,1) * a(1,2) - a(1,1) * a(3,2)
      b(3,3) = a(1,1) * a(2,2) - a(2,1) * a(1,2)

      det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)

      do i = 1,3
         do j = 1,3
            b(i,j) = b(i,j) / det
         end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Solve simultaneous equations using LU decomposition (Crout's method)
c  Result is stored in b(i)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine simeq(n,a,b)

      include 'ABA_PARAM.INC'

      dimension a(n,n), b(n), index(n)

      call kLU_Decomp(n,a,index)
      call kLU_BackSub(n,a,index,b)

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a matrix using 
c  LU decomposition (Crout's method)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine inverse(n,a,b)

      include 'ABA_PARAM.INC'

      dimension a(n,n), b(n,n), c(n,n), index(n)

      do i = 1,n
         do j = 1,n
            c(i,j) = a(i,j)
         end do
      end do

      do i = 1,n
         do j = 1,n
            b(i,j) = 0.0
         end do
         b(i,i) = 1.0
      end do

      call kLU_Decomp(n,c,index)
      do j = 1,n
         call kLU_BackSub(n,c,index,b(1,j))
      end do

      return
      end

c====================================================================
c====================================================================
c
c  This sub performs an LU Decomposition (Crout's method) on the 
c  matrix "a". It uses partial pivoting for stability. The index()
c  vector is used for the partial pivoting.  The v() vector is 
c  a dummy work area.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine kLU_Decomp(n,a,index)

      include 'ABA_PARAM.INC'

      dimension a(n,n), index(n), v(n)

      tiny = 1.0e-20

c--------------------------------------------------------------------
c  Loop over the rows to get the implicit scaling info.
c--------------------------------------------------------------------

      do i = 1,n
         a_max = 0.0
         do j = 1,n
            a_max = max(a_max,abs(a(i,j)))
         end do !j
         v(i) = 1.0 / a_max
      end do !i

c--------------------------------------------------------------------
c  Begin big loop over all the columns.
c--------------------------------------------------------------------

      do j = 1,n

         do i = 1,j-1
            sum = a(i,j)
            do k = 1,i-1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
         end do

         a_max = 0.0
         do i = j,n
            sum = a(i,j)
            do k = 1,j-1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
            dummy = v(i) * abs(sum)
            if ( dummy .gt. a_max ) then
               imax = i
               a_max = dummy
            end if
         end do

c--------------------------------------------------------------------
c  Pivot rows if necessary.
c--------------------------------------------------------------------

         if ( j .ne. imax ) then
            do k = 1,n
               dummy = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dummy
            end do
            v(imax) = v(j)
         end if
         index(j) = imax

c--------------------------------------------------------------------
c  Divide by the pivot element.
c--------------------------------------------------------------------

         if ( a(j,j) .eq. 0.0 ) a(j,j) = tiny
         if ( j .ne. n ) then
            dummy = 1.0 / a(j,j)
            do i = j+1,n
               a(i,j) = a(i,j) * dummy
            end do
         end if

      end do !j

      return
      end

c====================================================================
c====================================================================
c
c  Solves a set of simultaneous equations by doing back substitution.
c  The answer in returned in the b() vector.  The a(,) matrix
c  must have already been "LU Decomposed" by the above subroutine.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine kLU_BackSub(n,a,index,b)

      include 'ABA_PARAM.INC'

      dimension a(n,n), index(n), b(n)

      ii = 0

c--------------------------------------------------------------------
c  Do the forward substitution.
c--------------------------------------------------------------------

      do i = 1,n
         m = index(i)
         sum = b(m)
         b(m) = b(i)
         if ( ii .ne. 0 ) then
            do j = ii,i-1
               sum = sum - a(i,j) * b(j)
            end do
         else if ( sum .ne. 0.0 ) then
            ii = i
         end if
         b(i) = sum
      end do

c--------------------------------------------------------------------
c  Do the back substitution.
c--------------------------------------------------------------------

      do i = n,1,-1
         sum = b(i)
         if ( i .lt. n ) then
            do j = i+1,n
               sum = sum - a(i,j) * b(j)
            end do
         end if
         b(i) = sum / a(i,i)
      end do

      return
      end
      
c====================================================================
c====================================================================
c
c  Restore a symmetric 4th rank tensor stored in Voigt notation 
c  back to its 4th rank form.
c
c--------------------------------------------------------------------

      subroutine Voigt_to_forth(b,a)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = 1,3
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = 1,3
          ib = k
          if (k.ne.l) ib=9-k-l
          a(i,j,k,l) = b(ia,ib)
          if (ia.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
          if (ib.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
         end do
        end do
       end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Store a SYMMETRIC 4th rank tensor in Voigt notation.
c
c--------------------------------------------------------------------

      subroutine forth_to_Voigt(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=9-k-l
          b(ia,ib) = a(i,j,k,l)
         end do
        end do
       end do
      end do

      return
      end



c====================================================================
c====================================================================
c
c  Perform x**y but while retaining the sign of x.
c
c--------------------------------------------------------------------

      function power(x,y)

      include 'ABA_PARAM.INC'

      if (x.eq.0.0) then
        if (y.gt.0.0) then
          power = 0.0
        else if (y .lt. 0.0) then
          power = 1.0d+300
        else
          power = 1.0
        end if
      else
         power = y * log10(abs(x))
         if (power .gt. 300.) then
           power = 1.d+300
         else
           power = 10.d0 ** power
         end if
         if (x .lt. 0.0) power = -power
      end if

      return
      end

c====================================================================
c====================================================================
c
c  Return the sign of a number.
c
c--------------------------------------------------------------------

      function sgn(a)

      include 'ABA_PARAM.INC'

      sgn = 1.0
      if (a .lt. 0.0) sgn = -1.0

      return
      end 

c====================================================================
c====================================================================
c
c  Calculate the determinant of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      function determinant(a)

      include 'ABA_PARAM.INC'

      dimension a(3,3)

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

      return
      end

c===================================================================
c===================================================================
c
c  Print out an array.
c
c-------------------------------------------------------------------

      subroutine kprint_array(n,a)

      include 'ABA_PARAM.INC'
      dimension a(n,n)

      do i = 1,n
         write(7,'(10f12.5)')(a(i,j),j=1,n)
      end do
      print*,' '

      return
      end




c===================================================================
c===================================================================
c
c  Print out Euler angles in Kocks notation.
c
c-------------------------------------------------------------------

      subroutine kocks_angles(npt,time,array1)

      include 'ABA_PARAM.INC'

      dimension array1(3,3)

      pi = 4 * atan(1.0)

      if (abs(array1(3,3)) .gt. 0.999999) then
        psi   = atan2(array1(2,1),array1(1,1))
        theta = 0.0
        phi   = 0.0
      else
        psi   = atan2(array1(2,3),array1(1,3))
        theta = acos(array1(3,3))
        phi   = atan2(array1(3,2),-array1(3,1))
      end if

      psi   = 180 * psi   / pi
      theta = 180 * theta / pi
      phi   = 180 * phi   / pi

      print*,time,psi,theta,phi

      return
      end
c=======================================================================
c=======================================================================
c
c  Evaluate function to be minimized to zero.  gamma_dot()'s are
c  input and several things are output.
c
c-----------------------------------------------------------------------


      subroutine eval_func( xs0,	dtime,		gamma_dot,	
     &			    xm0,	F1,		num_slip_sys,
     &			    C,		F_p_inv,	F_p_inv_0,
     &			    d0,		Hdir,		Hdyn,
     &			    a0,		Adir,		Adyn,
     &			    g_sat,	F_el,		flow_exp,
     &			    sig,	tau,		gamma_dot_zero,
     &			    d,		func,		xL_p,
     &			    xLatent,	sse,	a, t, t0, t1, a1,
     &                Spk2,gamma,text)

      include 'ABA_PARAM.INC'
      
      dimension
     &	xs0(3,num_slip_sys),	xm0(3,num_slip_sys),  F_p_inv_0(3,3),
     &	F1(3,3),		C(3,3,3,3),		F_p_inv(3,3),
     &	d0(num_slip_sys),	xL_p_inter(3,3),	F_el(3,3),
     &	E_el(3,3),		Spk2(3,3),		sig(3,3),
     &	tau(num_slip_sys),	d(num_slip_sys),	array1(3,3),
     &	func(num_slip_sys),	gamma_dot(num_slip_sys),array2(3,3),
     &	F_el_inv(3,3),		xL_p(3,3),
     &	xs(3,num_slip_sys),	xm(3,num_slip_sys),
     &	a0(num_slip_sys),	a(num_slip_sys), t0(num_slip_sys),
     &	t(num_slip_sys), t1(num_slip_sys), gamma(num_slip_sys),text(3,3)

c*** Note that xs0 and xm0 are in INTERMEDIATE configuration!!!

c--------------------------------------------------------------------
c  Calculate the plastic part of the
c  velocity gradient in the intermediate configuration.
c--------------------------------------------------------------------
c      print*, 'step 1'
      do i = 1,3
        do j = 1,3
          xL_p_inter(i,j) = 0.0
          do k = 1,num_slip_sys
            xL_p_inter(i,j) = xL_p_inter(i,j) + 
     &			xs0(i,k) * xm0(j,k) * gamma_dot(k)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Begin calculation process of F_p_n+1 = exp(xL_p_inter*dt).F_p_n
c--------------------------------------------------------------------
c      print*, 'step 2'
      do i = 1,3
        do j = 1,3
          array1(i,j) = xL_p_inter(i,j) * dtime
        end do
      end do

c--------------------------------------------------------------------
c  Calculate omega.
c--------------------------------------------------------------------
c      print*, 'step 3'
      sum = 0
      do i = 1,3
        do j = 1,3
          sum = sum + array1(i,j) * array1(i,j)
        end do
      end do
      omega = sqrt(0.5 * sum)

c--------------------------------------------------------------------
c  Continue calculating intermediate stuff needed for F_p_n+1
c--------------------------------------------------------------------
c      print*, 'step 4'
      call aa_dot_bb(3,array1,array1,array2)

      do i = 1,3
        if (omega .ne. 0.0) then   ! if omega=0 then no need.
          do j = 1,3
            array1(i,j) = array1(i,j) * sin(omega) / omega +
     &			array2(i,j) * (1-cos(omega)) / omega**2
          end do
        end if
        array1(i,i) = 1 + array1(i,i)
      end do

c--------------------------------------------------------------------
c   Finally multiply arrays to get F_p_inv at end of time step.
c--------------------------------------------------------------------
c      print*, 'step 5'
      call inverse_3x3(array1,array2)
      call aa_dot_bb(3,F_p_inv_0,array2,F_p_inv)


c--------------------------------------------------------------------
c  Multiply F() by F_p_inv() to get F_el()
c--------------------------------------------------------------------
c      print*, 'step 6'
      call aa_dot_bb(3,F1,F_p_inv,F_el)
c      call inverse_3x3(F_el,F_el_inv)

c--------------------------------------------------------------------
c  Rotate director vectors from intermediate configuration to
c  the current configuration.
c--------------------------------------------------------------------

c      do n = 1,num_slip_sys
c        do i = 1,3
c          xs(i,n) = 0.0
c          xm(i,n) = 0.0
c          do j = 1,3
c            xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
c            xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
c          end do
c        end do
c      end do

c--------------------------------------------------------------------
c  Calculate elastic Green Strain
c--------------------------------------------------------------------
c      print*, 'step 7'
      call transpose(3,F_el,array1)
      call aa_dot_bb(3,array1,F_el,E_el)
c        print*,'green strain'
      do i = 1,3
        E_el(i,i) = E_el(i,i) - 1
        do j = 1,3
          E_el(i,j) = E_el(i,j) / 2
        end do 
      end do 

c--------------------------------------------------------------------
c  Multiply the stiffness tensor by the Green strain to get
c  the 2nd Piola Kirkhhoff stress
c--------------------------------------------------------------------
c      print*, 'step 8'
      call aaaa_dot_dot_bb(3,C,E_el,Spk2)

c--------------------------------------------------------------------
c  Convert from PK2 stress to Cauchy stress
c--------------------------------------------------------------------

c      det = determinant(F_el)
c      call transpose(3,F_el,array2)
c      call aa_dot_bb(3,F_el,Spk2,array1)
c      call aa_dot_bb(3,array1,array2,sig)
c        if (det.eq.0) then
c      do i = 1,3
c        do j = 1,3
c          sig(i,j) = 0.
c        end do
c      end do
c        endif

c        if (det.ne.0) then
c      do i = 1,3
c        do j = 1,3
c          sig(i,j) = sig(i,j) / det
c        end do
c      end do
c         endif

c--------------------------------------------------------------------
c  Calculate resolved shear stress for each slip system.
c--------------------------------------------------------------------
c      print*, 'step 9'
      do k = 1,num_slip_sys
        tau(k) = 0.0
        do j = 1,3
          do i = 1,3
c            tau(k) = tau(k) + xs(i,k) * xm(j,k) * sig(i,j)
	      tau(k) = tau(k) + xs0(i,k) * xm0(j,k) * Spk2(i,j)
          end do 
        end do
      end do
      
c--------------------------------------------------------------------
c  Calculate hardening law for drag stress
c--------------------------------------------------------------------
c      print*, 'step 10'
      sum = 0
      do i = 1,num_slip_sys
        sum = sum + abs(gamma_dot(i))
      end do

      do k = 1,num_slip_sys
        sum00 = xLatent * sum - (xLatent - 1) * abs(gamma_dot(k))
        d(k) = (d0(k) + Hdir*sum00*dtime) / (1 + Hdyn*sum*dtime)
      end do


c      do k = 1,num_slip_sys
c        d(k) = d0(k)
c      end do

c--------------------------------------------------------------------
c  Calculate hardening law for back stress
c--------------------------------------------------------------------
c      print*, 'step 11'
      do k = 1,num_slip_sys
        a(k) = (a0(k) + Adir * dtime * gamma_dot(k))
     &			 / (1 + Adyn * dtime * abs(gamma_dot(k)))
c       Print*,'within eval, ai is', a0(k),a(k),Adir,dtime
c       Print*,'within eval, ai p2 is',gamma_dot(k),Adyn
      end do
c--------------------------------------------------------------------
c  Calculate hardening law for threshold stress
c--------------------------------------------------------------------
c      print*, 'step 12'
      do k = 1,num_slip_sys
        t(k) = t0(k)
      end do

c--------------------------------------------------------------------
c  Modify threshold stress to include non-Schmid terms
c--------------------------------------------------------------------
c      print*, 'step 13'
      do k = 1,num_slip_sys
        t1(k) =t(k)
      end do

	t1(4) = t(4) + a1*tau(7) - a1*tau(8)
	t1(5) = t(5) + a1*tau(9) - a1*tau(10)
	t1(6) = t(6) + a1*tau(11) - a1*tau(12)


c--------------------------------------------------------------------
c  Calculate function values.
c--------------------------------------------------------------------
c      print*, 'step 14'
c      do k = 1,num_slip_sys
c        func(k) = tau(k) - a(k) + t(k) - d(k) * power( (gamma_dot(k)/
c     &			gamma_dot_zero), (1./flow_exp) )
c      end do

      do k = 1,num_slip_sys
	  taud = abs(tau(k)-a(k))-t1(k)
	  if (taud.gt.0.) then
	    func(k) = -power((taud/d(k)),flow_exp) + (gamma_dot(k)/
     &    gamma_dot_zero)*sgn(tau(k)-a(k))
	  else
	    func(k) = 0.
	  end if
	end do
c--------------------------------------------------------------------
c  Calculate root mean square error that is to be minimized to zero.
c--------------------------------------------------------------------
c      print*, 'step 15'
      sse = 0.0
      gamma_dot_max = abs(gamma_dot(1))
      do k = 1,num_slip_sys
        sse = sse + abs(gamma_dot(k)) * func(k) ** 2
c        sse = sse + (abs(gamma_dot(k)) * func(k)) ** 2
        gamma_dot_max = max(abs(gamma_dot(k)),gamma_dot_max)
      end do
c      print*, 'end do loop'       
c      if (gamma_dot_max .gt. 0.0) sse = sse / gamma_dot_max
       if (gamma_dot_max .gt. 0.0) sse = sqrt(sse / gamma_dot_max) / 
     &                                                  num_slip_sys


      return
      end
c====================================================================
c  Abaqus subroutine used for writing out data to files from the UMAT
c lop = 0 -  called at the start of the analysis
c lop = 1 -  called at the start of an analysis increment
c lop = 2 -  called at the end of an analysis increment
c lop = 3 -  called at the end of the analysis
c lop = 4 -  called at teh beginning of a restart analysis
c====================================================================
      subroutine uexternaldb (lop,lrestart,time,dtime,kstep,kinc)

      include 'ABA_PARAM.INC'

c      common/flag1/storsdv(9,4),jnoel(9),xcord(9,3)
      common/flag1/storsdv(9,4),jnoel(9),xcord(9,3),kincp,
     &icount
      
      dimension time(2)        

      if (lop.eq.0) then

       open(unit=101,
     & file='/crunch/jason/ms_thesis/external.txt',
     & status='unknown',form='formatted')
       close(101)

      else if ((lop.eq.2).and.(kinc.eq.1).and.(kstep.eq.1)) then

       open(unit=102,
     & file='/crunch/jason/ms_thesis/centroid.txt',
     & status='unknown',form='formatted')
       do i = 1,9
        write(102,116),jnoel(i),xcord(i,1),xcord(i,2),xcord(i,3)
       end do
       close(102)

      else if ((lop.eq.2).and.(time(2).ge.5.0)) then

       open(unit=101,
     & file='/crunch/jason/ms_thesis/external.txt',
     & status='unknown',form='formatted',access='append')
        do i = 1,9
         write(101,115),jnoel(i),time(2),storsdv(i,1),storsdv(i,2),
     &   storsdv(i,3),storsdv(i,4) 
        end do
       close(101)       
  115  format(i5,2x,e13.4e2,2x,e13.4e2,2x,e13.4e2,2x,e13.4e2,
     & 2x,e13.4e2)

      else if (lop.eq.3) then

        close(101)
        close(102)
      
      end if

  116 format(i5,2x,f10.4,2x,f10.4,2x,f10.4)   


      return
      end
c====================================================================
c====================================================================
c====================================================================
c
c  Generate uniform random numbers between 0 and 1.
c  This is a fairly coarse generator in that it returns only
c  714,025 different values.
c
c  ***IMPORTANT***
c  Set idum equal to any positive number between 1 and 714,024
c  to initialize or reinitialize the sequence!!!
c
c  call with:   x = ran0(idum)   for example.
c
c  Reference: "Numerical Recipes" Section 7.1  p. 195
c
c--------------------------------------------------------------------

c      function ran0(idum)

c      implicit double precision (a-h,o-z)

c      Parameter (M=714025, IA=1366, IC=150889, RM=1./M)

c      idum = mod(ia*idum+ic,m)
c      ran0 = float(idum)/float(m)

c      return
c      end

c -----------------



