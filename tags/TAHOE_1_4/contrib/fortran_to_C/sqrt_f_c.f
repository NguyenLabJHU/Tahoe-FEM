      subroutine sqrt_f_c(a_double)
      implicit double precision (a-h,o-z)
	  
	  external sqrt_c
	  call sqrt_c(a_double)
	  
	  return
	  end
	 
