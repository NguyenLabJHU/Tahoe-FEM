      subroutine print_double(a_double)
      implicit double precision (a-h,o-z)
	  
	  write(*,998) a_double
 998  format('double = ',e18.13)

	  return
	  end
c
c	 
      subroutine print_int(i)
      implicit double precision (a-h,o-z)
	  
	  write(*,999) i
 999  format('integer = ',i8)

	  return
	  end
	 
