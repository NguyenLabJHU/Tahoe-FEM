      program test
      implicit double precision (a-h,o-z)
	  
	  a_double = 1.1E+00
	  i_integer = 2

      write(*,998) i_integer
      call double_int(i_integer)
      write(*,998) i_integer

	  write(*,999) a_double, i_integer
	  call double_it(a_double, i_integer)
	  write(*,999) a_double, i_integer
	  
998   format(i8)
999   format(e18.13,i8)
	  
	  end
