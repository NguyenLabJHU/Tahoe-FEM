      subroutine double_it(a_double, i_integer)
      implicit double precision (a-h,o-z)
	  
	  a_double = a_double*2.0
	  i_integer = i_integer*2
	  
	  return
	  end
c
c
      subroutine double_int(i_integer)
      implicit double precision (a-h,o-z)
       
          i_integer = i_integer*2

          return
          end
