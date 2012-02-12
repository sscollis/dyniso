!========================================================================
!  Miscellaneous auxiliary routines
!========================================================================
real*4 function second()
  real*4 t1
  real*4 tarray( 2 )
#ifndef __GFORTRAN__
  real*4 etime
  external etime
#endif
  t1 = etime( tarray )
  second = tarray( 1 )
  return 
end function second
