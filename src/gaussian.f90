module randommod
implicit none
contains
real(c_float) function gaussian() bind(c,name="gaussian")
  use,intrinsic :: iso_c_binding
  implicit none
  !	real rd

  real(c_float) x1,x2,w,y1,y2,n1,n2
  integer(c_int) use_last

  use_last=0
  y2=0
  w=2.0
  if(use_last.ne.0) then
    y1=y2
    use_last=0
  else
    do while (w.ge.1.0)
      call random_number(n1)
      call random_number(n2) 
      x1=2.0*n1-1.0
      x2=2.0*n2-1.0
      w = x1 * x1 + x2 * x2
    enddo
    w=((-2.0*log(w))/w)**0.5
    y1=x1*w
    y2=x2*w
    use_last=1
  endif
  gaussian=y1
end function

end module
