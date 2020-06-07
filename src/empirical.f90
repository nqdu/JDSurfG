module empirical
contains
subroutine empirical_relation(vsz,vpz,rhoz)bind(c,name="empirical_relation")
    use,intrinsic :: iso_c_binding
  
    real(c_float), intent(in) ::  vsz
    real(c_float),intent(out) ::vpz,rhoz
  
    vpz=0.9409 + 2.0947*vsz - 0.8206*vsz**2+ &
            0.2683*vsz**3 - 0.0251*vsz**4
    rhoz=1.6612*vpz- 0.4721*vpz**2 + &
            0.0671*vpz**3 - 0.0043*vpz**4 + & 
            0.000106*vpz**5
  
end subroutine empirical_relation

! compute tmp1=drho/dvp and tmp2=dvp/dvs, by Brocher's relation
subroutine empirical_deriv(vp,vs,tmp1,tmp2)bind(c,name="empirical_deriv")
    use,intrinsic :: iso_c_binding
    real(c_float),value,INTENT(IN) :: vp,vs
    real(c_float),intent(out) ::tmp1,tmp2

    tmp1 = 1.6612 - 0.4721*2*vp + 0.0671*3*vp**2 - 0.0043*4*vp**3 + 0.000106*5*vp**4
    tmp2 = 2.0947 - 0.8206*2*vs + 0.2683*3 * vs**2 - 0.0251*4*vs**3
    
end subroutine empirical_deriv
end