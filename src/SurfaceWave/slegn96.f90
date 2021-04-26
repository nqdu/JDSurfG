!c----------------------------------------------------------------------c
!c                                                                      c
!c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
!c      VOLUME III                                                      c
!c                                                                      c
!c      PROGRAM: SLEGN96                                                c
!c                                                                      c
!c      COPYRIGHT 1996, 2010                                            c
!c      R. B. Herrmann                                                  c
!c      Department of Earth and Atmospheric Sciences                    c
!c      Saint Louis University                                          c
!c      221 North Grand Boulevard                                       c
!c      St. Louis, Missouri 63103                                       c
!c      U. S. A.                                                        c
!c                                                                      c
!c----------------------------------------------------------------------c
!c
!c
!c     This program calculates the group velocity and partial
!c     derivatives of Love waves for any plane multi-layered
!c     model.  The propagator-matrix, instead of numerical-
!c     integration method is used, in which the Haskell rather
!c     than Harkrider formalisms are concerned.
!c
!c     Developed by C. Y. Wang and R. B. Herrmann, St. Louis
!c     University, Oct. 10, 1981.  Modified for use in surface
!c     wave inversion, with addition of spherical earth flattening
!c     transformation and numerical calculation of group velocity
!c     partial derivatives by David R. Russell, St. Louis
!c     University, Jan. 1984.
!c
!c     Rewrite of theory to agree more with hspec96 wavenumber
!c     integration code
!c
!c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!c Revision history:
!c       07 AUG 2002 - make string lengths 120 characters from 80 
!c       13 OCT 2006 - verbose output of energy integrals
!c       26 SEP 2008 - fixed undefined LOT in subroutine up
!c       14 JUN 2009 - reformulate in general terms as per book
!c       01 AUG 2010 - change code to agree with book  for
!c                     migration to TI. Also clean up code using
!c                     implicit none
!c       20 SEP 2012 - cleaned up some compiler warnings
!c       25 FEB 2017 - fixed error in subroutine varsv. The
!c           computation of sin (nu z)/nu yielded NaN when nu = 0
!c           The corrected code gives lim  sin(nu z) / nu = z
!c                                         nu->0
!c       29 JUL 2017 - modify insert if depth is in halfspace
!c       22 FEB 2018 - verbose also gives LAGR/(omega^2 I0)
!        21 Apr 2021 -- Modified by Nanqiao Du, Fortran 90 version
!c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!===================================
! MODULE : Love wave Model
!===================================
module LoveWaveModel
use,intrinsic         :: iso_c_binding
implicit none

! global variables
!============================================
! model param
integer(c_int)              :: mmax
real(c_double),ALLOCATABLE  :: zd(:),zb(:),zrho(:),za(:),xmu(:)
integer(c_int),ALLOCATABLE  :: iwat(:)

! eigen function
real(c_double),ALLOCATABLE  :: uu(:),tt(:),dcdb(:),dcdh(:),dcdr(:)
real(c_double)              :: uu0(4)

! spherical earth
real(c_double),ALLOCATABLE  :: vtp(:),dtp(:),rtp(:)

! save
real(c_double),ALLOCATABLE  :: exl(:)

! love wave var
real(c_double)              :: cosq,sinq,yl,zl,mu

! sumi
real(c_double)              :: sumi0,sumi1,sumi2,flagr,ale,ugr

! openmp wrapper
!$omp threadprivate(mmax,zd,zb,zrho,za,uu,tt,dcdb,dcdh,dcdr,uu0)
!$omp threadprivate(vtp,dtp,rtp,xmu,exl,iwat)
!$omp threadprivate(cosq,sinq,yl,zl,mu)
!$omp threadprivate(sumi0,sumi1,sumi2,flagr,ale,ugr)

contains 
subroutine Alloc_Model(nlayer)
  implicit none
  integer(c_int)        ::nlayer 
  mmax = nlayer 
  ALLOCATE(zd(mmax),zb(mmax),za(mmax),uu(mmax),tt(mmax))
  ALLOCATE(dcdb(mmax),dcdh(mmax),dcdr(mmax))
  ALLOCATE(iwat(mmax))
  ALLOCATE(xmu(mmax),exl(mmax),zrho(mmax))
end subroutine Alloc_Model

subroutine Dealloc_Model()
  implicit none
  DEALLOCATE(zb,zd,za,uu,tt,dcdb,dcdh,dcdr,zrho)
  DEALLOCATE(iwat,xmu,exl)
end subroutine Dealloc_Model

subroutine bldsph()
  !c-----
  !c       Transform spherical earth to flat earth
  !c
  !c       Schwab, F. A., and L. Knopoff (1972). 
  !c            Fast surface wave and free
  !c       mode computations, 
  !c            in  Methods in Computational Physics, Volume 11,
  !c       Seismology: Surface Waves and Earth Oscillations,  
  !c            B. A. Bolt (ed),
  !c       Academic Press, New York
  !c
  !c       Love Wave Equations  44, 45 , 41 pp 112-113
  !c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
  !c
  !c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
  !c     slowness instead of using average velocity for the layer
  !c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
  !c
  !c       This is for Love Waves
  !c
  use,intrinsic                 :: iso_c_binding
  implicit none
  integer(c_int)                :: i
  real(c_double)                :: ar, dr, r0, r1, z0, z1,tmp
  !c-----
  !c       vtp is the factor used to convert 
  !c            spherical velocity to flat velocity
  !c       rtp is the factor used to convert 
  !c            spherical density  to flat density 
  !c       dtp is the factor used to convert 
  !c            spherical boundary to flat boundary
  !c-----
  !c-----
  !c       duplicate computations of srwvds(IV)
  !c-----
  ar=6371.0d0
  dr=0.0d0
  r0=ar
  zd(mmax)=1.0 
  do i=1,mmax
    dr= dr + zd(i)
    r1=ar-dr
    z0=ar*dlog(ar/r0)
    z1=ar*dlog(ar/r1)
    dtp(i) = ar/r0 !+ ar/r1

    !use layer midpoint
    TMP=(2.0d0 *ar)/(r0+r1)
    vtp(i) = tmp
    rtp(i) = tmp**(-5)

    ! earth flattening transform 
    zb(i) = zb(i) * tmp 
    zrho(i) = zrho(i) * rtp(i)
    zd(i) = z1 - z0
    r0=r1
  ENDDO
  zd(mmax) = 0.0
  return
end subroutine bldsph

end module LoveWaveModel

!===================================
! MODULE : Love wave Kernel
!===================================
module LoveWaveKernel
use,intrinsic               :: iso_c_binding
implicit none

contains
subroutine shfunc(omega,wvno)
  !c-----
  !c       This routine evaluates the eigenfunctions by calling sub up.
  !c-----
  !implicit double precision (a-h,o-z)
  use,intrinsic   :: iso_c_binding
  use LoveWaveModel,only : mmax,iwat,uu,tt,uu0,exl
  implicit none
  real(c_double)        :: fl,omega,wvno

  ! local
  real(c_double)        :: ext,umax,fact
  integer(c_int)        :: k
  !c-----
  !c       get the eigenfunctions in an extended floating point form
  !c       by propagating from the bottom upward. Note since we move in this
  !c       direction the propagator matrix A(d) is evaluated as A(-d)
  !c-----
  call up(omega,wvno,fl)
  uu0(1)=1.0
  !c-----
  !c       uu0(2)=stress0 is actually the value of period equation.
  !c       uu0(3) is used to print out the period equation value before
  !c       the root is refined.
  !c-----
  uu0(2)=fl
  uu0(3)=0.0
  uu0(4)=0.0
  !c-----
  !c       convert to actual displacement from the extended floating point
  !c       working top to down, where the amplitude should be small
  !c       Also normalize so that the V(0) = 1
  !c-----
  ext=0.0
  umax = uu(1)
  tt(1) = 0.0
  do k=2,mmax
    if(iwat(k).eq.0)then
        ext=ext+exl(k-1)
        fact=0.0
        if(ext.lt.80.0) fact=1./dexp(ext)
        uu(k)=uu(k)*fact
        tt(k)=tt(k)*fact
    else
        uu(k) = 0.0
        tt(k) = 0.0
    endif
    if(abs(uu(k)).gt.abs(umax))then
        umax = uu(k)
    endif
  enddo
  if(uu(1).ne.0.0)then
    umax = uu(1)
  endif
  if(abs(umax).gt.0.0)then
    do k=1,mmax
      if(iwat(k).eq.0)then
          uu(k) = uu(k) / umax
          tt(k) = tt(k) / umax
      endif
    enddo
  endif
  !c-----
  !c       force boundary condition at the free surface
  !c       free surface is stress free
  !c-----
  return
end subroutine shfunc

subroutine varl(m,rb,omega,wvno,xkb,dpth,eexl)
  !c-----
  !c     find variables cosQ,  sinQ, etc.
  !c-----
  !c     To handle the hyperbolic functions correctly for large
  !c     arguments, we use an extended precision procedure,
  !c     keeping in mind that the maximum precision in double
  !c     precision is on the order of 16 decimal places.
  !c
  !c     So  cosq = 0.5 ( exp(+q) + exp(-q))
  !c              = exp(q) * 0.5 * ( 1.0 + exp(-2q) )
  !c     becomes
  !c         cosq = 0.5 * (1.0 + exp(-2q) ) with an exponent q
  !c     In performing matrix multiplication, we multiply the modified
  !c     cosp terms and add the exponents. At the last step
  !c     when it is necessary to obtain a true amplitude,
  !c     we then form exp(p). For normalized amplitudes at any depth,
  !c     we carry an exponent for the numerator and the denominator, and
  !c     scale the resulting ratio by exp(NUMexp - DENexp)
  !c
  !!!c
  !c     HSKA        cosq  sinq
  !c
  !c     When the extended floating point is used, we use the
  !c     largest exponent for each, which is  the following:
  !c
  !c     Let sex = s exponent > 0 for evanescent waves = 0 otherwise
  !c
  !c     Then the modified matrix elements are as follow:
  !c
  !c     Haskell:  cosq -> 0.5 ( 1 + exp(-2q) ) exponent = pex
  !c-----
  !implicit double precision (a-h,o-z)
  use LoveWaveModel,only    : zb,zrho,cosq,sinq,&
                              yl,zl,mu
  implicit none 
  integer m
  real(c_double)            ::  rb, xkb, omega, wvno, dpth, eexl,q

  real(c_double)            :: wvnop, wvnom , fac

  !c-----
  !c     define the  horizontal wavenumber for the S-wave
  !c-----
  xkb = omega/zb(m)
  !c-----
  !c     define the vertical wavenumber for the given wvno
  !c-----
  wvnop = wvno + xkb
  wvnom = dabs(wvno - xkb)
  fac = wvnop*wvnom
  !C      WrItE(*,*)'fac:',m,fac
  rb=dsqrt(wvnop*wvnom)
  q = rb*dpth
  mu = zrho(m)*zb(m)*zb(m)
  !c-----
  !c     examine S-wave eigenfunctions
  !c        checking whether c > vs, c = vs, c < vs
  !c-----
  eexl = 0.0
  if(wvno.lt.xkb)then
    sinq=dsin(q)
    yl=sinq/rb
    zl=-rb*sinq
    cosq=dcos(q)
  else if(wvno.eq.xkb)then
    cosq=1.0d+00
    yl=dpth
    zl=0.0d+00
  else if(wvno.gt.xkb)then
    eexl = q
    fac = 0.0d+00
    if(q.lt.18.0d+00)fac = dexp(-2.0d+0*q)
    cosq = ( 1.0d+00 + fac ) * 0.5d+00
    sinq = ( 1.0d+00 - fac ) * 0.5d+00
    yl = sinq/rb
    zl = rb*sinq
  endif

  return
end subroutine varl

subroutine hskl(hl,iwat)
  use LoveWaveModel, only : cosq,yl,zl,mu
  implicit none
  !c-----
  !c      procedure arguments
  !c-----
  real(c_double)          :: hl(2,2)
  integer(c_int)          :: iwat

  if(iwat.eq.0)then
    hl(1,1) = cosq
    hl(1,2) = yl/mu
    hl(2,1) = zl*mu
    hl(2,2) = cosq
  else
    hl(1,1) = 1.0d+00
    hl(1,2) = 0.0d+00
    hl(2,1) = 0.0d+00
    hl(2,2) = 1.0d+00
  endif
  return
end subroutine hskl

subroutine emat(xnub,wvno,el,einvl)
  implicit none
  !c-----
  !c     compute the E and E^-1 matrices
  !c-----
  real(c_double)              :: wvno
  complex(c_double_complex)   :: xnub, el(2,2), einvl(2,2)
        
  einvl(1,1) =  0.5/wvno
  einvl(1,2) =  0.5/(wvno*xnub)
  einvl(2,1) =  0.5/wvno
  einvl(2,2) = -0.5/(wvno*xnub)
  el(1,1) =  wvno
  el(1,2) =  wvno
  el(2,1) =  wvno*xnub
  el(2,2) = -wvno*xnub
  return
end subroutine emat

subroutine up(omega,wvno,fl)
  use LoveWaveModel, only  : mmax,zb,zd,iwat,uu,tt,xmu,exl
  implicit none
  !c-----
  !c       This routine calculates the elements of Haskell matrix,
  !c       and finds the eigenfunctions by analytic solution.
  !c-----
  !c       History:
  !c       
  !c       16 JUN 2009 - generalized the routine so that  the
  !c         nature of the propagator matrices is buried in the
  !c         subroutine hskl
  !c         This will permit simpler extension to TI
  !c-----
  !implicit double precision (a-h,o-z)

  real(c_double)        :: fl,omega,wvno
  real(c_double)        :: a11, a12, a21, a22,rb, xkb, dpth,hl(2,2)
  real(c_double)        :: eexl,amp0,str0,rr,ss,ttlast
  integer(c_int)        :: mmx1,k,k1

  !c
  !c-----
  !c       apply the boundary conditions to the halfspace
  !c-----
  !c-----
  !c     kludge for fluid core
  !c-----
  if(zb(mmax).gt.0.01)then
    dpth = 0.0
    call varl(mmax,rb,omega,wvno,xkb,dpth,eexl)
    if(wvno.lt.xkb)then
        !write(LOT,*) ' imaginary nub derivl'
        !write(LOT,*)'omega,wvno,b(mmax)',omega,wvno,zb(mmax)
    endif
    uu(mmax)=1.0
    tt(mmax)=-xmu(mmax)*rb
  else
    uu(mmax)=1.0
    tt(mmax)=0.0
  endif
  exl(mmax) = 0.0
  mmx1=mmax-1
  ttlast = 0.0
  do k=mmx1,1,-1
    if(iwat(k).eq.0)then
      dpth=zd(k)
      call varl(k,rb,omega,wvno,xkb,dpth,eexl)
      call hskl(hl,iwat(k))
      k1=k+1
      !c We actually use A^-1 since we go from the bottom to top
      a11 = hl(1,1)
      a22 = hl(2,2)
      a12 = - hl(1,2)
      a21 = - hl(2,1)
      !C            WRITE(6,*)'eexl:',eexl
      !C            WRITE(6,*)'hl:',hl

      amp0 = a11 * uu(k1) + a12*(tt(k1))
      str0 = a21 * uu(k1) + a22*(tt(k1))

      !c now normalize and save the normalization
      rr=dabs(amp0)
      ss=dabs(str0)
      if(ss.gt.rr) rr=ss
      if(rr.lt.1.d-30) rr=1.d+00
      exl(k)=dlog(rr)+eexl
      uu(k)=amp0/rr
      tt(k)=str0/rr
      ttlast = tt(k)
    endif
  enddo
  fl=ttlast
end subroutine up

subroutine energy(omega,wvno,&
                Eut,Edut,Ed2ut,Eut0,Ett0,&
                lss,lrr)
  !!c-----
  !!c       This routine calculates the values of integrals I0, I1,
  !!c       and I2 using analytic solutions. It is found
  !!c       that such a formulation is more efficient and practical.
  !!c
  !!c       This is based on a suggestion by David Harkrider in 1980.
  !!c       Given the eigenfunctions, which are properly defined,
  !!c       define the potential coefficients at the bottom and
  !!c       top of each layer. If the V(z) = A exp (nu z) + B exp (-nu z)
  !!c       We want the B at the top of the layer for z=0, and the A at the
  !!c       bottom of the layer, from the relation K=[Kup Kdn]^T
  !!c       = [A B]^T = E^-1 U
  !!c-----
  !!c       History:
  !!c       
  !!c       16 JUN 2009 - generalized the routine so that  the
  !!c         nature of the propagator matrices is buried in the
  !!c         This will permit simpler extension to TI
  !!c-----
  !implicit double precision (a-h,o-z)
  use LoveWaveModel, only        :mmax,zd,zb,zrho,xmu,iwat,&
                                  uu,tt,dcdb,dcdh,dcdr,sumi0,&
                                  sumi1,sumi2,flagr,ale,ugr
  implicit none
  complex(c_double_complex)     :: nub,xnub,exqq,f1,f2,f3
  real(c_double)                :: omega,wvno
  COMPLEX(c_double_complex)     :: CDEXP,kmup, km1dn
  complex(c_double_complex)     :: einvl(2,2), el(2,2)
  complex(c_double_complex)     :: f, g
  real(c_double)                :: TN, TL,Eut,Edut,Ed2ut,Eut0,Ett0,upup,dupdup
  real(c_double)                :: DCDBH, DCDBV, DCDRSH,dcr,dcb,fac,&
                                    c2,llflag,dvdz,dfac
  real(c_double)                :: VSHH, VSHV,c,omega2,wvno2,drho,dpth,&
                                    dmu,rb,xkb,eexl
  integer(c_int)                :: k,k1,LSS,lrr

  c=omega/wvno
  omega2=omega*omega
  wvno2=wvno*wvno
  sumi0=0.0d+00
  sumi1=0.0d+00
  sumi2=0.0d+00
  !c-----
  !c RBH beware when nub = 0
  !c-----
  do k=1,mmax
    if(iwat(k).eq.0)then
      TN = zrho(k)*zb(k)*zb(k)
      TL = zrho(k)*zb(k)*zb(k)
      VSHH = zb(k)
      VSHV = zb(k)
      k1=k+1
      drho=zrho(k)
      dpth=zd(k)
      call varl(k,rb,omega,wvno,xkb,dpth,eexl)
      dmu=xmu(k)
      if(rb .lt. 1.0d-10)rb = 1.0d-10

      !c-----
      !c      for safety do not permit rb = 0
      !c-----
      if(k.eq.mmax) then
        upup  =(0.5d+00/rb)*uu(mmax)*uu(mmax)
        dupdup=(0.5d+00*rb)*uu(mmax)*uu(mmax)
      else
        !c-----
        !c               use Einv to get the potential coefficients
        !c               from the U and T
        !c-----
        nub=dcmplx(rb,0.0d+00)
        if(wvno.lt.xkb) nub=dcmplx(0.0d+00,rb)
        xnub=dmu*nub
        call emat(xnub,wvno,el,einvl)
        km1dn = einvl(2,1)*uu(k)  + einvl(2,2)*tt(k)
        kmup  = einvl(1,1)*uu(k1) + einvl(1,2)*tt(k1)

        ! c               get the f function
        f3=nub*dpth
        exqq=0.0d+00
        if(dreal(f3).lt.40.0) exqq=cdexp(-2.d+00*f3)
        f=(1.d+00-exqq)/(2.d+00*nub)

        !c get the g function
        exqq=0.0d+00
        if(dreal(f3).lt.75.0) exqq=cdexp(-f3)
        g=dpth*exqq

        f1 = f*(el(1,1)*el(1,1)*kmup*kmup  &
        + el(1,2)*el(1,2)*km1dn*km1dn)
        f2 = g*(el(1,1)*el(1,2)+el(1,1)*el(1,2))*kmup*km1dn

        !c cast to a real  upup
        upup = real(f1 + f2)

        !c cast to a real dupdup
        dupdup = real(nub*nub*(f1 - f2))
      endif

      sumi0=sumi0+drho*upup
      sumi1=sumi1+TN*upup
      sumi2=sumi2+TL*dupdup
      dcr=-0.5d+00*c*c*c*upup
      dcb=0.5d+00*c*(upup+dupdup/wvno2)
      DCDBH = c*drho*VSHH*upup
      DCDBV = c*drho*VSHV*dupdup/wvno2
      dcdb(k)=DCDBH + DCDBV
      DCDRSH= 0.5*c*(-c*c*upup + VSHH*VSHH*upup + &
      VSHV*VSHV*dupdup/wvno2)
      dcdr(k)=DCDRSH
    else
      dcdb(k)=0.0
      dcdr(k)=0.0
    endif
  enddo
  do k=1,mmax
    if(iwat(k).eq.0)then
      dcdb(k)=dcdb(k)/sumi1
      dcdr(k)=dcdr(k)/sumi1
    else
      dcdb(k)=0.0
      dcdr(k)=0.0
    endif
  enddo
  flagr=omega2*sumi0-wvno2*sumi1-sumi2
  ugr=sumi1/(c*sumi0)
  !print*,ugr
  ale=0.5d+00/sumi1

  !c-----
  !c       define partial with respect to layer thickness
  !c-----
  !c       fac = 0.5d+00*c**3/(omega2*sumi1)
  fac = ale*c/wvno2
  c2 = c*c
  llflag = 0
  do k=1,mmax
    if(iwat(k).eq.0)then
      if(llflag.eq.0)then
        drho = zrho(k)
        dmu  = xmu(k)
        dvdz = 0.0
      else 
        drho = zrho(k) - zrho(k-1)
        dmu  = xmu(k) - xmu(k-1)
        dvdz = tt(k)*tt(k)*(1.0/xmu(k) - 1.0/xmu(k-1))
      endif
      dfac = fac * ( uu(k)*uu(k)* &
      (omega2*drho - wvno2*dmu) + dvdz)
      if(dabs(dfac).lt.1.0d-38)then
        dcdh(k) = 0.0
      else
        dcdh(k) = dfac
      endif
        llflag = llflag + 1
    else
        dcdh(k) = 0.0
    endif
  enddo
  !c-----
  !c       compute the eigenfuntions and depth derivatives
  !c       at the source depth
  !c-----
  if(iwat(lss).eq.0)then
    Eut = uu(lss)
    Edut = tt(lss)/xmu(lss)
    Ed2ut = ( - (omega/zb(lss))**2 + wvno2)*Eut
  else
    Eut = 0.0d+00
    Edut = 0.0d+00
    Ed2ut = 0.0d+00
  endif
  if(iwat(lrr).eq.0)then
    Eut0 = uu(lrr)
    Ett0 = tt(lrr)
  else
    Eut0 = 0.0d+00
    Ett0 = 0.0d+00
  endif
  return
end subroutine energy

subroutine splove(om,c,csph,usph,ugr)
  !!c-----
  !!c       Transform spherical earth to flat earth
  !c       and relate the corresponding flat earth dispersion to spherical
  !!c
  !!c       Schwab, F. A., and L. Knopoff (1972). Fast surface wave 
  !!c            and free
  !!c       mode computations, in  
  !!c            Methods in Computational Physics, Volume 11,
  !!c       Seismology: Surface Waves and Earth Oscillations,  
  !!c            B. A. Bolt (ed),
  !!c       Academic Press, New York
  !!c
  !!c       Rayleigh Wave Equations 111, 114 p 144
  !!c
  !!c       Partial with respect to parameter uses the relation
  !!c       For phase velocity, for example,
  !!c
  !!c       dcs/dps = dcs/dpf * dpf/dps, c is phase velocity, p is
  !!c       parameter, and f is flat model
  !!c
  !!c       om      R*8     angular frequency
  !!c       c       R*8     phase velocity
  !!c       mmax    I*4     number of layers 
  !!c-----
  use LoveWaveModel, only        : mmax,vtp,dtp,rtp,dcdb,dcdh,dcdr
  real(c_double)                :: om,tm,c,csph,usph
  real(c_double)                :: ugr,a
  INTEGER(c_int)                :: i 
  a = 6371.0d0
  tm=sqrt(1.+(3.0*c/(2.*a*om))**2)
  do i=1,mmax
    dcdb(i)=dcdb(i)*  vtp(i)/(tm**3)
    dcdh(i)=dcdh(i)*  dtp(i)/ (tm**3)
    dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
  enddo
  csph = c / tm
  usph = ugr * tm
  return
  end

subroutine slegn96(thk,vs,rhom,nlayer,&
                  t,cp,cg,disp,stress,dc2db,dc2dh,dc2dr,&
                  iflsph)bind(c,name='slegn96_')
  use LoveWaveModel
  implicit none
  integer(c_int),value,INTENT(IN)   ::nlayer,iflsph
  real(c_float),INTENT(IN)          :: thk(nlayer),vs(nlayer),rhom(nlayer)
  real(c_double),INTENT(IN)         :: t 
  real(c_double),INTENT(INOUT)      :: cp,cg 
  real(c_double),INTENT(INOUT)      :: disp(nlayer),stress(nlayer),&
                                       dc2db(nlayer),dc2dh(nlayer),&
                                       dc2dr(nlayer)
  
  ! local variables
  integer(c_int)           :: i,nsph
  real(c_double)           :: c, omega, wvno, gamma, csph, usph,twopi
  real(c_double)           :: Eut, Edut, Eut0, Ett0, Ed2ut
  real(c_double),PARAMETER :: pi = 3.1415926535898

  ! velocity model
  real(c_double)           :: sums
  !real(c_float)            :: sut,sdut,sd2ut,suz,sduz,sd2uz,sale,wvnsrc,sur0

  ! allocate arrays

  ! copy model
  call Alloc_Model(nlayer)
  zb(1:mmax) = vs(1:mmax);
  zrho(1:mmax) = rhom(1:mmax); zd(1:mmax) = thk(1:mmax)      

  ! check water layer
  do i=1,mmax 
    if(zb(i) == 0.0)then
      iwat(i) = 1
    else
      iwat(i) = 0
    endif
  enddo

  ! earth flattening transformation
  nsph = iflsph
  if(nsph > 0) then 
    ALLOCATE(vtp(mmax),dtp(mmax),rtp(mmax))
    call bldsph()
  endif
  
  ! compute mu
  xmu(1:mmax) = zrho(1:mmax) * zb(1:mmax)**2   

  ! main part
  twopi = 2.0 * pi
  omega = twopi/t
  c=cp
  wvno=omega/c
  call shfunc(omega,wvno)
  call energy(omega,wvno,Eut,Edut,Ed2ut,Eut0,Ett0,1,1)
  
  ! no attenuation
  gamma = 0.0

  ! earth flattening transform for kernel
  if(nsph > 0) then 
    call splove(omega,c,csph,usph,ugr)
    wvno = omega/csph
  else
    csph = c 
    usph = ugr
  endif

  !check for underflow
  if(dabs(Eut).lt. 1.0d-36)Eut = 0.0d+00 
  if(dabs(Edut).lt. 1.0d-36)Edut = 0.0d+00 
  if(dabs(Ed2ut).lt. 1.0d-36)Ed2ut = 0.0d+00 
  if(dabs(Eut0).lt. 1.0d-36)Eut0 = 0.0d+00 
  if(dabs(Ett0).lt. 1.0d-36)Ett0 = 0.0d+00 
  if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
  if(dabs(sumi0).lt.1.0d-36)sumi0=0.0d+00
  if(dabs(sumi1).lt.1.0d-36)sumi1=0.0d+00
  if(dabs(sumi2).lt.1.0d-36)sumi2=0.0d+00
  if(dabs(ale).lt.1.0d-36)ale=0.0d+00
  if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
  if(dabs(gamma).lt.1.0d-36)gamma=0.0d+00

  ! get the derivatives of the eigenfunctions required
  ! for source excitation from the definition of stress. For
  ! completeness get the second derivative from the first
  ! derivatives and the equation of motion for the medium
  dc2dh(:) = dcdh(:)
  do i=1,mmax-1
    sums = sum(dc2dh(i+1:mmax)) 
    !print*,sum(dc2dh(i+1:mmax)),-sum(dc2dh(1:i))
    !sums = sums * 0.5
    dcdh(i) = sums
  enddo
  dcdh(mmax) = 0.0

  ! copy to input variables
  disp(:) = uu(:)
  stress(:) = uu(:)
  dc2db(:) = dcdb(:)
  dc2dr(:) = dcdr(:)
  dc2dh(:) = dcdh(:)
  cp = csph
  cg = usph

  call Dealloc_Model()

  if(nsph == 1) then 
    DEALLOCATE(vtp,dtp,rtp)
  endif

end subroutine slegn96

subroutine slegnpu(thk,vs,rhom,nlayer,&
                  t,cp,cg,disp,stress,&
                  t1,cp1,t2,cp2,dc2db,dc2dh,dc2dr,&
                  du2db,du2dh,du2dr,iflsph)bind(c,name='slegnpu_')
  use LoveWaveModel
  implicit none
  integer(c_int),value,INTENT(IN) ::nlayer,iflsph
  real(c_float),INTENT(IN)        :: thk(nlayer),vs(nlayer),rhom(nlayer)
  real(c_double),INTENT(IN)       :: t 
  real(c_double),INTENT(INOUT)    :: cp,cg 
  real(c_double),INTENT(INOUT)    :: disp(nlayer),stress(nlayer),&
                                    dc2db(nlayer),dc2dh(nlayer),dc2dr(nlayer)
  real(c_double),intent(inout)    :: du2db(nlayer),du2dh(nlayer),du2dr(nlayer)
  real(c_double),intent(inout)    :: t1,t2,cp1,cp2 ! period and phasev at slightly differentperiod

  !local
  integer(c_int)              :: i
  integer(c_int)              :: nsph
  real(c_double),PARAMETER    :: pi = atan(1.0) * 4.0
  real(c_double)              :: c, omega, wvno, twopi,sums
  real(c_double)              :: cg1,cg2,uc1
  real(c_double),ALLOCATABLE  :: dc2db1(:),dc2db2(:),&
                                dc2dh1(:),dc2dh2(:),dc2dr1(:),dc2dr2(:)
  real(c_double)              :: ar,tm,tm1
  real(c_double)              :: Eut, Edut, Eut0, Ett0, Ed2ut

  ! copy model
  call Alloc_Model(nlayer)
  zb(1:mmax) = vs(1:mmax)
  zrho(1:mmax) = rhom(1:mmax); zd(1:mmax) = thk(1:mmax)
  ALLOCATE(dc2db1(mmax),dc2db2(mmax),dc2dh1(mmax),&
          dc2dh2(mmax),dc2dr1(mmax),dc2dr2(mmax))

  ! check water layer
  do i=1,mmax
    if(zb(i).gt.0.0)then
      iwat(i) = 0
    else 
      iwat(i) = 1
    endif
  enddo

  ! earth flattening transformation
  nsph = iflsph
  if(nsph > 0) then 
    ALLOCATE(vtp(mmax),dtp(mmax),rtp(mmax))
    call bldsph()
  endif

  ! compute mu and lambda
  xmu(:) = zrho(:) * zb(:)**2

  ! compute group velocity for t and phase kernel
  twopi = 2.0 * pi
  omega = twopi/t
  c=cp
  wvno=omega/c
  call shfunc(omega,wvno)
  call energy(omega,wvno,Eut,Edut,Ed2ut,Eut0,Ett0,1,1)
  cg = ugr 
  dc2db(:) = dcdb(:); stress(:) = tt(:);
  dc2dr(:) = dcdr(:); disp(:) = uu(:)
  dc2dh(:) = dcdh(:);

  ! compute group velocity for t1 and phase kernel
  twopi = 2.0 * pi
  omega = twopi/t1
  c=cp1
  wvno=omega/c
  call shfunc(omega,wvno)
  call energy(omega,wvno,Eut,Edut,Ed2ut,Eut0,Ett0,1,1)
  cg1 = ugr 
  dc2db1(:) = dcdb(:)
  dc2dr1(:) = dcdr(:)
  dc2dh1(:) = dcdh(:)

  ! compute group velocity for t2 and phase kernel
  twopi = 2.0 * pi
  omega = twopi/t2
  c=cp2
  wvno=omega/c
  call shfunc(omega,wvno)
  call energy(omega,wvno,Eut,Edut,Ed2ut,Eut0,Ett0,1,1)
  cg2 = ugr 
  dc2db2(:) = dcdb(:)
  dc2dr2(:) = dcdr(:)
  dc2dh2(:) = dcdh(:)

  ! compute duf/dm
  uc1 = cg / cp
  du2db = uc1 * (2.0 - uc1) * dcdb - uc1**2 * t * (dc2db2 - dc2db1) / (t2-t1) 
  du2dr = uc1 * (2.0 - uc1) * dcdr - uc1**2 * t * (dc2dr2 - dc2dr1) / (t2-t1)
  du2dh = uc1 * (2.0 - uc1) * dcdh - uc1**2 * t * (dc2dh2 - dc2dh1) / (t2-t1)  

  ! earth flattening transform for kernel
  if(nsph >0)then
    ar = 6371.0 
    omega = twopi / t 
    tm = sqrt(1.+(3.0 * cp/(2.*ar*omega))**2)
    tm1 = (1.5/(ar * omega))**2 / tm

    ! convert duf/dm to dus/dm 
    du2db = (tm * du2db + cg * cp * dc2db * tm1) * vtp
    du2dr = (tm * du2dr + cg * cp * dc2dr * tm1) * rtp
    du2dh = (tm * du2dh + cg * cp * dc2dh * tm1) * dtp

    ! convert dcf/dm to dcs/dm
    dc2db = dc2db / tm**3 * vtp 
    dc2dr = dc2dr / tm**3 * rtp 
    dc2dh = dc2dh / tm**3 * dtp

    ! convert cp and cg to spherical one
    cp = cp  / tm 
    cg = cg * tm
  endif

  ! change dc/dz,du/dz to dc/dthk,du/dthk
  do i=1,mmax-1
    sums = sum(dc2dh(i+1:mmax))
    dc2dh(i) = sums
    sums = sum(du2dh(i+1:mmax))
    du2dh(i) = sums
  enddo
  dc2dh(mmax) = 0.0
  du2dh(mmax) = 0.0

  call Dealloc_Model()
  DEALLOCATE(dc2db1,dc2db2,dc2dh1,dc2dh2,dc2dr1,dc2dr2)

  if(nsph == 1) then 
    DEALLOCATE(vtp,dtp,rtp)
  endif

end subroutine slegnpu

end module LoveWaveKernel
