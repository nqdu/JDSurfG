!c----------------------------------------------------------------------c
!c                                                                      c
!c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
!c      VOLUME III                                                      c
!c                                                                      c
!c      PROGRAM: SREGN96                                                c
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

module RayleighWaveModel
use,intrinsic           :: iso_c_binding
implicit none

!! global variables
!============================================
! model param
integer(c_int)                          :: mmax
real(c_double),dimension(:),ALLOCATABLE :: zd,zrho,za,zb,xmu,xlam
real(c_double)                          :: uu0(4)

! eigen function
real(c_double),dimension(:),ALLOCATABLE ::ur,uz,tz,tr, &
                                           dcda,dcdb,dcdr,dcdh
integer(c_int),DIMENSION(:),ALLOCATABLE :: iwat

! sumi
real(c_double)                          ::sumi0,sumi1,sumi2,sumi3,&
                                          flagr,are,ugr
! spherical earth
real(c_double),ALLOCATABLE              :: vtp(:),dtp(:),rtp(:)

! haskel matrix
real(c_double),ALLOCATABLE              :: exe(:),exa(:),cd(:,:),vv(:,:)

! allfluid
logical(c_bool)                         :: allfluid

! hwat
real(c_double),ALLOCATABLE              :: wh(:,:),hex(:)

! save
real(c_double)                          :: dd(5,5),aa(4,4),ex1,ex2
real(c_double),ALLOCATABLE              :: ww(:),xx(:),yy(:),zz(:),&
                                           cospp(:),cosqq(:)

! engerw wra,wd,wba
real(c_double)                          :: wra,wd,wba

!common/emat/e,einv,ra,rb
complex(c_double_complex)               :: ra,rb
complex(c_double_complex)               :: e(4,4), einv(4,4)

! openmp wrapper
!$omp threadprivate(mmax,zd,zrho,za,zb,xmu,xlam,uu0)
!$omp threadprivate(ur,uz,tz,tr,dcda,dcdb,dcdr,dcdh)
!$omp threadprivate(iwat,sumi0,sumi1,sumi2,sumi3,flagr,are,ugr)
!$omp threadprivate(vtp,dtp,rtp,exe,exa,cd,vv)
!$omp threadprivate(allfluid,wh,hex,dd,aa,ex1,ex2)
!$omp threadprivate(ww,xx,yy,zz,cospp,cosqq,wra,wd,wba)
!$omp threadprivate(ra,rb,e,einv)

contains 
subroutine Alloc_Model(nlayer)
  implicit none
  integer(c_int)          :: nlayer 
  mmax = nlayer
  ALLOCATE(zd(mmax),zb(mmax),za(mmax),zrho(mmax),xmu(mmax),xlam(mmax))
  ALLOCATE(ur(mmax),uz(mmax),tz(mmax),tr(mmax))
  ALLOCATE(dcda(mmax),dcdb(mmax),dcdh(mmax),dcdr(mmax))
  ALLOCATE(iwat(mmax))
  ALLOCATE(exe(mmax),exa(mmax))
  ALLOCATE(cd(mmax,5),vv(mmax,4),wh(mmax,2),hex(mmax))
  ALLOCATE(ww(mmax),xx(mmax),yy(mmax),zz(mmax),cospp(mmax),&
          cosqq(mmax))
end subroutine Alloc_Model

subroutine Dealloc_Model()
  implicit none
  DEALLOCATE(zd,zrho,za,zb,xmu,xlam)
  DEALLOCATE(ur,uz,tz,tr,dcda,dcdb,dcdr,dcdh)
  DEALLOCATE(iwat)
  DEALLOCATE(exe,exa,cd,vv,wh,hex,ww,xx,yy,zz)
  DEALLOCATE(cospp,cosqq)
end subroutine Dealloc_Model

subroutine bldsph()
  !! Earth Flattening Transformation
  !! convert spherical model to flatten model
  !! and compute transfer coefficients
  !!c       Schwab, F. A., and L. Knopoff (1972). 
  !!c          Fast surface wave and free
  !!c       mode computations, in  
  !!c          Methods in Computational Physics, Volume 11,
  !!c       Seismology: Surface Waves and Earth Oscillations,  
  !!c          B. A. Bolt (ed),
  !!c       Academic Press, New York
  !!c
  !!c       Love Wave Equations  44, 45 , 41 pp 112-113
  !!c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
  !!c
  !!c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
  !!c     slowness instead of using average velocity for the layer
  !!c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
  !!c
  !!c       This is for Rayleigh Waves
  implicit none 
  real(c_double)        ::  ar, dr, r0, r1, z0, z1
  real(c_double)        :: tmp
  integer(c_int)        :: i
  ar=6371.0d0
  dr=0.0d0
  r0=ar
  dtp(:) = 0.0
  do i=1,mmax
    if(i==mmax) then 
      dr = dr + 1.0
    else 
      dr = dr + zd(i)
    endif
    r1 = ar - dr
    z0 = ar * log(ar/r0)
    z1 = ar * log(ar/r1)

    !use layer midpoint
    TMP=(2.0d0 *ar)/(r0+r1)
    vtp(i) = tmp
    rtp(i) = tmp**(-2.275)

    ! compute dtp
    dtp(i) = ar / r0

    ! earth flattening transformation
    za(i) = za(i) * tmp 
    zb(i) = zb(i) * tmp 
    zrho(i) = zrho(i) * rtp(i)
    zd(i) = z1 - z0
    r0=r1
  enddo
  zd(mmax) = 0.0
end subroutine bldsph 

end module RayleighWaveModel

module RayleighWaveKernel
use,intrinsic               :: iso_c_binding
implicit none

contains
subroutine svfunc(omega,wvno)
  !!c
  !!c     This combines the Haskell vector from sub down and
  !!c     Dunkin vector from sub up to form the eigenfunctions.
  !!c
  use RayleighWaveModel, only :  ur,uz,tz,tr,uu0,mmax,allfluid,&
                                 cd,vv,exe,exa,iwat
  implicit none
  real(c_double),INTENT(IN)   :: omega, wvno
  ! internal variables

  real(c_double)              :: fr,uu1, uu2, uu3, uu4, ext, f1213, fact
  real(c_double)              :: cd1, cd2, cd3, cd4, cd5, cd6
  real(c_double)              :: tz1, tz2, tz3, tz4
  integer(c_int)              :: i, jwat

  ! get compound matrix from bottom to top
  call up(omega,wvno,fr)
  call down(omega,wvno)
  
  ! get propagator from top to bottom 
  f1213 = -cd(1,2)
  ur(1) = cd(1,3)/cd(1,2)
  uz(1) = 1.0d+00
  tz(1) = 0.0d+00
  tr(1) = 0.0d+00
  uu0(1) = ur(1)
  uu0(2) = 1.0
  uu0(3) = fr
  uu0(4) = fr
  !c------
  !c      following
  !c                12       -1
  !c       Bk =    X|       Z
  !c                5-l,5-k  4l
  !c             -----------------
  !c                  12
  !c              - R |
  !c                  13
  !c-----
  !c       Compound Matrix (6x6)   Here (5x5)
  !c
  !c       X|12/12 ==      cd(1)
  !c       X|12/13 ==      cd(2)
  !!c       X|12/23 ==     -cd(3)
  !c       X|12/24 ==      cd(4)
  !c       X|12/34 ==      cd(5)
  !c       X|12/21 =  -X|12/12
  !c       X|12/ij =  -X|12/ji
  !c
  !c       where
  !c
  !c        |12
  !c       X|  = X|12/ij == X(1,i)*X(2,j) - X(1,j)*X(2,i)
  !c!        |ij
  
  !c   and using the shorthand
  !c       True [11 12 13 14 15 16 ] = Reduced [ 11 12 13 -13 14 15 ]
  !c     we have
  !c-----
  !c   In long hand with full notation and using the complex matrix short hand, e.g.,
  !C   that the denominator is -R12 and recalling R11 == 0
  !c
  !c
  !c    |    |         |                               | |      |
  !c    |    |         |           12      12      12  | |   -1 |
  !c    | B1 |         |   0      X|     -X|      X|   | |  Z   |
  !c    |    |         |           34      24      23  | |   41 |
  !c    |    |         |                               | |      |
  !c    |    |         |   12              12      12  | |   -1 |
  !c    | B2 |         | -X|       0      X|     -X|   | |  Z   |
  !c    |    |     1   |   34              14      13  | |   42 |
  !c    |    | = _____ |                               | |      |
  !c    |    |      12 |   12      12              12  | |   -1 |
  !c    | B3 |   -R|   |  X|     -X|       0      X|   | |  Z   |
  !c    |    |      13 |   24      14              12  | |   43 |
  !c    |    |         |   12      12      12          | |   -1 |
  !c    | B4 |         | -X|      X|     -X|       0   | |  Z   |
  !c    |    |         |   23      13      12          | |   44 |
  !c
  !c
  !c-----
  !c    Here the vv are the first column of the  Haskell down, and the
  !c    cd are the first row if the X up. At the surface R=X 
  !c------
  ext = 0.0
  do i=2,mmax
    cd1=  cd(i,1)
    cd2=  cd(i,2)
    cd3=  cd(i,3)
    cd4= -cd(i,3)
    cd5=  cd(i,4)
    cd6=  cd(i,5)
    tz1 = -vv(i,4)
    tz2 = -vv(i,3)
    tz3 =  vv(i,2)
    tz4 =  vv(i,1)

    !this is correct since the convention here is
    !that Uz = positive up and the z is positive up.
    !This also changes Tr
    uu1 =tz2*cd6 - tz3*cd5 + tz4*cd4
    uu2 = -tz1*cd6 + tz3*cd3 - tz4*cd2
    uu3 = tz1*cd5 - tz2*cd3 + tz4*cd1
    uu4 = -tz1*cd4 + tz2*cd2 - tz3*cd1
  
    ext=exa(i)+exe(i) -exe(1)
    if(ext.gt.-80.0 .and. ext.lt.80.0 ) then
      fact=dexp(ext)
      ur(i)=uu1*fact/f1213
      uz(i)=uu2*fact/f1213
      tz(i)=uu3*fact/f1213
      tr(i)=uu4*fact/f1213
    else
      ur(i) = 0.0
      uz(i) = 0.0
      tz(i) = 0.0
      tr(i) = 0.0
    endif
  enddo

  ! correction for fluid layers on top if not
  ! all fluid
  if(.not. allfluid)then
    jwat = 0
    do i=1,mmax
      if(iwat(i).gt.0)then
        jwat = i
      else
        exit
      endif
    enddo
    if(jwat.gt.0)then
      do i=1,jwat
        ur(i) = 0.0
        tr(i) = 0.0
      enddo
    endif
  endif
  
  !CRBHc-----
  !CRBHc       Continue KLUDGE if top layers are water
  !CRBHc       The problem here was that even though the 
  !CRBHc          phase velocity is determined,
  !CRBHc       the computed eigenfunctions do not match perfectly,
  !CRBHc       and we do not match the free surface requirement that Tz=0 at
  !CRBHc       top of fluid. So RBH used the fluid propagator from top
  !CRBHc       down to the solid-fluid interface, the Haskell from
  !CRBHc       base up to the same interface, and then adjust the amplitudes
  !CRBHc       to match. The purpose is to get a better estimate of the
  !CRBHc       eigenfunction in the fluid.
  !CRBHc
  !CRBHc       Ultimately use Kennett
  !CRBHc-----
  !CRBHc       first get the largest exponent
  !CRBHc-----
  !CRBH        if(iwat(1).eq.1)then
  !CRBHc-----
  !CRBHc           OK first layer is fluid
  !CRBHc-----
  !CRBH            hmax = 0.0
  !CRBH            jwat=0
  !CRBH            do 300 i=1,mmax-1
  !CRBH                if(iwat(i).eq.1)then
  !C!RBH                    jwat = i
  !CRBH                    hmax = hex(i+1)
  !CRBH                endif
  !CRBH  300       continue
  !CRBHc-----
  !CRBHc           now get correction factor
  !CRBHc-----
  !CRBH            corrd = sqrt(uz(jwat+1)**2 + tz(jwat+1)**2 )
  !CRBH            corrw = sqrt(wh(jwat+1,1)**2 + wh(jwat+1,2)**2 )
  !CRBH            if(corrd.eq.0.0)then
  !CRBHc-----
  !CRBHc       give precedence to surface values
  !CRBHc-----
  !CRBH                fac = 0.0
  !CRBH            else
  !CRBH                fac = corrw/corrd
  !CRBH            endif
  !CRBH            fac = fac * exp(hmax - hmax)
  !CRBH            fac = fac * 2.0
  !CRBHc-----
  !CRBHc           eventually normalize
  !C!RBHc-----
  !CRBH        IF(.not. ALLFLUID)then
  !CRBH            do 400 i=1,mmax
  !CRBH                if(i.le.jwat)then
  !CRBH                    efac = exp(hex(i) - hmax)
  !CRBH                    efac = efac * 2.0
  !CRBH                    ur(i) = 0.0
  !CRBH                    uz(i) = wh(i,1)*efac
  !CRBH                    tz(i) = wh(i,2)*efac
  !CRBH                    tr(i) = 0.0
  !CRBH                else
  !CRBH                    ur(i) = ur(i) * fac
  !CRBH                    uz(i) = uz(i) * fac
  !CRBH                    tr(i) = tr(i) * fac
  !CRBH                    tz(i) = tz(i) * fac
  !CRBH                endif
  !CRBH  400       continue
  !CRBH        ENDIF
  !CRBH        endif
          
  return;
end subroutine svfunc

subroutine up(omega,wvno,fr)
  !!c-----
  !!c       This finds the values of the Dunkin vectors at
  !1c       each layer boundaries from bottom layer upward.
  !!c!
  !!c       Note that these are the reduced vectors, e.g.,
  !!c       The true compound (G An-1 ... A1 H) is
  !!c       True [11 12 13 14 15 16 ] = Reduced [ 11 12 13 -13 14 15 ]
  !!c-----
  use RayleighWaveModel,        only : mmax,za,zb,zd,cd,exe,iwat,&
                                        zrho
  implicit none
  real(c_double),INTENT(IN)         :: omega, wvno
  real(c_double),intent(inout)      :: fr
  
  !The save labeled common is used to pass the
  !Dunkin 5x5 and Haskell 4x4 matrices to the main program
  !exp(ex1) is the scaling for the Dunkin-Thrower compound matrix
  !exp(ex2) is the scaling for the Haskell matrices

  ! local variables
  real(c_double)                    :: wvno2, om2,cr, ee(5), exn,&
                                      ca(5,5),xka, xkb, pex, svex
  real(c_double)                    :: cosp,rsinp,sinpr,cossv,& 
                                        rsinsv,sinsvr,exsum
  complex(c_double_complex)         ::rp, rsv, p, q,gbr(2,5)
  integer(c_int)                    :: i,j,m,mmm1,nmat

  ! initialize base vector
  ! set up starting values for bottom halfspace
  !c-----
  wvno2=wvno*wvno
  om2 = omega*omega
  call evalg(0,mmax,mmax-1,gbr,1,&
            wvno,omega,om2,wvno2)
  
  !note for surface waves, the ra, rb are real for the halfspace as long
  !as the phase velocity is less than the respective halfspace velocity
  
  cd(mmax,1)= dreal(gbr(1,1))
  cd(mmax,2)= dreal(gbr(1,2))
  cd(mmax,3)= dreal(gbr(1,3))
  cd(mmax,4)= dreal(gbr(1,4))
  cd(mmax,5)= dreal(gbr(1,5))
  exe(mmax) = 0.0d+00

  !matrix multiplication from bottom layer upward
  mmm1 = mmax-1
  exsum = 0.0
  do m = mmm1,1,-1
    xka = omega/za(m)
    if(zb(m).gt.0.0)then
      xkb = omega/zb(m)
    else
      xkb = 0.0
    endif
    rp =cdsqrt(dcmplx(wvno2-xka*xka,0.0d+00))
    rsv=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
    p = rp  * zd(m)
    q = rsv * zd(m)
    call varsv(p,q,rp, rsv,&
              cosp, cossv, rsinp, rsinsv,&
              sinpr, sinsvr, pex,svex,iwat(m),zd(m))
    call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,&
              sngl(zrho(m)),sngl(zb(m)),iwat(m),pex,pex+svex,&
              wvno,wvno2,om2)
  
    nmat = 5
    do i=1,nmat
      cr=0.0d+00
      do j=1,nmat
        cr=cr+cd(m+1,j)*ca(j,i)
      enddo
      ee(i)=cr
    enddo
    exn= 0.0d+00
    call normc(ee,exn,nmat)
    exsum = exsum + pex + svex + exn
    exe(m) = exsum
    do i = 1,nmat
      cd(m,i)=ee(i)
    enddo
  enddo
  !c-----
  ! define period equation
  !----
  fr=cd(1,1)
  return
end subroutine up

subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,&
                rho,b,iwat,ex,exa,wvno,wvno2,om2)
  implicit none
  real(c_double)        ::  ca(5,5), wvno, wvno2, om2
  real(c_float)         :: rho,b
  real(c_double)        :: cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
  integer(c_int)        :: iwat

  ! local
  real(c_float)         :: rho2
  real(c_double)        :: gam,ex, exa,a0,cpcq,cpy,&
                            cpz,cqw,cqx,xy,xz,wy,wz
  real(c_double)        :: gam2,gamm1,gamm2,a0c,xz2,wy2,temp,&
                            dfac,cqww2, cqxw2, g1wy2, &
                            gxz2, g2wy2, g2xz2,gg1, a0cgg1
  integer(c_int)        :: i,j
  !c-----
  !c        A11     A12     A13    -A13     A15     A16
  !c        A21     A22     A23    -A23     A25     A15
  !c        A31     A32     A33    1-A33   -A23    -A13
  !c       -A31    -A32    1-A33    A33     A23     A13
  !c        A51     A52    -A32     A32     A22     A12
  !c        A61     A51    -A31     A31     A21     A11
  !c-----
  !c       this will be multipled on the left by the G matrix
  !c
  !c       [ G11   G12 G13 -G13    G15 G16 ]
  !c
  !c-----
  !c       or on the right by
  !c
  !c       [ H11   H21 H31 -H31    H51 H61  ] ^T
  !c-----
  !c       the number of multiplications can be reduced from 36 to 25 
  !c       if we define a new matrices
  !c       related to the original matrices by
  !c-----
  !c         A11     A12     A13         A15     A16
  !c         A21     A22     A23         A25     A15
  !c        2 A31   2 A32   2 A33 -1   -2 A23  -2 A13
  !c         A51     A52    -A32         A22     A12
  !c         A61     A51    -A31         A21     A11
  !c-----
  !c
  !c       [ G11   G12  G13    G15 G16  ]
  !c       [ H11   H21 2 H31   H51 H61  ] ^T
  !c
  !c-----
  !c       this means that some of the original definitions of the 
  !c       Aij elements must be changed for the
  !c       definition of the modified 5x5 compount A matrix
  !
  !c       old 6x6                 new 5x5
  !c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
  !c       A53 = -A32              A43 = - (1/2) new A32
  !c       A63 = -A31              A53 = - (1/2) new A31
  !c-----
  !c       To recover the needed elements, 
  !c          we note that the old G14 = -old G14 = new G13
  !c-----
  !c-----
  if(iwat.eq.1)then
    !fluid layer
    do j=1,5
      do i=1,5
          ca(i,j) = 0.0d+00
      enddo
    enddo
    if(ex.gt.35.0d+00)then
        dfac = 0.0d+00
    else
        dfac = dexp(-ex)
    endif
    ca(3,3) = dfac
    ca(1,1) = cosp
    ca(5,5) = cosp
    ca(1,2) = - rsinp/(rho*om2)
    ca(2,1) = - rho*sinpr*om2
    ca(2,2) = cosp
    ca(4,4) = cosp
    ca(4,5) = ca(1,2)
    ca(5,4) = ca(2,1)
  else ! elastic layer
    if ( exa.lt. 60.0)then
        a0=dexp(-exa)
    else
        a0 = 0.0d+00
    endif
    cpcq = cosp * cossv
    cpy  = cosp*sinsvr
    cpz  = cosp*rsinsv
    cqw  = cossv*sinpr
    cqx  = cossv*rsinp
    xy   = rsinp*sinsvr
    xz   = rsinp*rsinsv
    wy   = sinpr*sinsvr  
    wz   = sinpr*rsinsv
    rho2= rho*rho
    gam = 2.0*b*b*wvno2/om2
    gam2  = gam*gam
    gamm1 = gam-1.
    gamm2 = gamm1*gamm1
    cqww2 = cqw * wvno2
    cqxw2 = cqx / wvno2
    gg1 = gam*gamm1
    a0c  = real(dcmplx(2.0d+00,0.0d+00)*&
          (dcmplx(a0,0.0d+00)-cpcq))
    xz2  = xz/wvno2
    gxz2 = gam*xz2
    g2xz2 = gam2 * xz2
    a0cgg1 = a0c*(gam+gamm1)
    wy2  = wy*wvno2
    g2wy2 = gamm2 * wy2
    g1wy2 = gamm1 * wy2

    !OK by symmetry
    temp = a0c*gg1 + g2xz2 + g2wy2
    ca(3,3) = a0 + temp + temp
    ca(1,1) = cpcq-temp
    ca(1,2) = (-cqx + wvno2*cpy)/(rho*om2)
    temp = dreal(dcmplx(0.5d+00,0.0d+00)*a0cgg1 + gxz2 + g1wy2)
    ca(1,3) = wvno*temp/(rho*om2)

    ca(1,4) = (-cqww2+cpz)/(rho*om2)
    temp = wvno2*(a0c + wy2) + xz
    ca(1,5) = -temp/(rho2*om2*om2)

    ca(2,1) = (-gamm2*cqw + gam2*cpz/wvno2)*rho*om2
    ca(2,2) = cpcq
    ca(2,3) = (gamm1*cqww2 - gam*cpz)/wvno
    ca(2,4) = -wz
    ca(2,5)=ca(1,4)


    temp =dreal(dcmplx(0.5d+00,0.0d+00)*a0cgg1*gg1 &
      + gam2*gxz2 + gamm2*g1wy2)
    ca(3,1) = real(-dcmplx(2.0d+00,0.0d+00)*temp*rho*om2/wvno)
    ca(3,2) = real(-wvno*(gam*cqxw2 - gamm1*cpy)* &
                  dcmplx(2.0d+00,0.0d+00))
    ca(3,4)=-2.0d+00*ca(2,3)
    ca(3,5)=-2.0d+00*ca(1,3)

    ca(4,1) = (-gam2*cqxw2 + gamm2*cpy)*rho*om2
    ca(4,2) = -xy
    ca(4,3)= -ca(3,2)/2.0d+00
    ca(4,4)=ca(2,2)
    ca(4,5)=ca(1,2)

    temp = gamm2*(a0c*gam2 + g2wy2) + gam2*g2xz2
    ca(5,1) = -rho2*om2*om2*temp/wvno2
    ca(5,2)=ca(4,1)
    ca(5,3)=-ca(3,1)/2.0d+00
    ca(5,4)=ca(2,1)
    ca(5,5)=ca(1,1)
  endif
  return
end subroutine dnka 

subroutine evalg(jbdry,m,m1,gbr,inp,&
                wvno,om,om2,wvno2)
  !
  !!c       this is from hspec96, but not everything is required
  !!c       tot he layered halfspace prroblem
  use RayleighWaveModel, only: iwat,za,zb,ra,rb,allfluid,&
                                e,einv,zrho
  implicit none
  integer(c_int)            ::jbdry, m, m1, inp
  complex(c_double_complex) :: gbr(2,5) 
  real(c_double)            :: wvno, wvno2, om, om2

  !internal parameters
  real(c_double)            :: xka,xkb,gam,gamm1
     
  ! set up halfspace conditions
  m1 = m1 * 1
  xka = om/za(m)
  if(zb(m).gt.0.0)then
    xkb = om/zb(m)
  else
    xkb = 0.0
  endif
  ra=CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
  rb=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
  gam = dble(zb(m))*wvno/om
  gam = 2.0d+00 * (gam * gam)
  gamm1 = real(gam - dcmplx(1.0d+00,0.0d+00))

  !set up halfspace boundary conditions
  !c
  !c       jbdry   = -1  RIGID
  !c           =  0  ELASTIC
  !c           = +1  FREE SURFACE
  !c
  if(jbdry.lt.0)then
    !RIGID - check properties of layer above
    if(zb(m) .gt. 0.0)then
      !ELASTIC ABOVE - RIGID
      gbr(inp,1) = dcmplx(1.0d+00,0.0d+00)
      gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
    else
      !FLUID ABOVE - RIGID
      gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
      if(allfluid)then
        gbr(inp,1) = dcmplx(1.0d+00,0.0d+00)
      else
        gbr(inp,4) = dcmplx(1.0d+00,0.0d+00)
      endif
      !(pseudo SH)
    endif
  else if(jbdry.eq.0)then

    !HALFSPACE
    if(iwat(m).eq.0)then
      !ELASTIC HALFSPACE
      !multiply G of Herrmann 2001 by - rho^2 om^4 k^2 ra rb
      !should have no effect since it is in both numerator and
      !denominator -- however will not give the correct potential
      !coefficient -- so rethink?
      E(1,1)= wvno
      E(1,2)= rb
      E(1,3)= wvno
      E(1,4)= -rb

      E(2,1)= ra
      E(2,2)= wvno
      E(2,3)= -ra
      E(2,4)= wvno

      E(3,1)= zrho(m)*om2*gamm1
      E(3,2)= zrho(m)*om2*gam*rb/wvno
      E(3,3)= zrho(m)*om2*gamm1
      E(3,4)= -zrho(m)*om2*gam*rb/wvno

      E(4,1)= zrho(m)*om2*gam*ra/wvno
      E(4,2)= zrho(m)*om2*gamm1
      E(4,3)= -zrho(m)*om2*gam*ra/wvno
      E(4,4)= zrho(m)*om2*gamm1

      EINV(1,1)= 0.5*gam/wvno
      EINV(1,2)= -0.5*gamm1/ra
      EINV(1,3)= -0.5/(zrho(m)*om2)
      EINV(1,4)= 0.5*wvno/(zrho(m)*om2*ra)

      EINV(2,1)= -0.5*gamm1/rb
      EINV(2,2)= 0.5*gam/wvno
      EINV(2,3)= 0.5*wvno/(zrho(m)*om2*rb)
      EINV(2,4)= -0.5/(zrho(m)*om2)

      EINV(3,1)= 0.5*gam/wvno
      EINV(3,2)=  0.5*gamm1/ra
      EINV(3,3)= -0.5/(zrho(m)*om2)
      EINV(3,4)= -0.5*wvno/(zrho(m)*om2*ra)

      EINV(4,1)= 0.5*gamm1/rb
      EINV(4,2)= 0.5*gam/wvno
      EINV(4,3)= -0.5*wvno/(zrho(m)*om2*rb)
      EINV(4,4)= -0.5/(zrho(m)*om2)

      gbr(inp,1)=dble(zrho(m)*zrho(m))*om2*om2*&
        (-gam*gam*ra*rb+wvno2*gamm1*gamm1)
      gbr(inp,2)=-dble(zrho(m))*(wvno2*ra)*om2
      gbr(inp,3)=-dble(zrho(m))*(-gam*ra*rb+wvno2*gamm1) &
        *om2*wvno
      gbr(inp,4)=dble(zrho(m))*(wvno2*rb)*om2
      gbr(inp,5)=wvno2*(wvno2-ra*rb)
      gbr(inp,1)=0.25*gbr(inp,1)/ &
                (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
      gbr(inp,2)=0.25*gbr(inp,2)/ &
                (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
      gbr(inp,3)=0.25*gbr(inp,3)/&
                (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
      gbr(inp,4)=0.25*gbr(inp,4)/&
                  (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
      gbr(inp,5)=0.25*gbr(inp,5)/&
                (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)

    else if(iwat(m).eq.1)then
      !FLUID HALFSPACE
      if(allfluid)then
        gbr(inp,1) = dble(0.5) / ra
        gbr(inp,2) = dcmplx(0.5d+00,0.0d+00)/ &
          (-dble(zrho(m))*om2)
        gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
        gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
        gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
      else
        gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
        gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
        gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
        gbr(inp,4) = dble(0.5*zrho(m)*om2) / ra
        gbr(inp,5) = dcmplx(-0.5d+00,0.0d+00)
      endif
      !for safety null the matrices and then fill with a 2x2
      E(1:4,1:4) = dcmplx(0.0,0.0)
      Einv(1:4,1:4) = dcmplx(0.0,0.0)
      E(1,1)=  ra
      E(1,2)= -ra
      E(2,1)= -zrho(m)*om2
      E(2,2)= -zrho(m)*om2

      EINV(1,1)= 0.5/ra
      EINV(1,2)= -0.5/(zrho(m)*om2)
      EINV(2,1)= -0.5/ra
      EINV(2,2)= -0.5/(zrho(m)*om2)
    endif
  else if(jbdry.eq.1)then
    !FREE - check properties of layer above
    if(zb(m) .gt. 0.0)then
      gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,5) = dcmplx(1.0d+00,0.0d+00)
    else
      gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
      gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
      if(allfluid)then
        gbr(inp,2) = dcmplx(1.0d+00,0.0d+00)
      else
        gbr(inp,5) = dcmplx(1.0d+00,0.0d+00)
      endif
    endif
  endif
  return
end subroutine evalg

subroutine varsv(p,q, rp, rsv, & 
                cosp, cosq, rsinp, rsinq, &
                sinpr, sinqr, pex,svex,iwat,zd)
  !!c       p = rp  * h
  !!c       q = rsv * h
  !!c       rp  vertical wave number for P
  !!c       rsv vertical wave number for SV
  !!c       cosp=cosh(p)  rsinp =rp *sinh(p)  sinpr = sinh(p)/rp
  !!c       cosq=cosh(q)  rsinsv=rsv*sinh(p)  sinpq = sinh(p)/rsv
  !!c
  !!c       for pure real or imaginary rp, and rsv, the
  !!c       functions returned are pure real
  implicit none
  complex(c_double_complex)   :: p,q,rp,rsv
  real(c_double)              :: cosp, cosq, rsinp, rsinq, &
                                 sinpr,sinqr,pex,svex,zd
  integer(c_int)              :: iwat
  real(c_double)              :: pr, pi, qr, qi,PFAC, SVFAC
  complex(c_double_complex)   :: epp, epm, eqp, eqm,sinp, sinq
     
  pex  = 0.0d+00
  svex = 0.0d+00
  pr = dreal(p)
  pi = dimag(p)
  qr = dreal(q)
  qi = dimag(q)
  pex   = pr
  if(iwat.eq.1)then
    !fluid layer
    epp = dcmplx(dcos(pi), dsin(pi))/2.0
    epm = dconjg(epp)
    if(pr.lt.30.) then
      pfac=dexp(-2.*pr)
    else
      pfac  = 0.0d+00
    endif
    cosp = dreal(epp + pfac*epm)
    sinp = epp - pfac*epm
    rsinp = dreal(rp *sinp)
    if(dabs(pr) .lt. 1.0e-5 .and. cdabs(rp).lt.1.0e-5)then
      sinpr = zd
    else
      sinpr = real(sinp/rp)
    endif
    cosq  = 1.0d+00
    rsinq = 0.0d+00
    sinqr = 0.0d+00
  else
    !elastic layer
    svex = qr
    epp = dcmplx(dcos(pi), dsin(pi))/2.0
    epm = dconjg(epp)
    eqp = dcmplx(dcos(qi), dsin(qi))/2.0
    eqm = dconjg(eqp)
    if(pr.lt.30.) then
        pfac=dexp(-2.*pr)
    else
        pfac  = 0.0d+00
    endif
    cosp = dreal(epp + pfac*epm)
    sinp = epp - pfac*epm
    rsinp = dreal(rp *sinp)
    if(dabs(pr) .lt. 1.0e-5 .and. cdabs(rp).lt.1.0e-5)then
      sinpr = zd
    else
      sinpr = real(sinp/rp)
    endif

    if(qr.lt.30.) then
      svfac=dexp(-2.*qr)
    else
      svfac  = 0.0d+00
    endif
    cosq = dreal(eqp + svfac*eqm)
    sinq = eqp - svfac*eqm
    rsinq = dreal(rsv*sinq)
    if(dabs(qr) .lt. 1.0e-5 .and. cdabs(rsv).lt.1.0e-5)then
      sinqr = zd
    else
      sinqr = real(sinq/rsv)
    endif

  endif
  return
end subroutine varsv

subroutine hska(AA,cosp,rsinp,sinpr,tcossv,trsinsv,tsinsvr, &
              rho,b,iwat,pex,svex,wvno,wvno2,om2)
  implicit none

  !subroutine variables
  real(c_double)        :: AA(4,4),cosp , tcossv, rsinp, trsinsv,&
                            sinpr, tsinsvr
  real(c_float)         :: rho,b
  integer(c_int)        :: iwat
  real(c_double)        :: pex, svex,om2, wvno, wvno2

  !local variables
  real(c_double)        :: gam, gamm1,cossv,rsinsv,sinsvr,dfac

  if(iwat.eq.1)then
    !fluid layer
    AA(1:4,1:4) = 0.0
    if(pex.gt.35.0d+00)then
        dfac = 0.0d+00
    else
        dfac = dexp(-pex)
    endif
    AA(1,1) = dfac
    AA(4,4) = dfac
    AA(2,2) = cosp
    AA(3,3) = cosp
    AA(2,3) = -rsinp/(rho*om2)
    AA(3,2) = - rho*om2*sinpr
  else
    !elastic layer
    !adjust the definitions of the SV trig functions
    !Note that the exp(pex) and exp(svex) have
    !been factored from the these functions, so that
    !cosp     = exp(pex) ( 1 + exp(-2 pex)/2 = exp(pex) cosp
    !True    local
    !cossv     = exp(svex) ( 1 + exp(-2 svex)/2 = exp(svex) cossv
    !True   local
    !This separation is fine for the combound matrices which have
    !terms such as cosp cosq, however for the Haskell propagator we need
    !cosp - cossv, or since pex >= svex
    !(cosp - cossv)    = exp(pex) [ cosp      - exp(svex-pex) cossv    ]
    !true                  local                       local
    !It is for this reason that we must modify the tcossv trsinsv rsinsvr in
    !the command line 
    if( (pex-svex) .gt. 70.0)then
        dfac = 0.0d+00
    else
        dfac = dexp(svex-pex)
    endif
    cossv = dfac * tcossv
    rsinsv = dfac * trsinsv
    sinsvr = dfac * tsinsvr
        
    gam = 2.0*b*b*wvno2/om2
    gamm1 = gam -1.0
    ! AA(1,1) =  gam*cosp - gamm1*cossv
    AA(1,1) =  cossv + gam*(cosp - cossv)
    AA(1,2) =  -wvno*gamm1*sinpr + gam*rsinsv/wvno
    AA(1,3) =  -wvno*(cosp-cossv)/(rho*om2)
    AA(1,4) =   (wvno2*sinpr - rsinsv)/(rho*om2)
    AA(2,1) =   gam*rsinp/wvno - wvno*gamm1*sinsvr
    AA(2,2) =   cosp - gam*(cosp- cossv)
    AA(2,3) =   ( -rsinp + wvno2*sinsvr)/(rho*om2)
    AA(2,4) = - AA(1,3)
    AA(3,1) =   rho*om2*gam*gamm1*(cosp-cossv)/wvno
    AA(3,2) = rho*om2*(-gamm1*gamm1*sinpr+gam*gam*rsinsv/wvno2)
    AA(3,3) =   AA(2,2)
    AA(3,4) = - AA(1,2)
    AA(4,1) =   rho*om2*(gam*gam*rsinp/wvno2-gamm1*gamm1*sinsvr)
    AA(4,2) = - AA(3,1)
    AA(4,3) = - AA(2,1)
    AA(4,4) =   AA(1,1)
  endif
  return
end subroutine hska

subroutine down(omega,wvno)
  !!This finds the values of the Haskell vectors at
  !!each layer boundaries from top layer downward.
  use RayleighWaveModel, only   : mmax,zd,za,zb,zrho,iwat,vv,&
                                  exa
  implicit none
  real(c_double)                :: omega,wvno

  !internal variables
  real(c_double)                :: wvno2, om2,aa(4,4),xka,xkb,pex,svex,&
                                    ex2,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,&
                                    cc,aa0(4),exsum
  complex(c_double_complex)     :: rp,rsv,p,q
  integer(c_int)                :: m,i,j
  om2 = omega*omega
  wvno2 = wvno*wvno
  !propagate down to get the first column of the
  !Am .. A2 A1
  !
  !initialize the top surface for the
  !first column of the Haskell propagator
  do i=1,4
    if(i.eq.1)then
      vv(1,i) = 1.0d+00
    else
      vv(1,i) = 0.0d+00
    endif
  enddo
  exa(1) = 0.0
        
  !matrix multiplication from top layer downward
  exsum  = 0.0
  do m= 1, mmax -1
    xka = omega/za(m)
    if(zb(m).gt.0.0)then
      xkb = omega/zb(m)
    else
      xkb = 0.0
    endif
    rp =CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
    rsv=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
    p = rp  * zd(m)
    q = rsv * zd(m)
    call varsv(p,q,rp, rsv,&
              cosp, cossv, rsinp, rsinsv,&
              sinpr, sinsvr, pex,svex,iwat(m),zd(m))

    call hska(AA,cosp,rsinp,sinpr,&
              cossv,rsinsv,sinsvr,&
              sngl(zrho(m)),sngl(zb(m)),iwat(m),&
              pex,svex,wvno,wvno2,om2)
    do i=1,4
      cc=0.0d+00
      do j=1,4
        cc=cc+aa(i,j)*vv(m,j)
      enddo
      aa0(i)=cc
    enddo
    ex2 = 0.0
    call normc(aa0,ex2,4)
    exsum = exsum + pex + ex2
    exa(m+1)= exsum
      
    do i=1,4
      vv(m+1,i)=aa0(i)
    enddo
  enddo

  !vv is  column 1 of the Haskell propagator at the top of layer m
  return
end subroutine down

subroutine energy(om,wvno,mmax)
  !!determine the energy integrals and also the
  !!phase velocity partials
  use RayleighWaveModel,only : sumi0,sumi1,sumi2,sumi3,iwat,&
                               dcda,dcdr,ur,tz,dcdb,flagr,are,&
                               ugr
  implicit none
  real(c_double)            :: om,wvno
  integer(c_int)            :: mmax

  !internal variables
  real(c_double)            :: om2, wvno2,TA, TC, TF, TL, TN
  real(c_double)            :: INT11, INT13, INT22, INT44, INT33, INT24
  real(c_double)            :: URUR, UZUZ, DURDUR, DUZDUZ, URDUZ, UZDUR
  real(c_double)            :: c,fac,facah, facav, facbh, facbv,facr
  
  integer(c_int)            :: m
  !coefficients of the ODE
  real(c_double)            :: a12, a14, a21, a23
  real(c_double)            :: ah, av, bh, bv, eta, rho
  
  integer(c_int)            :: TYPELYR
  
  !TYPELYR = -1  layer represents upper halfspace
  !           0  layer is true internal layer
  !           +1  layer represents lower halfspace
  !c-----

  !initialize
  sumi0 = 0.0
  sumi1 = 0.0
  sumi2 = 0.0
  sumi3 = 0.0
  c = om/wvno
  om2 = om*om
  wvno2 = wvno*wvno
  do m=1,mmax
    call getmat(m,wvno,om,a12, a14, a21, a23,&
                ah,av,bh,bv,eta,rho,&
                TA,TC,TF,TL,TN,iwat(m))

    ! get the integrals for layers over a halfspace
    ! this is here is we ever adopt the code to the coal 
    ! seam problem
    if(m.eq.mmax)then
      typelyr = 1
    else
      typelyr = 0
    endif
    INT11 = intijr(1,1,m,typelyr,om,om2,wvno,wvno2)
    INT13 = intijr(1,3,m,typelyr,om,om2,wvno,wvno2)
    INT22 = intijr(2,2,m,typelyr,om,om2,wvno,wvno2)
    INT24 = intijr(2,4,m,typelyr,om,om2,wvno,wvno2)
    INT33 = intijr(3,3,m,typelyr,om,om2,wvno,wvno2)
    INT44 = intijr(4,4,m,typelyr,om,om2,wvno,wvno2)
  
  
    if(iwat(m).eq.1)then
      !fluid - note these are for the 4x4 formulation
      URUR   = INT22*(wvno/(rho*om2))**2
      UZUZ   = INT11
      URDUZ  =  - (wvno/(rho*om2))*a12*INT22
      DUZDUZ = a12*a12*INT22
      sumi0  = sumi0 + rho*(URUR + UZUZ)
      sumi1  = sumi1 + TA*URUR
      sumi2  = sumi2 - TF*URDUZ
      sumi3  = sumi3 + TC*DUZDUZ
      facah = rho*ah*(URUR -2.*eta*URDUZ/wvno)
      facav = rho*av*DUZDUZ / wvno2
      dcda(m) = facah + facav
      facr = -0.5*c*c*(URUR + UZUZ)
      dcdr(m) = 0.5*(av*facav + ah*facah )/rho + facr
  
      ! define the ur in the fliud from the Tz - this will be at
      ! the top of the layer
      ur(m) = - wvno*tz(m)/(rho*om2)
  
    else
      !solid
      URUR   = INT11
      UZUZ   = INT22
      DURDUR = a12*a12*INT22 + 2.*a12*a14*INT24 + a14*a14*INT44
      DUZDUZ = a21*a21*INT11 + 2.*a21*a23*INT13 + a23*a23*INT33
      URDUZ  = a21*INT11 + a23*INT13
      UZDUR  = a12*INT22 + a14*INT24
      sumi0  = sumi0 + rho*(URUR + UZUZ)
      sumi1  = sumi1 + TL*UZUZ + TA*URUR
      sumi2  = sumi2 + TL*UZDUR - TF*URDUZ
      sumi3  = sumi3 + TL*DURDUR + TC*DUZDUZ
  
      !partial derivatives of phase velocity with
      !respect to medium parameters. Note that these are
      !later divided by (U sumi0) when these are finalized
      !note that the code distainguishes between
      !ah and av, bh and bv and uses eta. These are
      !actually TI parameters, and for isotropic media
      !ah = av, bh = bv aned eta = 1. The code is written this
      !way for easier conversion to a TI case
      facah = rho*ah*(URUR -2.*eta*URDUZ/wvno)
      facav = rho*av*DUZDUZ / wvno2
      facbh = 0.0
      facbv = rho*bv*(UZUZ + 2.*UZDUR/wvno + DURDUR/wvno2 + &
                 4.*eta*URDUZ/wvno )
      dcda(m) = facah + facav
      dcdb(m) = facbv + facbh
      !face = - TF*URDUZ/(wvno*eta)
      !dcdah(m) = facah
      !dcdav(m) = facav
      !dcdbh(m) = facbh
      !cdbv(m) = facbv
      !dcdeta(m) = faceta
      facr = -0.5*c*c*(URUR + UZUZ)
      !this is correct for TI
      dcdr(m) = 0.5*(av*facav + ah*facah + bv*facbv)/rho + facr
      !partial with layer thickness needs the value at the layer, e.g.,
      !dcdh(1) is top surface
    endif
  enddo

  ! determine final parameters
  flagr=om2*sumi0-wvno2*sumi1-2.d+00*wvno*sumi2-sumi3
  ugr=(wvno*sumi1+sumi2)/(om*sumi0)
  are=wvno/(2.d+00*om*ugr*sumi0)
  fac = are*c/wvno2

  ! use the final factors to apply the 1/(U I0) factor
  do m=1,mmax
    dcda(m) = dcda(m) /(ugr*sumi0)
    dcdb(m) = dcdb(m) /(ugr*sumi0)
    dcdr(m) = dcdr(m) /(ugr*sumi0)
  enddo

  !determine the dcdh
  call getdcdh(om2,wvno,wvno2,fac)
  return
end subroutine energy

function intijr(i,j,m,typelyr,om,om2,wvno,wvno2)
  use RayleighWaveModel, only: iwat,tz,ra,e,tr,zd,einv,&
                                uz,ur,rb
  implicit none
  real(c_double)            :: intijr
  integer(c_int)            :: i,j,m,TYPELYR
  real(c_double)            :: om, om2, wvno,wvno2
  !TYPELYR = -1  layer represents upper halfspace
  !           0  layer is true intenal layer
  !           +1  layer represents lower halfspace
  !c NOTE DO NOT PERMIT RA RB to be ZERO add a small number TAKEN CARE OF
  !c IN FUNC CALLS
  !c beware of E matrix for fluid - internally I use a 4x4
  !c but the book uses a 2x2 in the 8-10 problem
  !c need typelyr, wvno, om, om2, wvno2, m, mmax, do not need medium
  !c-----
  !c       Potential coefficients
  !c       kmpu, kmsu are the upward coefficients at bottom of layer
  !c       km1pd, km1sd are the downward coefficients at top of layer
  !c-----
  !internal variables
  complex(c_double_complex)       :: cintijr, FA, GA, FB, GB, H1, H2
  complex(c_double_complex)       :: gbr(2,5),kmpu, kmsu, km1pd, km1sd

  !call evalg to get the E and EINV matrices
  !evalg knows about water - we ignore everything else

  call evalg(0,m,m-1,gbr,1,&
            wvno,om,om2,wvno2)
  !for an elastic solid
  !T        PU  SvU  PD  SvD T
  ! [Ur, Uz, Tr, Tz]  = E [ K  ,K   ,K  ,K   ]
  ! and
  !PU  SvU  PD  SvD T   -1                 T
  ![ K  ,K   ,K  ,K   ] = E   [Ur, Uz, Tr, Tz]
  !The text used the notation K m-1 to represent the potentials
  !at the top of a layer and K m those at the bottom of the layer
  !Because of Fortran indexing UZ(1) is the Z displacement at the
  !top of layer 1 and UZ(2) is the value at the bottom.  This is
  !reason for the slight incondistency in presentation below
  !    

  if(iwat(m).eq.1)then
    !fluid
    if(typelyr .lt. 0)then
        kmpu = einv(1,1)*uz(m+1) &
              + einv(1,2)*tz(m+1) 
        cintijr = e(i,1)*e(j,1)*kmpu*kmpu/(2.0*ra)
    else if(typelyr .eq. 0)then
      km1pd = einv(2,1)*uz(m) &
              + einv(2,2)*tz(m) 
      kmpu = einv(1,1)*uz(m+1) &
              + einv(1,2)*tz(m+1)  
      FA=ffunc(ra,zd(m))
      GA=gfunc(ra,zd(m))
      cintijr = e(i,1)*e(j,1)*kmpu*kmpu*FA &
              +(e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*km1pd*GA &  
              + e(i,2)*e(j,2)*km1pd*km1pd*FA
    else
      km1pd = einv(2,1)*uz(m) &
              + einv(2,2)*tz(m)
      cintijr = e(i,2)*e(j,2)*km1pd*km1pd/(2.0*ra)
    endif
  else
    !solid
    if(typelyr .lt. 0)then
      !potential coefficients for upper halfspace
      !based on the displacement, stress at top
      !Basically m = 1
      kmpu = einv(1,1)*ur(m) + einv(1,2)*uz(m) &
             + einv(1,3)*tz(m) + einv(1,4)*tr(m)
      kmsu = einv(2,1)*ur(m) + einv(2,2)*uz(m) &
             + einv(2,3)*tz(m) + einv(2,4)*tr(m)
      cintijr = e(i,1)*e(j,1)*kmpu*kmpu/(2.0*ra) &
                + (e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*kmsu/(ra+rb) &
                + e(i,2)*e(j,2)*kmsu*kmsu/(2.0*rb)
    else if(typelyr .eq. 0)then
      !downward potentials coefficients at top of layer m
      km1pd = einv(3,1)*ur(m) + einv(3,2)*uz(m) &
              + einv(3,3)*tz(m) + einv(3,4)*tr(m)
      km1sd = einv(4,1)*ur(m) + einv(4,2)*uz(m) &
              + einv(4,3)*tz(m) + einv(4,4)*tr(m)
     !upward potentials coefficients at bottom of layer m
      kmpu = einv(1,1)*ur(m+1) + einv(1,2)*uz(m+1) &
             + einv(1,3)*tz(m+1) + einv(1,4)*tr(m+1)
      kmsu = einv(2,1)*ur(m+1) + einv(2,2)*uz(m+1) &
             + einv(2,3)*tz(m+1) + einv(2,4)*tr(m+1)
      FA=ffunc(ra,zd(m))
      GA=gfunc(ra,zd(m))
      FB=ffunc(rb,zd(m))
      GB=gfunc(rb,zd(m))
      H1=h1func(ra,rb,zd(m))
      H2=h2func(ra,rb,zd(m))
      cintijr = e(i,1)*e(j,1)*kmpu*kmpu*FA &
                + e(i,3)*e(j,3)*km1pd*km1pd*FA &
                + e(i,2)*e(j,2)*kmsu*kmsu*FB &
                + e(i,4)*e(j,4)*km1sd*km1sd*FB &
                + H1*((e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*kmsu + &
                (e(i,3)*e(j,4)+e(i,4)*e(j,3))*km1pd*km1sd) &
                + H2*((e(i,1)*e(j,4)+e(i,4)*e(j,1))*kmpu*km1sd + &
                (e(i,2)*e(j,3)+e(i,3)*e(j,2))*km1pd*kmsu) &
                + GA*(e(i,1)*e(j,3)+e(i,3)*e(j,1))*kmpu*km1pd &
                + GB*(e(i,2)*e(j,4)+e(i,4)*e(j,2))*kmsu*km1sd 
    else
      !downward potential coefficients for lower halfspace
      !based on the displacement, stress at bottom
      !Basically m = mmax
      km1pd = einv(3,1)*ur(m)+einv(3,2)*uz(m) &
             +einv(3,3)*tz(m)+einv(3,4)*tr(m)
      km1sd = einv(4,1)*ur(m)+einv(4,2)*uz(m) &
              +einv(4,3)*tz(m)+einv(4,4)*tr(m)
      cintijr = e(i,3)*e(j,3)*km1pd*km1pd/(2.0*ra)  &
                +(e(i,3)*e(j,4)+e(i,4)*e(j,3))*km1pd*km1sd/(ra+rb) &
                +e(i,4)*e(j,4)*km1sd*km1sd/(2.0*rb)
    endif
  endif

  intijr = dreal(cintijr)
  return
end function intijr

function ffunc(nub,dm)
  implicit none
  complex(c_double_complex)       :: ffunc,nub
  real(c_double)                  :: dm
  complex(c_double_complex)       :: exqq,argcd

  !get the f function
  if(cdabs(nub).lt. 1.0d-08)then
    ffunc = dm
  else
    argcd = nub*dm
    if(dreal(argcd).lt.40.0)then
          exqq = cdexp(-2.0d+00*argcd)
    else
          exqq=0.0d+00
    endif
    ffunc = (1.d+00-exqq)/(2.d+00*nub)
  endif
  return
end function ffunc

function gfunc(nub,dm)
  implicit none
  complex(c_double_complex)   ::  gfunc,nub,argcd
  real(c_double)              :: dm
  argcd = nub*dm
  if(dreal(argcd).lt.75)then
    gfunc = cdexp(-argcd)*dm
  else
    gfunc = dcmplx(0.0d+00, 0.0d+00)
  endif
  return
end function gfunc

function h1func(nua,nub,dm)
  implicit none
  complex(c_double_complex)       :: h1func,nua, nub
  real(c_double)                  :: dm
  complex(c_double_complex)       :: argcd,exqq
  if(cdabs(nub+nua).lt. 1.0d-08)then
    h1func = dm
  else
    argcd = (nua+nub)*dm
    if(dreal(argcd).lt.40.0)then
      exqq = cdexp(-argcd)
    else
      exqq=0.0d+00
    endif
    h1func = (1.d+00-exqq)/(nub+nua)
  endif
  return
end function h1func

function h2func(nua,nub,dm)
  implicit none
  complex(c_double_complex)       :: h2func,nua,nub
  real(c_double)                  :: dm 
  complex(c_double_complex)       :: argcd,exqq,exqp

  if(cdabs(nub-nua).lt. 1.0d-08)then
    !this should never occur for surface waves
    h2func =  dm
  else
    argcd = nua*dm
    if(dreal(argcd).lt.40.0)then
          exqp = cdexp(-argcd)
    else
          exqp=0.0d+00
    endif
    argcd = nub*dm
    if(dreal(argcd).lt.40.0)then
          exqq = cdexp(-argcd)
    else
          exqq=0.0d+00
    endif
    h2func = (exqq - exqp)/(nua-nub)
  endif
  return
end function h2func

subroutine normc(ee,ex,nmat)
  !! This routine is an important step to control over- or
  !! underflow.
  !! The Haskell or Dunkin vectors are normalized before
  !! the layer matrix stacking.
  !! Note that some precision will be lost during normalization.
  !!c
  implicit none
  integer(c_int)          :: nmat
  real(c_double)          :: ex,ee(nmat)
  
  !local
  real(c_double)          :: t1,t2
  INTEGER(c_int)          :: i
  ex = 0.0d+00
  t1 = 0.0d+00
  do i = 1,nmat
    if(dabs(ee(i)).gt.t1) t1 = dabs(ee(i))
  enddo
  if(t1.lt.1.d-40) t1=1.d+00
  do i =1,nmat
    t2=ee(i)
    t2=t2/t1
    ee(i)=t2
  enddo
  
  ! store the normalization factor in exponential form.
  ex=dlog(t1)
  return
end subroutine normc

subroutine getdcdh(om2,wvno,wvno2,fac)
  use RayleighWaveModel,only     :mmax,xmu,ur,uz,tr,tz,&
                                  iwat,zrho,xlam,dcdh
  implicit none
  real(c_double)                ::om2,wvno,wvno2,fac

  !internal variables
  real(c_double)                :: tur, tuz, ttz, ttr
  real(c_double)                :: dfac, gfac1, gfac2, gfac3, &
                                    gfac4, gfac5, gfac6
  real(c_double)                :: drho, dmu, dlm
  real(c_double)                :: duzdzp, dlur2, xl2mp, xl2mm
  real(c_double)                :: dl2mu, duzdzm, drur2
  real(c_double)                :: URB, DURDZM, DURDZP
  integer(c_int)                :: m

  do m=1,mmax
    tuz = uz(m)
    ttz = tz(m)
    ttr = tr(m)
    if(iwat(m).eq.1)then
      tur = -wvno*ttz/(zrho(m)*om2)
    else
      tur = ur(m)
    endif

    !this assumes that the top is a halfspace
    if(m.eq.1)then
      drho = zrho(1) - 0.0
      dmu  = xmu(1) - 0.0
      dlm = xlam(1) - 0.0
      dl2mu = dlm + dmu + dmu
      xl2mp = xlam(m) + xmu(m) + xmu(m)
      duzdzp = (ttz + wvno*xlam(m)*tur)/xl2mp
      if(iwat(m) .eq.1)then
        durdzp = wvno*tuz
      else
        durdzp = (ttr/xmu(m)) - wvno*tuz
      endif
      drur2 = tur*tur*drho
      dlur2 = tur*tur*dl2mu

      gfac1 =  om2*drho*tuz**2
      gfac2 =  om2*drur2
      gfac3 = -wvno2*dmu*tuz**2
      gfac4 = -wvno2*dlur2
      gfac5 =  (xl2mp*duzdzp**2)
      gfac6 =  (xmu(m)*durdzp**2 )
    else
      drho = zrho(m) - zrho(m-1)
      dmu = xmu(m) - xmu(m-1)
      dlm = xlam(m) - xlam(m-1)
      dl2mu = dlm + dmu + dmu
      xl2mp = xlam(m)   + xmu(m)   + xmu(m)
      xl2mm = xlam(m-1) + xmu(m-1) + xmu(m-1)
      duzdzp = (ttz + wvno*xlam(m)*tur)/xl2mp
      if(xmu(m).eq.0.0)then
          durdzp = wvno*tuz
      else
          durdzp = (ttr/xmu(m)) - wvno*tuz
      endif
      if(xmu(m-1).eq.0.0 )then
          durdzm = wvno*tuz
      else
          durdzm = (ttr/xmu(m-1)) - wvno*tuz
      endif
      !attempt to fix for water layer, since Ur is not continuous
      !across fluid - solid interface or fluid-fluid interface
      if(iwat(m-1).eq.1 .and. iwat(m).eq.0 )then
          URB = -wvno*tz(m)/(zrho(m-1)*om2)
          drur2 = tur*tur*zrho(m)-URB*URB*zrho(m-1)
          dlur2 = tur*tur*xl2mp -URB*URB*xl2mm
          duzdzm = (ttz + wvno*xlam(m-1)*URB)/xlam(m-1)
      else if(iwat(m-1).eq.1 .and. iwat(m).eq.1 )then
          URB = -wvno*tz(m)/(zrho(m-1)*om2)
          drur2 = tur*tur*zrho(m)- URB*URB*zrho(m-1)
          dlur2 = tur*tur*xl2mp  - URB*URB*xl2mm
          duzdzm = (ttz + wvno*xlam(m-1)*URB)/xl2mm
      else
          drur2 = tur*tur*drho
          dlur2 = tur*tur*dl2mu
          duzdzm = (ttz + wvno*xlam(m-1)*tur)/xl2mm
      endif
      gfac1 =  om2*drho*tuz**2
      gfac2 =  om2*drur2
      gfac3 = -wvno2*dmu*tuz**2
      gfac4 = -wvno2*dlur2
      gfac5 =  (xl2mp*duzdzp**2-xl2mm*duzdzm**2)
      gfac6 =  (xmu(m)*durdzp**2 - xmu(m-1)*durdzm**2)
    endif
    dfac = fac * (gfac1 + gfac2 + gfac3 + gfac4 &
                  + gfac5 + gfac6 )
    if(dabs(dfac).lt.1.0d-38)then
        dfac = 0.0d+00
    endif

    dcdh(m) = dfac
  enddo
  return
end subroutine getdcdh

subroutine getmat(m,wvno,om,a12, a14, a21, a23, &
                   ah,av,bh,bv,eta,rho,&
                  TA,TC,TF,TL,TN,iwat)
     
  !!get matrix elements of the ODE
  !! We do not require all eight non-zero elements
  !!       also return the model characterization parameters
  use RayleighWaveModel,only : za,zrho,zb
  implicit none
  real(c_double)        :: wvno,om,a12, a14, a21, a23
  real(c_double)        :: ah,av,bh,bv,eta,rho
  real(c_double)        :: TA,TC,TF,TL,TN
  integer(c_int)        :: m, iwat

  if(iwat.eq.1)then
    !non-gravitating fluid - using 2x2 ODE
    ah= za(m)
    av= za(m)
    bh= 0.0
    bv= 0.0
    rho=zrho(m)
    eta = 1.0

    TL = 0.0
    TN = 0.0
    TC = zrho(m)*za(m)*za(m)
    TA = zrho(m)*za(m)*za(m)
    TF = TA - 2.*TN
    a12 = - ( wvno*wvno - om*om/(ah*ah))/(rho*om*om)
  else
    !elastic - using 4x4 ODE
    ah= za(m)
    av= za(m)
    bh= zb(m)
    bv= zb(m)
    rho=zrho(m)
    eta = 1.0


    TL = zrho(m)*zb(m)*zb(m)
    TN = zrho(m)*zb(m)*zb(m)
    TC = zrho(m)*za(m)*za(m)
    TA = zrho(m)*za(m)*za(m)
    TF = TA - 2.*TN

    a12 = -wvno
    a14 = 1.0/TL
    a21 = wvno * TF/TC
    a23 = 1.0/TC
  endif

  return
end subroutine getmat

subroutine sprayl(om,c,mmax,csph,usph,ugr)
  !!Transform spherical earth to flat earth
  !!and relate the corresponding flat earth dispersion to spherical
  
  !!Schwab, F. A., and L. Knopoff (1972). 
  !!Fast surface wave and free mode computations, in Methods in Computational Physics, 
  !!Volume 11,  Seismology: Surface Waves and Earth Oscillations,  
  !!B. A. Bolt (ed),
  !!Academic Press, New York
  !!Rayleigh Wave Equations 111, 114 p 144
  !!Partial with respect to parameter uses the relation
  !!For phase velocity, for example,
  !!c
  !!dcs/dps = dcs/dpf * dpf/dps, c is phase velocity, p is
  !!parameter, and f is flat model
  !!
  !!om      R*8     angular frequency
  !!c       R*8     phase velocity
  !!mmax    I*4     number of layers 
  !!
  use RayleighWaveModel,only:vtp,dtp,dcda,dcdb,&
                              dcdh,dcdr,rtp
  implicit none
  real(c_double)           :: om,c,csph,usph,ugr
  integer(c_int)           :: mmax

  !local arguments
  real(c_double)           :: ar,tm
  integer(c_int)           :: i
  
  ar = 6371.0d0
  tm=sqrt(1.+(c/(2.*ar*om))**2)
  do i=1,mmax
    dcda(i)=dcda(i)*  vtp(i) / (tm**3)
    dcdb(i)=dcdb(i)*  vtp(i) / (tm**3)
    dcdh(i)=dcdh(i) * dtp(i) / (tm**3)
    dcdr(i)=dcdr(i)*  rtp(i) / (tm**3)
  enddo

  ! convert cf and uf to cs and us
  usph = ugr*tm
  csph = c/tm

  return
end subroutine sprayl

subroutine sregn96(thk,vp,vs,rhom,nlayer,&
                  t,cp,cg,dispu,dispw,stressu,stressw,&
                  dc2da,dc2db,dc2dh,dc2dr,iflsph)bind(c,name='sregn96_')
  use RayleighWaveModel
  implicit none
  integer(c_int),value        :: nlayer,iflsph
  real(c_float),INTENT(IN)    :: thk(nlayer),vp(nlayer),rhom(nlayer),&
                                  vs(nlayer)
  real(c_double),intent(inout):: dispu(nlayer),dispw(nlayer),&
                                stressu(nlayer),stressw(nlayer),&
                                  dc2da(nlayer),dc2db(nlayer),dc2dh(nlayer),&
                                  dc2dr(nlayer)
  real(c_double)              :: t,cp,cg

  !local
  integer(c_int)              :: i
  integer(c_int)              :: nsph
  real(c_double),PARAMETER    :: pi = 3.1415926535898
  real(c_double)              :: c, omega, wvno, gammar, csph, usph,twopi,sums

  ! copy model
  call Alloc_Model(nlayer)
  zb(1:mmax) = vs(1:mmax); za(1:mmax) = vp(1:mmax)
  zrho(1:mmax) = rhom(1:mmax); zd(1:mmax) = thk(1:mmax)    

  ! check water layer
  allfluid = .true.
  do i=1,mmax
    if(zb(i).gt.0.0)then
      allfluid = .false.
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
  xlam(:) = zrho(:) * za(:)**2 -2 * xmu(:)

  ! main part
  twopi = 2.0 * pi
  omega = twopi/t
  c=cp; ugr = c 
  wvno=omega/c
  call svfunc(omega,wvno)
  call energy(omega,wvno,mmax)

  ! no attenuation
  gammar = 0.0d+00

  !also check for possible conversion errors in IEEE
  !conversion from double precision to single precision
  if(dabs(uu0(1)).lt.1.0d-36)uu0(1)=0.0d+00
  if(dabs(uu0(3)).lt.1.0d-36)uu0(3)=0.0d+00
  if(dabs(c).lt.1.0d-36)c=0.0d+00
  if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
  if(dabs(are).lt.1.0d-36)are=0.0d+00
  if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
  if(dabs(gammar).lt.1.0d-36)gammar=0.0d+00

  !earth flattening transform for kernel
  if(nsph >0)then
    !sphericity correction for partial derivatives 
    !of the original model
    call sprayl(omega,c,mmax,csph,usph,ugr)
    wvno = omega / csph
  else 
    csph = c 
    usph = ugr 
  endif

  !up to this point the dcdh are changes to phase velocity if
  !if the layer boundary changes. Here we change this to mean
  !the dc/dh for a change in layer thickness
  !A layer becomes thicker if the base increases and the top
  !decreases its position. The dcdh to this point indicates 
  !the effect of moving a boundary down. Now we convert to
  !the effect of changing a layer thickness.
  do i=1,mmax-1
    sums = sum(dcdh(i+1:mmax))
    dcdh(i) = sums
  enddo
  dcdh(mmax) = 0.0

  ! copy to input variables
  dc2da(:) = dcda(:);dc2db(:) = dcdb(:);
  dc2dr(:) = dcdr(:); dc2dh(:) = dcdh(:)
  dispu(:) = ur(:); dispw(:) = uz(:)
  stressu(:) = tr(:); stressw(:) = tz(:)
  cp = csph; cg = usph
  call Dealloc_Model()

  if(nsph == 1) then 
    DEALLOCATE(vtp,dtp,rtp)
  endif

end subroutine sregn96

subroutine sregnpu(thk,vp,vs,rhom,nlayer,&
                  t,cp,cg,dispu,dispw,stressu,stressw,&
                  t1,cp1,t2,cp2,&
                  dc2da,dc2db,dc2dh,dc2dr,&
                  du2da,du2db,du2dh,du2dr,&
                  iflsph)bind(c,name='sregnpu_')
  use RayleighWaveModel
  implicit none
  integer(c_int),value        :: nlayer,iflsph
  real(c_float),INTENT(IN)    :: thk(nlayer),vp(nlayer),rhom(nlayer),&
                                  vs(nlayer)
  real(c_double),intent(inout):: dispu(nlayer),dispw(nlayer),&
                                  stressu(nlayer),stressw(nlayer)
  real(c_double),intent(inout):: dc2da(nlayer),dc2db(nlayer),dc2dh(nlayer),dc2dr(nlayer)
  real(c_double),intent(inout):: du2da(nlayer),du2db(nlayer),du2dh(nlayer),du2dr(nlayer)
  real(c_double)              :: t,cp,cg
  real(c_double),intent(inout):: t1,t2,cp1,cp2 ! period and phasev at slightly differentperiod

  !local
  integer(c_int)              :: i
  integer(c_int)              :: nsph
  real(c_double),PARAMETER    :: pi = atan(1.0) * 4.0
  real(c_double)              :: c, omega, wvno,twopi,sums
  real(c_double)              :: cg1,cg2,uc1
  real(c_double),ALLOCATABLE  :: dc2da1(:),dc2da2(:),dc2db1(:),dc2db2(:),&
                                dc2dh1(:),dc2dh2(:),dc2dr1(:),dc2dr2(:)
  real(c_double)              :: ar,tm,tm1
  
  ! copy model
  call Alloc_Model(nlayer)
  zb(1:mmax) = vs(1:mmax); za(1:mmax) = vp(1:mmax)
  zrho(1:mmax) = rhom(1:mmax); zd(1:mmax) = thk(1:mmax)
  ALLOCATE(dc2da1(mmax),dc2da2(mmax),dc2db1(mmax),dc2db2(mmax),&
            dc2dh1(mmax),dc2dh2(mmax),dc2dr1(mmax),dc2dr2(mmax))

  ! check water layer
  allfluid = .true.
  do i=1,mmax
    if(zb(i).gt.0.0)then
      allfluid = .false.
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
  xlam(:) = zrho(:) * za(:)**2 -2 * xmu(:)

  ! compute group velocity for t and phase kernel
  twopi = 2.0 * pi
  omega = twopi/t
  c=cp
  wvno=omega/c
  call svfunc(omega,wvno)
  call energy(omega,wvno,mmax)
  cg = ugr 
  dc2da(:) = dcda(:);dc2db(:) = dcdb(:);
  dc2dr(:) = dcdr(:); dc2dh(:) = dcdh(:)
  dispu(:) = ur(:); dispw(:) = uz(:);
  stressu(:) = tr(:); stressw(:) = tz(:)

  ! compute group velocity for t1 and phase kernel
  twopi = 2.0 * pi
  omega = twopi/t1
  c=cp1
  wvno=omega/c
  call svfunc(omega,wvno)
  call energy(omega,wvno,mmax)
  cg1 = ugr 
  dc2da1(:) = dcda(:);dc2db1(:) = dcdb(:);
  dc2dr1(:) = dcdr(:); dc2dh1(:) = dcdh(:)

  ! compute group velocity for t2 and phase kernel
  twopi = 2.0 * pi
  omega = twopi/t2
  c=cp2
  wvno=omega/c
  call svfunc(omega,wvno)
  call energy(omega,wvno,mmax)
  cg2 = ugr
  dc2da2(:) = dcda(:);dc2db2(:) = dcdb(:);
  dc2dr2(:) = dcdr(:); dc2dh2(:) = dcdh(:)

  ! compute duf/dm
  uc1 = cg / cp
  du2da = uc1 * (2.0 - uc1) * dcda - uc1**2 * t * (dc2da2 - dc2da1) / (t2-t1) 
  du2db = uc1 * (2.0 - uc1) * dcdb - uc1**2 * t * (dc2db2 - dc2db1) / (t2-t1) 
  du2dr = uc1 * (2.0 - uc1) * dcdr - uc1**2 * t * (dc2dr2 - dc2dr1) / (t2-t1)
  du2dh = uc1 * (2.0 - uc1) * dcdh - uc1**2 * t * (dc2dh2 - dc2dh1) / (t2-t1)  

  ! earth flattening transform for kernel
  if(nsph >0)then
    ar = 6371.0 
    omega = twopi / t 
    tm = sqrt(1.+(cp/(2.*ar*omega))**2)
    tm1 = (0.5/(ar * omega))**2 / tm

    ! convert duf/dm to dus/dm 
    du2da = (tm * du2da + cg * cp * dc2da * tm1) * vtp
    du2db = (tm * du2db + cg * cp * dc2db * tm1) * vtp
    du2dr = (tm * du2dr + cg * cp * dc2dr * tm1) * rtp
    du2dh = (tm * du2dh + cg * cp * dc2dh * tm1) * dtp

    ! convert dcf/dm to dcs/dm
    dc2da = dc2da / tm**3 * vtp 
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
  DEALLOCATE(dc2da1,dc2da2,dc2db1,dc2db2,&
            dc2dh1,dc2dh2,dc2dr1,dc2dr2)

  if(nsph == 1) then 
    DEALLOCATE(vtp,dtp,rtp)
  endif

end subroutine sregnpu

end module RayleighWaveKernel