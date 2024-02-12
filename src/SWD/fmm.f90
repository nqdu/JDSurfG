!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module. 
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
  use,intrinsic :: iso_c_binding
  use omp_lib
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp=c_float
  INTEGER, PARAMETER :: dp=c_double
  INTEGER(c_int), SAVE :: nvx,nvz,nnx,nnz,fom,gdx,gdz
  INTEGER(c_int), SAVE :: vnl,vnr,vnt,vnb,nrnx,nrnz,sgdl,rbint
  INTEGER(c_int), SAVE :: nnxr,nnzr,asgr
  INTEGER(c_int), DIMENSION (:,:), ALLOCATABLE :: nsts,nstsr,srs
  REAL(KIND=sp), SAVE :: gox,goz,dnx,dnz,dvx,dvz,snb,earth
  REAL(KIND=sp), SAVE :: goxd,gozd,dvxd,dvzd,dnxd,dnzd
  REAL(KIND=sp), SAVE :: drnx,drnz,gorx,gorz
  REAL(KIND=sp), SAVE :: dnxr,dnzr,goxr,gozr
  REAL(KIND=sp), DIMENSION (:,:), ALLOCATABLE, SAVE :: velv,veln,velnb
  REAL(KIND=sp), DIMENSION (:,:), ALLOCATABLE, SAVE :: ttn,ttnr
  !REAL(KIND=sp), DIMENSION (:), ALLOCATABLE, SAVE :: rcx,rcz
  REAL(KIND=sp), PARAMETER :: pi=3.1415926535898 
  !!!--------------------------------------------------------------
  !!	modified by Hongjian Fang @ USTC
  !	real,dimension(:),allocatable,save::rw
  !	integer,dimension(:),allocatable,save::iw,col
  !	real,dimension(:,:,:),allocatable::vpf,vsf
  !	real,dimension(:),allocatable,save::obst,cbst,wt,dtres
  !!	integer,dimension(:),allocatable,save::cbst_stat
  !	real,dimension(:,:,:),allocatable,save::sen_vs,sen_vp,sen_rho
  !!!	real,dimension(:,:,:),allocatable,save::sen_vsRc,sen_vpRc,sen_rhoRc
  !!!	real,dimension(:,:,:),allocatable,save::sen_vsRg,sen_vpRg,sen_rhoRg
  !!!	real,dimension(:,:,:),allocatable,save::sen_vsLc,sen_vpLc,sen_rhoLc
  !!!	real,dimension(:,:,:),allocatable,save::sen_vsLg,sen_vpLg,sen_rhoLg
  !!!	integer,save:: count1,count2
  !	integer*8,save:: nar
  !	integer,save:: iter,maxiter
  !!!--------------------------------------------------------------
  !
  ! nvx,nvz = B-spline vertex values
  ! dvx,dvz = B-spline vertex separation
  ! velv(i,j) = velocity values at control points
  ! nnx,nnz = Number of nodes of grid in x and z
  ! nnxr,nnzr = Number of nodes of refined grid in x and z
  ! gox,goz = Origin of grid (theta,phi)
  ! goxr, gozr = Origin of refined grid (theta,phi)
  ! dnx,dnz = Node separation of grid in  x and z
  ! dnxr,dnzr = Node separation of refined grid in x and z
  ! veln(i,j) = velocity values on a refined grid of nodes
  ! velnb(i,j) = Backup of veln required for source grid refinement
  ! ttn(i,j) = traveltime field on the refined grid of nodes
  ! ttnr(i,j) = ttn for refined grid
  ! nsts(i,j) = node status (-1=far,0=alive,>0=close)
  ! nstsr(i,j) = nsts for refined grid
  ! checkstat = check status of memory allocation
  ! fom = use first-order(0) or mixed-order(1) scheme
  ! snb = Maximum size of narrow band as fraction of nnx*nnz
  ! nrc = number of receivers
  ! rcx(i),rcz(i) = (x,z) coordinates of receivers
  ! earth = radius of Earth (in km)
  ! goxd,gozd = gox,goz in degrees
  ! dvxd,dvzd = dvx,dvz in degrees
  ! dnzd,dnzd = dnx,dnz in degrees
  ! gdx,gdz = grid dicing in x and z
  ! vnl,vnr,vnb,vnt = Bounds of refined grid
  ! nrnx,nrnz = Number of nodes in x and z for refined grid
  ! gorx,gorz = Grid origin of refined grid
  ! sgdl = Source grid dicing level
  ! rbint = Ray-boundary intersection (0=no, 1=yes).
  ! asgr = Apply source grid refinement (0=no,1=yes)
  ! srs = Source-receiver status (0=no path, 1=path exists)
  !

  ! modified by Nanqiao Du for parallelization
  !$omp threadprivate(nvx,nvz,nnx,nnz,fom,gdx,gdz)
  !$omp threadprivate(vnl,vnr,vnt,vnb,nrnx,nrnz,sgdl,rbint)
  !$omp threadprivate(nnxr,nnzr,asgr)
  !$omp threadprivate(nsts,nstsr,srs)
  !$omp threadprivate(gox,goz,dnx,dnz,dvx,dvz,snb,earth)
  !$omp threadprivate(goxd,gozd,dvxd,dvzd,dnxd,dnzd)
  !$omp threadprivate(drnx,drnz,gorx,gorz)
  !$omp threadprivate(dnxr,dnzr,goxr,gozr)
  !$omp threadprivate(velv,veln,velnb)
  !$omp threadprivate(ttn,ttnr)
  
END MODULE globalp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module contains all the subroutines used to calculate
! the first-arrival traveltime field through the grid.
! Subroutines are:
! (1) travel
! (2) fouds1
! (3) fouds2
! (4) addtree
! (5) downtree
! (6) updtree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE traveltime
  USE globalp
  use omp_lib
  IMPLICIT NONE
  INTEGER ntr
  TYPE backpointer
    !INTEGER(KIND=2) :: px,pz
  INTEGER(c_int) :: px,pz
  END TYPE backpointer
  TYPE(backpointer), DIMENSION (:), ALLOCATABLE :: btg
  ! modified by Nanqiao Du
  !$omp threadprivate(ntr,btg)

  !
  ! btg = backpointer to relate grid nodes to binary tree entries
  ! px = grid-point in x
  ! pz = grid-point in z
  ! ntr = number of entries in binary tree
  !
  public::travel
  CONTAINS
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TYPE: SUBROUTINE
    ! CODE: FORTRAN 90
    ! This subroutine is passed the location of a source, and from
    ! this point the first-arrival traveltime field through the
    ! velocity grid is determined.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE travel(scx,scz,urg)
      IMPLICIT NONE
      INTEGER :: isx,isz,sw,i,j,ix,iz,urg,swrg
      REAL(KIND=sp) :: scx,scz,vsrc,dsx,dsz,ds
      REAL(KIND=sp), DIMENSION (2,2) :: vss
      ! isx,isz = grid cell indices (i,j,k) which contains source
      ! scx,scz = (r,x,y) location of source
      ! sw = a switch (0=off,1=on)
      ! ix,iz = j,k position of "close" point with minimum traveltime
      ! maxbt = maximum size of narrow band binary tree
      ! rd2,rd3 = substitution variables
      ! vsrc = velocity at source
      ! vss = velocity at nodes surrounding source
      ! dsx, dsz = distance from source to cell boundary in x and z
      ! ds = distance from source to nearby node
      ! urg = use refined grid (0=no,1=yes,2=previously used)
      ! swrg = switch to end refined source grid computation
      !
      ! The first step is to find out where the source resides
      ! in the grid of nodes. The cell in which it resides is
      ! identified by the "north-west" node of the cell. If the
      ! source lies on the edge or corner (a node) of the cell, then
      ! this scheme still applies.
      !
      isx=INT((scx-gox)/dnx)+1
      isz=INT((scz-goz)/dnz)+1
      sw=0
      IF(isx.lt.1.or.isx.gt.nnx)sw=1
      IF(isz.lt.1.or.isz.gt.nnz)sw=1
      IF(sw.eq.1)then
        scx=90.0-scx*180.0/pi
        scz=scz*180.0/pi
        WRITE(6,*)"Source lies outside bounds of model (lat,long)= ",scx,scz
        WRITE(6,*)"TERMINATING PROGRAM!!!"
        STOP
      ENDIF
      IF(isx.eq.nnx)isx=isx-1
      IF(isz.eq.nnz)isz=isz-1
      !
      ! Set all values of nsts to -1 if beginning from a source
      ! point.
      !
      IF(urg.NE.2)nsts=-1
      !
      ! set initial size of binary tree to zero
      !
      ntr=0
      IF(urg.EQ.2)THEN
        !
        !  In this case, source grid refinement has been applied, so
        !  the initial narrow band will come from resampling the
        !  refined grid.
        !
        DO i=1,nnx
          DO j=1,nnz
            IF(nsts(j,i).GT.0)THEN
              CALL addtree(j,i)
            ENDIF
          ENDDO
        ENDDO
      ELSE 
        !
        !  In general, the source point need not lie on a grid point.
        !  Bi-linear interpolation is used to find velocity at the
        !  source point.
        !
        nsts=-1
        DO i=1,2
          DO j=1,2
            vss(i,j)=veln(isz-1+j,isx-1+i)
          ENDDO
        ENDDO
        dsx=(scx-gox)-(isx-1)*dnx
        dsz=(scz-goz)-(isz-1)*dnz
        CALL bilinear(vss,dsx,dsz,vsrc)
        !
        !  Now find the traveltime at the four surrounding grid points. This
        !  is calculated approximately by assuming the traveltime from the
        !  source point to each node is equal to the the distance between
        !  the two points divided by the average velocity of the points
        !
        DO i=1,2
          DO j=1,2
            ds=SQRT((dsx-(i-1)*dnx)**2+(dsz-(j-1)*dnz)**2)
            ttn(isz-1+j,isx-1+i)=2.0*ds/(vss(i,j)+vsrc)
            CALL addtree(isz-1+j,isx-1+i)
          ENDDO
        ENDDO
      ENDIF
      !
      ! Now calculate the first-arrival traveltimes at the
      ! remaining grid points. This is done via a loop which
      ! repeats the procedure of finding the first-arrival
      ! of all "close" points, adding it to the set of "alive"
      ! points and updating the points surrounding the new "alive"
      ! point. The process ceases when the binary tree is empty,
      ! in which case all grid points are "alive".
      !
      DO WHILE(ntr.gt.0)
        !
        ! First, check whether source grid refinement is
        ! being applied; if so, then there is a special
        ! exit condition.
        !
        IF(urg.EQ.1)THEN
          ix=btg(1)%px
          iz=btg(1)%pz
          swrg=0
          IF(ix.EQ.1)THEN
            IF(vnl.NE.1)swrg=1
          ENDIF
          IF(ix.EQ.nnx)THEN
            IF(vnr.NE.nnx)swrg=1
          ENDIF
          IF(iz.EQ.1)THEN
            IF(vnt.NE.1)swrg=1
          ENDIF
          IF(iz.EQ.nnz)THEN
            IF(vnb.NE.nnz)swrg=1
          ENDIF
          IF(swrg.EQ.1)THEN
            nsts(iz,ix)=0
            EXIT
          ENDIF
        ENDIF
        !
        ! Set the "close" point with minimum traveltime
        ! to "alive"
        !
        ix=btg(1)%px
        iz=btg(1)%pz
        nsts(iz,ix)=0
        !
        ! Update the binary tree by removing the root and
        ! sweeping down the tree.
        !
        CALL downtree
        !
        ! Now update or find values of up to four grid points
        ! that surround the new "alive" point.
        !
        ! Test points that vary in x
        !
        DO i=ix-1,ix+1,2
          IF(i.ge.1.and.i.le.nnx)THEN
            IF(nsts(iz,i).eq.-1)THEN
              !
              ! This option occurs when a far point is added to the list
              ! of "close" points
              !
              IF(fom.eq.0)THEN
                CALL fouds1(iz,i)
              ELSE
                CALL fouds2(iz,i)
              ENDIF
              CALL addtree(iz,i)
            ELSE IF(nsts(iz,i).gt.0)THEN
              !
              ! This happens when a "close" point is updated
              !
              IF(fom.eq.0)THEN
                CALL fouds1(iz,i)
              ELSE
                CALL fouds2(iz,i)
              ENDIF
              CALL updtree(iz,i)
            ENDIF
          ENDIF
        ENDDO
        !
        ! Test points that vary in z
        !
        DO i=iz-1,iz+1,2
          IF(i.ge.1.and.i.le.nnz)THEN
            IF(nsts(i,ix).eq.-1)THEN
              !
              ! This option occurs when a far point is added to the list
              ! of "close" points
              !
              IF(fom.eq.0)THEN
                CALL fouds1(i,ix)
              ELSE
                CALL fouds2(i,ix)
              ENDIF
              CALL addtree(i,ix)
            ELSE IF(nsts(i,ix).gt.0)THEN
              !
              ! This happens when a "close" point is updated
              !
              IF(fom.eq.0)THEN
                CALL fouds1(i,ix)
              ELSE
                CALL fouds2(i,ix)
              ENDIF
              CALL updtree(i,ix)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    END SUBROUTINE travel
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TYPE: SUBROUTINE
    ! CODE: FORTRAN 90
    ! This subroutine calculates a trial first-arrival traveltime
    ! at a given node from surrounding nodes using the
    ! First-Order Upwind Difference Scheme (FOUDS) of
    ! Sethian and Popovici (1999).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE fouds1(iz,ix)
      IMPLICIT NONE
      INTEGER :: j,k,ix,iz,tsw1,swsol
      REAL(KIND=sp) :: trav,travm,slown,tdsh,tref
      REAL(KIND=sp) :: a,b,c,u,v,em,ri,risti
      REAL(KIND=sp) :: rd1
      !
      ! ix = NS position of node coordinate for determination
      ! iz = EW vertical position of node coordinate for determination
      ! trav = traveltime calculated for trial node
      ! travm = minimum traveltime calculated for trial node
      ! slown = slowness at (iz,ix)
      ! tsw1 = traveltime switch (0=first time,1=previously)
      ! a,b,c,u,v,em = Convenience variables for solving quadratic
      ! tdsh = local traveltime from neighbouring node
      ! tref = reference traveltime at neighbouring node
      ! ri = Radial distance
      ! risti = ri*sin(theta) at point (iz,ix)
      ! rd1 = dummy variable
      ! swsol = switch for solution (0=no solution, 1=solution)
      !
      ! Inspect each of the four quadrants for the minimum time
      ! solution.
      !
      tsw1=0
      slown=1.0/veln(iz,ix)
      ri=earth
      risti=ri*sin(gox+(ix-1)*dnx)
      DO j=ix-1,ix+1,2
        DO k=iz-1,iz+1,2 
          IF(j.GE.1.AND.j.LE.nnx)THEN
            IF(k.GE.1.AND.k.LE.nnz)THEN
              !
              !           There are seven solution options in
              !           each quadrant.
              !
              swsol=0
              IF(nsts(iz,j).EQ.0)THEN
                swsol=1
                IF(nsts(k,ix).EQ.0)THEN
                  u=ri*dnx
                  v=risti*dnz
                  em=ttn(k,ix)-ttn(iz,j)
                  a=u**2+v**2
                  b=-2.0*u**2*em
                  c=u**2*(em**2-v**2*slown**2)
                  tref=ttn(iz,j)
                ELSE
                  a=1.0
                  b=0.0
                  c=-slown**2*ri**2*dnx**2
                  tref=ttn(iz,j)
                ENDIF
              ELSE IF(nsts(k,ix).EQ.0)THEN
                swsol=1
                a=1.0
                b=0.0
                c=-(slown*risti*dnz)**2
                tref=ttn(k,ix)
              ENDIF
              !
              !           Now find the solution of the quadratic equation
              !
              IF(swsol.EQ.1)THEN
                rd1=b**2-4.0*a*c
                IF(rd1.LT.0.0)rd1=0.0
                tdsh=(-b+sqrt(rd1))/(2.0*a)
                trav=tref+tdsh
                IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
                ELSE
                  travm=trav
                  tsw1=1
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      ttn(iz,ix)=travm
    END SUBROUTINE fouds1
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TYPE: SUBROUTINE
    ! CODE: FORTRAN 90
    ! This subroutine calculates a trial first-arrival traveltime
    ! at a given node from surrounding nodes using the
    ! Mixed-Order (2nd) Upwind Difference Scheme (FOUDS) of
    ! Popovici and Sethian (2002).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE fouds2(iz,ix)
      IMPLICIT NONE
      INTEGER :: j,k,j2,k2,ix,iz,tsw1
      INTEGER :: swj,swk,swsol
      REAL(KIND=sp) :: trav,travm,slown,tdsh,tref,tdiv
      REAL(KIND=sp) :: a,b,c,u,v,em,ri,risti,rd1
      !
      ! ix = NS position of node coordinate for determination
      ! iz = EW vertical position of node coordinate for determination
      ! trav = traveltime calculated for trial node
      ! travm = minimum traveltime calculated for trial node
      ! slown = slowness at (iz,ix)
      ! tsw1 = traveltime switch (0=first time,1=previously)
      ! a,b,c,u,v,em = Convenience variables for solving quadratic
      ! tdsh = local traveltime from neighbouring node
      ! tref = reference traveltime at neighbouring node
      ! ri = Radial distance
      ! risti = ri*sin(theta) at point (iz,ix)
      ! swj,swk = switches for second order operators
      ! tdiv = term to divide tref by depending on operator order
      ! swsol = switch for solution (0=no solution, 1=solution)
      !
      ! Inspect each of the four quadrants for the minimum time
      ! solution.
      !
      tsw1=0
      slown=1.0/veln(iz,ix)
      ri=earth
      risti=ri*sin(gox+(ix-1)*dnx)
      DO j=ix-1,ix+1,2
        IF(j.GE.1.AND.j.LE.nnx)THEN
          swj=-1
          IF(j.eq.ix-1)THEN
            j2=j-1
            IF(j2.GE.1)THEN
              IF(nsts(iz,j2).EQ.0)swj=0
            ENDIF
          ELSE
            j2=j+1
            IF(j2.LE.nnx)THEN
              IF(nsts(iz,j2).EQ.0)swj=0
            ENDIF
          ENDIF
          IF(nsts(iz,j).EQ.0.AND.swj.EQ.0)THEN
            swj=-1
            IF(ttn(iz,j).GT.ttn(iz,j2))THEN
              swj=0
            ENDIF
          ELSE
            swj=-1
          ENDIF
          DO k=iz-1,iz+1,2
            IF(k.GE.1.AND.k.LE.nnz)THEN
              swk=-1
              IF(k.eq.iz-1)THEN
                k2=k-1
                IF(k2.GE.1)THEN
                  IF(nsts(k2,ix).EQ.0)swk=0
                ENDIF
              ELSE
                k2=k+1
                IF(k2.LE.nnz)THEN
                  IF(nsts(k2,ix).EQ.0)swk=0
                ENDIF
              ENDIF
              IF(nsts(k,ix).EQ.0.AND.swk.EQ.0)THEN
                swk=-1
                IF(ttn(k,ix).GT.ttn(k2,ix))THEN
                  swk=0
                ENDIF
              ELSE
                swk=-1
              ENDIF
              !
              !           There are 8 solution options in
              !           each quadrant.
              !
              swsol=0
              IF(swj.EQ.0)THEN
                swsol=1
                IF(swk.EQ.0)THEN
                  u=2.0*ri*dnx
                  v=2.0*risti*dnz
                  em=4.0*ttn(iz,j)-ttn(iz,j2)-4.0*ttn(k,ix)
                  em=em+ttn(k2,ix)
                  a=v**2+u**2
                  b=2.0*em*u**2
                  c=u**2*(em**2-slown**2*v**2)
                  tref=4.0*ttn(iz,j)-ttn(iz,j2)
                  tdiv=3.0
                ELSE IF(nsts(k,ix).EQ.0)THEN
                  u=risti*dnz
                  v=2.0*ri*dnx
                  em=3.0*ttn(k,ix)-4.0*ttn(iz,j)+ttn(iz,j2)
                  a=v**2+9.0*u**2
                  b=6.0*em*u**2
                  c=u**2*(em**2-slown**2*v**2)
                  tref=ttn(k,ix)
                  tdiv=1.0
                ELSE
                  u=2.0*ri*dnx
                  a=1.0
                  b=0.0
                  c=-u**2*slown**2
                  tref=4.0*ttn(iz,j)-ttn(iz,j2)
                  tdiv=3.0
                ENDIF
              ELSE IF(nsts(iz,j).EQ.0)THEN
                swsol=1
                IF(swk.EQ.0)THEN
                  u=ri*dnx
                  v=2.0*risti*dnz
                  em=3.0*ttn(iz,j)-4.0*ttn(k,ix)+ttn(k2,ix)
                  a=v**2+9.0*u**2
                  b=6.0*em*u**2
                  c=u**2*(em**2-v**2*slown**2)
                  tref=ttn(iz,j)
                  tdiv=1.0
                ELSE IF(nsts(k,ix).EQ.0)THEN
                  u=ri*dnx
                  v=risti*dnz
                  em=ttn(k,ix)-ttn(iz,j)
                  a=u**2+v**2
                  b=-2.0*u**2*em
                  c=u**2*(em**2-v**2*slown**2)
                  tref=ttn(iz,j)
                  tdiv=1.0
                ELSE
                  a=1.0
                  b=0.0
                  c=-slown**2*ri**2*dnx**2
                  tref=ttn(iz,j)
                  tdiv=1.0
                ENDIF
              ELSE
                IF(swk.EQ.0)THEN
                  swsol=1
                  u=2.0*risti*dnz
                  a=1.0
                  b=0.0
                  c=-u**2*slown**2
                  tref=4.0*ttn(k,ix)-ttn(k2,ix)
                  tdiv=3.0
                ELSE IF(nsts(k,ix).EQ.0)THEN
                  swsol=1
                  a=1.0
                  b=0.0
                  c=-slown**2*risti**2*dnz**2
                  tref=ttn(k,ix)
                  tdiv=1.0
                ENDIF
              ENDIF
              !
              !           Now find the solution of the quadratic equation
              !
              IF(swsol.EQ.1)THEN
                rd1=b**2-4.0*a*c
                IF(rd1.LT.0.0)rd1=0.0
                tdsh=(-b+sqrt(rd1))/(2.0*a)
                trav=(tref+tdsh)/tdiv
                IF(tsw1.EQ.1)THEN
                  travm=MIN(trav,travm)
                ELSE
                  travm=trav
                  tsw1=1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      ttn(iz,ix)=travm
    END SUBROUTINE fouds2
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TYPE: SUBROUTINE
    ! CODE: FORTRAN 90
    ! This subroutine adds a value to the binary tree by
    ! placing a value at the bottom and pushing it up
    ! to its correct position.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE addtree(iz,ix)
      IMPLICIT NONE
      INTEGER :: ix,iz,tpp,tpc
      TYPE(backpointer) :: exch
      !
      ! ix,iz = grid position of new addition to tree
      ! tpp = tree position of parent
      ! tpc = tree position of child
      ! exch = dummy to exchange btg values
      !
      ! First, increase the size of the tree by one.
      !
      ntr=ntr+1
      !
      ! Put new value at base of tree
      !
      nsts(iz,ix)=ntr
      btg(ntr)%px=ix
      btg(ntr)%pz=iz
      !
      ! Now filter the new value up to its correct position
      !
      tpc=ntr
      tpp=tpc/2
      DO WHILE(tpp.gt.0)
        IF(ttn(iz,ix).lt.ttn(btg(tpp)%pz,btg(tpp)%px))THEN
          nsts(iz,ix)=tpp
          nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
          exch=btg(tpc)
          btg(tpc)=btg(tpp)
          btg(tpp)=exch
          tpc=tpp
          tpp=tpc/2
        ELSE
          tpp=0
        ENDIF
      ENDDO
    END SUBROUTINE addtree
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TYPE: SUBROUTINE
    ! CODE: FORTRAN 90
    ! This subroutine updates the binary tree after the root
    ! value has been used. The root is replaced by the value
    ! at the bottom of the tree, which is then filtered down
    ! to its correct position. This ensures that the tree remains
    ! balanced.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE downtree
      IMPLICIT NONE
      INTEGER :: tpp,tpc
      REAL(KIND=sp) :: rd1,rd2
      TYPE(backpointer) :: exch
      !
      ! tpp = tree position of parent
      ! tpc = tree position of child
      ! exch = dummy to exchange btg values
      ! rd1,rd2 = substitution variables
      !
      ! Replace root of tree with its last value
      !
      IF(ntr.EQ.1)THEN
        ntr=ntr-1
        RETURN
      ENDIF
      nsts(btg(ntr)%pz,btg(ntr)%px)=1
      btg(1)=btg(ntr)
      !
      ! Reduce size of tree by one
      !
      ntr=ntr-1
      !
      ! Now filter new root down to its correct position
      !
      tpp=1
      tpc=2*tpp
      DO WHILE(tpc.lt.ntr)
        !
        ! Check which of the two children is smallest - use the smallest
        !
        rd1=ttn(btg(tpc)%pz,btg(tpc)%px)
        rd2=ttn(btg(tpc+1)%pz,btg(tpc+1)%px)
        IF(rd1.gt.rd2)THEN
          tpc=tpc+1
        ENDIF
        !
        !  Check whether the child is smaller than the parent; if so, then swap,
        !  if not, then we are done
        !
        rd1=ttn(btg(tpc)%pz,btg(tpc)%px)
        rd2=ttn(btg(tpp)%pz,btg(tpp)%px)
        IF(rd1.lt.rd2)THEN
          nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
          nsts(btg(tpc)%pz,btg(tpc)%px)=tpp
          exch=btg(tpc)
          btg(tpc)=btg(tpp)
          btg(tpp)=exch
          tpp=tpc
          tpc=2*tpp
        ELSE
          tpc=ntr+1
        ENDIF
      ENDDO
      !
      ! If ntr is an even number, then we still have one more test to do
      !
      IF(tpc.eq.ntr)THEN
        rd1=ttn(btg(tpc)%pz,btg(tpc)%px)
        rd2=ttn(btg(tpp)%pz,btg(tpp)%px)
        IF(rd1.lt.rd2)THEN
          nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
          nsts(btg(tpc)%pz,btg(tpc)%px)=tpp
          exch=btg(tpc)
          btg(tpc)=btg(tpp)
          btg(tpp)=exch
        ENDIF
      ENDIF
    END SUBROUTINE downtree
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TYPE: SUBROUTINE
    ! CODE: FORTRAN 90
    ! This subroutine updates a value on the binary tree. The FMM
    ! should only produce updated values that are less than their
    ! prior values.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updtree(iz,ix)
      IMPLICIT NONE
      INTEGER :: ix,iz,tpp,tpc
      TYPE(backpointer) :: exch
      !
      ! ix,iz = grid position of new addition to tree
      ! tpp = tree position of parent
      ! tpc = tree position of child
      ! exch = dummy to exchange btg values
      !
      ! Filter the updated value to its correct position
      !
      tpc=nsts(iz,ix)
      tpp=tpc/2
      DO WHILE(tpp.gt.0)
        IF(ttn(iz,ix).lt.ttn(btg(tpp)%pz,btg(tpp)%px))THEN
          nsts(iz,ix)=tpp
          nsts(btg(tpp)%pz,btg(tpp)%px)=tpc
          exch=btg(tpc)
          btg(tpc)=btg(tpp)
          btg(tpp)=exch
          tpc=tpp
          tpp=tpc/2
        ELSE
          tpp=0
        ENDIF
      ENDDO
    END SUBROUTINE updtree
  
  END MODULE traveltime

  SUBROUTINE gridder(pv)
    !subroutine gridder(pv)
    !subroutine gridder()
    USE globalp
    IMPLICIT NONE
    INTEGER :: i,j,l,m,i1,j1,conx,conz,stx,stz,checkstat
    REAL(KIND=sp) :: u,sumi,sumj
    REAL(KIND=sp), DIMENSION(:,:), ALLOCATABLE :: ui,vi
    !real(kind=sp):: ui(gdx+1,4),vi(gdz+1,4)
    !CHARACTER (LEN=30) :: grid
    !
    ! u = independent parameter for b-spline
    ! ui,vi = bspline basis functions
    ! conx,conz = variables for edge of B-spline grid
    ! stx,stz = counters for veln grid points
    ! sumi,sumj = summation variables for computing b-spline
    !
    !C---------------------------------------------------------------
    real(dp),INTENT(IN) ::  pv(*)
    !integer count1
    !C---------------------------------------------------------------
    ! Open the grid file and read in the velocity grid.
    !
    !OPEN(UNIT=10,FILE=grid,STATUS='old')
    !READ(10,*)nvx,nvz
    !READ(10,*)goxd,gozd
    !READ(10,*)dvxd,dvzd
    !count1=0
    DO i=0,nvz+1
      DO j=0,nvx+1
        !	count1=count1+1
        !      READ(10,*)velv(i,j)
        !	velv(i,j)=real(pv(count1))
        velv(i,j)=real(pv(i*(nvx+2)+j+1))
      ENDDO
    ENDDO
    !CLOSE(10)
    !
    ! Convert from degrees to radians
    !
    !
    ! Now dice up the grid
    !
    ALLOCATE(ui(gdx+1,4), STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: Subroutine gridder: REAL ui'
    ENDIF
    DO i=1,gdx+1
      u=gdx
      u=(i-1)/u
      ui(i,1)=(1.0-u)**3/6.0
      ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
      ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
      ui(i,4)=u**3/6.0
    ENDDO
    ALLOCATE(vi(gdz+1,4), STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: Subroutine gridder: REAL vi'
    ENDIF
    DO i=1,gdz+1
      u=gdz
      u=(i-1)/u
      vi(i,1)=(1.0-u)**3/6.0
      vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
      vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
      vi(i,4)=u**3/6.0
    ENDDO
    DO i=1,nvz-1
      conz=gdz
      IF(i==nvz-1)conz=gdz+1
      DO j=1,nvx-1
        conx=gdx
        IF(j==nvx-1)conx=gdx+1
        DO l=1,conz
          stz=gdz*(i-1)+l
          DO m=1,conx
            stx=gdx*(j-1)+m
            sumi=0.0
            DO i1=1,4
              sumj=0.0
              DO j1=1,4
                sumj=sumj+ui(m,j1)*velv(i-2+i1,j-2+j1)
              ENDDO
              sumi=sumi+vi(l,i1)*sumj
            ENDDO
            veln(stz,stx)=sumi
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(ui,vi, STAT=checkstat)
    IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE gridder: REAL ui,vi'
    ENDIF
  END SUBROUTINE gridder
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TYPE: SUBROUTINE
  ! CODE: FORTRAN 90
  ! This subroutine is similar to bsplreg except that it has been
  ! modified to deal with source grid refinement
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE bsplrefine
    USE globalp
    INTEGER :: i,j,k,l,i1,j1,st1,st2,nrzr,nrxr
    INTEGER :: origx,origz,conx,conz,idm1,idm2
    REAL(KIND=sp) :: u,v
    REAL(KIND=sp), DIMENSION (4) :: sum
    REAL(KIND=sp), DIMENSION(gdx*sgdl+1,gdz*sgdl+1,4) :: ui,vi
    !
    ! nrxr,nrzr = grid refinement level for source grid in x,z
    ! origx,origz = local origin of refined source grid
    !
    ! Begin by calculating the values of the basis functions
    !
    nrxr=gdx*sgdl
    nrzr=gdz*sgdl
    DO i=1,nrzr+1
      v=nrzr
      v=(i-1)/v
      DO j=1,nrxr+1
        u=nrxr
        u=(j-1)/u
        ui(j,i,1)=(1.0-u)**3/6.0
        ui(j,i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
        ui(j,i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
        ui(j,i,4)=u**3/6.0
        vi(j,i,1)=(1.0-v)**3/6.0
        vi(j,i,2)=(4.0-6.0*v**2+3.0*v**3)/6.0
        vi(j,i,3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
        vi(j,i,4)=v**3/6.0
      ENDDO
    ENDDO
    !
    ! Calculate the velocity values.
    !
    origx=(vnl-1)*sgdl+1
    origz=(vnt-1)*sgdl+1
    DO i=1,nvz-1
      conz=nrzr
      IF(i==nvz-1)conz=nrzr+1
      DO j=1,nvx-1
        conx=nrxr
        IF(j==nvx-1)conx=nrxr+1
        DO k=1,conz
          st1=gdz*(i-1)+(k-1)/sgdl+1
          IF(st1.LT.vnt.OR.st1.GT.vnb)CYCLE
          st1=nrzr*(i-1)+k
          DO l=1,conx
            st2=gdx*(j-1)+(l-1)/sgdl+1
            IF(st2.LT.vnl.OR.st2.GT.vnr)CYCLE
            st2=nrxr*(j-1)+l
            DO i1=1,4
              sum(i1)=0.0
              DO j1=1,4
                sum(i1)=sum(i1)+ui(l,k,j1)*velv(i-2+i1,j-2+j1)
              ENDDO
              sum(i1)=vi(l,k,i1)*sum(i1)
            ENDDO
            idm1=st1-origz+1
            idm2=st2-origx+1
            IF(idm1.LT.1.OR.idm1.GT.nnz)CYCLE
            IF(idm2.LT.1.OR.idm2.GT.nnx)CYCLE
            veln(idm1,idm2)=sum(1)+sum(2)+sum(3)+sum(4)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE bsplrefine 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TYPE: SUBROUTINE
  ! CODE: FORTRAN 90
  ! This subroutine calculates all receiver traveltimes for
  ! a given source and writes the results to file.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !SUBROUTINE srtimes(scx,scz,rcx1,rcz1,cbst1)
  SUBROUTINE srtimes(scx,scz,rcx1,rcz1,cbst1)
    USE globalp
    IMPLICIT NONE
    real(sp)  :: scx,scz,rcx1,rcz1
    real(sp),intent(inout) :: cbst1
    INTEGER :: k,l,irx,irz,sw,isx,isz
    INTEGER, PARAMETER :: noray=0,yesray=1
    !INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(6)
    REAL(sp) :: trr
    REAL(KIND=sp), PARAMETER :: norayt=0.0
    REAL(KIND=sp) :: drx,drz,produ
    REAL(KIND=sp) :: sred,dpl,rd1,vels,velr
    REAL(KIND=sp), DIMENSION (2,2) :: vss
    !!------------------------------------------------------
    !	modified by Hongjian Fang @ USTC
    !integer no_p,nsrc
    !real dist
    !	real cbst(*) !note that the type difference(kind=i5 vs real)
    !	integer cbst_stat(*)
    !!------------------------------------------------------
    !
    ! irx,irz = Coordinates of cell containing receiver
    ! trr = traveltime value at receiver
    ! produ = dummy multiplier
    ! drx,drz = receiver distance from (i,j,k) grid node
    ! scx,scz = source coordinates
    ! isx,isz = source cell location
    ! sred = Distance from source to receiver
    ! dpl = Minimum path length in source neighbourhood.
    ! vels,velr = velocity at source and receiver
    ! vss = velocity at four grid points about source or receiver.
    ! csid = current source ID
    ! noray = switch to indicate no ray present
    ! norayt = default value given to null ray
    ! yesray = switch to indicate that ray is present
    !
    ! Determine source-receiver traveltimes one at a time.
    !
    !0605DO i=1,nrc
    !0605   IF(srs(i,csid).EQ.0)THEN
    !0605!      WRITE(10,*)noray,norayt
    !0605      CYCLE
    !0605   ENDIF
    ! 
    !  The first step is to locate the receiver in the grid.
    !
    irx=INT((rcx1-gox)/dnx)+1
    irz=INT((rcz1-goz)/dnz)+1
    sw=0
    IF(irx.lt.1.or.irx.gt.nnx)sw=1
    IF(irz.lt.1.or.irz.gt.nnz)sw=1
    IF(sw.eq.1)then
      rcx1=90.0-rcx1*180.0/pi
      rcz1=rcz1*180.0/pi
      WRITE(6,*)"Receiver lies outside model (lat,long)= ",rcx1,rcz1
      WRITE(6,*)"TERMINATING PROGRAM!!!!"
      STOP
    ENDIF
    IF(irx.eq.nnx)irx=irx-1
    IF(irz.eq.nnz)irz=irz-1
    !
    !  Location of receiver successfully found within the grid. Now approximate
    !  traveltime at receiver using bilinear interpolation from four
    !  surrounding grid points. Note that bilinear interpolation is a poor
    !  approximation when traveltime gradient varies significantly across a cell,
    !  particularly near the source. Thus, we use an improved approximation in this
    !  case. First, locate current source cell.
    !
    isx=INT((scx-gox)/dnx)+1
    isz=INT((scz-goz)/dnz)+1
    dpl=dnx*earth
    rd1=dnz*earth*SIN(gox)
    IF(rd1.LT.dpl)dpl=rd1
    rd1=dnz*earth*SIN(gox+(nnx-1)*dnx)
    IF(rd1.LT.dpl)dpl=rd1
    sred=((scx-rcx1)*earth)**2
    sred=sred+((scz-rcz1)*earth*SIN(rcx1))**2
    sred=SQRT(sred)
    IF(sred.LT.dpl)sw=1
    IF(isx.EQ.irx)THEN
      IF(isz.EQ.irz)sw=1
    ENDIF
    IF(sw.EQ.1)THEN
      !
      !     Compute velocity at source and receiver
      !
      DO k=1,2
        DO l=1,2
          vss(k,l)=veln(isz-1+l,isx-1+k)
        ENDDO
      ENDDO
      drx=(scx-gox)-(isx-1)*dnx
      drz=(scz-goz)-(isz-1)*dnz
      CALL bilinear(vss,drx,drz,vels)
      DO k=1,2
        DO l=1,2
          vss(k,l)=veln(irz-1+l,irx-1+k)
        ENDDO
      ENDDO
      drx=(rcx1-gox)-(irx-1)*dnx
      drz=(rcz1-goz)-(irz-1)*dnz
      CALL bilinear(vss,drx,drz,velr)
      trr=2.0*sred/(vels+velr)
    ELSE
      drx=(rcx1-gox)-(irx-1)*dnx
      drz=(rcz1-goz)-(irz-1)*dnz
      trr=0.0
      DO k=1,2
        DO l=1,2
          produ=(1.0-ABS(((l-1)*dnz-drz)/dnz))*(1.0-ABS(((k-1)*dnx-drx)/dnx))
          trr=trr+ttn(irz-1+l,irx-1+k)*produ
        ENDDO
      ENDDO
    ENDIF
    !   WRITE(10,*)yesray,trr
    !!-----------------------------------------------------------------
    !	modified bu Hongjian Fang @ USTC
    !	count2=count2+1
    !	cbst((no_p-1)*nsrc*nrc+(csid-1)*nrc+i)=trr
    cbst1=trr
    !	call delsph(scx,scz,rcx(i),rcz(i),dist)
    !	travel_path(count2)=dist
    !cbst_stat((no_p-1)*nsrc*nrc+(csid-1)*nrc+i)=yesray
    !0605ENDDO
  END SUBROUTINE srtimes
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TYPE: SUBROUTINE
  ! CODE: FORTRAN 90
  ! This subroutine calculates ray path geometries for each
  ! source-receiver combination. It will also compute
  ! Frechet derivatives using these ray paths if required.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !SUBROUTINE rpaths(wrgf,csid,cfd,scx,scz)
  !SUBROUTINE rpaths()
  SUBROUTINE rpaths(scx,scz,fdm,surfrcx,surfrcz,igr,Tid,writepath)
    USE globalp
    IMPLICIT NONE
    !INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
    INTEGER, PARAMETER :: nopath=0
    INTEGER :: j,k,l,m,n,ipx,ipz,ipxr,ipzr,nrp,sw,checkstat
    !fang!INTEGER :: wrgf,cfd,csid,ipxo,ipzo,isx,isz
    INTEGER :: ipxo,ipzo,isx,isz
    INTEGER :: ivx,ivz,ivxo,ivzo,nhp,maxrp
    INTEGER :: ivxt,ivzt,ipxt,ipzt,igref
    INTEGER, DIMENSION (4) :: chp
    REAL(KIND=sp) :: rayx,rayz
    REAL(KIND=sp) :: dpl,rd1,rd2,xi,zi,vel,velo
    REAL(KIND=sp) :: v,w,rigz,rigx,dinc,scx,scz
    REAL(KIND=sp) :: dtx,dtz,drx,drz,produ,sred
    REAL(KIND=sp), DIMENSION (:), ALLOCATABLE :: rgx,rgz
    !fang!REAL(KIND=i5), DIMENSION (:,:), ALLOCATABLE :: fdm
    REAL(KIND=sp), DIMENSION (4) :: vrat,vi,wi,vio,wio
    !fang!------------------------------------------------
    real(sp),INTENT(inout)  :: fdm(0:nvz+1,0:nvx+1)
    REAL(KIND=sp) surfrcx,surfrcz
    
    !modified by nqdu
    character*100 outfile,outpath
    integer Tid,igr,writepath
    !fang!------------------------------------------------
    !
    ! ipx,ipz = Coordinates of cell containing current point
    ! ipxr,ipzr = Same as ipx,apz except for refined grid
    ! ipxo,ipzo = Coordinates of previous point
    ! rgx,rgz = (x,z) coordinates of ray geometry
    ! ivx,ivz = Coordinates of B-spline vertex containing current point
    ! ivxo,ivzo = Coordinates of previous point
    ! maxrp = maximum number of ray points
    ! nrp = number of points to describe ray
    ! dpl = incremental path length of ray
    ! xi,zi = edge of model coordinates
    ! dtx,dtz = components of gradT
    ! wrgf = Write out raypaths? (<0=all,0=no,>0=souce id)
    ! cfd = calculate Frechet derivatives? (0=no,1=yes)
    ! csid = current source id
    ! fdm = Frechet derivative matrix
    ! nhp = Number of ray segment-B-spline cell hit points
    ! vrat = length ratio of ray sub-segment
    ! chp = pointer to incremental change in x or z cell
    ! drx,drz = distance from reference node of cell
    ! produ = variable for trilinear interpolation
    ! vel = velocity at current point
    ! velo = velocity at previous point
    ! v,w = local variables of x,z
    ! vi,wi = B-spline basis functions at current point
    ! vio,wio = vi,wi for previous point
    ! ivxt,ivzt = temporary ivr,ivx,ivz values
    ! rigx,rigz = end point of sub-segment of ray path
    ! ipxt,ipzt = temporary ipx,ipz values
    ! dinc = path length of ray sub-segment
    ! rayr,rayx,rayz = ray path coordinates in single precision
    ! isx,isz = current source cell location
    ! scx,scz = current source coordinates
    ! sred = source to ray endpoint distance
    ! igref = ray endpoint lies in refined grid? (0=no,1=yes)
    ! nopath = switch to indicate that no path is present
    !
    ! Allocate memory to arrays for storing ray path geometry
    !
    maxrp=nnx*nnz
    ALLOCATE(rgx(maxrp+1), STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgx'
    ENDIF
    ALLOCATE(rgz(maxrp+1), STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgz'
    ENDIF
    !
    ! Allocate memory to partial derivative array
    !
    !fang!IF(cfd.EQ.1)THEN
    !fang!   ALLOCATE(fdm(0:nvz+1,0:nvx+1), STAT=checkstat)
    !fang!   IF(checkstat > 0)THEN
    !fang!      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL fdm'
    !fang!   ENDIF
    !fang!ENDIF
    !
    ! Locate current source cell
    !
    IF(asgr.EQ.1)THEN
      isx=INT((scx-goxr)/dnxr)+1
      isz=INT((scz-gozr)/dnzr)+1
    ELSE
      isx=INT((scx-gox)/dnx)+1
      isz=INT((scz-goz)/dnz)+1
    ENDIF
    !
    ! Set ray incremental path length equal to half width
    ! of cell
    !
    dpl=dnx*earth
    rd1=dnz*earth*SIN(gox)
    IF(rd1.LT.dpl)dpl=rd1
    rd1=dnz*earth*SIN(gox+(nnx-1)*dnx)
    IF(rd1.LT.dpl)dpl=rd1
    dpl=0.5*dpl
    !
    ! Loop through all the receivers
    !
    !fang!DO i=1,nrc
    !
    !  If path does not exist, then cycle the loop
    !
    fdm=0.0
    !fang!   IF(cfd.EQ.1)THEN
    !fang!      fdm=0.0
    !fang!   ENDIF
    !fang!   IF(srs(i,csid).EQ.0)THEN
    !fang!      IF(wrgf.EQ.csid.OR.wrgf.LT.0)THEN
    !fang!         WRITE(40)nopath
    !fang!      ENDIF
    !fang!      IF(cfd.EQ.1)THEN
    !fang!         WRITE(50)nopath
    !fang!      ENDIF
    !fang!      CYCLE
    !fang!   ENDIF
    !
    !  The first step is to locate the receiver in the grid.
    !
    ipx=INT((surfrcx-gox)/dnx)+1
    ipz=INT((surfrcz-goz)/dnz)+1
    sw=0
    IF(ipx.lt.1.or.ipx.ge.nnx)sw=1
    IF(ipz.lt.1.or.ipz.ge.nnz)sw=1
    IF(sw.eq.1)then
      surfrcx=90.0-surfrcx*180.0/pi
      surfrcz=surfrcz*180.0/pi
      WRITE(6,*)"rpath Receiver lies outside model (lat,long)= ",surfrcx,surfrcz
      WRITE(6,*)"TERMINATING PROGRAM!!!"
      STOP
    ENDIF
    IF(ipx.eq.nnx)ipx=ipx-1
    IF(ipz.eq.nnz)ipz=ipz-1
    !
    !  First point of the ray path is the receiver
    !
    rgx(1)=surfrcx
    rgz(1)=surfrcz
    !
    !  Test to see if receiver is in source neighbourhood
    !
    sred=((scx-rgx(1))*earth)**2
    sred=sred+((scz-rgz(1))*earth*SIN(rgx(1)))**2
    sred=SQRT(sred)
    IF(sred.LT.2.0*dpl)THEN
      rgx(2)=scx
      rgz(2)=scz
      nrp=2
      sw=1
    ENDIF
    !
    !  If required, see if receiver lies within refined grid
    !
    IF(asgr.EQ.1)THEN
      ipxr=INT((surfrcx-goxr)/dnxr)+1
      ipzr=INT((surfrcz-gozr)/dnzr)+1
      igref=1
      IF(ipxr.LT.1.OR.ipxr.GE.nnxr)igref=0
      IF(ipzr.LT.1.OR.ipzr.GE.nnzr)igref=0
      IF(igref.EQ.1)THEN
        IF(nstsr(ipzr,ipxr).NE.0.OR.nstsr(ipzr+1,ipxr).NE.0)igref=0
        IF(nstsr(ipzr,ipxr+1).NE.0.OR.nstsr(ipzr+1,ipxr+1).NE.0)igref=0
      ENDIF
    ELSE
      igref=0
    ENDIF
    !
    !  Due to the method for calculating traveltime gradient, if the
    !  the ray end point lies in the source cell, then we are also done.
    !
    IF(sw.EQ.0)THEN
      IF(asgr.EQ.1)THEN
        IF(igref.EQ.1)THEN
          IF(ipxr.EQ.isx)THEN
            IF(ipzr.EQ.isz)THEN
              rgx(2)=scx
              rgz(2)=scz
              nrp=2
              sw=1
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IF(ipx.EQ.isx)THEN
          IF(ipz.EQ.isz)THEN
            rgx(2)=scx
            rgz(2)=scz
            nrp=2
            sw=1
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    !
    !  Now trace ray from receiver to "source"
    !
    DO j=1,maxrp
      IF(sw.EQ.1)EXIT
      !
      !     Calculate traveltime gradient vector for current cell using
      !     a first-order or second-order scheme.
      !
      IF(igref.EQ.1)THEN
        !
        !        In this case, we are in the refined grid.
        !
        !        First order scheme applied here.
        !
        dtx=ttnr(ipzr,ipxr+1)-ttnr(ipzr,ipxr)
        dtx=dtx+ttnr(ipzr+1,ipxr+1)-ttnr(ipzr+1,ipxr)
        dtx=dtx/(2.0*earth*dnxr)
        dtz=ttnr(ipzr+1,ipxr)-ttnr(ipzr,ipxr)
        dtz=dtz+ttnr(ipzr+1,ipxr+1)-ttnr(ipzr,ipxr+1)
        dtz=dtz/(2.0*earth*SIN(rgx(j))*dnzr)
      ELSE
        !
        !        Here, we are in the coarse grid.
        !
        !        First order scheme applied here.
        !
        dtx=ttn(ipz,ipx+1)-ttn(ipz,ipx)
        dtx=dtx+ttn(ipz+1,ipx+1)-ttn(ipz+1,ipx)
        dtx=dtx/(2.0*earth*dnx)
        dtz=ttn(ipz+1,ipx)-ttn(ipz,ipx)
        dtz=dtz+ttn(ipz+1,ipx+1)-ttn(ipz,ipx+1)
        dtz=dtz/(2.0*earth*SIN(rgx(j))*dnz)
      ENDIF
      !
      !     Calculate the next ray path point
      !
      rd1=SQRT(dtx**2+dtz**2)
      rgx(j+1)=rgx(j)-dpl*dtx/(earth*rd1)
      rgz(j+1)=rgz(j)-dpl*dtz/(earth*SIN(rgx(j))*rd1)
      !
      !     Determine which cell the new ray endpoint
      !     lies in.
      !
      ipxo=ipx
      ipzo=ipz
      IF(asgr.EQ.1)THEN
        !
        !        Here, we test to see whether the ray endpoint lies
        !        within a cell of the refined grid
        !
        ipxr=INT((rgx(j+1)-goxr)/dnxr)+1
        ipzr=INT((rgz(j+1)-gozr)/dnzr)+1
        igref=1
        IF(ipxr.LT.1.OR.ipxr.GE.nnxr)igref=0
        IF(ipzr.LT.1.OR.ipzr.GE.nnzr)igref=0
        IF(igref.EQ.1)THEN
          IF(nstsr(ipzr,ipxr).NE.0.OR.nstsr(ipzr+1,ipxr).NE.0)igref=0
          IF(nstsr(ipzr,ipxr+1).NE.0.OR.nstsr(ipzr+1,ipxr+1).NE.0)igref=0
        ENDIF
        ipx=INT((rgx(j+1)-gox)/dnx)+1
        ipz=INT((rgz(j+1)-goz)/dnz)+1
      ELSE
        ipx=INT((rgx(j+1)-gox)/dnx)+1
        ipz=INT((rgz(j+1)-goz)/dnz)+1
        igref=0
      ENDIF
      !
      !     Test the proximity of the source to the ray end point.
      !     If it is less than dpl then we are done
      !
      sred=((scx-rgx(j+1))*earth)**2
      sred=sred+((scz-rgz(j+1))*earth*SIN(rgx(j+1)))**2
      sred=SQRT(sred)
      sw=0
      IF(sred.LT.2.0*dpl)THEN
        rgx(j+2)=scx
        rgz(j+2)=scz
        nrp=j+2
        sw=1
        !fang!         IF(cfd.NE.1)EXIT
      ENDIF
      !
      !     Due to the method for calculating traveltime gradient, if the
      !     the ray end point lies in the source cell, then we are also done.
      !
      IF(sw.EQ.0)THEN
        IF(asgr.EQ.1)THEN
          IF(igref.EQ.1)THEN
            IF(ipxr.EQ.isx)THEN
              IF(ipzr.EQ.isz)THEN
                rgx(j+2)=scx
                rgz(j+2)=scz
                nrp=j+2
                sw=1
                !fang!                    IF(cfd.NE.1)EXIT
              ENDIF
            ENDIF
          ENDIF
        ELSE
          IF(ipx.EQ.isx)THEN
            IF(ipz.EQ.isz)THEN
              rgx(j+2)=scx
              rgz(j+2)=scz
              nrp=j+2
              sw=1
              !fang!                 IF(cfd.NE.1)EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      !
      !     Test whether ray path segment extends beyond
      !     box boundaries
      !
      IF(ipx.LT.1)THEN
        rgx(j+1)=gox
        ipx=1
        rbint=1
      ENDIF
      IF(ipx.GE.nnx)THEN
        rgx(j+1)=gox+(nnx-1)*dnx
        ipx=nnx-1
        rbint=1
      ENDIF
      IF(ipz.LT.1)THEN
        rgz(j+1)=goz
        ipz=1
        rbint=1
      ENDIF
      IF(ipz.GE.nnz)THEN
        rgz(j+1)=goz+(nnz-1)*dnz
        ipz=nnz-1
        rbint=1
      ENDIF
      !
      !     Calculate the Frechet derivatives if required.
      !
      !fang!     IF(cfd.EQ.1)THEN
      !
      !        First determine which B-spline cell the refined cells
      !        containing the ray path segment lies in. If they lie
      !        in more than one, then we need to divide the problem
      !        into separate parts (up to three).
      !
      ivx=INT((ipx-1)/gdx)+1
      ivz=INT((ipz-1)/gdz)+1
      ivxo=INT((ipxo-1)/gdx)+1
      ivzo=INT((ipzo-1)/gdz)+1
      !
      !        Calculate up to two hit points between straight
      !        ray segment and cell faces.
      !
      nhp=0
      IF(ivx.NE.ivxo)THEN
        nhp=nhp+1
        IF(ivx.GT.ivxo)THEN
          xi=gox+(ivx-1)*dvx
        ELSE
          xi=gox+ivx*dvx
        ENDIF
        vrat(nhp)=(xi-rgx(j))/(rgx(j+1)-rgx(j))
        chp(nhp)=1
      ENDIF
      IF(ivz.NE.ivzo)THEN
        nhp=nhp+1
        IF(ivz.GT.ivzo)THEN
          zi=goz+(ivz-1)*dvz
        ELSE
          zi=goz+ivz*dvz
        ENDIF
        rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
        IF(nhp.EQ.1)THEN
          vrat(nhp)=rd1
          chp(nhp)=2
        ELSE
          IF(rd1.GE.vrat(nhp-1))THEN
            vrat(nhp)=rd1
            chp(nhp)=2
          ELSE
            vrat(nhp)=vrat(nhp-1)
            chp(nhp)=chp(nhp-1)
            vrat(nhp-1)=rd1
            chp(nhp-1)=2
          ENDIF
        ENDIF
      ENDIF
      nhp=nhp+1
      vrat(nhp)=1.0
      chp(nhp)=0
      !
      !        Calculate the velocity, v and w values of the
      !        first point
      !
      drx=(rgx(j)-gox)-(ipxo-1)*dnx
      drz=(rgz(j)-goz)-(ipzo-1)*dnz
      vel=0.0
      DO l=1,2
        DO m=1,2
          produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
          produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
          IF(ipzo-1+m.LE.nnz.AND.ipxo-1+l.LE.nnx)THEN
            vel=vel+veln(ipzo-1+m,ipxo-1+l)*produ
          ENDIF
        ENDDO
      ENDDO
      drx=(rgx(j)-gox)-(ivxo-1)*dvx
      drz=(rgz(j)-goz)-(ivzo-1)*dvz
      v=drx/dvx
      w=drz/dvz
      !
      !        Calculate the 12 basis values at the point
      !
      vi(1)=(1.0-v)**3/6.0
      vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
      vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
      vi(4)=v**3/6.0
      wi(1)=(1.0-w)**3/6.0
      wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
      wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
      wi(4)=w**3/6.0
      ivxt=ivxo
      ivzt=ivzo
      !
      !        Now loop through the one or more sub-segments of the
      !        ray path segment and calculate partial derivatives
      !
      DO k=1,nhp
        velo=vel
        vio=vi
        wio=wi
        IF(k.GT.1)THEN
          IF(chp(k-1).EQ.1)THEN
            ivxt=ivx
          ELSE IF(chp(k-1).EQ.2)THEN
            ivzt=ivz
          ENDIF
        ENDIF
        !
        !           Calculate the velocity, v and w values of the
        !           new point
        !
        rigz=rgz(j)+vrat(k)*(rgz(j+1)-rgz(j))
        rigx=rgx(j)+vrat(k)*(rgx(j+1)-rgx(j))
        ipxt=INT((rigx-gox)/dnx)+1
        ipzt=INT((rigz-goz)/dnz)+1
        drx=(rigx-gox)-(ipxt-1)*dnx
        drz=(rigz-goz)-(ipzt-1)*dnz
        vel=0.0
        DO m=1,2
          DO n=1,2
            produ=(1.0-ABS(((n-1)*dnz-drz)/dnz))
            produ=produ*(1.0-ABS(((m-1)*dnx-drx)/dnx))
            IF(ipzt-1+n.LE.nnz.AND.ipxt-1+m.LE.nnx)THEN
              vel=vel+veln(ipzt-1+n,ipxt-1+m)*produ
            ENDIF
          ENDDO
        ENDDO
        drx=(rigx-gox)-(ivxt-1)*dvx
        drz=(rigz-goz)-(ivzt-1)*dvz
        v=drx/dvx
        w=drz/dvz
        !
        !           Calculate the 8 basis values at the new point
        !
        vi(1)=(1.0-v)**3/6.0
        vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
        vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
        vi(4)=v**3/6.0
        wi(1)=(1.0-w)**3/6.0
        wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
        wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
        wi(4)=w**3/6.0
        !
        !           Calculate the incremental path length
        !
        IF(k.EQ.1)THEN
          dinc=vrat(k)*dpl
        ELSE
          dinc=(vrat(k)-vrat(k-1))*dpl
        ENDIF
        !
        !           Now compute the 16 contributions to the partial
        !           derivatives.
        !
        DO l=1,4
          DO m=1,4
            rd1=vi(m)*wi(l)/vel**2
            rd2=vio(m)*wio(l)/velo**2
            rd1=-(rd1+rd2)*dinc/2.0
            !fang!                 rd1=vi(m)*wi(l)
            !fang!                 rd2=vio(m)*wio(l)
            !fang!                 rd1=(rd1+rd2)*dinc/2.0
            rd2=fdm(ivzt-2+l,ivxt-2+m)
            fdm(ivzt-2+l,ivxt-2+m)=rd1+rd2
          ENDDO
        ENDDO
      ENDDO
      !fang!     ENDIF
      !fang!      IF(j.EQ.maxrp.AND.sw.EQ.0)THEN
      !fang!         WRITE(6,*)'Error with ray path detected!!!'
      !fang!         WRITE(6,*)'Source id: ',csid
      !fang!         WRITE(6,*)'Receiver id: ',i
      !fang!      ENDIF
    ENDDO
    !
    !  Write ray paths to output file
    !
    !fang!   IF(wrgf.EQ.csid.OR.wrgf.LT.0)THEN
    !if(writepath == 1) then
  
    ! nqdu
    ! save raypath in a file
    if(writepath.eq.1) then
      write(outfile,*)Tid
        if(igr==0) then
            outfile='raypath/'//Trim(adjustl(outpath))//'/phase/'//trim(adjustl(outfile))//'.dat'
        else
            outfile='raypath/'//Trim(adjustl(outpath))//'/group/'//trim(adjustl(outfile))//'.dat'
      endif
      open(40,position='append',file=outfile)
      WRITE(40,*)'#',nrp
          DO j=1,nrp
          rayx=(pi/2-rgx(j))*180.0/pi
          rayz=rgz(j)*180.0/pi
          WRITE(40,*)rayx,rayz
        ENDDO
      close(40)
    endif
    !endif
    !fang!   ENDIF
    !
    !  Write partial derivatives to output file
    !
    !fang!   IF(cfd.EQ.1)THEN
    !fang!!
    !fang!!     Determine the number of non-zero elements.
    !fang!!
    !fang!      isum=0
    !fang!      DO j=0,nvz+1
    !fang!         DO k=0,nvx+1
    !fang!            IF(ABS(fdm(j,k)).GE.ftol)isum=isum+1
    !fang!         ENDDO
    !fang!      ENDDO
    !fang!      WRITE(50)isum
    !fang!      isum=0
    !fang!      DO j=0,nvz+1
    !fang!         DO k=0,nvx+1
    !fang!            isum=isum+1
    !fang!            IF(ABS(fdm(j,k)).GE.ftol)WRITE(50)isum,fdm(j,k)
    !fang!         ENDDO
    !fang!      ENDDO
    !fang!   ENDIF
    !fang!ENDDO
    !fang!IF(cfd.EQ.1)THEN
    !fang!   DEALLOCATE(fdm, STAT=checkstat)
    !fang!   IF(checkstat > 0)THEN
    !fang!      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: fdm'
    !fang!   ENDIF
    !fang!ENDIF
    DEALLOCATE(rgx,rgz, STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: rgx,rgz'
    ENDIF
  END SUBROUTINE rpaths
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TYPE: SUBROUTINE
  ! CODE: FORTRAN 90
  ! This subroutine is passed four node values which lie on
  ! the corners of a rectangle and the coordinates of a point
  ! lying within the rectangle. It calculates the value at
  ! the internal point by using bilinear interpolation.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE bilinear(nv,dsx,dsz,biv)
    USE globalp
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(KIND=sp) :: dsx,dsz,biv
    REAL(KIND=sp), DIMENSION(2,2) :: nv
    REAL(KIND=sp) :: produ
    !
    ! nv = four node vertex values
    ! dsx,dsz = distance between internal point and top left node
    ! dnx,dnz = width and height of node rectangle
    ! biv = value at internal point calculated by bilinear interpolation
    ! produ = product variable
    !
    biv=0.0
    DO i=1,2
      DO j=1,2
        produ=(1.0-ABS(((i-1)*dnx-dsx)/dnx))*(1.0-ABS(((j-1)*dnz-dsz)/dnz))
        biv=biv+nv(i,j)*produ
      ENDDO
    ENDDO
  END SUBROUTINE bilinear