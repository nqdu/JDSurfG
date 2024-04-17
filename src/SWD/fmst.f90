module fmst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!public :: fmst_init,fmst_free,fmst_traveltime,fmst_syn,fmst_
real,allocatable :: velnew(:,:)
!$omp threadprivate(velnew)
contains

subroutine fmst_init(nx,ny,goxdf,gozdf,dvxdf,dvzdf) bind(C,name='fmst_init')
  !! initialize fmst
  !! Input Parameters:
  !!    nx,ny         : nx,ny grid dimension for lat/lon direction
  !!    goxdf,gozdf   : upper left grid coordinates in colat/lon, in degree 
  !!    dvxdf,dvzdf   : grid interval, in degree 
  use iso_c_binding, only          : c_int,c_float
  use globalp
  use traveltime
  implicit none

  ! input vars
  integer(c_int),value,INTENT(IN) :: nx,ny ! no. of receivers, grid dimensions

  ! local variables
  integer(c_int)                  :: checkstat
  real(sp),value,INTENT(IN)       :: goxdf,gozdf,dvxdf,dvzdf ! grid information,degree

    ! pass information to global variables
  gdx = 8
  gdz = 8
  asgr=1
  sgdl=8
  earth=6371.0
  fom=1
  snb=0.5
  nvx=nx-2
  nvz=ny-2
  goxd=goxdf - dvxdf
  gozd=gozdf + dvzdf
  dvxd=dvxdf
  dvzd=dvzdf

  ! allocate space
  ALLOCATE(velv(0:nvz+1,0:nvx+1), STAT=checkstat)
  IF(checkstat > 0)THEN
    WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL velv'
  ENDIF

  !
  ! Convert from degrees to radians
  !
  dvx=dvxd*pi/180.0
  dvz=dvzd*pi/180.0
  gox=(90.0-goxd)*pi/180.0
  goz=gozd*pi/180.0
  !
  ! Compute corresponding values for propagation grid.
  !
  nnx=(nvx-1)*gdx+1
  nnz=(nvz-1)*gdz+1
  dnx=dvx/gdx
  dnz=dvz/gdz
  dnxd=dvxd/gdx
  dnzd=dvzd/gdz
  ALLOCATE(veln(nnz,nnx), STAT=checkstat)
  IF(checkstat > 0)THEN
    WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL veln'
  ENDIF

  ALLOCATE(ttn(nnz,nnx), STAT=checkstat)
  IF(checkstat > 0)THEN
    WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL ttn'
  ENDIF
  rbint=0
  !
  ! Allocate memory for node status and binary trees
  !
  ALLOCATE(nsts(nnz,nnx))
  ALLOCATE(btg(NINT(snb*nnx*nnz)))
  if(asgr.EQ.1) then 
    ALLOCATE(velnb(nnz,nnx))
  endif  
end subroutine fmst_init

subroutine fmst_finalize() bind(c,name='fmst_finalize')
  !! deallocate memory used in fmm 
  use globalp
  use traveltime
  implicit none

  !local 
  integer     :: checkstat

  ! tackle problems and deallocate space
  IF(asgr.EQ.1)DEALLOCATE(ttnr,nstsr)
  IF(asgr.EQ.1)THEN
    DEALLOCATE (velnb, STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
    ENDIF
  ENDIF
  deallocate(velv,veln,ttn,nsts,btg,STAT=checkstat)
  if (allocated(velnew)) deallocate(velnew)
  IF(checkstat > 0)THEN
    WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
  ENDIF
end subroutine fmst_finalize

subroutine fmst_travel(velf,srcx,srcz,rcx,rcz,nr,ttime)bind(c,name='fmst_run')
  !! compute source-receivers pair traveltime
  !! Input Parameters:
  !!    nx,ny         : no. of grid nodes in lat and lon direction
  !!    velf(nx*ny)   : dispersion map
  !!    goxdf,gozdf   : upper left grid coordinates, in degree 
  !!    dvxdf,dvzdf   : grid interval, in degree 
  !!    srcx,srcz     : source coordinates(colat,lon), in rad
  !!    nr            : no of receivers
  !!    rcx(nr),rcz(nr)     : receiver coordinates(colat,lon),in rad 
  !!    ttime(nr)   : travel times
  use,intrinsic :: iso_c_binding
  use globalp
  use traveltime
  implicit none

  ! input variables defined
  integer(c_int),value,INTENT(IN) :: nr! no. of receivers
  real(sp),value,INTENT(IN)       :: srcx,srcz ! source coordinates :: colat,lon
  real(sp),INTENT(IN)             :: rcx(nr),rcz(nr) ! receiver coordinates :: colat,lon
  real(dp),INTENT(IN)             :: velf(nvx+2,nvz+2) ! dispersion map
  real(sp),INTENT(INOUT)          :: ttime(nr) ! travel time for this source-receivers pair

  ! local variables
  integer(c_int)                  :: i,j,k,l,urg
  integer(c_int)                  :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
  INTEGER(c_int)                  :: ogx,ogz,grdfx,grdfz,maxbt
  REAL(sp)                        :: x,z,goxb,gozb,dnxb,dnzb

  ! prepare grid for fmm
  sgs = 8
  call gridder(velf)
  x = srcx
  z = srcz
  !
  !  Begin by computing refined source grid if required
  !
  urg=0
  IF(asgr.EQ.1)THEN
    !
    !     Back up coarse velocity grid to a holding matrix
    !
    !ALLOCATE(velnb(nnz,nnx))
    ! MODIFIEDY BY HONGJIAN FANG @ USTC 2014/04/17
    velnb(1:nnz,1:nnx)=veln(1:nnz,1:nnx)
    nnxb=nnx
    nnzb=nnz
    dnxb=dnx
    dnzb=dnz
    goxb=gox
    gozb=goz
    !
    !     Identify nearest neighbouring node to source
    !
    isx=INT((x-gox)/dnx)+1
    isz=INT((z-goz)/dnz)+1
    sw=0
    IF(isx.lt.1.or.isx.gt.nnx)sw=1
    IF(isz.lt.1.or.isz.gt.nnz)sw=1
    IF(sw.eq.1)then
      x=90.0-x*180.0/pi
      z=z*180.0/pi
      WRITE(6,*)"Source lies outside bounds of model (lat,long)= ",x,z
      WRITE(6,*)"TERMINATING PROGRAM!!!"
      STOP
    ENDIF
    IF(isx.eq.nnx)isx=isx-1
    IF(isz.eq.nnz)isz=isz-1
    !
    !     Now find rectangular box that extends outward from the nearest source node
    !     to "sgs" nodes away.
    !
    vnl=isx-sgs
    IF(vnl.lt.1)vnl=1
    vnr=isx+sgs
    IF(vnr.gt.nnx)vnr=nnx
    vnt=isz-sgs
    IF(vnt.lt.1)vnt=1
    vnb=isz+sgs
    IF(vnb.gt.nnz)vnb=nnz
    nrnx=(vnr-vnl)*sgdl+1
    nrnz=(vnb-vnt)*sgdl+1
    drnx=dvx/REAL(gdx*sgdl)
    drnz=dvz/REAL(gdz*sgdl)
    gorx=gox+dnx*(vnl-1)
    gorz=goz+dnz*(vnt-1)
    nnx=nrnx
    nnz=nrnz
    dnx=drnx
    dnz=drnz
    gox=gorx
    goz=gorz
    !
    !     Reallocate velocity and traveltime arrays if nnx>nnxb or
    !     nnz<nnzb.
    !
    IF(nnx.GT.nnxb.OR.nnz.GT.nnzb)THEN
      idm1=nnx
      IF(nnxb.GT.idm1)idm1=nnxb
      idm2=nnz
      IF(nnzb.GT.idm2)idm2=nnzb
      DEALLOCATE(veln,ttn,nsts,btg)
      ALLOCATE(veln(idm2,idm1))
      ALLOCATE(ttn(idm2,idm1))
      ALLOCATE(nsts(idm2,idm1))
      maxbt=NINT(snb*idm1*idm2)
      ALLOCATE(btg(maxbt))
    ENDIF
    !
    !     Call a subroutine to compute values of refined velocity nodes
    !
    CALL bsplrefine
    !
    !     Compute first-arrival traveltime field through refined grid.
    !
    urg=1
    CALL travel(x,z,urg)
    !
    !     Now map refined grid onto coarse grid.
    !
    ALLOCATE(ttnr(nnzb,nnxb))
    ALLOCATE(nstsr(nnzb,nnxb))
    IF(nnx.GT.nnxb.OR.nnz.GT.nnzb)THEN
      idm1=nnx
      IF(nnxb.GT.idm1)idm1=nnxb
      idm2=nnz
      IF(nnzb.GT.idm2)idm2=nnzb
      DEALLOCATE(ttnr,nstsr)
      ALLOCATE(ttnr(idm2,idm1))
      ALLOCATE(nstsr(idm2,idm1))
    ENDIF
    ttnr=ttn
    nstsr=nsts
    ogx=vnl
    ogz=vnt
    grdfx=sgdl
    grdfz=sgdl
    nsts=-1
    DO k=1,nnz,grdfz
      idm1=ogz+(k-1)/grdfz
      DO l=1,nnx,grdfx
        idm2=ogx+(l-1)/grdfx
        nsts(idm1,idm2)=nstsr(k,l)
        IF(nsts(idm1,idm2).GE.0)THEN
          ttn(idm1,idm2)=ttnr(k,l)
        ENDIF
      ENDDO
    ENDDO
    !
    !     Backup refined grid information
    !
    nnxr=nnx
    nnzr=nnz
    goxr=gox
    gozr=goz
    dnxr=dnx
    dnzr=dnz
    !
    !     Restore remaining values.
    !
    nnx=nnxb
    nnz=nnzb
    dnx=dnxb
    dnz=dnzb
    gox=goxb
    goz=gozb
    DO j=1,nnx
      DO k=1,nnz 
        veln(k,j)=velnb(k,j)
      ENDDO
    ENDDO
    !
    !     Ensure that the narrow band is complete; if
    !     not, then some alive points will need to be
    !     made close.
    !
    DO k=1,nnx
      DO l=1,nnz
        IF(nsts(l,k).EQ.0)THEN
          IF(l-1.GE.1)THEN
            IF(nsts(l-1,k).EQ.-1)nsts(l,k)=1
          ENDIF
          IF(l+1.LE.nnz)THEN
            IF(nsts(l+1,k).EQ.-1)nsts(l,k)=1
          ENDIF
          IF(k-1.GE.1)THEN
            IF(nsts(l,k-1).EQ.-1)nsts(l,k)=1
          ENDIF
          IF(k+1.LE.nnx)THEN
            IF(nsts(l,k+1).EQ.-1)nsts(l,k)=1
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !
    !     Finally, call routine for computing traveltimes once
    !     again.
    !
    urg=2
    CALL travel(x,z,urg)
  ELSE
    !
    !     Call a subroutine that works out the first-arrival traveltime
    !     field.
    !
    CALL travel(x,z,urg)
  ENDIF

  !
  !  Find source-receiver traveltimes if required
  !
  do i=1,nr
    call srtimes(x,z,rcx(i),rcz(i),ttime(i))
  enddo

  IF(rbint.EQ.1)THEN
    WRITE(6,*)'Note that at least one two-point ray path'
    WRITE(6,*)'tracked along the boundary of the model.'
    WRITE(6,*)'This class of path is unlikely to be'
    WRITE(6,*)'a true path, and it is STRONGLY RECOMMENDED'
    WRITE(6,*)'that you adjust the dimensions of your grid'
    WRITE(6,*)'to prevent this from occurring.'
  ENDIF
  
end subroutine fmst_travel

subroutine fmst_reset(velin) bind(c,name='fmst_reset')
  use globalp,only : nvx,nvz,sp,nnx,nnz,veln
  use iso_c_binding

  real(kind=c_double),intent(in) :: velin(nvx+2,nvz+2)
  real(sp)                       :: temp(nnz,nnx)

  ! allocate space 
  if (.not. allocated(velnew)) allocate(velnew(nnz,nnx))

  temp = veln
  call gridder(velin)
  velnew = veln; veln = temp 

end subroutine fmst_reset

function fmst_raypath(scx,scz,surfrcx,surfrcz,fdm) result(tsyn) bind(c,name='fmst_raypath')
  USE globalp
  use iso_c_binding
  
  IMPLICIT NONE

  interface 
    function delsph(colatrad1,lonrad1,colatrad2,lonrad2) result(s) bind(c,name="delsph")
      use iso_c_binding,only : c_float
      implicit none
      real(c_float),value :: colatrad1,lonrad1,colatrad2,lonrad2
      real(c_float) :: s
    end function delsph
  end interface 
  
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
  REAL(KIND=sp) :: v,w,rigz,rigx,dinc
  REAL(KIND=sp) :: dtx,dtz,drx,drz,produ,sred
  REAL(KIND=sp), DIMENSION (:), ALLOCATABLE :: rgx,rgz
  !fang!REAL(KIND=i5), DIMENSION (:,:), ALLOCATABLE :: fdm
  REAL(KIND=sp), DIMENSION (4) :: vrat,vi,wi,vio,wio

  !nqdu
  REAL(KIND=sp),value,intent(in) ::  surfrcx,surfrcz,scx,scz
  real(kind=sp),INTENT(inout)  :: fdm(0:nvx+1,0:nvz+1)

  ! nqdu local
  real(kind=sp)       :: tsyn
  real(sp)            :: fdm0(0:nvz+1,0:nvx+1)
  

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
  fdm0(:,:)=0.0
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
    rayx=90.0-surfrcx*180.0/pi
    rayz=surfrcz*180.0/pi
    WRITE(6,*)"rpath Receiver lies outside model (lat,long)= ",rayx,rayz
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
          vel=vel+velnew(ipzo-1+m,ipxo-1+l)*produ
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
            vel=vel+velnew(ipzt-1+n,ipxt-1+m)*produ
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
          rd2=fdm0(ivzt-2+l,ivxt-2+m)
          fdm0(ivzt-2+l,ivxt-2+m)=rd1+rd2
        ENDDO
      ENDDO
      
      ! travel time
      ! rayx=(pi/2-rgx(j))*180.0/pi
      ! rayz=rgz(j)*180.0/pi
      ! print*,rayz,rayx 
      ! rayx=(pi/2-rgx(j+1))*180.0/pi
      ! rayz=rgz(j+1)*180.0/pi
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

  ! compute travel time
  tsyn = 0
  do j=1,nrp-1
    ! distance 
    dinc = delsph(rgx(j),rgz(j),rgx(j+1),rgz(j+1))

    ! old points
    ipxo=INT((rgx(j)-gox)/dnx)+1
    ipzo=INT((rgz(j)-goz)/dnz)+1
    
    drx = (rgx(j)-gox)-(ipxo-1)*dnx
    drz = (rgz(j)-goz)-(ipzo-1)*dnz
    vel=0.0
    DO l=1,2
      DO m=1,2
        produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
        produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
        IF(ipzo-1+m.LE.nnz.AND.ipxo-1+l.LE.nnx)THEN
          vel=vel+velnew(ipzo-1+m,ipxo-1+l)*produ
        ENDIF
      ENDDO
    ENDDO
    velo = vel
    
    ipxo=INT((rgx(j+1)-gox)/dnx)+1
    ipzo=INT((rgz(j+1)-goz)/dnz)+1
    drx = (rgx(j+1)-gox)-(ipxo-1)*dnx
    drz = (rgz(j+1)-goz)-(ipzo-1)*dnz
    vel=0.0
    DO l=1,2
      DO m=1,2
        produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
        produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
        IF(ipzo-1+m.LE.nnz.AND.ipxo-1+l.LE.nnx)THEN
          vel=vel+velnew(ipzo-1+m,ipxo-1+l)*produ
        ENDIF
      ENDDO
    ENDDO
    !print*,velo,vel
    tsyn = tsyn + dinc * 0.5 * (1. / vel + 1. / velo)
  enddo

  ! ! write raypath
  ! open(40,file='raypath.txt',position='append')
  ! !WRITE(40,*)'#',nrp
  ! DO j=1,nrp
  !   !tsyn = tsyn + delsph(rgx(j),rgz(j),rgx(j+1),rgz(j+1)) / 3.
  !   rayx=(pi/2-rgx(j))*180.0/pi
  !   rayz=rgz(j)*180.0/pi
  !   WRITE(40,*)rayz,rayx
  ! ENDDO
  ! write(40,'(a)') '>'
  ! close(40)

  DEALLOCATE(rgx,rgz, STAT=checkstat)
  IF(checkstat > 0)THEN
    WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: rgx,rgz'
  ENDIF

  ! copy fdm0 to fdm
  fdm = transpose(fdm0)
END function fmst_raypath

end module fmst 