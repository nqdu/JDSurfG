module fmmsph
implicit none
public :: synthetic,CalFrechet
contains

subroutine synthetic(nx,ny,velf,goxdf,gozdf,dvxdf,dvzdf,&
                      srcx,srcz,rcx,rcz,nr,ttime)bind(c,name="synthetic")
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
  integer(c_int),value,INTENT(IN) :: nr,nx,ny ! no. of receivers, grid dimensions
  real(sp),value,INTENT(IN)       :: goxdf,gozdf,dvxdf,dvzdf ! grid information,degree
  real(sp),value,INTENT(IN)       :: srcx,srcz ! source coordinates :: colat,lon
  real(sp),INTENT(IN)             :: rcx(nr),rcz(nr) ! receiver coordinates :: colat,lon
  real(dp),INTENT(IN)             :: velf(nx*ny) ! dispersion map
  real(sp),INTENT(INOUT)          :: ttime(nr) ! travel time for this source-receivers pair

  ! local variables
  integer(c_int)                  :: i,j,k,l,urg,checkstat
  integer(c_int)                  :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
  INTEGER(c_int)                  :: ogx,ogz,grdfx,grdfz,maxbt
  REAL(sp)                        :: x,z,goxb,gozb,dnxb,dnzb

  ! pass information to global variables
  gdx=8
  gdz=8
  asgr=1
  sgdl=8
  sgs=8
  earth=6371.0
  fom=1
  snb=0.5
  nvx=nx-2
  nvz=ny-2
  goxd=goxdf
  gozd=gozdf
  dvxd=dvxdf
  dvzd=dvzdf

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
  maxbt=NINT(snb*nnx*nnz)
  ALLOCATE(btg(maxbt))

  ! prepare grid for fmm
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
    ALLOCATE(velnb(nnz,nnx))
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

  ! tackle problems and deallocate space
  IF(asgr.EQ.1)DEALLOCATE(ttnr,nstsr)
  IF(rbint.EQ.1)THEN
    WRITE(6,*)'Note that at least one two-point ray path'
    WRITE(6,*)'tracked along the boundary of the model.'
    WRITE(6,*)'This class of path is unlikely to be'
    WRITE(6,*)'a true path, and it is STRONGLY RECOMMENDED'
    WRITE(6,*)'that you adjust the dimensions of your grid'
    WRITE(6,*)'to prevent this from occurring.'
  ENDIF
  IF(asgr.EQ.1)THEN
    DEALLOCATE (velnb, STAT=checkstat)
    IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
    ENDIF
  ENDIF
  deallocate(velv,veln,ttn,nsts,btg)
end subroutine

subroutine CalFrechet(nx,ny,nz,velf,goxdf,gozdf,dvxdf,dvzdf,&
                      srcx,srcz,rcx,rcz,nr,ttime,kernel,frechet)&
                      bind(c,name="CalFrechet")
  !! compute source-receivers pair traveltime and frechet kernel
  !! Input Parameters:
  !!    nx,ny         	: no. of grid nodes in lat and lon direction
  !!	nz 			  	: no. of grid points in z direction
  !!    velf(nx*ny)   	: dispersion map
  !!    goxdf,gozdf   	: upper left grid coordinates, in degree 
  !!    dvxdf,dvzdf   	: grid interval, in degree 
  !!    srcx,srcz     	: source coordinates(colat,lon), in rad
  !!    nr            	: no of receivers
  !!    rcx(nr),rcz(nr)	: receiver coordinates(colat,lon),in rad 
  !!    ttime(nr)   	: travel times
  !!	kernel(nx*ny,nz): 1-D surface wave Frechet kernel
  !!	frechet(n,nr)	: pseudo 3-D frechet kernel for each station pair 
  use,intrinsic :: iso_c_binding
  use globalp
  use traveltime
  implicit none

  ! input variables defined
  integer(c_int),value,INTENT(IN):: nr,nx,ny,nz ! no. of receivers, grid dimensions
  real(sp),value,INTENT(IN)      :: goxdf,gozdf,dvxdf,dvzdf ! grid information,degree
  real(sp),value,INTENT(IN)      :: srcx,srcz ! source coordinates :: colat,lon
  real(sp),INTENT(IN)            :: rcx(nr),rcz(nr) ! receiver coordinates :: colat,lon
  real(dp),INTENT(IN)            :: velf(nx*ny) ! dispersion map
  real(sp),INTENT(INOUT)         :: ttime(nr) ! travel time for this source-receivers pair
  real(dp),INTENT(IN)            :: kernel(nx*ny,nz) ! kernel for surface wave
  real(sp),INTENT(INOUT)         :: frechet((nx-2)*(ny-2)*(nz-1),nr)

  ! local variables
  INTEGER(c_int)                 :: i,j,k,l,urg,checkstat
  INTEGER(c_int)                 :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
  INTEGER(c_int)                 :: ogx,ogz,grdfx,grdfz,maxbt
  REAL(sp)                       :: x,z,goxb,gozb,dnxb,dnzb,cbst(nz-1)
  !real(dp)                       :: velf0(ny*nx)
  integer(c_int)                 :: nparpi
  real(sp),allocatable           :: fdm(:,:)
  integer(c_int)                 :: ii,jj,kk
  real(sp),parameter             :: ftol=1.0e-4

  ! pass information to global variables
  gdx=8
  gdz=8
  asgr=1
  sgdl=8
  sgs=8
  earth=6371.0
  fom=1
  snb=0.5
  goxd=goxdf
  gozd=gozdf
  dvxd=dvxdf
  dvzd=dvzdf
  nvx=nx-2
  nvz=ny-2
  nparpi = (nx-2)*(ny-2)*(nz-1)
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
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE CalFrechet: REAL veln'
  ENDIF

  !
  ! Call a subroutine which reads in the velocity grid
  !
  !CALL gridder(grid)
  !
  ! Read in all source coordinates.
  !
  !
  ! Now work out, source by source, the first-arrival traveltime
  ! field plus source-receiver traveltimes
  ! and ray paths if required. First, allocate memory to the
  ! traveltime field array
  !
  ALLOCATE(ttn(nnz,nnx), STAT=checkstat)
  IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE CalFrechet: REAL ttn'
  ENDIF
  rbint=0
  !
  ! Allocate memory for node status and binary trees
  !
  ALLOCATE(nsts(nnz,nnx))
  maxbt=NINT(snb*nnx*nnz)
  ALLOCATE(btg(maxbt))

  allocate(fdm(0:nvz+1,0:nvx+1),STAT=checkstat)
  IF(checkstat > 0)THEN
    WRITE(6,*)'Error with ALLOCATE: SUBROUTINE CalFrechet :real fdm'
  ENDIF

  ! prepare grid for fmm
  call gridder(velf)
  x = srcx
  z = srcz
  !
  !  Begin by computing refined source grid if required
  !
  urg=0
  IF(asgr.EQ.1)THEN
      !
      !Back up coarse velocity grid to a holding matrix
      !
      ALLOCATE(velnb(nnz,nnx))
      ! MODIFIEDY BY HONGJIAN FANG @ USTC 2014/04/17
      velnb(1:nnz,1:nnx)=veln(1:nnz,1:nnx)
      nnxb=nnx
      nnzb=nnz
      dnxb=dnx
      dnzb=dnz
      goxb=gox
      gozb=goz
      !
      ! Identify nearest neighbouring node to source
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
      ! Now find rectangular box that extends outward from the nearest source node
      ! to "sgs" nodes away.
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
      ! Reallocate velocity and traveltime arrays if nnx>nnxb or
      ! nnz<nnzb.
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
      ! Now map refined grid onto coarse grid.
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
      !ttnr(1:nnzb,1:nnxb)=ttn(1:nnzb,1:nnxb)
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
      !Backup refined grid information
      !
      nnxr=nnx
      nnzr=nnz
      goxr=gox
      gozr=goz
      dnxr=dnx
      dnzr=dnz
      !
      !Restore remaining values.
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
      ! Ensure that the narrow band is complete; if
      ! not, then some alive points will need to be
      ! made close.
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
      ! Finally, call routine for computing traveltimes once
      ! again.
      !
      urg=2
      CALL travel(x,z,urg)
  ELSE
      !
      !Call a subroutine that works out the first-arrival traveltime field.
      !
      CALL travel(x,z,urg)
  ENDIF
  !
  !  Find source-receiver traveltimes and frechet kernel if required
  !
  frechet(1:nparpi,1:nr) = 0.0
  do i=1,nr
    call srtimes(x,z,rcx(i),rcz(i),ttime(i))
    call rpaths(x,z,fdm,rcx(i),rcz(i),0,0,0)
    
    ! compute 3-D pseudo sensitivity kernel, save it in frechet
    do jj=1,nvz
      do kk=1,nvx
        if(abs(fdm(jj,kk)) > ftol) then
          cbst(1:nz-1) = real(kernel(jj*(nvx+2)+kk+1,1:nz-1) * fdm(jj,kk))
          do ii=1,nz-1
            if(abs(cbst(ii)) > ftol) then
             frechet((ii-1)*nvz*nvx + (jj-1) * nvx + kk,i) = cbst(ii)
           endif
          enddo
        endif
        !cbst(1:nz-1) = real(kernel(jj*(nvx+2)+kk+1,1:nz-1) * fdm(jj,kk))
        !do ii = 1,nz-1
        ! if(abs(cbst(ii)) > ftol) then
        !    frechet((ii-1)*nvz*nvx + (jj-1) * nvx + kk,i) = cbst(ii)
        !  endif
        !enddo
      enddo
    enddo
  enddo
    
  ! deallocate space
  IF(asgr.EQ.1)THEN
    DEALLOCATE (velnb, STAT=checkstat)
    IF(checkstat > 0)THEN
    WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
    ENDIF
  ENDIF
  IF(asgr.EQ.1)DEALLOCATE(ttnr,nstsr)

  IF(rbint.EQ.1)THEN
    WRITE(6,*)'Note that at least one two-point ray path'
    WRITE(6,*)'tracked along the boundary of the model.'
    WRITE(6,*)'This class of path is unlikely to be'
    WRITE(6,*)'a true path, and it is STRONGLY RECOMMENDED'
    WRITE(6,*)'that you adjust the dimensions of your grid'
    WRITE(6,*)'to prevent this from occurring.'
  ENDIF
  deallocate(fdm)
  deallocate(velv,veln,ttn,nsts,btg)
end

end