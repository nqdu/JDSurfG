module calsurf
use empirical
contains
subroutine refineGrid2LayerMdl(minthk0,mmax,dep,vp,vs,rho,&
    rmax,rvp,rvs,rrho,rthk)
  !!--------------------------------------------------------------------c
  !!refine grid based model to layerd based model
  !!INPUT:   minthk0: no. of grids inserted in each raw layer
  !!         mmax: number of depth grid points in the model
  !!         dep, vp, vs, rho: the depth-grid model parameters
  !!         rmax: number of layers in the fined layered model
  !!         rdep, rvp, rvs, rrho, rthk: the refined layered velocity model : change in subroutine
  !!         
  implicit none

  integer,intent(in) :: mmax,rmax,minthk0
  real,intent(in) :: dep(*),vp(*),vs(*),rho(*)
  real rvp(rmax),rvs(rmax),rrho(rmax),rthk(rmax)

  real thk
  integer i,j,n,k

  n = minthk0 + 1
  do i = 1, mmax-1
    thk = (dep(i+1)-dep(i)) / n
    do j = 1, n
        k = (i-1)*n + j
        rthk(k) = thk
        rvp(k) = vp(i)+(2*j-1)*(vp(i+1)-vp(i))/(2*n)
        rvs(k) = vs(i)+(2*j-1)*(vs(i+1)-vs(i))/(2*n)
        rrho(k) = rho(i)+(2*j-1)*(rho(i+1)-rho(i))/(2*n)
    enddo
  enddo

  !! half space model
  k = rmax;
  rthk(k) = 0.0
  rvp(k) = vp(mmax)
  rvs(k) = vs(mmax)
  rrho(k) = rho(mmax)        

end subroutine refineGrid2LayerMdl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
! This subroutine is used to compute Fundamental Mode Dispersion Map for a spherical layer
! Here every node in the free surface  is viewed as a layermodel
! Thompson-Haskell Matrix Method is used
! openmpi is also useded
! Input:
!   nx,ny,nz : dimension of velocity grid
!   kmaxRc   : length of tRc
!   vel      : velocity array, shape(nx,ny,nz)
!   iwave    : wavetype : 0 for love and 1 for Rayleigh
!   igr      : velotype : 0 for phase and 1 for group
!   tRc      : period vector
!   minthk   : insert minthk layers in grid-based model
! Output:
!   pvRc     : dispersion map, shape(nx*ny,kmaxRc)
subroutine caldispersion(nx,ny,nz,vel,pvRc, &
    iwave,igr,kmaxRc,tRc,depz,minthk)
    use omp_lib
    implicit none

    integer nx,ny,nz
    real vel(nx,ny,nz)
    integer iwave,igr,iflsph,mode
    integer minthk
    real depz(nz)
    integer kmaxRc
    real*8 tRc(kmaxRc)
    real*8 pvRc(nx*ny,kmaxRc)

    integer rmax
    
    ! private variables
    integer ii,jj,kk,nn
    real vsz(nz),vpz(nz),rhoz(nz)
    real rthk(nz + (nz-1) * minthk),rrho(nz + (nz-1) * minthk)
    real rvp(nz + (nz-1) * minthk),rvs(nz + (nz-1) * minthk)
    real*8 cgRc(kmaxRc)

    iflsph = 1
    mode = 1
    rmax = nz + (nz-1) * minthk
    call omp_set_num_threads( 4 )

    !$omp parallel & 
    !$omp default(shared) &
    !$omp private(jj,ii,kk,nn,vsz,vpz,rhoz,rthk,rrho,rvp,rvs,cgRc)
    !$omp do
    do nn=1,nx*ny
        jj = (nn-1) / nx + 1
        ii = mod(nn-1,nx) + 1
        vsz(1:nz) = vel(ii,jj,1:nz)
        do kk=1,nz
            call empirical_relation(vsz(kk),vpz(kk),rhoz(kk))
        enddo
        call refineGrid2LayerMdl(minthk,nz,depz,vpz,vsz,rhoz,&
                rmax,rvp,rvs,rrho,rthk)
        call surfdisp96(rthk,rvp,rvs,rrho,rmax,&
                    iflsph,iwave,mode,igr,kmaxRc,&
                    tRc,cgRc)
        pvRc((jj-1)*nx+ii,1:kmaxRc)=cgRc(1:kmaxRc)
    enddo
    !$omp end do
    !$omp end parallel
end subroutine caldispersion

! compute dispersion map and sensitivity kernel 
subroutine depthkernel(nx,ny,nz,vel,pvRc,sen_Rc,&
  iwave,igr,kmaxRc,tRc,depz,minthk)
  implicit none

  integer nx,ny,nz
  real vel(nx,ny,nz)
  real*8 sen_Rc(ny*nx,kmaxRc,nz)
  integer iwave,igr
  integer minthk
  real depz(nz)
  integer kmaxRc
  real*8 tRc(kmaxRc)
  real*8 pvRc(nx*ny,kmaxRc)

  integer iflsph,mode,rmax
  real dvs,v0

  ! private variables
  integer ii,jj,k,kk,nn
  real vsz(nz),vpz(nz),rhoz(nz)
  real rthk(nz + (nz-1) * minthk),rrho(nz + (nz-1) * minthk)
  real rvp(nz + (nz-1) * minthk),rvs(nz + (nz-1) * minthk)
  real*8 cgRc(kmaxRc),c1(kmaxRc),c2(kmaxRc)

  iflsph = 1
  mode = 1
  rmax = nz + (nz-1) * minthk
  dvs  = 0.01
  call omp_set_num_threads( 4 )

  !$omp parallel & 
  !$omp default(shared) &
  !$omp private(jj,ii,k,kk,nn,vsz,v0,vpz,rhoz,rthk,rrho,rvp,rvs,cgRc,c1,c2)
  !$omp do
  do nn=1,ny*nx
    ii = mod(nn-1,nx) + 1
    jj = (nn-1) /nx + 1
    ! compute dispersion curve
    vsz(1:nz) = vel(ii,jj,1:nz)
    do kk=1,nz
        call empirical_relation(vsz(kk),vpz(kk),rhoz(kk))
    enddo
    call refineGrid2LayerMdl(minthk,nz,depz,vpz,vsz,rhoz,&
      rmax,rvp,rvs,rrho,rthk)
    call surfdisp96(rthk,rvp,rvs,rrho,rmax,&
          iflsph,iwave,mode,igr,kmaxRc,&
          tRc,cgRc)
    pvRc((jj-1)*nx+ii,1:kmaxRc)=cgRc(1:kmaxRc)

    do k=1,nz
      v0 = vsz(k)
      vsz(k) = v0 * (1 + 0.5 *dvs)
      !do kk=1,nz
      !    call empirical_relation(vsz(kk),vpz(kk),rhoz(kk))
      !enddo
      call empirical_relation(vsz(k),vpz(k),rhoz(k))
      call refineGrid2LayerMdl(minthk,nz,depz,vpz,vsz,rhoz,&
        rmax,rvp,rvs,rrho,rthk)
      call surfdisp96(rthk,rvp,rvs,rrho,rmax,&
        iflsph,iwave,mode,igr,kmaxRc,&
        tRc,c2)

      vsz(k) = v0 * (1 - 0.5 *dvs)
      !do kk=1,nz
      !    call empirical_relation(vsz(kk),vpz(kk),rhoz(kk))
      !enddo
      call empirical_relation(vsz(k),vpz(k),rhoz(k))
      call refineGrid2LayerMdl(minthk,nz,depz,vpz,vsz,rhoz,&
        rmax,rvp,rvs,rrho,rthk)
      call surfdisp96(rthk,rvp,rvs,rrho,rmax,&
        iflsph,iwave,mode,igr,kmaxRc,&
        tRc,c1)
      vsz(k) = v0
      call empirical_relation(v0,vpz(k),rhoz(k))

      sen_Rc((jj-1)*nx+ii,1:kmaxRc,k) = (c2-c1) /(dvs*v0)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

end subroutine depthkernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! synthetic surface wave traveltime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine synthetic(nx,ny,nz,vels,obst, &
    goxdf,gozdf,dvxdf,dvzdf,kmaxRc,kmaxRg,kmaxLc,kmaxLg, &
    tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk, &
    scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf)bind(c,name="synthetic") 

  USE globalp
  USE traveltime
  use,intrinsic :: iso_c_binding

  IMPLICIT NONE
  !CHARACTER (LEN=30) ::grid,frechet
  !CHARACTER (LEN=40) :: sources,receivers,otimes
  !CHARACTER (LEN=30) :: travelt,rtravel,wrays,cdum
  INTEGER :: i,j,k,l,urg
  INTEGER :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
  INTEGER :: ogx,ogz,grdfx,grdfz,maxbt
  REAL(KIND=i10) :: x,z,goxb,gozb,dnxb,dnzb
  !REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE :: scxf,sczf
  !REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE :: rcxf,rczf
  !
  ! sources = File containing source locations
  ! receivers = File containing receiver locations
  ! grid = File containing grid of velocity vertices for
  !        resampling on a finer grid with cubic B-splines
  ! frechet = output file containing matrix of frechet derivatives
  ! travelt = File name for storage of traveltime field
  ! wttf = Write traveltimes to file? (0=no,>0=source id)
  ! fom = Use first-order(0) or mixed-order(1) scheme
  ! nsrc = number of sources
  ! scx,scz = source location in r,x,z
  ! scx,scz = source location in r,x,z
  ! x,z = temporary variables for source location
  ! fsrt = find source-receiver traveltimes? (0=no,1=yes)
  ! rtravel = output file for source-receiver traveltimes
  ! cdum = dummy character variable ! wrgf = write ray geometries to file? (<0=all,0=no,>0=source id.)
  ! wrays = file containing raypath geometries
  ! cfd = calculate Frechet derivatives? (0=no, 1=yes)
  ! tnr = total number of receivers
  ! sgs = Extent of refined source grid
  ! isx,isz = cell containing source
  ! nnxb,nnzb = Backup for nnz,nnx
  ! goxb,gozb = Backup for gox,goz
  ! dnxb,dnzb = Backup for dnx,dnz
  ! ogx,ogz = Location of refined grid origin
  ! gridfx,grdfz = Number of refined nodes per cell
  ! urg = use refined grid (0=no,1=yes,2=previously used)
  ! maxbt = maximum size of narrow band binary tree
  ! otimes = file containing source-receiver association information
  !c-----------------------------------------------------------------
  ! variables defined by Hongjian Fang
  integer(c_int),value,intent(in) :: nx,ny,nz
  integer(c_int),value,intent(in) :: kmax,nsrcsurf,nrcf
  real(c_float),intent(in) :: vels(nx,ny,nz)
  real(c_float),intent(inout) :: obst(*)
  real(c_float),value,intent(in) :: goxdf,gozdf,dvxdf,dvzdf
  integer(c_int),value,intent(in) :: kmaxRc,kmaxRg,kmaxLc,kmaxLg
  real(c_double),intent(in) :: tRc(*),tRg(*),tLc(*),tLg(*)
  integer(c_int),intent(in) :: wavetype(nsrcsurf,kmax),periods(nsrcsurf,kmax),&
                            nrc1(nsrcsurf,kmax),nsrcsurf1(kmax),igrt(nsrcsurf,kmax)
  real(c_float),intent(in) :: scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax),&
                            rcxf(nrcf,nsrcsurf,kmax),rczf(nrcf,nsrcsurf,kmax)
  integer(c_int),value,intent(in) :: minthk
  !integer nparpi

  real(c_float),intent(in) :: depz(nz)
  real*8 pvRc(nx*ny,kmaxRc),pvRg(nx*ny,kmaxRg),pvLc(nx*ny,kmaxLc),pvLg(nx*ny,kmaxLg)
  real*8 velf(ny*nx)
  integer kmax1,kmax2,kmax3,count1
  integer igr
  integer iwave
  integer knumi,srcnum
  real cbst1
  !real gaussian
  !external gaussian
  integer istep
  gdx=5
  gdz=5
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
    WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL ttn'
  ENDIF
  rbint=0
  !
  ! Allocate memory for node status and binary trees
  !
  ALLOCATE(nsts(nnz,nnx))
  maxbt=NINT(snb*nnx*nnz)
  ALLOCATE(btg(maxbt))

  !allocate(fdm(0:nvz+1,0:nvx+1))

  if(kmaxRc.gt.0) then
    iwave=2
    igr=0
    call caldispersion(nx,ny,nz,vels,pvRc, &
      iwave,igr,kmaxRc,tRc,depz,minthk)

      open(62,file='velmap2dRc.dat')
      do k = 1,kmaxRc
      do j=1,ny-2
      do i=1,nx-2
      !write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tRc(k),pvRc((j+1)*nx+i+1,k)
      write(62,'(5f8.4)') gozd+j*dvzd,goxd-i*dvxd,tRc(k),pvRc((j+1)*nx+i+1,k)
      enddo
      enddo
      enddo
      close(62)
  endif

  if(kmaxRg.gt.0) then
    iwave=2
    igr=1
    call caldispersion(nx,ny,nz,vels,pvRg, &
      iwave,igr,kmaxRg,tRg,depz,minthk)
      open(62,file='velmap2dRg.dat')
      do k = 1,kmaxRg
      do j=1,ny-2
      do i=1,nx-2
      !write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tRg(k),pvRg((j+1)*nx+i+1,k)
      write(62,'(5f8.4)') gozd+j*dvzd,goxd-i*dvxd,tRg(k),pvRg((j+1)*nx+i+1,k)
      enddo
      enddo
      enddo
      close(62)
  endif

  if(kmaxLc.gt.0) then
    iwave=1
    igr=0
    call caldispersion(nx,ny,nz,vels,pvLc, &
      iwave,igr,kmaxLc,tLc,depz,minthk)

      open(62,file='velmap2dLc.dat')
      do k = 1,kmaxLc
      do j=1,ny-2
      do i=1,nx-2
      write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tLc(k),pvLc((j+1)*nx+i+1,k)
      enddo
      enddo
      enddo
      close(62)
  endif

  if(kmaxLg.gt.0) then
    iwave=1
    igr=1
    call caldispersion(nx,ny,nz,vels,pvLg, &
      iwave,igr,kmaxLg,tLg,depz,minthk)

      open(62,file='velmap2dLg.dat')
      do k = 1,kmaxLg
      do j=1,ny-2
      do i=1,nx-2
      write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tLg(k),pvLg((j+1)*nx+i+1,k)
      enddo
      enddo
      enddo
      close(62)
  endif


  !nar=0
  count1=0

  !sen_vs=0
  !sen_vp=0
  !sen_rho=0
  kmax1=kmaxRc
  kmax2=kmaxRc+kmaxRg
  kmax3=kmaxRc+kmaxRg+kmaxLc
  do knumi=1,kmax
    do srcnum=1,nsrcsurf1(knumi)
      if(wavetype(srcnum,knumi)==2.and.igrt(srcnum,knumi)==0) then
        velf(1:nx*ny)=pvRc(1:nx*ny,periods(srcnum,knumi))
        !        sen_vs(:,1:kmax1,:)=sen_vsRc(:,1:kmaxRc,:)!(:,nt(istep),:)
        !        sen_vp(:,1:kmax1,:)=sen_vpRc(:,1:kmaxRc,:)!(:,nt(istep),:)
        !        sen_rho(:,1:kmax1,:)=sen_rhoRc(:,1:kmaxRc,:)!(:,nt(istep),:)
      endif
      if(wavetype(srcnum,knumi)==2.and.igrt(srcnum,knumi)==1) then
        velf(1:nx*ny)=pvRg(1:nx*ny,periods(srcnum,knumi))
        !        sen_vs(:,kmax1+1:kmax2,:)=sen_vsRg(:,1:kmaxRg,:)!(:,nt,:)
        !        sen_vp(:,kmax1+1:kmax2,:)=sen_vpRg(:,1:kmaxRg,:)!(:,nt,:)
        !        sen_rho(:,kmax1+1:kmax2,:)=sen_rhoRg(:,1:kmaxRg,:)!(:,nt,:)
      endif
      if(wavetype(srcnum,knumi)==1.and.igrt(srcnum,knumi)==0) then
        velf(1:nx*ny)=pvLc(1:nx*ny,periods(srcnum,knumi))
        !        sen_vs(:,kmax2+1:kmax3,:)=sen_vsLc(:,1:kmaxLc,:)!(:,nt,:)
        !        sen_vp(:,kmax2+1:kmax3,:)=sen_vpLc(:,1:kmaxLc,:)!(:,nt,:)
        !        sen_rho(:,kmax2+1:kmax3,:)=sen_rhoLc(:,1:kmaxLc,:)!(:,nt,:)
      endif
      if(wavetype(srcnum,knumi)==1.and.igrt(srcnum,knumi)==1) then
        velf(1:nx*ny)=pvLg(1:nx*ny,periods(srcnum,knumi))
        !        sen_vs(:,kmax3+1:kmax,:)=sen_vsLg(:,1:kmaxLg,:)!(:,nt,:)
        !        sen_vp(:,kmax3+1:kmax,:)=sen_vpLg(:,1:kmaxLg,:)!(:,nt,:)
        !        sen_rho(:,kmax3+1:kmax,:)=sen_rhoLg(:,1:kmaxLg,:)!(:,nt,:)
      endif

      call gridder(velf)
      x=scxf(srcnum,knumi)
      z=sczf(srcnum,knumi)
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
      !  
      do istep=1,nrc1(srcnum,knumi)
        CALL srtimes(x,z,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),cbst1)
        count1=count1+1
        obst(count1)=cbst1
      enddo



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
    enddo
  enddo
  !deallocate(fdm)
  deallocate(velv,veln,ttn,nsts,btg)
end subroutine synthetic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! synthetic surface wave traveltime and compute Frechet Kernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalSurfG(nx,ny,nz,nonzeros,vels,rw,col,val,dsurf, &
  goxdf,gozdf,dvxdf,dvzdf,kmaxRc,kmaxRg,kmaxLc,kmaxLg, &
  tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk, &
  scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf, &
  nar) bind(c,name="CalSurfG")
  USE globalp
  USE traveltime
  use,intrinsic :: iso_c_binding
  IMPLICIT NONE

  ! variables defined by Hongjian Fang
  ! add c interoperability by Nanqiao Du
  integer(c_int),value,intent(in) :: nx,ny,nz,nonzeros
  integer(c_int),value,intent(in) :: kmax,nsrcsurf,nrcf,minthk
  real(c_float),intent(in) :: vels(nx,ny,nz),depz(nz)
  real(c_float),intent(inout) :: val(*)
  integer(c_int),intent(inout) :: rw(*),col(*)
  real(c_float),intent(inout) :: dsurf(*)
  real(c_float),value,intent(in) ::goxdf,gozdf,dvxdf,dvzdf
  integer(c_int),value,intent(in) :: kmaxRc,kmaxRg,kmaxLc,kmaxLg
  real(c_double),intent(in) ::tRc(*),tRg(*),tLc(*),tLg(*)
  integer(c_int),intent(in) :: wavetype(nsrcsurf,kmax)
  integer(c_int),intent(in) :: periods(nsrcsurf,kmax),nrc1(nsrcsurf,kmax),nsrcsurf1(kmax)
  integer(c_int),intent(in) :: igrt(nsrcsurf,kmax)
  real(c_float),intent(in) ::scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax),&
                            rcxf(nrcf,nsrcsurf,kmax),rczf(nrcf,nsrcsurf,kmax)
  integer(c_int),intent(inout) :: nar

  !CHARACTER (LEN=30) ::grid,frechet
  !CHARACTER (LEN=40) :: sources,receivers,otimes
  !CHARACTER (LEN=30) :: travelt,rtravel,wrays,cdum
  INTEGER :: j,k,l,urg
  INTEGER :: sgs,isx,isz,sw,idm1,idm2,nnxb,nnzb
  INTEGER :: ogx,ogz,grdfx,grdfz,maxbt
  REAL(KIND=i10) :: x,z,goxb,gozb,dnxb,dnzb
  !REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE :: scxf,sczf
  !REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE :: rcxf,rczf
  !
  ! sources = File containing source locations
  ! receivers = File containing receiver locations
  ! grid = File containing grid of velocity vertices for
  !        resampling on a finer grid with cubic B-splines
  ! frechet = output file containing matrix of frechet derivatives
  ! travelt = File name for storage of traveltime field
  ! wttf = Write traveltimes to file? (0=no,>0=source id)
  ! fom = Use first-order(0) or mixed-order(1) scheme
  ! nsrc = number of sources
  ! scx,scz = source location in r,x,z
  ! scx,scz = source location in r,x,z
  ! x,z = temporary variables for source location
  ! fsrt = find source-receiver traveltimes? (0=no,1=yes)
  ! rtravel = output file for source-receiver traveltimes
  ! cdum = dummy character variable ! wrgf = write ray geometries to file? (<0=all,0=no,>0=source id.)
  ! wrays = file containing raypath geometries
  ! cfd = calculate Frechet derivatives? (0=no, 1=yes)
  ! tnr = total number of receivers
  ! sgs = Extent of refined source grid
  ! isx,isz = cell containing source
  ! nnxb,nnzb = Backup for nnz,nnx
  ! goxb,gozb = Backup for gox,goz
  ! dnxb,dnzb = Backup for dnx,dnz
  ! ogx,ogz = Location of refined grid origin
  ! gridfx,grdfz = Number of refined nodes per cell
  ! urg = use refined grid (0=no,1=yes,2=previously used)
  ! maxbt = maximum size of narrow band binary tree
  ! otimes = file containing source-receiver association information
  !c-----------------------------------------------------------------


  real(c_double) :: pvRc(nx*ny,kmax),pvRg(nx*ny,kmaxRg),pvLc(nx*ny,kmax),pvLg(nx*ny,kmaxLg)
  real(c_double) :: sen_Rc(nx*ny,kmaxRc,nz),sen_Rg(nx*ny,kmaxRc,nz)
  real(c_double) :: sen_Lc(nx*ny,kmaxRc,nz),sen_Lg(nx*ny,kmaxRc,nz)
  real(c_double) :: kernel(nx*ny,kmax,nz)
  !nqdu
  !real*8 sen_vsRc(nx*ny,kmaxRc,nz),sen_vpRc(nx*ny,kmaxRc,nz)
  !real*8 sen_rhoRc(nx*ny,kmaxRc,nz)
  !real*8 sen_vsRg(nx*ny,kmaxRg,nz),sen_vpRg(nx*ny,kmaxRg,nz)
  !real*8 sen_rhoRg(nx*ny,kmaxRg,nz)
  !real*8 sen_vsLc(nx*ny,kmaxLc,nz),sen_vpLc(nx*ny,kmaxLc,nz)
  !real*8 sen_rhoLc(nx*ny,kmaxLc,nz)
  !real*8 sen_vsLg(nx*ny,kmaxLg,nz),sen_vpLg(nx*ny,kmaxLg,nz)
  !real*8 sen_rhoLg(nx*ny,kmaxLg,nz)
  !real*8 sen_vs(nx*ny,kmax,nz),sen_vp(nx*ny,kmax,nz)
  !real*8 sen_rho(nx*ny,kmax,nz)
  !real coe_rho(nz-1),coe_a(nz-1)
  real*8 velf(ny*nx),velf0(ny*nx)
  integer(c_int) :: nparpi
  integer kmax1,kmax2,kmax3,count1,count11
  integer igr
  integer iwave
  integer knumi,srcnum
  real(c_float),dimension(:,:),allocatable:: fdm
  real(c_float),allocatable :: row(:)
  real cbst1
  integer jj,kk,nn,istep,i1,j1,k1
  real(c_float),parameter::ftol=1e-4
  integer ig, igroup
  !modified by nqdu
  !character(len=512) ::cfilealpha
  !character(len=512) ::cfilebeta
  !character(len=512) ::cfilerho

  nparpi = (nx -2) * (ny - 2) * (nz -1)
  allocate(row(nparpi))
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
      WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL ttn'
  ENDIF
  rbint=0
  !
  ! Allocate memory for node status and binary trees
  !
  ALLOCATE(nsts(nnz,nnx))
  maxbt=NINT(snb*nnx*nnz)
  ALLOCATE(btg(maxbt))

  allocate(fdm(0:nvz+1,0:nvx+1))
  if(kmaxRc.gt.0) then
      write(*,*)'computing Rayleigh phase Frechet Kernel ...'
      iwave=2
      igr=0
      !call depthkernel(nx,ny,nz,vels,pvRc,sen_vsRc,sen_vpRc, &
      !sen_rhoRc,iwave,igr,kmaxRc,tRc,depz,minthk)
      call depthkernel(nx,ny,nz,vels,pvRc,sen_Rc,iwave,&
                      igr,kmaxRc,tRc,depz,minthk)

      ! save phase velocity
      open(62,file='velmap2dRc.dat')
      do k1 = 1,kmaxRc
          do j1=1,ny-2
              do i1=1,nx-2
                  !write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tRg(k),pvRg((j+1)*nx+i+1,k)
                  write(62,'(5f8.4)') gozd+j1*dvzd,goxd-i1*dvxd,tRc(k1),pvRc((j1+1)*nx+i1+1,k1)
              enddo
          enddo
      enddo
      close(62)
  endif
  !write(*,*)'this is the end of phase velocity depthkernel...'    

  if(kmaxRg.gt.0) then
      write(*,*)'computing Rayleigh Group Frechet Kernel ...'
      iwave=2
      igr=1
      call depthkernel(nx,ny,nz,vels,pvRg,sen_Rg,iwave,&
                      igr,kmaxRg,tRg,depz,minthk)

      ! save group velocity
      open(62,file='velmap2dRg.dat')
      do k1 = 1,kmaxRg
          do j1=1,ny-2
              do i1=1,nx-2
                  !write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tRg(k),pvRg((j+1)*nx+i+1,k)
                  write(62,'(5f8.4)') gozd+j1*dvzd,goxd-i1*dvxd,tRg(k1),pvRg((j1+1)*nx+i1+1,k1)
              enddo
          enddo
      enddo
      close(62)
  endif

  !write(*,*)'this is the end of group velocity depthkernel...'

  if(kmaxLc.gt.0) then
      write(*,*)'computing Love phase Frechet Kernel ...'
      iwave=1
      igr=0
      call depthkernel(nx,ny,nz,vels,pvLc,sen_Lc,iwave,&
                      igr,kmaxLc,tLc,depz,minthk)
      ! save phase velocity
      open(62,file='velmap2dLc.dat')
      do k1 = 1,kmaxRg
          do j1=1,ny-2
              do i1=1,nx-2
                  !write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tRg(k),pvRg((j+1)*nx+i+1,k)
                  write(62,'(5f8.4)') gozd+j1*dvzd,goxd-i1*dvxd,tLc(k1),pvLc((j1+1)*nx+i1+1,k1)
              enddo
          enddo
      enddo
      close(62)
  endif

  if(kmaxLg.gt.0) then
      write(*,*) 'computing Love Group Frechet Kernel ...'
      iwave=1
      igr=0
      !call caldispersion(nx,ny,nz,vels,pvLc, &
      !    iwave,igr,kmax,tLg,depz,minthk)
      igr=1
      call depthkernel(nx,ny,nz,vels,pvLg,sen_Lg,iwave,&
                      igr,kmaxLg,tLg,depz,minthk)
      ! save group velocity
      open(62,file='velmap2dLg.dat')
      do k1 = 1,kmaxRg
          do j1=1,ny-2
              do i1=1,nx-2
                  !write(62,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,tRg(k),pvRg((j+1)*nx+i+1,k)
                  write(62,'(5f8.4)') gozd+j1*dvzd,goxd-i1*dvxd,tLg(k1),pvLg((j1+1)*nx+i1+1,k1)
              enddo
          enddo
      enddo
      close(62)
  endif

  nar=0
  count1=0
  
  write(*,*)'computing surface wave traveltimes and Frechet Kernel ...'
  kernel = 0.0
  kmax1=kmaxRc
  kmax2=kmaxRc+kmaxRg
  kmax3=kmaxRc+kmaxRg+kmaxLc
  !Modified by nqdu
  !for every period knumi and every source index srcnum:
  do knumi=1,kmax 
      do srcnum=1,nsrcsurf1(knumi)
      !store sensitivity matrix in sen_*(ny*nx,kmax,nz)
      !store phase/group velocity
      if(wavetype(srcnum,knumi)==2.and.igrt(srcnum,knumi)==0) then
          velf(1:nx*ny)=pvRc(1:nx*ny,periods(srcnum,knumi))
          kernel(:,1:kmax1,:)=sen_Rc(:,1:kmaxRc,:)!(:,nt(istep),:)
          !sen_vp(:,1:kmax1,:)=sen_vpRc(:,1:kmaxRc,:)!(:,nt(istep),:)
          !sen_rho(:,1:kmax1,:)=sen_rhoRc(:,1:kmaxRc,:)!(:,nt(istep),:)
      endif
      if(wavetype(srcnum,knumi)==2.and.igrt(srcnum,knumi)==1) then
          velf(1:nx*ny)=pvRg(1:nx*ny,periods(srcnum,knumi))
          kernel(:,kmax1+1:kmax2,:)=sen_Rg(:,1:kmaxRg,:)!(:,nt,:)
          !sen_vp(:,kmax1+1:kmax2,:)=sen_vpRg(:,1:kmaxRg,:)!(:,nt,:)
          !sen_rho(:,kmax1+1:kmax2,:)=sen_rhoRg(:,1:kmaxRg,:)!(:,nt,:)
      endif
      if(wavetype(srcnum,knumi)==1.and.igrt(srcnum,knumi)==0) then
          velf(1:nx*ny)=pvLc(1:nx*ny,periods(srcnum,knumi))
          kernel(:,kmax2+1:kmax3,:)=sen_Lc(:,1:kmaxLc,:)!(:,nt,:)
          !sen_vp(:,kmax2+1:kmax3,:)=sen_vpLc(:,1:kmaxLc,:)!(:,nt,:)
          !sen_rho(:,kmax2+1:kmax3,:)=sen_rhoLc(:,1:kmaxLc,:)!(:,nt,:)
      endif
      if(wavetype(srcnum,knumi)==1.and.igrt(srcnum,knumi)==1) then
          velf(1:nx*ny)=pvLg(1:nx*ny,periods(srcnum,knumi))
          kernel(:,kmax3+1:kmax,:)=sen_Lg(:,1:kmaxLg,:)!(:,nt,:)
          !sen_vp(:,kmax3+1:kmax,:)=sen_vpLg(:,1:kmaxLg,:)!(:,nt,:)
          !sen_rho(:,kmax3+1:kmax,:)=sen_rhoLg(:,1:kmaxLg,:)!(:,nt,:)
      endif

      ! only for Rayleigh wave group velocity, revise this later for Love wave group velocity
      if (igrt(srcnum,knumi)==1) then
          igroup = 2
      else
          igroup = 1
      endif
      velf0 = velf
      count11 = count1
      do ig = 1,igroup
      if (ig ==2 .and. wavetype(srcnum,knumi) == 2) then
          velf(1:nx*ny) = pvRc(1:nx*ny,periods(srcnum,knumi))
      endif
      if (ig ==2 .and. wavetype(srcnum,knumi) == 1) then
          velf(1:nx*ny) = pvLc(1:nx*ny,periods(srcnum,knumi))
      endif
      call gridder(velf)
      x=scxf(srcnum,knumi)
      z=sczf(srcnum,knumi)
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
      !  
      do istep=1,nrc1(srcnum,knumi)
          if (ig == 1) then
          CALL srtimes(x,z,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),cbst1)
          count1=count1+1
          dsurf(count1)=cbst1
          endif
          !!-------------------------------------------------------------
          !   ENDIF
          !
          !  Calculate raypath geometries and write to file if required.
          !  Calculate Frechet derivatives with the same subroutine
          !  if required.
          !
          if (igrt(srcnum,knumi) == 0 .or. (ig == 2 .and. igrt(srcnum,knumi) == 1)) then
          ! a little stupid, remember to change latter
          if (igrt(srcnum,knumi) == 1) then
          call gridder(velf0)
          endif
          count11=count11+1
          CALL rpaths(x,z,fdm,rcxf(istep,srcnum,knumi),rczf(istep,srcnum,knumi),&
                      igrt(srcnum,knumi),periods(srcnum,knumi),0)
          row(1:nparpi)=0.0
          do jj=1,nvz
              do kk=1,nvx
                  if(abs(fdm(jj,kk)).ge.ftol) then
                      row((jj-1)*nvx+kk:(nz-2)*nvz*nvx+(jj-1)*nvx+kk:nvx*nvz)=&
                              real(kernel(jj*(nvx+2)+kk+1,knumi,1:nz-1) * fdm(jj,kk))
                          !(sen_vp(jj*(nvx+2)+kk+1,knumi,1:nz-1)*coe_a+&
                          !sen_rho(jj*(nvx+2)+kk+1,knumi,1:nz-1)*coe_rho+&
                          !sen_vs(jj*(nvx+2)+kk+1,knumi,1:nz-1))*fdm(jj,kk)
                  endif
              enddo
          enddo
          do nn=1,nparpi
          if(abs(row(nn)).gt.ftol) then
              nar=nar+1
              if(nar > nonzeros) THEN
                write(*,'(a)')'Warnings! please increase the sparse ratio'
                stop;
              endif
              !rw(nar)=real(row(nn))
              val(nar) = row(nn)
              col(nar) = nn - 1 
              rw(nar) = count11 - 1 
              !iw(nar+1)= count11
              !col(nar)=nn
          endif
          enddo
      endif ! 'if' before rpath
      enddo
      IF(asgr.EQ.1)THEN
          DEALLOCATE (velnb, STAT=checkstat)
          IF(checkstat > 0)THEN
          WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin2d: velnb'
          ENDIF
      ENDIF
      IF(asgr.EQ.1)DEALLOCATE(ttnr,nstsr)
      enddo ! 'do' before gridder 


      IF(rbint.EQ.1)THEN
          WRITE(6,*)'Note that at least one two-point ray path'
          WRITE(6,*)'tracked along the boundary of the model.'
          WRITE(6,*)'This class of path is unlikely to be'
          WRITE(6,*)'a true path, and it is STRONGLY RECOMMENDED'
          WRITE(6,*)'that you adjust the dimensions of your grid'
          WRITE(6,*)'to prevent this from occurring.'
      ENDIF
      enddo
  enddo
  deallocate(fdm)
  deallocate(row)
  deallocate(velv,veln,ttn,nsts,btg)
end subroutine CalSurfG

end module
