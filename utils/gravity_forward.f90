module empirical
implicit none
public:: empirical_relation

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
end module empirical

module sparsemat
use,intrinsic                 :: iso_c_binding
implicit none

real(c_float),ALLOCATABLE     :: val(:)
INTEGER(c_int),ALLOCATABLE    :: rw(:),col(:)
INTEGER(c_int)                :: nonzeros
  
contains
subroutine aprod(x,y,m,n) ! y = y + A * x
  implicit none
  INTEGER,INTENT(IN)          :: m,n
  real(c_float),INTENT(INOUT) :: y(m),x(n)

  INTEGER                     ::i

  y(1:m) = 0.0
  do i=1,nonzeros
    y(rw(i)) = y(rw(i)) + val(i) * x(col(i))
  enddo
end subroutine aprod

subroutine forward(nx,ny,nz,vr,vt,m,dsyn)
  use empirical,only                     : empirical_relation
  implicit none
  integer(c_int),INTENT(IN)              :: nx,ny,nz,m
  real(c_float),INTENT(IN)               :: vr(nx,ny,nz),vt(nx,ny,nz)
  real(c_float),INTENT(INOUT)            :: dsyn(m)

  ! local variables
  real(c_float)                          :: rho,vp,vs
  real(c_float),ALLOCATABLE              :: drho(:)
  INTEGER(c_int)                         :: i,j,k,nn,n

  n = (nx -2) * (ny -2) * (nz -1)
  ALLOCATE(drho(n))

  ! compute density difference
  do k=1,nz-1
    do j=1,ny-2
      do i=1,nx-2
        nn = (k-1) * (nx-2) * (ny -2) + (j-1) * (nx-2) + i
        vs = vr(i+1,j+1,k)
        call empirical_relation(vs,vp,rho)
        drho(nn) = rho
        vs = vt(i+1,j+1,k)
        call empirical_relation(vs,vp,rho)
        drho(nn) = -drho(nn) + rho
      enddo
    enddo
  enddo
  
  ! compute gravity
  dsyn(1:m) = 0.0
  call aprod(drho,dsyn,m,n)

  DEALLOCATE(drho)
  return;
end subroutine forward

end module sparsemat

program main
  use sparsemat
  implicit none

  CHARACTER(50)               :: refmod,truemod,paramfile,obsgrav,gravmat
  CHARACTER(50)               :: outfile,flag_remove_mean
  real(c_float),ALLOCATABLE   :: vr(:,:,:),vt(:,:,:),dsyn(:),lat(:),lon(:)
  INTEGER(c_int)              :: nx,ny,nz,m,eof
  real(c_float)               :: a,b,c,mean
  INTEGER(c_int)              :: i,j,k,dataid,num,rm_mean
  CHARACTER(256)              :: line
  CHARACTER                   :: dummy

  ! extract input args
  if(iargc() /= 7) then
    print*,"please run ./this paramfile refmod truemod obsgrav gravmat &
            remove_mean(false or true) outfile"
    stop
  else
    call getarg(1,paramfile)
    call getarg(2,refmod)
    call getarg(3,truemod)
    call getarg(4,obsgrav)
    call getarg(5,gravmat)
    call getarg(6,flag_remove_mean)
    call getarg(7,outfile)

    if(flag_remove_mean == "true") then 
        rm_mean = 1
    else if (flag_remove_mean == "false")then
        rm_mean = 0
    else
        print*, "remove_mean should be one of [true,false]!"
        stop
    endif
  endif

  ! read no. of nodes in each direction
  open(11,file=paramfile)
  read(11,*)nx,ny,nz
  close(11)
  ALLOCATE(vr(nx,ny,nz),vt(nx,ny,nz))

  ! read reference model,format is same to MOD
  open(11,file=refmod)
  read(11,*)(vr(1,1,i),i=1,nz)
  do k=1,nz
    do j=1,ny
      read(11,*)(vr(i,j,k),i=1,nx)
    enddo
    ! remove mean if required
    if(rm_mean == 1) then
      mean = sum(vr(2:nx-1,2:ny-1,k)) / real((nx-2)*(ny-2))
      vr(:,:,k) = mean
    endif
  enddo
  close(11)

  ! read true model( read from results)
  open(11,file=truemod)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        read(11,*)a,b,c,vt(i,j,k)
      enddo
    enddo
  enddo
  close(11)
  print*,'finish reading reference and final model'

  ! read gravity data
  m = 0
  open(11,file=obsgrav)
  do
    read(11,*,iostat=eof)
    if(eof==0) then
      m = m + 1
    else
      exit
    endif
  enddo
  close(11)
  ALLOCATE(lon(m),lat(m),dsyn(m))
  open(11,file=obsgrav)
  do i=1,m
    read(11,*)lon(i),lat(i)
  enddo
  close(11)
  print*,'finish reading gravity data'

  ! read gravity matrix
  nonzeros = 0
  open(11,file=gravmat)
  do
    read(11,'(a)',iostat=eof) line
    if(eof==0) then
      if(line(1:1) /= '#') nonzeros = nonzeros + 1
    else
      exit
    endif
  enddo
  close(11)
  ALLOCATE(rw(nonzeros),val(nonzeros),col(nonzeros))
  i = 1
  open(11,file=gravmat)
  do
    read(11,'(a)',iostat=eof) line
    if(eof/=0) exit 
    read(line,*) dummy,dataid,num
    rw(i:i+num-1) = dataid + 1
    do j=0,num - 1
      read(11,*)col(i + j),val(i + j)
      col(i+j) = col(i+j) + 1
    enddo
    i = i + num
  enddo
  close(11)
  print*,'finish reading gravity matrix, now begin forward computing'

  ! compute gravity
  call forward(nx,ny,nz,vr,vt,m,dsyn)
  mean = sum(dsyn(1:m)) / real(m)
  dsyn(1:m) = dsyn(1:m) - mean

  ! write
  open(11,file=outfile)
  do i=1,m
    write(11,'(3f10.4)')lon(i),lat(i),dsyn(i)
  enddo
  close(11)

  DEALLOCATE(rw,col,val,lat,lon,vr,vt,dsyn)
end program main
