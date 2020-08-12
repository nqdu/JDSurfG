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
end

module sparsemat
use,intrinsic                 :: iso_c_binding
implicit none

real(c_float),ALLOCATABLE     :: val(:)
INTEGER(c_int),ALLOCATABLE    :: rw(:),col(:)
INTEGER                       :: nonzeros
  
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
  use empirical
  implicit none
  integer(c_int),INTENT(INOUT)           :: nx,ny,nz,m
  real(c_float),INTENT(IN)               :: vr(nx,ny,nz),vt(nx,ny,nz)
  real(c_float),INTENT(INOUT)            :: dsyn(m)

  ! local variables
  real(c_float)                          :: rho,vp,vs
  real(c_float),ALLOCATABLE              :: drho(:)
  INTEGER(c_int)                         :: i,j,k,nn,n

  n = (nx -2) * (ny -2) * (nz -1)
  ALLOCATE(drho(n))

  do k=1,nz-1
    do j=1,ny-2
      do i=1,nx-2
        nn = k * (nx-2) * (ny -2) + j * (nx-2) + i
        vs = vr(i+1,j+1,k)
        call empirical_relation(vs,vp,rho)
        drho(nn) = rho
        vs = vt(i+1,j+1,k)
        call empirical_relation(vs,vp,rho)
        drho(nn) = drho(nn) - rho
      enddo
    enddo
  enddo

  dsyn(1:m) = 0.0

  call aprod(drho,dsyn,m,n)
  DEALLOCATE(drho)
  return;
end subroutine forward

end


program main
  use sparsemat
  implicit none

  CHARACTER(50)               :: refmod,truemod
  real,ALLOCATABLE            :: vr(:,:,:),vt(:,:,:),dsyn(:),lat(:),lon(:)
  INTEGER                     :: nx,ny,nz,m,eof
  real(c_float)               :: a,b,c,mean
  INTEGER                     :: i,j,k

  ! read no. of nodes in each direction
  open(11,file='DSurfTomo.in')
  read(11,*)nx,ny,nz
  close(11)
  ALLOCATE(vr(nx,ny,nz),vt(nx,ny,nz))

  ! read reference model,format is same to MOD
  ! you could change it to your own reference file
  refmod='MOD'
  open(11,file=refmod)
  read(11,*)(vr(0,0,i),i=1,nz)
  do k=1,nz
    do j=1,ny
      read(11,*)(vr(i,j,k),i=1,nx)
    enddo
    ! remove mean if required
    !mean = sum(vr(:,:,k)) / real(nx*ny)
    !vr(:,:,k) = vr(:,:,k) - mean
  enddo
  close(11)

  ! read true model( read from results)
  truemod='results/mod_iter_10.dat'
  open(11,file=truemod)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        read(11,*)a,b,c,vt(i,j,k)
      enddo
    enddo
  enddo
  close(11)

  ! read gravity data
  m = 0
  open(11,file='obsgrav.dat')
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
  open(11,file='obsgrav.dat')
  do i=1,m
    read(11,*)lon(i),lat(i)
  enddo
  close(11)

  ! read gravity matrix
  nonzeros = 0
  open(11,file='gravmat.dat')
  do
    read(11,*,iostat=eof)
    if(eof==0) then
      nonzeros = nonzeros + 1
    else
      exit
    endif
  enddo
  close(11)
  ALLOCATE(rw(nonzeros),val(nonzeros),col(nonzeros))
  open(11,file='gravmat.dat')
  do i=1,nonzeros
    read(11,*)rw(i),col(i),val(i)
  enddo
  close(11)

  ! compute
  call forward(nx,ny,nz,vr,vt,m,dsyn)
  mean = sum(dsyn(1:m)) / real(m)
  dsyn(1:m) = dsyn(1:m) - mean

  ! write
  open(11,file='out.dat')
  do i=1,m
    write(11,'(3f10.4)')lon(i),lat(i),dsyn(i)
  enddo
  close(11)

  DEALLOCATE(rw,col,val,lat,lon,vr,vt,dsyn)
end program main
