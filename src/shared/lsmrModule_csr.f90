!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrDataModule.f90
!
! Defines real(dp) and a few constants for use in other modules.
!
! 24 Oct 2007: Allows floating-point precision dp to be defined
!              in exactly one place (here).  Note that we need
!                 use lsmrDataModule
!              at the beginning of modules AND inside interfaces.
!              zero and one are not currently used by LSMR,
!              but this shows how they should be declared
!              by a user routine that does need them.
! 16 Jul 2010: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsmrDataModule
   use,intrinsic :: iso_c_binding
   implicit none
 
   intrinsic                   ::      selected_real_kind
   !nqdu
   integer(c_int),  parameter, public :: dp = c_float
   real(dp), parameter, public :: zero = 0.0_dp, one = 1.0_dp
 
end module lsmrDataModule

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrblasInterface.f90
!
!    BLAS1 Interfaces:   ddot    dnrm2    dscal
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 19 Dec 2008: lsqrblasInterface module implemented.
!              Metcalf and Reid recommend putting interfaces in a module.
! 16 Jul 2010: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsmrblasInterface
   implicit none
   public   :: ddot, dnrm2, dscal
 
   interface                              ! Level 1 BLAS
      function ddot  (n,dx,incx,dy,incy)
        use lsmrDataModule, only : dp
        use,intrinsic :: iso_c_binding

        integer(c_int),  intent(in)    :: n,incx,incy
        real(dp), intent(in)    :: dx(*),dy(*)
        real(dp)                :: ddot
      end function ddot
 
      function dnrm2 (n,dx,incx)
        use lsmrDataModule, only : dp
        use,intrinsic :: iso_c_binding
        integer(c_int),  intent(in)    :: n,incx
        real(dp), intent(in)    :: dx(*)
        real(dp)                :: dnrm2
      end function dnrm2
 
      subroutine dscal (n,sa,x,incx)
        use lsmrDataModule, only : dp
        use,intrinsic :: iso_c_binding
        integer(c_int),  intent(in)    :: n,incx
        real(dp), intent(in)    :: sa
        real(dp), intent(inout) :: x(*)
      end subroutine dscal
   end interface
 
 end module lsmrblasInterface
 
!nqdu
module lsmr_csr_matrix
use,intrinsic :: iso_c_binding
implicit none

interface aprod_funcs
subroutine myaprod_csr(mode,x,y,indices,indptr,row,col,val,nprocs) bind(c,name="myaprod_csr")
   use,intrinsic :: iso_c_binding
   implicit none
   integer,PARAMETER                :: dp = c_float
   integer(c_int),value,intent(in)  :: mode,row,col,nprocs
   integer(c_int),INTENT(IN)        :: indices(*),indptr(*)
   real(dp),INTENT(IN)              :: val(*)
   real(dp),INTENT(INOUT)           :: x(col),y(row)
end subroutine myaprod_csr
end interface 

contains

subroutine aprod_parallel(mode, m, n, x, y,val,indices,indptr,nprocs)
   !!compute y = y + A *x or x = x + A.T * y 
   !! parallel version
   use,intrinsic :: iso_c_binding 
   implicit none
   integer,PARAMETER         :: dp = c_float
   integer(c_int),intent(in) :: mode,m,n,nprocs
   integer(c_int),INTENT(IN) :: indices(*),indptr(*)
   real(dp),INTENT(IN)       :: val(*)
   real(dp),INTENT(INOUT)    :: x(n),y(m)

   call myaprod_csr(mode,x,y,indices,indptr,m,n,val,nprocs)
   
end subroutine aprod_parallel

end module lsmr_csr_matrix

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File lsmrblas.f90   (double precision)
!
!     This file contains the following BLAS routines
!        dcopy, ddot, dnrm2, dscal
!     required by subroutines LSMR and Acheck.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!! DCOPY copies a vector X to a vector Y.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!    The routine uses unrolled loops for increments equal to one.
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries of DX.
!
!    Output, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries of DY.


      subroutine  dcopy(n,dx,incx,dy,incy)
         use lsmrDataModule, only : dp
         use,intrinsic :: iso_c_binding
      implicit none
!      double precision dx(*),dy(*)
      real(dp) dx(*),dy(*)
      integer(c_int) i,incx,incy,ix,iy,m,n

      if ( n <= 0 ) then
         return
      end if

      if ( incx == 1 .and. incy == 1 ) then

         m = mod ( n, 7 )

         if ( m /= 0 ) then
            dy(1:m) = dx(1:m)
         end if

         do i = m+1, n, 7
            dy(i) = dx(i)
            dy(i + 1) = dx(i + 1)
            dy(i + 2) = dx(i + 2)
            dy(i + 3) = dx(i + 3)
            dy(i + 4) = dx(i + 4)
            dy(i + 5) = dx(i + 5)
            dy(i + 6) = dx(i + 6)
         end do

        else

           if ( 0 <= incx ) then
              ix = 1
           else
              ix = ( -n + 1 ) * incx + 1
           end if

           if ( 0 <= incy ) then
              iy = 1
           else
              iy = ( -n + 1 ) * incy + 1
           end if

           do i = 1, n
              dy(iy) = dx(ix)
              ix = ix + incx
              iy = iy + incy
           end do
        end if
        return
end subroutine dcopy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!    This routine uses unrolled loops for increments equal to one.
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries in DX.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries in DY.
!
!    Output, real ( kind = 8 ) DDOT, the sum of the product of the 
!    corresponding entries of DX and DY.


     ! double precision function ddot(n,dx,incx,dy,incy)
      real(dp) function ddot(n,dx,incx,dy,incy)
      use lsmrDataModule, only : dp
      use,intrinsic :: iso_c_binding
      implicit         none
    ! double precision dx(*),dy(*),dtemp
      integer(c_int),  intent(in)    :: n,incx,incy
      real(dp), intent(in)    :: dx(*),dy(*)
      real(dp) dtemp
      integer(c_int)          i,ix,iy,m

      ddot = 0.0d0
      dtemp = 0.0d0
      if ( n <= 0 ) then
         return
      end if

!  Code for unequal increments or equal increments
!  not equal to 1.

      if ( incx /= 1 .or. incy /= 1 ) then

         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         if ( 0 <= incy ) then
            iy = 1
         else
            iy = ( - n + 1 ) * incy + 1
         end if

         do i = 1, n
            dtemp = dtemp + dx(ix) * dy(iy)
            ix = ix + incx
            iy = iy + incy
         end do

!  Code for both increments equal to 1.

        else

           m = mod ( n, 5 )

           do i = 1, m
              dtemp = dtemp + dx(i) * dy(i)
           end do

           do i = m+1, n, 5
              dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) &
                                          + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
           end do

        end if

        ddot = dtemp
        return
end function ddot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Sven Hammarling
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!

    !  double precision function dnrm2 ( n, x, incx)
      real(dp) function dnrm2 ( n, x, incx)
      use lsmrDataModule, only : dp
      use,intrinsic :: iso_c_binding
      implicit         none
      integer(c_int),  intent(in)    :: n,incx
      real(dp), intent(in)    :: x(*)
      
      integer(c_int)          ix
     ! double precision x(*), ssq,absxi,norm,scale
      real(dp) ssq,absxi,norm,scale

      if ( n < 1 .or. incx < 1 ) then
         norm  = 0.d0
      else if ( n == 1 ) then
         norm  = abs ( x(1) )
      else
         scale = 0.d0
         ssq = 1.d0

         do ix = 1, 1 + ( n - 1 )*incx, incx
            if ( x(ix) /= 0.d0 ) then
               absxi = abs ( x(ix) )
               if ( scale < absxi ) then
                  ssq = 1.0 + ssq * ( scale / absxi )**2
                  scale = absxi
               else
                  ssq = ssq + ( absxi / scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt ( ssq )
      end if

      dnrm2 = norm
      return
end function dnrm2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DSCAL scales a vector by a constant.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!
!  Modified:
!    08 April 1999
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries of X.
!

      subroutine  dscal(n,sa,x,incx)
         use lsmrDataModule, only : dp
         use,intrinsic :: iso_c_binding
      implicit none

      integer(c_int),  intent(in)    :: n,incx
      real(dp), intent(in)    :: sa
      real(dp), intent(inout) :: x(*)

      integer(c_int) i
      integer(c_int) ix
      integer(c_int) m
      !double precision sa
      !double precision x(*)
      if ( n <= 0 ) then
         return
      else if ( incx == 1 ) then
         m = mod ( n, 5 )
         x(1:m) = sa * x(1:m)

         do i = m+1, n, 5
            x(i)   = sa * x(i)
            x(i+1) = sa * x(i+1)
            x(i+2) = sa * x(i+2)
            x(i+3) = sa * x(i+3)
            x(i+4) = sa * x(i+4)
         end do
      else
         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         do i = 1, n
            x(ix) = sa * x(ix)
            ix = ix + incx
         end do

      end if

      return
end subroutine dscal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrModule.f90
!
!     LSMR
!
! LSMR   solves Ax = b or min ||Ax - b|| with or without damping,
! using the iterative algorithm of David Fong and Michael Saunders:
!     http://www.stanford.edu/group/SOL/software/lsmr.html
!
! Maintained by
!     David Fong       <clfong@stanford.edu>
!     Michael Saunders <saunders@stanford.edu>
!     Systems Optimization Laboratory (SOL)
!     Stanford University
!     Stanford, CA 94305-4026, USA
!
! 17 Jul 2010: F90 LSMR derived from F90 LSQR and lsqr.m.
! 07 Sep 2010: Local reorthogonalization now works (localSize > 0).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsmrModule

  use  lsmrDataModule,    only : dp
  use  lsmrblasInterface, only : dnrm2, dscal
  use,intrinsic :: iso_c_binding
  implicit none
  private
  public   :: LSMR

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!  subroutine LSMR  ( m, n, Aprod1, Aprod2, b, damp,                    &
!                     atol, btol, conlim, itnlim, localSize, nout,      &
!                     x, istop, itn, normA, condA, normr, normAr, normx )

subroutine LSMR(m,n,val,indices,indptr,b, damp,atol, btol, conlim, itnlim, &
               localSize,x, istop, itn, normA,&
               condA, normr, normAr, normx,verbose,num_threads) bind(C,name="LSMR_csr")
!subroutine LSMR  ( m, n, leniw, lenrw,iw,rw, b, damp,               &
!                     atol, btol, conlim, itnlim, localSize, nout,      &
!                     x, istop, itn, normA, condA, normr, normAr, normx )
   use lsmr_csr_matrix,only      : aprod_parallel
   integer(c_int),value, intent(in) :: m,n,verbose,num_threads
   integer(c_int), INTENT(IN) :: indices(*),indptr(*)
   real(dp),INTENT(IN) :: val(*)

   integer(c_int),value,intent(in)  :: itnlim, localSize
   integer(c_int),intent(out) :: istop, itn
   real(dp), intent(in)  :: b(m)
   real(dp), intent(out) :: x(n)
   real(dp), value,intent(in)  :: atol, btol, conlim, damp
   real(dp), intent(out) :: normA, condA, normr, normAr, normx

    !-------------------------------------------------------------------
    ! LSMR  finds a solution x to the following problems:
    !
    ! 1. Unsymmetric equations:    Solve  A*x = b
    !
    ! 2. Linear least squares:     Solve  A*x = b
    !                              in the least-squares sense
    !
    ! 3. Damped least squares:     Solve  (   A    )*x = ( b )
    !                                     ( damp*I )     ( 0 )
    !                              in the least-squares sense
    !
    ! where A is a matrix with m rows and n columns, b is an m-vector,
    ! and damp is a scalar.  (All quantities are real.)
    ! The matrix A is treated as a linear operator.  It is accessed
    ! by means of subroutine calls with the following purpose:
    !
    ! call Aprod1(m,n,x,y)  must compute y = y + A*x  without altering x.
    ! call Aprod2(m,n,x,y)  must compute x = x + A'*y without altering y.
    !
    ! LSMR uses an iterative method to approximate the solution.
    ! The number of iterations required to reach a certain accuracy
    ! depends strongly on the scaling of the problem.  Poor scaling of
    ! the rows or columns of A should therefore be avoided where
    ! possible.
    !
    ! For example, in problem 1 the solution is unaltered by
    ! row-scaling.  If a row of A is very small or large compared to
    ! the other rows of A, the corresponding row of ( A  b ) should be
    ! scaled up or down.
    !
    ! In problems 1 and 2, the solution x is easily recovered
    ! following column-scaling.  Unless better information is known,
    ! the nonzero columns of A should be scaled so that they all have
    ! the same Euclidean norm (e.g., 1.0).
    !
    ! In problem 3, there is no freedom to re-scale if damp is
    ! nonzero.  However, the value of damp should be assigned only
    ! after attention has been paid to the scaling of A.
    !
    ! The parameter damp is intended to help regularize
    ! ill-conditioned systems, by preventing the true solution from
    ! being very large.  Another aid to regularization is provided by
    ! the parameter condA, which may be used to terminate iterations
    ! before the computed solution becomes very large.
    !
    ! Note that x is not an input parameter.
    ! If some initial estimate x0 is known and if damp = 0,
    ! one could proceed as follows:
    !
    ! 1. Compute a residual vector     r0 = b - A*x0.
    ! 2. Use LSMR to solve the system  A*dx = r0.
    ! 3. Add the correction dx to obtain a final solution x = x0 + dx.
    !
    ! This requires that x0 be available before and after the call
    ! to LSMR.  To judge the benefits, suppose LSMR takes k1 iterations
    ! to solve A*x = b and k2 iterations to solve A*dx = r0.
    ! If x0 is "good", norm(r0) will be smaller than norm(b).
    ! If the same stopping tolerances atol and btol are used for each
    ! system, k1 and k2 will be similar, but the final solution x0 + dx
    ! should be more accurate.  The only way to reduce the total work
    ! is to use a larger stopping tolerance for the second system.
    ! If some value btol is suitable for A*x = b, the larger value
    ! btol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
    !
    ! Preconditioning is another way to reduce the number of iterations.
    ! If it is possible to solve a related system M*x = b efficiently,
    ! where M approximates A in some helpful way
    ! (e.g. M - A has low rank or its elements are small relative to
    ! those of A), LSMR may converge more rapidly on the system
    !       A*M(inverse)*z = b,
    ! after which x can be recovered by solving M*x = z.
    !
    ! NOTE: If A is symmetric, LSMR should not be used!
    ! Alternatives are the symmetric conjugate-gradient method (CG)
    ! and/or SYMMLQ.
    ! SYMMLQ is an implementation of symmetric CG that applies to
    ! any symmetric A and will converge more rapidly than LSMR.
    ! If A is positive definite, there are other implementations of
    ! symmetric CG that require slightly less work per iteration
    ! than SYMMLQ (but will take the same number of iterations).
    !
    !
    ! Notation
    ! --------
    ! The following quantities are used in discussing the subroutine
    ! parameters:
    !
    ! Abar   =  (  A   ),        bbar  =  (b)
    !           (damp*I)                  (0)
    !
    ! r      =  b - A*x,         rbar  =  bbar - Abar*x
    !
    ! normr  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
    !        =  norm( rbar )
    !
    ! eps    =  the relative precision of floating-point arithmetic.
    !           On most machines, eps is about 1.0e-7 and 1.0e-16
    !           in single and double precision respectively.
    !           We expect eps to be about 1e-16 always.
    !
    ! LSMR  minimizes the function normr with respect to x.
    !
    !
    ! Parameters
    ! ----------
    ! m       input      m, the number of rows in A.
    !
    ! n       input      n, the number of columns in A.
    !
    ! Aprod1, Aprod2     See above.
    !
    ! damp    input      The damping parameter for problem 3 above.
    !                    (damp should be 0.0 for problems 1 and 2.)
    !                    If the system A*x = b is incompatible, values
    !                    of damp in the range 0 to sqrt(eps)*norm(A)
    !                    will probably have a negligible effect.
    !                    Larger values of damp will tend to decrease
    !                    the norm of x and reduce the number of 
    !                    iterations required by LSMR.
    !
    !                    The work per iteration and the storage needed
    !                    by LSMR are the same for all values of damp.
    !
    ! b(m)    input      The rhs vector b.
    !
    ! x(n)    output     Returns the computed solution x.
    !
    ! atol    input      An estimate of the relative error in the data
    !                    defining the matrix A.  For example, if A is
    !                    accurate to about 6 digits, set atol = 1.0e-6.
    !
    ! btol    input      An estimate of the relative error in the data
    !                    defining the rhs b.  For example, if b is
    !                    accurate to about 6 digits, set btol = 1.0e-6.
    !
    ! conlim  input      An upper limit on cond(Abar), the apparent
    !                    condition number of the matrix Abar.
    !                    Iterations will be terminated if a computed
    !                    estimate of cond(Abar) exceeds conlim.
    !                    This is intended to prevent certain small or
    !                    zero singular values of A or Abar from
    !                    coming into effect and causing unwanted growth
    !                    in the computed solution.
    !
    !                    conlim and damp may be used separately or
    !                    together to regularize ill-conditioned systems.
    !
    !                    Normally, conlim should be in the range
    !                    1000 to 1/eps.
    !                    Suggested value:
    !                    conlim = 1/(100*eps)  for compatible systems,
    !                    conlim = 1/(10*sqrt(eps)) for least squares.
    !
    !         Note: Any or all of atol, btol, conlim may be set to zero.
    !         The effect will be the same as the values eps, eps, 1/eps.
    !
    ! itnlim  input      An upper limit on the number of iterations.
    !                    Suggested value:
    !                    itnlim = n/2   for well-conditioned systems
    !                                   with clustered singular values,
    !                    itnlim = 4*n   otherwise.
    !
    ! localSize input    No. of vectors for local reorthogonalization.
    !            0       No reorthogonalization is performed.
    !           >0       This many n-vectors "v" (the most recent ones)
    !                    are saved for reorthogonalizing the next v.
    !                    localSize need not be more than min(m,n).
    !                    At most min(m,n) vectors will be allocated.
    !
    ! nout    input      File number for printed output.  If positive,
    !                    a summary will be printed on file nout.
    !
    ! istop   output     An integer giving the reason for termination:
    !
    !            0       x = 0  is the exact solution.
    !                    No iterations were performed.
    !
    !            1       The equations A*x = b are probably compatible.
    !                    Norm(A*x - b) is sufficiently small, given the
    !                    values of atol and btol.
    !
    !            2       damp is zero.  The system A*x = b is probably
    !                    not compatible.  A least-squares solution has
    !                    been obtained that is sufficiently accurate,
    !                    given the value of atol.
    !
    !            3       damp is nonzero.  A damped least-squares
    !                    solution has been obtained that is sufficiently
    !                    accurate, given the value of atol.
    !
    !            4       An estimate of cond(Abar) has exceeded conlim.
    !                    The system A*x = b appears to be ill-conditioned,
    !                    or there could be an error in Aprod1 or Aprod2.
    !
    !            5       The iteration limit itnlim was reached.
    !
    ! itn     output     The number of iterations performed.
    !
    ! normA   output     An estimate of the Frobenius norm of Abar.
    !                    This is the square-root of the sum of squares
    !                    of the elements of Abar.
    !                    If damp is small and the columns of A
    !                    have all been scaled to have length 1.0,
    !                    normA should increase to roughly sqrt(n).
    !                    A radically different value for normA may
    !                    indicate an error in Aprod1 or Aprod2.
    !
    ! condA   output     An estimate of cond(Abar), the condition
    !                    number of Abar.  A very high value of condA
    !                    may again indicate an error in Aprod1 or Aprod2.
    !
    ! normr   output     An estimate of the final value of norm(rbar),
    !                    the function being minimized (see notation
    !                    above).  This will be small if A*x = b has
    !                    a solution.
    !
    ! normAr  output     An estimate of the final value of
    !                    norm( Abar'*rbar ), the norm of
    !                    the residual for the normal equations.
    !                    This should be small in all cases.  (normAr
    !                    will often be smaller than the true value
    !                    computed from the output vector x.)
    !
    ! normx   output     An estimate of norm(x) for the final solution x.
    !
    ! Subroutines and functions used              
    ! ------------------------------
    ! BLAS               dscal, dnrm2
    ! USER               Aprod1, Aprod2
    !
    ! Precision
    ! ---------
    ! The number of iterations required by LSMR will decrease
    ! if the computation is performed in higher precision.
    ! At least 15-digit arithmetic should normally be used.
    ! "real(dp)" declarations should normally be 8-byte words.
    ! If this ever changes, the BLAS routines  dnrm2, dscal
    ! (Lawson, et al., 1979) will also need to be changed.
    !
    !
    ! Reference
    ! ---------
    ! http://www.stanford.edu/group/SOL/software/lsmr.html
    ! ------------------------------------------------------------------
    !
    ! LSMR development:
    ! 21 Sep 2007: Fortran 90 version of LSQR implemented.
    !              Aprod1, Aprod2 implemented via f90 interface.
    ! 17 Jul 2010: LSMR derived from LSQR and lsmr.m.
    ! 07 Sep 2010: Local reorthogonalization now working.
    !-------------------------------------------------------------------

    intrinsic :: abs, dot_product, min, max, sqrt
    ! Local arrays and variables
    real(dp)  :: h(n), hbar(n), u(m), v(n), w(n), localV(n,min(localSize,m,n))
    logical   :: damped, localOrtho, localVQueueFull, prnt, show
    integer(c_int)   :: localOrthoCount, localOrthoLimit, localPointer, localVecs, &
                 pcount, pfreq
    real(dp)  :: alpha, alphabar, alphahat, &
                 beta, betaacute, betacheck, betad, betadd, betahat, &
                 normb, c, cbar, chat, ctildeold, ctol,    &
                 d, maxrbar, minrbar, normA2, &
                 rho, rhobar, rhobarold, rhodold, rhoold, rhotemp, &
                 rhotildeold, rtol, s, sbar, shat, stildeold, &
                 t1, taud, tautildeold, test1, test2, test3, &
                 thetabar, thetanew, thetatilde, thetatildeold, &
                 zeta, zetabar, zetaold

    ! Local constants
    real(dp),         parameter :: zero = 0.0_dp, one = 1.0_dp
    character(len=*), parameter :: enter = ' Enter LSMR.  '
    character(len=*), parameter :: exitt = ' Exit  LSMR.  '
    character(len=*), parameter :: msg(0:7) =                     &
      (/ 'The exact solution is  x = 0                         ', &
         'Ax - b is small enough, given atol, btol             ', &
         'The least-squares solution is good enough, given atol', &
         'The estimate of cond(Abar) has exceeded conlim       ', &
         'Ax - b is small enough for this machine              ', &
         'The LS solution is good enough for this machine      ', &
         'Cond(Abar) seems to be too large for this machine    ', &
         'The iteration limit has been reached                 ' /)
    !-------------------------------------------------------------------


    ! Initialize.
   !open(nout,file="lsmr.log")
    !call set_threads(num_threads)
    localVecs = min(localSize,m,n)
    show = verbose > 0
    if (show) then
       write(*, 1000) enter,m,n,damp,atol,conlim,btol,itnlim,localVecs
    end if
    
    pfreq  = 20           ! print frequency (for repeating the heading)
    pcount = 0            ! print counter
    damped = damp > zero  !

    !-------------------------------------------------------------------
    ! Set up the first vectors u and v for the bidiagonalization.
    ! These satisfy  beta*u = b,  alpha*v = A(transpose)*u.
    !-------------------------------------------------------------------
    u(1:m) = b(1:m)
    v(1:n) = zero
    x(1:n) = zero

    alpha  = zero
    beta   = dnrm2 (m, u, 1)

    if (beta > zero) then
       call dscal (m, (one/beta), u, 1)
    !   call Aprod2(m, n, v, u)          ! v = A'*u
        !call aprod(2,m,n,v,u,leniw,lenrw,iw,rw)
       !call aprod(2,m,n,v,u,lenrw,rw,col,val)
       call aprod_parallel(2,m,n,v,u,val,indices,indptr,num_threads)
       !call aprod(2,spmat,v,u,m,n)
       alpha = dnrm2 (n, v, 1)
    end if

    if (alpha > zero) then
       call dscal (n, (one/alpha), v, 1)
       w = v
    end if

    normAr = alpha*beta
    if (normAr == zero) go to 800

    ! Initialization for local reorthogonalization.

    localOrtho = .false.
    if (localVecs > 0) then
       localPointer    = 1
       localOrtho      = .true.
       localVQueueFull = .false.
       localV(:,1)     = v
    end if

    ! Initialize variables for 1st iteration.

    itn      = 0
    zetabar  = alpha*beta
    alphabar = alpha
    rho      = 1
    rhobar   = 1
    cbar     = 1
    sbar     = 0

    h         = v
    hbar(1:n) = zero
    x(1:n)    = zero

    ! Initialize variables for estimation of ||r||.

    betadd      = beta
    betad       = 0
    rhodold     = 1
    tautildeold = 0
    thetatilde  = 0
    zeta        = 0
    d           = 0

    ! Initialize variables for estimation of ||A|| and cond(A).

    normA2  = alpha**2
    maxrbar = 0_dp
    minrbar = 1e+30_dp

    ! Items for use in stopping rules.
    normb  = beta
    istop  = 0
    ctol   = zero
    if (conlim > zero) ctol = one/conlim
    normr  = beta

    ! Exit if b=0 or A'b = 0.

    normAr = alpha * beta
    if (normAr == 0) then
       if (show) then
          write(*,'(a)') msg(1)
       end if
       return
    end if

    ! Heading for iteration log.

    if (show) then
       if (damped) then
          write(*,1300)
       else
          write(*,1200)
       end if
       test1 = one
       test2 = alpha/beta
       write(*, 1500) itn,x(1),normr,normAr,test1,test2
    end if

    !===================================================================
    ! Main iteration loop.
    !===================================================================
    do
       itn = itn + 1
	  
       !----------------------------------------------------------------
       ! Perform the next step of the bidiagonalization to obtain the
       ! next beta, u, alpha, v.  These satisfy
       !     beta*u = A*v  - alpha*u,
       !    alpha*v = A'*u -  beta*v.
       !----------------------------------------------------------------
       call dscal (m,(- alpha), u, 1)
      ! call Aprod1(m, n, v, u)             ! u = A*v
        !call aprod ( 1,m,n,v,u,leniw,lenrw,iw,rw )
       !call aprod(1,m,n,v,u,lenrw,rw,col,val)  
       call aprod_parallel(1,m,n,v,u,val,indices,indptr,num_threads)    
       beta   = dnrm2 (m, u, 1)

       if (beta > zero) then
          call dscal (m, (one/beta), u, 1)
          if (localOrtho) then    ! Store v into the circular buffer localV.
             call localVEnqueue   ! Store old v for local reorthog'n of new v.
          end if
          call dscal (n, (- beta), v, 1)

          !call Aprod2(m, n, v, u)          ! v = A'*u
          !call aprod ( 2,m,n,v,u,leniw,lenrw,iw,rw )
          !call aprod(2,m,n,v,u,lenrw,rw,col,val)
          call aprod_parallel(2,m,n,v,u,val,indices,indptr,num_threads)
          if (localOrtho) then    ! Perform local reorthogonalization of V.
             call localVOrtho     ! Local-reorthogonalization of new v.
          end if
          alpha  = dnrm2 (n, v, 1)
          if (alpha > zero) then
             call dscal (n, (one/alpha), v, 1)
          end if
       end if
    
       ! At this point, beta = beta_{k+1}, alpha = alpha_{k+1}.
    
       !----------------------------------------------------------------
       ! Construct rotation Qhat_{k,2k+1}.

       alphahat = d2norm(alphabar, damp)
       chat     = alphabar/alphahat
       shat     = damp/alphahat

       ! Use a plane rotation (Q_i) to turn B_i to R_i.

       rhoold   = rho
       rho      = d2norm(alphahat, beta)
       c        = alphahat/rho
       s        = beta/rho
       thetanew = s*alpha
       alphabar = c*alpha

       ! Use a plane rotation (Qbar_i) to turn R_i^T into R_i^bar.
    
       rhobarold = rhobar
       zetaold   = zeta
       thetabar  = sbar*rho
       rhotemp   = cbar*rho
       rhobar    = d2norm(cbar*rho, thetanew)
       cbar      = cbar*rho/rhobar
       sbar      = thetanew/rhobar
       zeta      =   cbar*zetabar
       zetabar   = - sbar*zetabar

       ! Update h, h_hat, x.

       hbar      = h - (thetabar*rho/(rhoold*rhobarold))*hbar
       x         = x + (zeta/(rho*rhobar))*hbar
       h         = v - (thetanew/rho)*h

       ! Estimate ||r||.
    
       ! Apply rotation Qhat_{k,2k+1}.
       betaacute =   chat* betadd
       betacheck = - shat* betadd

       ! Apply rotation Q_{k,k+1}.
       betahat   =   c*betaacute
       betadd    = - s*betaacute

       ! Apply rotation Qtilde_{k-1}.
       ! betad = betad_{k-1} here.

       thetatildeold = thetatilde
       rhotildeold   = d2norm(rhodold, thetabar)
       ctildeold     = rhodold/rhotildeold
       stildeold     = thetabar/rhotildeold
       thetatilde    = stildeold* rhobar
       rhodold       =   ctildeold* rhobar
       betad         = - stildeold*betad + ctildeold*betahat

       ! betad   = betad_k here.
       ! rhodold = rhod_k  here.

       tautildeold   = (zetaold - thetatildeold*tautildeold)/rhotildeold
       taud          = (zeta - thetatilde*tautildeold)/rhodold
       d             = d + betacheck**2
       normr         = sqrt(d + (betad - taud)**2 + betadd**2)
    
       ! Estimate ||A||.
       normA2        = normA2 + beta**2
       normA         = sqrt(normA2)
       normA2        = normA2 + alpha**2
    
       ! Estimate cond(A).
       maxrbar       = max(maxrbar,rhobarold)
       if (itn > 1) then 
          minrbar    = min(minrbar,rhobarold)
       end if
       condA         = max(maxrbar,rhotemp)/min(minrbar,rhotemp)

       !----------------------------------------------------------------
       ! Test for convergence.
       !----------------------------------------------------------------

       ! Compute norms for convergence testing.
       normAr  = abs(zetabar)
       normx   = dnrm2(n, x, 1)

       ! Now use these norms to estimate certain other quantities,
       ! some of which will be small near a solution.

       test1   = normr /normb
       test2   = normAr/(normA*normr)
       test3   =    one/condA
       t1      =  test1/(one + normA*normx/normb)
       rtol    =   btol + atol*normA*normx/normb

       ! The following tests guard against extremely small values of
       ! atol, btol or ctol.  (The user may have set any or all of
       ! the parameters atol, btol, conlim  to 0.)
       ! The effect is equivalent to the normAl tests using
       ! atol = eps,  btol = eps,  conlim = 1/eps.

       if (itn     >= itnlim) istop = 7
       if (one+test3 <=  one) istop = 6
       if (one+test2 <=  one) istop = 5
       if (one+t1    <=  one) istop = 4

       ! Allow for tolerances set by the user.

       if (  test3   <= ctol) istop = 3
       if (  test2   <= atol) istop = 2
       if (  test1   <= rtol) istop = 1

       !----------------------------------------------------------------
       ! See if it is time to print something.
       !----------------------------------------------------------------
       prnt = .false.
       if (show) then
          if (n     <=        40) prnt = .true.
          if (itn   <=        10) prnt = .true.
          if (itn   >= itnlim-10) prnt = .true.
          if (mod(itn,10)  ==  0) prnt = .true.
          if (test3 <=  1.1*ctol) prnt = .true.
          if (test2 <=  1.1*atol) prnt = .true.
          if (test1 <=  1.1*rtol) prnt = .true.
          if (istop /=         0) prnt = .true.

          if (prnt) then        ! Print a line for this iteration
             if (pcount >= pfreq) then  ! Print a heading first
                pcount = 0
                if (damped) then
                   write(*,1300)
                else
                   write(*,1200)
                end if
             end if
             pcount = pcount + 1
             write(*,1500) itn,x(1),normr,normAr,test1,test2,normA,condA
          end if
       end if

       if (istop /= 0) exit
    end do
    !===================================================================
    ! End of iteration loop.
    !===================================================================

    ! Come here if normAr = 0, or if normal exit.

800 if (damped .and. istop==2) istop=3  ! Decide if istop = 2 or 3.
    if (show) then                      ! Print the stopping condition.
       write(*, 2000)                &
            exitt,istop,itn,            &
            exitt,normA,condA,          &
            exitt,normb, normx,         &
            exitt,normr,normAr
       write(*, 3000)                &
            exitt, msg(istop)
    end if

    return

 1000 format(// a, '     Least-squares solution of  Ax = b'       &
      / ' The matrix  A  has', i7, ' rows   and', i7, ' columns'  &
      / ' damp   =', es22.14                                      &
      / ' atol   =', es10.2, 15x,        'conlim =', es10.2       &
      / ' btol   =', es10.2, 15x,        'itnlim =', i10          &
      / ' localSize (no. of vectors for local reorthogonalization) =', i7)
 1200 format(/ "   Itn       x(1)            norm r         A'r   ", &
      ' Compatible    LS      norm A    cond A')
 1300 format(/ "   Itn       x(1)           norm rbar    Abar'rbar", &
      ' Compatible    LS    norm Abar cond Abar')
 1500 format(i6, 2es17.9, 5es10.2)
 2000 format(/ a, 5x, 'istop  =', i2,   15x, 'itn    =', i8      &
      /      a, 5x, 'normA  =', es12.5, 5x, 'condA  =', es12.5   &
      /      a, 5x, 'normb  =', es12.5, 5x, 'normx  =', es12.5   &
      /      a, 5x, 'normr  =', es12.5, 5x, 'normAr =', es12.5)
 3000 format(a, 5x, a)
    
    !close(nout)

  contains

    function d2norm( a, b )

      real(dp)             :: d2norm
      real(dp), intent(in) :: a, b

      !-------------------------------------------------------------------
      ! d2norm returns sqrt( a**2 + b**2 )
      ! with precautions to avoid overflow.
      !
      ! 21 Mar 1990: First version.
      ! 17 Sep 2007: Fortran 90 version.
      ! 24 Oct 2007: User real(dp) instead of compiler option -r8.
      !-------------------------------------------------------------------

      intrinsic            :: abs, sqrt
      real(dp)             :: scale
      real(dp), parameter  :: zero = 0.0_dp

      scale = abs(a) + abs(b)
      if (scale == zero) then
         d2norm = zero
      else
         d2norm = scale*sqrt((a/scale)**2 + (b/scale)**2)
      end if

    end function d2norm

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine localVEnqueue

      ! Store v into the circular buffer localV.

      if (localPointer < localVecs) then
         localPointer = localPointer + 1
      else
         localPointer = 1
         localVQueueFull = .true.
      end if
      localV(:,localPointer) = v

    end subroutine localVEnqueue

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine localVOrtho

      ! Perform local reorthogonalization of current v.

      real(dp)  :: d

      if (localVQueueFull) then
         localOrthoLimit = localVecs
      else
         localOrthoLimit = localPointer
      end if

      do localOrthoCount = 1, localOrthoLimit
         d = dot_product(v,localV(:,localOrthoCount))
         v = v    -    d * localV(:,localOrthoCount)
      end do

    end subroutine localVOrtho

  end subroutine LSMR

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module LSMRmodule
