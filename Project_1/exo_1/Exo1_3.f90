      Subroutine ratint(xa,ya,n,x,y,dy)

!     !     xa : array of x-values of the grid
!     ya : corresponding function values
!     n  : number of points
!     x  : abscissa for interpolation
!     y  : interpolated value at x
!     dy : estimated interpolation error

      Implicit Real*8(a-h,o-z)
      PARAMETER (NMAX=50,TINY=1.d-25)
      Dimension xa(n),ya(n),c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0d0)stop 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

program Polynome_rationnel
  implicit none
  real*8 :: D, Rast, alpha, c1, c2, c3
  integer, parameter :: n = 8
  real*8,dimension(n) :: R = [1.8d0, 2.1d0, 2.4d0, 2.8d0, 3.3d0, 4.0d0, 5.0d0, 7.0d0]
  real*8,dimension(n) :: f
  integer :: i, j
  real*8 :: dy, x, y

  D = 3.8886d0
  Rast = 2.2818d0
  alpha = 3.3522498d0
  c1 = 3.6445906d0
  c2 = 3.9281238d0
  c3 = 2.0986689d0

  do i = 1, n
    f(i)= -D*(1.0+c1*(R(i)-Rast)+c2*(R(i)-Rast)**2+c3*(R(i)-Rast)**3)*exp(-alpha*(R(i)-Rast))
  end do

    ! Open the file txt to save the results
  open(unit = 10, file = 'potentiel.dat', status = 'replace')
  do j = 1, 100 !nombre de points
     x = 0.1d0*j
    call ratint(R,f,n,x,y,dy)
    write(10,*) x, y !dy :  estimation of error
  end do
  close(10)

 call system ('gnuplot -p potentiel.plt')

end program Polynome_rationnel
