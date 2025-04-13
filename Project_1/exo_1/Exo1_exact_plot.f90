program potentiel
  implicit none
  real*8 :: D, Rast, alpha, c1, c2, c3
  integer, parameter :: n = 8
  !real*8,dimension(n) :: R = [1.8d0, 2.1d0, 2.4d0, 2.8d0, 3.3d0, 4.0d0, 5.0d0, 7.0d0]
  real*8 :: f, x
  integer :: i

  D = 3.8886d0
  Rast = 2.2818d0
  alpha = 3.3522498d0
  c1 = 3.6445906d0
  c2 = 3.9281238d0
  c3 = 2.0986689d0

  open(unit = 10, file = 'potentiel.dat', status = 'replace')

  do i = 1, 100
     x = i*0.1d0
     f = -D*(1.0d0+c1*(x-Rast)+c2*(x-Rast)**2+c3*(x-Rast)**3)*exp(-alpha*(x-Rast))
     write(10,*) x, f
  end do

  close(10)
  call system ('gnuplot -p potentiel.plt')

end program potentiel