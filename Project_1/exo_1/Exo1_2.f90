program Polynome_en_puissances_de_1_divise_sur_R
  implicit none
  real*8 :: D, Rast, alpha, c1, c2, c3
  real*8 :: P, resultat, Rq
  integer, parameter :: n = 8
  real*8,dimension(n) :: R = [1.8d0, 2.1d0, 2.4d0, 2.8d0, 3.3d0, 4.0d0, 5.0d0, 7.0d0]
  real*8,dimension(n) :: f
  integer :: i, j

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
  do j = 1, n
    write(10,*) f(j), R(j)
  end do
  close(10)

!INTERPOLATION
  open(unit=11, file = 'myresults.dat', status = 'replace')
  Rq = 0.1d0
  do while  (Rq<=10)
  resultat = 0 !initialisation

  do i = 1, n
     P = 1.0d0
     do j = 1, n
        if (i/=j) then
        P = P*((1/Rq)-(1/R(j)))/((1/R(i))-(1/R(j))) 
        end if
     end do
     resultat = resultat + f(i) * P
  end do
  write(11,*) Rq, resultat
  Rq = Rq + 0.1
  end do

 call system ('gnuplot -p potentiel.plt')
 close(10)

end program Polynome_en_puissances_de_1_divise_sur_R
