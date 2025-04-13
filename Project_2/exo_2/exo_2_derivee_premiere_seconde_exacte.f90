program potentiel_premiere_exacte
    implicit none
    real*8 :: D, Rast, alpha, c1, c2, c3
    integer, parameter :: n = 8
    !real*8,dimension(n) :: R = [1.8d0, 2.1d0, 2.4d0, 2.8d0, 3.3d0, 4.0d0, 5.0d0, 7.0d0]
    real*8 :: f, x, dF1, dF2

    D = 3.8886d0
    Rast = 2.2818d0
    alpha = 3.3522498d0
    c1 = 3.6445906d0
    c2 = 3.9281238d0
    c3 = 2.0986689d0
    x = Rast

       f = -D*(1.0d0+c1*(x-Rast)+c2*(x-Rast)**2+c3*(x-Rast)**3)*exp(-alpha*(x-Rast))

      dF1=-D*((c1+2*c2*(x-Rast)+3*c3*(x-Rast)**2)*DEXP(-alpha*(x-Rast)) & 
       -alpha*(1+c1*(x-Rast)+c2*(x-Rast)**2+c3*(x-Rast)**3)*DEXP(-alpha*(x-Rast)))

      dF2 = -D*((2*c2 + 6*c3*(x-Rast))*DEXP(-alpha*(x-Rast)) &
      -2*alpha*(c1 + 2*c2*(x-Rast) + 3*c3*((x-Rast)**2))*DEXP(-alpha*(x-Rast)) &
      + (alpha**2)*(1 + c1*(x-Rast) + c2*((x-Rast)**2) + c3*((x-Rast)**3))*DEXP(-alpha*(x-Rast)))
      
  print *, dF1, dF2 

  end program potentiel_premiere_exacte
