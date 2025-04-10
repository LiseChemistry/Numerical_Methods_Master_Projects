function potentiel(x)
   real*8 :: D, Rast, alpha, c1, c2, c3
   integer, parameter :: n = 8 !Le nombre des points 
   real*8 :: potentiel
   real*8, intent (in) :: x
 
     D = 3.8886d0
     Rast = 2.2818d0
     alpha = 3.3522498d0
     c1 = 3.6445906d0
     c2 = 3.9281238d0
     c3 = 2.0986689d0
 
    potentiel = -D*(1.0d0+c1*(x-Rast)+c2*(x-Rast)**2+c3*(x-Rast)**3)*exp(-alpha*(x-Rast))
 
end function potentiel 

program potentiel_derivee_premiere
    implicit none
    real*8 :: delta_x = 0.15d0 !La largeur de la grille
    real*8 :: f_derivee_1_d_1, f_derivee_1_d_2, f_derivee_1_d_3, f_derivee_1_d_4, f_derivee_1_d_5
    real*8 :: Rast = 2.2818d0
    real*8 :: potentiel

    f_derivee_1_d_1 = (potentiel(Rast + delta_x) * (1.0d0/2.0d0) - (1.0d0/2.0d0) * potentiel(Rast - delta_x))* (1.0d0/delta_x)

    f_derivee_1_d_2 = (potentiel(Rast + delta_x) * (2.0d0/3.0d0) - (2.0d0/3.0d0) * potentiel(Rast - delta_x) &
    - (1.0d0/12.0d0) * potentiel(Rast + 2.0d0 * delta_x) + (1.0d0/12.0d0) * potentiel(Rast - 2.0d0* delta_x)) * (1.0d0/delta_x)

    f_derivee_1_d_3 = ((3.0d0/4.0d0) * potentiel(Rast + delta_x) - (3.0d0/4.0d0) * potentiel(Rast - delta_x) &
    - (3.0d0/20.0d0) * potentiel(Rast + 2.0d0 * delta_x) + (3.0d0/20.0d0) * potentiel(Rast - 2.0d0 * delta_x) &
    + (1.0d0/60.0d0) * potentiel(Rast + 3.0d0 * delta_x) - (1.0d0/60.0d0) * potentiel(Rast - 3.0d0 * delta_x)) * (1.0d0/delta_x)

    f_derivee_1_d_4 = ((4.0d0/5.0d0) * potentiel(Rast + delta_x) - (4.0d0/5.0d0) * potentiel(Rast - delta_x) &
    - (1.0d0/5.0d0) * potentiel(Rast + 2.0d0 * delta_x) + (1.0d0/5.0d0) * potentiel(Rast - 2.0d0 * delta_x) &
    + (4.0d0/105.0d0) * potentiel(Rast + 3.0d0 * delta_x) - (4.0d0/105.0d0) * potentiel(Rast - 3.0d0 * delta_x) &
    - (1.0d0/280.0d0) * potentiel(Rast + 4.0d0 * delta_x) + (1.0d0/280.0d0) * potentiel(Rast - 4.0d0 * delta_x)) * (1.0d0/delta_x)

    f_derivee_1_d_5 = ((5.0d0/6.0d0) * potentiel(Rast + delta_x) - (5.0d0/6.0d0) * potentiel(Rast - delta_x) &
    - (5.0d0/21.0d0) * potentiel(Rast + 2.0d0 * delta_x) +  (5.0d0/21.0d0) * potentiel(Rast - 2.0d0 * delta_x) &
     + (5.0d0/84.0d0) * potentiel(Rast + 3.0d0 * delta_x) - (5.0d0/84.0d0) * potentiel(Rast - 3.0d0 * delta_x) &
     - (5.0d0/504.0d0) * potentiel(Rast + 4.0d0 * delta_x) + (5.0d0/504.0d0) * potentiel(Rast - 4.0d0 * delta_x) &
     + (1.0d0/1260.0d0) * potentiel(Rast + 5.0d0*delta_x) - (1.0d0/1260.0d0) * potentiel(Rast - 5.0d0 * delta_x))* (1.0d0/delta_x)
    
   print *, "La dérivée première (d_1) :", f_derivee_1_d_1
   print *, "La dérivée première (d_2) :", f_derivee_1_d_2
   print *, "La dérivée première (d_3) :", f_derivee_1_d_3
   print *, "La dérivée première (d_4) :", f_derivee_1_d_4
   print *, "La dérivée première (d_5) :", f_derivee_1_d_5

  end program potentiel_derivee_premiere