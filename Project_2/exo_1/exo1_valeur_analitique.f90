program valeur_analytique
    implicit none
    real(8) :: sum, mean_distance
    real(8) :: l = 0.d0 ! orbital s
    real(8) :: n = 3.d0 ! ns = 3n
    real(8) :: a_o = 1.d0 ! atomic units
    real(8) :: Z = 1.d0 ! Charge
   
    sum = 1.d0 + ((1.d0/ 2.d0) * (1.d0 - ((l*(l+1.d0))/(n**2))))

    mean_distance = (((n**2) * a_o) / Z) * sum! Multiplication

    print *, "The mean  distance ⟨r⟩ 3s is :", mean_distance

end program valeur_analytique
