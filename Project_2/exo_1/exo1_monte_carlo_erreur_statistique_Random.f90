program monte_carlo_erreur_statistique
    implicit none
    real(8) :: r
    real(8) :: sum, sum_2
    real(8), parameter :: l = 0.0d0 ! Orbitale s
    real(8), parameter :: n = 3.0d0 ! ns = 3n
    real(8), parameter :: a_o = 1.0d0 ! Unité atomique
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: Z = 1.0d0 ! Charge du noyau
    integer, parameter :: M = 10000000 ! Nombre d'échantillons
    integer :: i
    real(8) :: V = 61.0d0
    real(8) :: phi_3s
    real(8) :: moy_r
    real(8) :: Random, moy_f, moy_f_2
    real(8) :: erreur_stat, variance
     
    sum = 0.0d0
    phi_3s = 0.0d0
    moy_f = 0.0d0
    moy_r = 0.0d0
    moy_f_2 = 0.d0
    sum_2 = 0.d0
    erreur_stat = 0.0d0
    variance = 0.0d0
    
    ! Génération de nombres aléatoires Random
    do i = 1, M
        Random = rand() * 61.0d0
        r = Random
        phi_3s = (1.0d0/ sqrt (4.0d0*pi)) * (2.0d0 / (81.0d0 * sqrt (3.0d0))) * &
        (27.0d0 - (18.0d0 * r) + 2.0d0*(r**2.0d0)) * &
        exp(-(r/3.0d0))
        moy_r = 4.0d0 * pi * (((phi_3s)**2.0d0) * (r**3.0d0)) !Calculs de la distance moyenne de r_3s
        !print *, moy_r , i
        sum = sum + moy_r
        sum_2 = sum_2 + (moy_r)**2
    end do
        !Calcul de la moyenne pondérée
        moy_f = (1.0d0/ M) * sum 
        moy_f_2 = (1.0d0/ M) * sum_2

        !L'erreur statistique et la variance
        erreur_stat = 61.0d0 * sqrt (((moy_f_2)-((moy_f)**2))/M)
        variance = (moy_f_2)-((moy_f)**2)

        print *, "Moyenne (moy_f) :", moy_f**2
        print *, "Moyenne des carrés (moy_f_2) :", moy_f_2
        print *, "Variance :", variance
        print *, "Erreur statistique :", erreur_stat

    end program monte_carlo_erreur_statistique