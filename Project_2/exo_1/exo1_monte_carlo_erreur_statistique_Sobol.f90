!######################################################################################################
!
!
! sobseq : routine numerical recipes pour la générateur de Sobol
!
!######################################################################################################

subroutine sobseq(n,x)
    integer n,MAXBIT,MAXDIM
    real(kind=8) x(*)
    parameter (MAXBIT=30,MAXDIM=6)
    
    !vector x(1..n) the next values from n of these sequences. (n must not be changed between
    !initializations.)
    integer i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT), iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
    real(kind=8) fac
    save ip,mdeg,ix,iv,in,fac
    equivalence (iv,iu) !To allow both 1D and 2D addressing.
    data ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
    data iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
        if (n.lt.0) then !Initialize, dont return a vector.
            do k=1,MAXDIM
                ix(k)=0
            enddo
            in=0
            if(iv(1).ne.1)return
            fac=1./2.**MAXBIT
            do  k=1,MAXDIM
                do  j=1,mdeg(k) !Stored values only require normalization.
                    iu(k,j)=iu(k,j)*2**(MAXBIT-j)
                enddo
                do j=mdeg(k)+1,MAXBIT !Use the recurrence to get other values.
                    ipp=ip(k)
                    i=iu(k,j-mdeg(k))
                    i=ieor(i,i/2**mdeg(k))
                    do l=mdeg(k)-1,1,-1
                        if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
                        ipp=ipp/2
                    enddo
                    iu(k,j)=i
                enddo
            enddo
        else !Calculate the next vector in the sequence.
            im=in
            do j=1,MAXBIT !Find the rightmost zero bit.
                if(iand(im,1).eq.0)goto 1
                im=im/2
            enddo
            !pause !MAXBIT too small in sobseq
           ! write(*,*) "MAXBIT too small in sobseq"
            stop
    1       im=(j-1)*MAXDIM
            do  k=1,min(n,MAXDIM) !XOR the appropriate direction number into each component of the vector and convert to a floating number.
                ix(k)=ieor(ix(k),iv(im+k))
                x(k)=ix(k)*fac
            enddo
            in=in+1
    endif
    return
    end

program monte_carlo_erreur_statistique
    implicit none
    real(8) :: r
    real(kind=8) :: sobol(1)
    integer(kind=8) :: i
    real(8) :: sum, sum_2
    real(8), parameter :: l = 0.0d0 ! Orbitale s
    real(8), parameter :: n = 3.0d0 ! ns = 3n
    real(8), parameter :: a_o = 1.0d0 ! Unité atomique
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: Z = 1.0d0 ! Charge du noyau
    integer, parameter :: M = 10000 ! Nombre d'échantillons
    real(8) :: upper_bound = 61.0d0
    real(8) :: phi_3s
    real(8) :: moy_r
    real(8) :: moy_f, moy_f_2
    real(8) :: erreur_stat, variance
     
    sum = 0.0d0
    phi_3s = 0.0d0
    moy_f = 0.0d0
    moy_r = 0.0d0
    moy_f_2 = 0.d0
    sum_2 = 0.d0
    erreur_stat = 0.0d0
    variance = 0.0d0
    
    call sobseq(-1, sobol)
     
    ! Génération de nombres aléatoires Random
    do i = 1, M
        call sobseq(1, sobol)
        r = sobol(1) * upper_bound
        phi_3s = (1.0d0/ sqrt (4.0d0*pi)) * (2.0d0 / (81.0d0 * sqrt (3.0d0))) * &
        (27.0d0 - (18.0d0 * r) + 2.0d0*(r**2.0d0)) * &
        exp(-(r/3.0d0))
        moy_r = 4.0d0 * pi * (((phi_3s)**2.0d0) * (r**3.0d0)) !Calculs de la distance moyenne de r_3s
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